#!/usr/bin/env python

import sys
import json
import math
import argparse
from collections import Counter

import numpy as np
from scipy.sparse import csr_matrix

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

# The target (v8) pan_graph_regions schema and the metric settings are declared
# inline on purpose. A migration must pin the schema and the computation exactly
# as they were at the time of the version bump and must NOT import anvio.tables or
# the graph code -- those always track the current head and would silently drift
# as they evolve in later versions, breaking this script retroactively.
#
# In v8 every region metric comes from its region alone; nothing is min-max
# scaled over the dataset any more:
#   * complexity is the scaled complexity (the former complexity_normalized);
#     its raw count is num_covering_syngcs.
#   * weight is the share of the genomes that carry the region (the former
#     weight_normalized); the raw count is num_genomes_present.
#   * diversity is the graded haplotype diversity: the expected Jaccard distance
#     between the sets of the region's SynGCs (alleles) that two random genomes
#     carry, a genome without the region carrying the empty set.
#     diversity_normalized is dropped.
#   * expansion_log is new: log(1 + max_expansion) / log(1 + 100), capped at 1.
#   * num_genomes is new: the number of genomes the region was scored against,
#     which the former weight / weight_normalized gives back.
#   * composite_variability_score is replaced by composite_variability_score_geometric
#     and composite_variability_score_arithmetic, the geometric and the arithmetic
#     mean of complexity, diversity, expansion_log and weight; the four *_mm_scaled
#     columns are dropped.
# The old diversity cannot be converted, and v7 databases made by earlier code
# hold earlier definitions of the other metrics (backbone expansions stored as 0,
# complexity counted one lower), so every metric of every region is recomputed
# from its nodes in pan_graph_nodes. The regions themselves (ids, types and the
# nodes they hold) are kept.
pan_graph_nodes_table_name = 'pan_graph_nodes'

pan_graph_regions_table_name      = 'pan_graph_regions'
pan_graph_regions_table_structure = ['component_id', 'region_id', 'region_type', 'x_min', 'x_max', 'num_synteny_gene_clusters', 'num_gene_clusters', 'num_gene_calls', 'num_genomes', 'max_expansion', 'min_expansion', 'expansion_log', 'num_covering_syngcs', 'complexity', 'diversity', 'num_genomes_present', 'weight', 'composite_variability_score_geometric', 'composite_variability_score_arithmetic']
pan_graph_regions_table_types     = [    'str'     ,    'str'   ,     'str'    , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']

LOG_SIZE_CAP = 100


def _column_names(pan_graph_db, table_name):
    """Return the column names of an existing SQLite table without reading rows."""
    return [d[0] for d in pan_graph_db._exec(
        f"SELECT * FROM {table_name} LIMIT 0").description]


def _nodes_per_region(pan_graph_db):
    """Return ``{(component_id, region_id): [(node_id, node_x, genomes)]}`` from
    pan_graph_nodes; the genomes of a node are the keys of its gene_calls_json."""
    regions = {}
    for node_id, component_id, region_id, node_x, gene_calls_json in pan_graph_db._exec(
            f"SELECT node_id, component_id, region_id, node_x, gene_calls_json FROM {pan_graph_nodes_table_name}").fetchall():
        if region_id is None or region_id == '':
            continue
        regions.setdefault((str(component_id), str(region_id)), []).append((node_id, node_x, list(json.loads(gene_calls_json))))
    return regions


def _summarize_region(nodes, num_genomes):
    """The v8 metrics of one region from its ``nodes``, scored against
    ``num_genomes`` genomes with the defaults of the time: the graded haplotype
    diversity with the empty allele, LOG_SIZE_CAP, and both the geometric and the
    arithmetic mean for the composite variability score."""
    node_ids, region_x, genomes_sets = zip(*nodes)

    # the allele of every genome: the set of the region's SynGCs it carries
    alleles_of_genomes = {}
    for node, genomes in zip(node_ids, genomes_sets):
        for genome in genomes:
            alleles_of_genomes.setdefault(genome, set()).add(node)
    genes_per_genome = [len(allele) for allele in alleles_of_genomes.values()]
    num_genomes_present = len(alleles_of_genomes)
    weight = num_genomes_present / num_genomes

    max_expansion = max(genes_per_genome)
    min_expansion = min(genes_per_genome) if num_genomes_present == num_genomes else 0
    expansion_log = min(math.log1p(max_expansion) / math.log1p(LOG_SIZE_CAP), 1.0)

    covered, num_covering_syngcs = set(), 0
    for genomes in sorted(map(set, genomes_sets), key=len):
        if not genomes <= covered:
            num_covering_syngcs += 1
            covered |= genomes
    complexity = (num_covering_syngcs - 1) / (num_genomes - 1) if num_genomes > 1 else 0

    # 1 - the mean Jaccard similarity of the alleles of all pairs of the N genomes;
    # two genomes without the region do not differ
    N = num_genomes
    diversity = 0.0
    if N > 1:
        genomes_per_allele = Counter(map(frozenset, alleles_of_genomes.values()))
        allele_counts = np.fromiter(genomes_per_allele.values(), dtype=float)
        alleles = list(genomes_per_allele)
        size = np.array([len(allele) for allele in alleles])
        node_index = {node: i for i, node in enumerate(node_ids)}
        rows = np.repeat(np.arange(len(alleles)), size)
        cols = [node_index[node] for allele in alleles for node in allele]
        A = csr_matrix((np.ones(len(cols)), (rows, cols)), shape=(len(alleles), len(node_ids)))
        shared = (A @ A.T).tocoo()
        pair = shared.row < shared.col
        a, b, n = shared.row[pair], shared.col[pair], shared.data[pair]
        # the pairs with the same allele (the empty one included), plus the Jaccard
        # index of every two different alleles
        similarity = ((allele_counts * (allele_counts - 1) / 2).sum()
                      + (N - num_genomes_present) * (N - num_genomes_present - 1) / 2)
        similarity += (allele_counts[a] * allele_counts[b] * n / (size[a] + size[b] - n)).sum()
        diversity = float(1 - similarity / (N * (N - 1) / 2))

    return {'x_min': min(region_x),
            'x_max': max(region_x),
            'num_synteny_gene_clusters': len(node_ids),
            'num_gene_clusters': len({node.rsplit('_', 1)[0] for node in node_ids}),
            'num_gene_calls': sum(genes_per_genome),
            'num_genomes': num_genomes,
            'max_expansion': max_expansion,
            'min_expansion': min_expansion,
            'expansion_log': expansion_log,
            'num_covering_syngcs': num_covering_syngcs,
            'complexity': complexity,
            'diversity': diversity,
            'num_genomes_present': num_genomes_present,
            'weight': weight,
            'composite_variability_score_geometric': (complexity * diversity * expansion_log * weight) ** (1 / 4),
            'composite_variability_score_arithmetic': (complexity + diversity + expansion_log + weight) / 4}


def _migrate_regions(pan_graph_db):
    """Recreate pan_graph_regions with the pinned v8 schema: keep the ids and
    types of the regions and recompute their metrics from their nodes. Returns
    the number of regions, and of regions whose nodes were not found (their
    metrics are left empty)."""
    old_columns = _column_names(pan_graph_db, pan_graph_regions_table_name)
    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_regions_table_name}").fetchall()
    nodes_per_region = _nodes_per_region(pan_graph_db)

    new_rows, num_without_nodes = [], 0
    for row in rows:
        old = dict(zip(old_columns, row))
        new = {col: old[col] for col in ('component_id', 'region_id', 'region_type')}

        nodes = nodes_per_region.get((str(old['component_id']), str(old['region_id'])))
        if nodes:
            # the former weight_normalized was weight / num_genomes
            new |= _summarize_region(nodes, round(old['weight'] / old['weight_normalized']))
        else:
            num_without_nodes += 1

        new_rows.append([new.get(col) for col in pan_graph_regions_table_structure])

    pan_graph_db._exec(f"DROP TABLE {pan_graph_regions_table_name}")
    pan_graph_db.create_table(pan_graph_regions_table_name, pan_graph_regions_table_structure, pan_graph_regions_table_types)
    placeholders = ','.join(['?'] * len(pan_graph_regions_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_regions_table_name} VALUES ({placeholders})", new_rows)

    return len(new_rows), num_without_nodes


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_graph_db_info = dbinfo.PanGraphDBInfo(db_path)
    if str(pan_graph_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided pan-graph database is {pan_graph_db_info.version}, "
            f"not the required version, {current_version}, so this script cannot upgrade the database.")

    pan_graph_db = pan_graph_db_info.load_db()

    progress.new("Migrating")

    progress.update("Recomputing the region metrics in pan_graph_regions ...")
    num_regions, num_without_nodes = _migrate_regions(pan_graph_db)

    progress.update("Updating version")
    pan_graph_db.remove_meta_key_value_pair('version')
    pan_graph_db.set_version(next_version)

    progress.update("Committing changes")
    pan_graph_db.disconnect()

    progress.end()

    if num_without_nodes:
        run.warning(f"Anvi'o could not find the nodes of {num_without_nodes} of the {num_regions} region(s) in "
                    f"pan_graph_nodes, so their metrics are left empty. "
                    f"Running anvi-pan-graph again will compute them.")

    run.info_single(
        f"Your pan-graph database is now version {next_version}. Anvi'o recomputed the metrics of the "
        f"{num_regions} region(s) from their nodes, and each now comes from its region alone, without "
        f"scaling over the dataset: complexity and weight are the scaled values (their raw counts are "
        f"num_covering_syngcs and num_genomes_present), diversity is the graded haplotype diversity, "
        f"expansion_log is new, and the composite variability score is stored twice, as the geometric "
        f"(composite_variability_score_geometric) and the arithmetic (composite_variability_score_arithmetic) "
        f"mean of the four. "
        f"The regions themselves and the rest of the graph did not change.",
        nl_before=1, nl_after=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f"A simple script to upgrade an anvi'o pan-graph database from version {current_version} to version {next_version}"
    )
    parser.add_argument('pan_graph_db', metavar='PAN_GRAPH_DB',
                        help=f"An anvi'o pan-graph database of version {current_version}")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_graph_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
