#!/usr/bin/env python

import sys
import json
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

# The target (v6) pan-graph table schemas are declared inline on purpose. A
# migration must pin the schema exactly as it was at the time of the version
# bump and must NOT import anvio.tables -- that module always tracks the current
# head and would silently drift as the schema evolves in later versions,
# breaking this script retroactively.
pan_graph_nodes_table_name      = 'pan_graph_nodes'
pan_graph_nodes_table_structure = ['node_id', 'node_type', 'region_id', 'gene_cluster_id', 'gene_calls_json', 'synteny_position_json', 'node_x', 'node_y', 'alignment_summary', 'component_id']
pan_graph_nodes_table_types     = [  'str'  ,    'str'   ,    'str'   ,       'str'      ,       'str'      ,          'str'         , 'numeric', 'numeric',       'str'       ,  'numeric'    ]

pan_graph_edges_table_name      = 'pan_graph_edges'
pan_graph_edges_table_structure = ['edge_id', 'source', 'target', 'weight' , 'genomes_json']
pan_graph_edges_table_types     = [  'str'  ,   'str' ,   'str' , 'numeric',     'str'      ]

pan_graph_regions_table_name      = 'pan_graph_regions'
pan_graph_regions_table_structure = ['component_id', 'region_id', 'region_type', 'x_min', 'x_max', 'num_synteny_gene_clusters', 'num_gene_clusters', 'num_gene_calls', 'max_expansion', 'min_expansion', 'complexity', 'complexity_normalized', 'diversity', 'diversity_normalized', 'weight', 'weight_normalized', 'composite_variability_score', 'complexity_mm_scaled', 'diversity_mm_scaled', 'expansion_mm_scaled', 'weight_mm_scaled']
pan_graph_regions_table_types     = [  'numeric'   ,    'str'   ,     'str'    , 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric']

# Every pan-graph produced by the master line was collapsed to its single
# largest connected component at generation time, so every v5 database is
# single-component -- component_id is 0 throughout.
SINGLE_COMPONENT_ID = 0


def _column_names(pan_graph_db, table_name):
    """Return the column names of an existing SQLite table without reading rows."""
    return [d[0] for d in pan_graph_db._exec(
        f"SELECT * FROM {table_name} LIMIT 0").description]


def _norm_region_id(rid):
    """v5 stored region_id as a number; v6 stores a plain per-component string
    id. Normalise node.region_id and region.region_id through the same function
    so the (component_id, region_id) pairs still join after migration."""
    if rid is None or rid == '':
        return rid
    try:
        return str(int(float(rid)))
    except (ValueError, TypeError):
        return str(rid)


def _migrate_nodes(pan_graph_db):
    """v6 adds a `component_id` column to pan_graph_nodes and stores region_id
    as a plain string. Every v5 db is single-component, so component_id is 0."""
    old_columns = _column_names(pan_graph_db, pan_graph_nodes_table_name)
    if 'component_id' in old_columns:
        # Already migrated (e.g. a manual repair) -- skip silently.
        return

    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_nodes_table_name}").fetchall()
    new_rows = []
    for row in rows:
        carryover = {c: v for c, v in zip(old_columns, row)}
        carryover['component_id'] = SINGLE_COMPONENT_ID
        carryover['region_id'] = _norm_region_id(carryover.get('region_id'))
        new_rows.append([carryover.get(col) for col in pan_graph_nodes_table_structure])

    pan_graph_db._exec(f"DROP TABLE {pan_graph_nodes_table_name}")
    pan_graph_db.create_table(
        pan_graph_nodes_table_name,
        pan_graph_nodes_table_structure,
        pan_graph_nodes_table_types,
    )
    placeholders = ','.join(['?'] * len(pan_graph_nodes_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_nodes_table_name} VALUES ({placeholders})",
        new_rows,
    )


def _edge_orientations(directions_raw):
    """The v5 `directions` column is a JSON object mapping each genome on the
    edge to its traversal orientation, 'L' (reverse) or 'R' (forward). Return
    the set of orientations present (upper-cased); an empty set if the payload
    is missing or unparseable."""
    if not directions_raw:
        return set()
    try:
        parsed = json.loads(directions_raw)
    except (ValueError, TypeError):
        return set()
    if isinstance(parsed, dict):
        values = parsed.values()
    elif isinstance(parsed, (list, tuple)):
        values = parsed
    else:
        values = [parsed]
    return {str(v).upper() for v in values}


def _migrate_edges(pan_graph_db):
    """v6 has no per-edge orientation: it replaces the v5 `directions` column
    with `genomes_json` (the set of genomes crossing the edge) and treats every
    edge as forward source->target.

    Reversed edges -- those the v5 renderer drew reversed, i.e. with no forward
    ('R') orientation -- cannot be represented in the v6 model and are DROPPED
    here. Surviving (forward / mixed) edges get an empty genomes_json because a
    v5 db never stored per-edge genome membership; edge colouring by genome will
    be blank until the pan-graph is regenerated.

    Returns the number of reversed edges removed."""
    old_columns = _column_names(pan_graph_db, pan_graph_edges_table_name)
    if 'genomes_json' in old_columns and 'directions' not in old_columns:
        # Already on the target shape.
        return 0

    dir_idx = old_columns.index('directions') if 'directions' in old_columns else None
    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_edges_table_name}").fetchall()

    new_rows = []
    removed = 0
    for row in rows:
        if dir_idx is not None:
            orientations = _edge_orientations(row[dir_idx])
            # A reversed edge has orientations but none of them forward ('R').
            if orientations and 'R' not in orientations:
                removed += 1
                continue
        carryover = {c: v for c, v in zip(old_columns, row)}
        carryover['genomes_json'] = '[]'
        # Only emit the v6 structure, dropping the legacy 'directions' column.
        new_rows.append([carryover.get(col) for col in pan_graph_edges_table_structure])

    pan_graph_db._exec(f"DROP TABLE {pan_graph_edges_table_name}")
    pan_graph_db.create_table(
        pan_graph_edges_table_name,
        pan_graph_edges_table_structure,
        pan_graph_edges_table_types,
    )
    placeholders = ','.join(['?'] * len(pan_graph_edges_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_edges_table_name} VALUES ({placeholders})",
        new_rows,
    )
    return removed


def _migrate_regions(pan_graph_db):
    """v6 prepends a `component_id` column to pan_graph_regions and stores
    region_id as a plain string. Single-component v5 dbs get component_id 0."""
    old_columns = _column_names(pan_graph_db, pan_graph_regions_table_name)
    if 'component_id' in old_columns:
        return

    rows = pan_graph_db._exec(f"SELECT * FROM {pan_graph_regions_table_name}").fetchall()
    new_rows = []
    for row in rows:
        carryover = {c: v for c, v in zip(old_columns, row)}
        carryover['component_id'] = SINGLE_COMPONENT_ID
        carryover['region_id'] = _norm_region_id(carryover.get('region_id'))
        new_rows.append([carryover.get(col) for col in pan_graph_regions_table_structure])

    pan_graph_db._exec(f"DROP TABLE {pan_graph_regions_table_name}")
    pan_graph_db.create_table(
        pan_graph_regions_table_name,
        pan_graph_regions_table_structure,
        pan_graph_regions_table_types,
    )
    placeholders = ','.join(['?'] * len(pan_graph_regions_table_structure))
    pan_graph_db._exec_many(
        f"INSERT INTO {pan_graph_regions_table_name} VALUES ({placeholders})",
        new_rows,
    )


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

    progress.update("Rewriting pan_graph_nodes (adding component_id = 0) ...")
    _migrate_nodes(pan_graph_db)

    progress.update("Rewriting pan_graph_edges (dropping reversed edges; directions -> genomes_json) ...")
    reversed_edges_removed = _migrate_edges(pan_graph_db)

    progress.update("Rewriting pan_graph_regions (adding component_id = 0) ...")
    _migrate_regions(pan_graph_db)

    progress.update("Updating version")
    pan_graph_db.remove_meta_key_value_pair('version')
    pan_graph_db.set_version(next_version)

    progress.update("Committing changes")
    pan_graph_db.disconnect()

    progress.end()

    run.warning(
        f"Your pan-graph database is now version {next_version}, but please read this carefully. "
        f"The whole pan-graph has been assigned to a single component (component_id = 0): every "
        f"pan-graph produced by the older (master) engine was already reduced to its single largest "
        f"connected component, so this is faithful. However, {reversed_edges_removed} reversed "
        f"edge(s) (edges the old engine traversed in the reverse orientation) were REMOVED during "
        f"migration, because the new engine has no notion of edge orientation and cannot represent "
        f"them. Surviving edges also lost their per-edge genome membership (it was never stored in "
        f"v5), so genome-based edge colouring will be blank. Because the two engine versions differ "
        f"substantially, this migrated database is a best-effort, structural upgrade only -- anvi'o "
        f"strongly recommends that you re-run the entire pan-graph generation on your data rather "
        f"than relying on this migrated database for analysis.",
        header=f"PAN-GRAPH DB MIGRATED TO v{next_version} -- PLEASE RE-RUN THE PAN-GRAPH",
        lc='yellow'
    )


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
