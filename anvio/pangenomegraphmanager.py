# pylint: disable=line-too-long
"""
Classes for pangenome graph operations.

This module provides PangenomeGraphManager, a directed graph (nx.DiGraph)
representation of a pangenome where nodes are synteny gene clusters and
edges capture gene order relationships across genomes. It supports graph
construction, connectivity checks, region classification (backbone vs.
variable), and graph-based genome distance calculations.
"""

import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
from collections import Counter

from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree

import anvio
import anvio.terminal as terminal
import anvio.clustering as clustering

from anvio.errors import ConfigError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alexander Henoch"
__email__ = "ahenoch@outlook.de"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

additional_param_sets_for_sequence_search = {'diamond'   : '--masking 0',
                                             'ncbi_blast': ''}


# ANCHOR - PangenomeGraphManager
class PangenomeGraphManager():
    """All in one pangenome graph object.

    The purpose of this class is to create a nx.DiGraph object with a fixed set of
    attributes per node and edge that is still
    useable with nx.DiGraph algorithms. This pangenome graph object inherits all
    functions from a nx.DiGraph object with an extra set of boundaries. Therefore it
    can be used with all existing networkx algorithms made for digraph objects. The
    visible changes to the inheriteted functions ensure that the pangenome graph object
    always has the same set of attributes to ensure smooth usage in the anvi'o pangenome
    graph suitcase. Most attributes are self explanatory. Group attributes connects
    single connection nodes e.g. operon or conserved regions. Active is necessary for
    UI changes, it describes wether the edge will be shown or was removed. Bended describes
    long edges over multiple positions. This whole class is created after the nx.DiGraph
    class in the networkx package.

    Definition of PangenomeGraphManager object nodes:
    id: syn_cluster ----> 'GC1_1'
    attributes:
        gene_cluster: --> 'GC1'
        position: ------> (x,y)
        gene_calls: ----> {G1: 0, G2: 5}
        group: ---------> 'GCG3'
        component_id: --> 0

    Definition of PangenomeGraphManager object edges:
    i: syn_cluster ----> 'GC1_1'
    j: syn_cluster ----> 'GC2_1'
    attributes:
        weight: --------> 2
        active: --------> True
        genomes: -------> ['G1', 'G2']
        route: --------> [(0,1), (0,2)]
    """
    def __init__(self, run=run, progress=progress):

        self.run = run
        self.progress = progress

        self.node_standard_attributes = {
            'gene_cluster': '',
            'position': (0,0),
            'gene_calls': {},
            'synteny': {},
            'type': '',
            'group': '',
            'layer': {},
            'alignment': '',
            'component_id': '',
        }
        self.edge_standard_attributes = {
            'name': '',
            'weight': 0.0,
            'active': True,
            'genomes': [],
            'route': [],
            'length': 0,
        }
        self.graph = nx.DiGraph()


    def sanity_check_node_attributes(self, node_attributes):
        if not set(node_attributes.keys()).issubset(self.node_standard_attributes.keys()):
            raise ConfigError("Node attributes not properly formated")

        for node_attribute in self.node_standard_attributes:
            if node_attribute in node_attributes:
                node_attribute_type = type(self.node_standard_attributes[node_attribute])
                if not isinstance(node_attributes[node_attribute], node_attribute_type):
                    raise ConfigError(f"Node attribute {node_attribute} is not a {node_attribute_type} type.")
            else:
                node_attributes[node_attribute] = self.node_standard_attributes[node_attribute]

        return(node_attributes)


    def sanity_check_edge_attributes(self, edge_attributes):
        if not set(edge_attributes.keys()).issubset(self.edge_standard_attributes.keys()):
            raise ConfigError("Edge attributes not properly formated")

        for edge_attribute in self.edge_standard_attributes:
            if edge_attribute in edge_attributes:
                edge_attribute_type = type(self.edge_standard_attributes[edge_attribute])
                if not isinstance(edge_attributes[edge_attribute], edge_attribute_type):
                    raise ConfigError(f"Edge attribute {edge_attribute} is not a {edge_attribute_type} type.")
            else:
                edge_attributes[edge_attribute] = self.edge_standard_attributes[edge_attribute]

        return(edge_attributes)


    def add_node_to_graph(self, syn_cluster, attributes):

        attributes = self.sanity_check_node_attributes(attributes)

        if not self.graph.has_node(syn_cluster):
            self.graph.add_node(syn_cluster, **attributes)
        else:
            self.graph.nodes[syn_cluster]['gene_calls'].update(attributes['gene_calls'])
            self.graph.nodes[syn_cluster]['synteny'].update(attributes['synteny'])


    def add_edge_to_graph(self, syn_cluster_i, syn_cluster_j, attributes):

        attributes = self.sanity_check_edge_attributes(attributes)

        if not self.graph.has_edge(*(syn_cluster_i, syn_cluster_j)):
            self.graph.add_edge(*(syn_cluster_i, syn_cluster_j), **attributes)
        else:
            merged_genomes = sorted(set(self.graph[syn_cluster_i][syn_cluster_j]['genomes'])
                                    | set(attributes['genomes']))
            self.graph[syn_cluster_i][syn_cluster_j]['genomes'] = merged_genomes
            self.graph[syn_cluster_i][syn_cluster_j]['weight'] = float(len(merged_genomes))


    def run_connectivity_check(self):
        """One of the main functions. It checks the graph for fragmented nature.
        This is very common if the underlying data is based on draft genomes. Some
        gene clusters might not be connected to the rest of the graph because e.g.
        the splitting of the gene clusters into syn clusters left them allone to
        prevent problems. Aside from only checking the graph, removes the small
        fragmented pieces that are not connected to the main graph.

        Parameters
        ==========
        self.graph

        Returns
        =======
        self.graph
        """

        self.run.warning("A fully connected graph is required for downstream processing. If the graph contains "
                         "disconnected components (common with draft genomes or fragmented assemblies), only the "
                         "largest connected component will be retained. Disconnected fragments typically arise from "
                         "incomplete assemblies or highly variable accessory genes.",
                         header="CHECKING GRAPH CONNECTIVITY", lc="green")

        connectivity = nx.is_connected(self.graph.to_undirected())

        if connectivity == False:
            weakly_components = list(nx.weakly_connected_components(self.graph))
            num_components = len(weakly_components)
            component_sizes = sorted([len(comp) for comp in weakly_components], reverse=True)

            self.run.info('Graph connectivity', 'Disconnected', mc='red')
            self.run.info('Number of components', num_components)
            self.run.info('Largest component size', f"{component_sizes[0]} nodes")
            if num_components > 1:
                self.run.info('Other component sizes', f"{component_sizes[1:10]}" + ("..." if num_components > 10 else ""))

            self.run.warning(f"The graph contains {num_components} disconnected component(s), likely due to "
                             f"fragmented genomes or highly variable accessory regions. Keeping only the largest component.")

            subgraph = self.graph.subgraph(max(weakly_components, key=len))
            nodes_removed = len(self.graph.nodes()) - len(subgraph.nodes())
            edges_removed = len(self.graph.edges()) - len(subgraph.edges())

            self.graph = nx.DiGraph(subgraph)

            self.run.info('Nodes removed', nodes_removed, mc='red')
            self.run.info('Edges removed', edges_removed, mc='red')
            self.run.info('Remaining nodes', len(self.graph.nodes()), mc='green')
            self.run.info('Remaining edges', len(self.graph.edges()), mc='green')

            connectivity = nx.is_connected(self.graph.to_undirected())
            if connectivity == True:
                self.run.info('Final connectivity status', 'Connected', mc='green')
            else:
                raise ConfigError("The graph is still fragmented after removing small components. This suggests "
                                  "severe inconsistencies in your data. Please verify genome quality and pangenome "
                                  "generation parameters.")
        else:
            self.run.info('Graph connectivity', 'Fully connected', mc='green')


    def walk_one_step(self, G, current, nodes_position_dict, visited):
        successors = [successor for successor in G.successors(current) if successor not in visited]
        predecessors = [predecessor for predecessor in G.predecessors(current) if predecessor not in visited]
        current_position_x, current_position_y = nodes_position_dict[current]

        if len(successors) > 0:
            nodes_position = [(nodes_position_dict[successor][0] - current_position_x, successor) for successor in successors]
            successor_position, successor = sorted(nodes_position)[0]
            successor_position_x, successor_position_y = nodes_position_dict[successor]

            return(successor_position_x, successor_position_y, successor)

        elif len(predecessors) > 0:
            nodes_position = [(current_position_x - nodes_position_dict[predecessor][0], predecessor) for predecessor in predecessors]
            predecessor_position, predecessor = sorted(nodes_position, reverse=True)[0]
            predecessor_position_x, predecessor_position_y = nodes_position_dict[predecessor]

            return(predecessor_position_x, predecessor_position_y, predecessor)

        else:
            return(-1, -1, '')


    def is_node_core(self, node, genome_names):
        if set(self.graph.nodes()[node]['gene_calls'].keys()) == set(genome_names):
            return(True)
        else:
            return(False)


    def merge_dict_list(self, dict1, dict2):

        merged_dict = {}

        for key in set(dict1) | set(dict2):
            merged_dict[key] = []

            if key in dict1:
                merged_dict[key].extend(dict1[key])

            if key in dict2:
                merged_dict[key].extend(dict2[key])
        return(merged_dict)


    def layout_all_components(self, gene_cluster_grouping_threshold, groupcompress):
        """Run TopologicalLayout on every weakly-connected component and
        apply the merged positions/edges/groups across the full graph.

        Each component is laid out in its own local coordinate frame so the
        viewer can render any one component cleanly. Cross-component
        positions overlap on purpose -- only one component is shown at a
        time (see dbops.get_json's component filter).

        Returns the number of components laid out.
        """
        from anvio.topologicallayout import TopologicalLayout

        n_components = sum(1 for _ in nx.weakly_connected_components(self.graph))

        merged_positions = {}
        merged_edges = {}
        merged_groups = {}
        for cid in range(n_components):
            pos, epos, grps = TopologicalLayout(r=self.run, p=self.progress).run_synteny_layout_algorithm(
                F=self.graph,
                gene_cluster_grouping_threshold=gene_cluster_grouping_threshold,
                groupcompress=groupcompress,
                component=cid,
            )
            merged_positions.update(pos)
            merged_edges.update(epos)
            merged_groups.update(grps)

        self.set_node_positions(merged_positions)
        self.set_edge_positions(merged_edges)
        self.set_node_groups(merged_groups)
        return n_components


    def summarize_all_components(self, scope='global', use_stored_backbone=True):
        """Run summarize() for every component and return the concatenated
        ``region_sides_df`` plus a ``{node -> backbone}`` mapping (1 if the
        node sits in a BR region, 0 if VR).

        Components are the ``component_id`` values present on the nodes,
        not ``CP_0001`` to ``CP_<number of components>``: a genome subset in
        the interface can dissolve components (and split others off under new
        ids), so the ids need not be contiguous.

        ``region_id`` values embed their component (``"CP_0001_1"``, see
        ``_assign_region_ids``), and each row also carries a ``component_id``
        column so consumers (the pan_graph_regions DB table, ``self.regions``,
        ``self.region_sides_info``) can key by ``(component_id, region_id)``.

        ``scope`` is forwarded to ``summarize`` and controls whether the
        genome-count denominator (used for BR/VR classification and the
        normalized metrics) spans the whole graph (``'global'``) or only
        the current weakly-connected component (``'component'``).
        ``use_stored_backbone`` is forwarded as well (see ``summarize``).
        """
        component_ids = sorted({d.get('component_id', 'CP_0001') for _, d in self.graph.nodes(data=True)},
                               key=lambda cid: int(cid.split('_')[1]))

        region_dfs = []
        backbone_by_node = {}
        for cid in component_ids:
            region_sides_df, nodes_df, _ = self.summarize(component_id=cid, scope=scope,
                                                          use_stored_backbone=use_stored_backbone)
            if region_sides_df.empty:
                continue
            region_dfs.append(region_sides_df)
            merged = pd.merge(
                region_sides_df.reset_index(),
                nodes_df.reset_index(),
                how="left", on="region_id",
            ).set_index('syn_cluster')
            for syn, row in merged.iterrows():
                backbone_by_node[syn] = 1 if row['region'] == 'BR' else 0

        all_region_sides_df = (pd.concat(region_dfs) if region_dfs
                               else pd.DataFrame())
        return all_region_sides_df, backbone_by_node


    def summarize(self, component_id='CP_0001', scope='global', use_stored_backbone=True):
        """Compute region-level summaries for one weakly connected component.

        The component is selected by the ``component_id`` attribute set on
        every node (see Task 3 in HANDOFF.md): ``component_id="CP_0001"`` is the
        largest. ``region_id`` embeds its component and a 1-based index
        (``"CP_0001_1"``, ``"CP_0001_2"``, ...), so it is globally unique on its
        own; the ``component_id`` column is still attached to ``region_sides_df``
        and stored alongside on each region row in the DB.

        ``use_stored_backbone`` decides whether a ``backbone`` flag already in a
        node's layer data overrides the coverage rule (see ``_node_is_core``).
        Pass False when the graph is a genome subset of the one the flags were
        computed for, so core-ness follows the genomes actually in the graph.

        Returns ``(region_sides_df, nodes_df, gene_calls_df)``. All three may
        be empty DataFrames if the requested component is empty or contains
        no nodes whose x-position falls inside any region span.
        """
        nodes_in_component = [n for n, d in self.graph.nodes(data=True)
                              if d.get('component_id', 'CP_0001') == component_id]
        if not nodes_in_component:
            self.run.info_single(f"Component {component_id} is empty — nothing to summarize.")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

        G = self.graph.subgraph(nodes_in_component)

        if not nx.is_directed_acyclic_graph(G):
            raise ConfigError("Cyclic graphs are not implemented.")

        # The genome-count denominator drives both BR/VR classification
        # (`len(genomes_involved) == num_genomes`) and the normalized
        # metrics. ``scope='global'`` counts genomes across the whole
        # pangenome (default; previously hardcoded). ``scope='component'``
        # counts only the genomes present in this WCC -- a region is BR
        # only relative to the genomes that actually reach this
        # component, which makes sense when each component is treated as
        # an independent sub-pangenome rather than a slice of the global
        # one. The original hardcoded global behavior was introduced to
        # avoid trivially classifying a 1-genome component as backbone
        # by virtue of its denominator collapsing to 1.
        if scope == 'component':
            genome_names = set(it.chain(*[list(d.keys())
                                          for _, d in G.nodes(data='gene_calls')]))
        elif scope == 'global':
            genome_names = set(it.chain(*[list(d.keys())
                                          for _, d in self.graph.nodes(data='gene_calls')]))
        else:
            raise ConfigError(f"Unknown region scope {scope!r}; expected 'global' or 'component' :/")
        num_genomes = len(genome_names)

        all_positions = sorted({d['position'][0] for _, d in G.nodes(data=True)})

        def _node_is_core(d):
            """Whether a node's x-column counts as 'core' -- a backbone column
            that partitions the region layout.

            By default this is coverage-driven (the node spans every genome in
            scope). But if the node already carries a per-node ``backbone`` flag
            (e.g. one computed by an earlier summary / an older engine and stored
            in its layer data), honour that flag instead: a node deliberately
            kept OFF the backbone must not re-split a region here just because
            its coverage happens to reach all genomes. The flag is only present
            once a prior summary has run, so a fresh generation still falls back
            to the coverage rule. With ``use_stored_backbone`` False the flag is
            ignored: it was computed for a different genome set."""
            backbone = (d.get('layer') or {}).get('backbone') if use_stored_backbone else None
            try:
                if backbone is not None:
                    return int(float(backbone)) == 1
            except (TypeError, ValueError):
                pass
            return len(d['gene_calls'].keys()) == num_genomes

        core_positions = sorted({d['position'][0]
                                 for _, d in G.nodes(data=True)
                                 if _node_is_core(d)})

        regions_dict, backbone_region_ids = self._assign_region_ids(
            core_positions, all_positions, component_id)

        nodes_df, regions_info_dict = self._attach_node_regions(G, regions_dict)

        regions_summary_rows = []
        for region_id, values_list in sorted(regions_info_dict.items()):
            regions_summary_rows.append(
                self._summarize_one_region(region_id, values_list, num_genomes,
                                           region_id in backbone_region_ids))

        region_sides_df = pd.DataFrame(regions_summary_rows).set_index('region_id') \
                          if regions_summary_rows else pd.DataFrame()

        if not region_sides_df.empty:
            # Tag every row with the component it belongs to so callers that
            # concatenate across components can disambiguate the (now plain)
            # per-component region_ids.
            region_sides_df['component_id'] = component_id

        gene_calls_df = self._build_gene_calls_df(G)

        return region_sides_df, nodes_df, gene_calls_df


    def _assign_region_ids(self, core_positions, all_positions, component_id):
        """Map every x position to a region id that embeds its component and a
        1-based index (``"CP_0001_1"``, ``"CP_0001_2"``, ...) and report which
        of those regions are backbone runs.

        Walks the integer x columns from ``all_min`` to ``all_max`` and
        partitions them into maximal runs of all-core vs all-non-core
        columns. Each run becomes one region; the index is sequential from 1
        within the component and the ``component_id`` is prefixed so the id is
        globally unique on its own.

        Returns ``(regions_dict, backbone_region_ids)`` where
        ``backbone_region_ids`` is the set of region_id strings whose columns
        are core -- i.e. the backbone (BR) regions. The BR/VR label of every
        region follows directly from this partition, so whatever fed the
        boundaries (the stored backbone flag when present, else coverage) also
        determines backbone/variable downstream.
        """
        if not all_positions:
            return {}, set()

        all_min, all_max = min(all_positions), max(all_positions)

        if not core_positions:
            # No core columns -> the whole component is a single variable region.
            return {pos: f"{component_id}_1" for pos in range(all_min, all_max + 1)}, set()

        core_set = set(core_positions)
        regions_dict = {}
        backbone_region_ids = set()
        rid = 1
        current_is_core = (all_min in core_set)
        if current_is_core:
            backbone_region_ids.add(f"{component_id}_{rid}")
        for x in range(all_min, all_max + 1):
            is_core = x in core_set
            if is_core != current_is_core:
                rid += 1
                current_is_core = is_core
                if is_core:
                    backbone_region_ids.add(f"{component_id}_{rid}")
            regions_dict[x] = f"{component_id}_{rid}"
        return regions_dict, backbone_region_ids


    def _attach_node_regions(self, G, regions_dict):
        """Set ``region_id`` on every node whose x falls in a known region
        and return ``(nodes_df, regions_info_dict)``.

        ``regions_info_dict`` is keyed by region_id and holds the
        ``(x, y, node, genomes)`` tuples used by ``_summarize_one_region``.
        """
        regions_info_dict = {}
        node_rows = []
        for node, data in G.nodes(data=True):
            x = data['position'][0]
            y = data['position'][1]
            if x not in regions_dict:
                continue
            region_id = regions_dict[x]
            genomes = list(data['gene_calls'].keys())

            node_rows.append({'syn_cluster': node, 'x': x, 'y': y, 'region_id': region_id})
            G.nodes[node]['region_id'] = region_id
            regions_info_dict.setdefault(region_id, []).append((x, y, node, genomes))

        nodes_df = (pd.DataFrame(node_rows).set_index('syn_cluster')
                    if node_rows else pd.DataFrame())
        return nodes_df, regions_info_dict


    def _summarize_one_region(self, region_id, values_list, num_genomes, is_backbone, log_size_cap = 100, haplotype_type = 'dissimilarity', count_empty_allel = True):
        """Compute one row of the region summary table (BR / VR + the
        complexity / diversity / expansion / weight metrics and the composite
        variability score). Every metric comes from this region alone, so a row
        does not depend on the other regions of the dataset.

        ``is_backbone`` comes straight from the region-boundary partition
        (``_assign_region_ids``): a region built from core columns is BR, any
        other region is VR. The variability metrics below are still computed for
        display, but no longer decide the BR/VR label -- that follows the
        boundaries, so whatever fed them (the stored backbone flag when present,
        else coverage) is what classifies the region.

        ``log_size_cap``       genes at which ``expansion_log`` reaches 1.
        ``haplotype_type``     the diversity metric, 'dissimilarity' or 'diversity'.
        ``count_empty_allel``  whether the genomes without the region count as an
                               allele (the empty one) in the diversity.

        The composite variability score is reported twice, as the geometric and
        as the arithmetic mean of the four 0-to-1 metrics."""
        if haplotype_type not in ('dissimilarity', 'diversity'):
            raise ConfigError(f"Unknown haplotype type {haplotype_type!r}; expected 'dissimilarity' or 'diversity' :/")

        region_x, _, nodes, genomes_sets = zip(*values_list)

        # Every genome has one allele in the region: the set of the region's SynGCs
        # it carries (empty if it carries none). A genome is in a SynGC at most
        # once, so the size of its allele is its number of genes in the region.
        alleles_of_genomes = {}
        for node, genomes in zip(nodes, genomes_sets):
            for genome in genomes:
                alleles_of_genomes.setdefault(genome, set()).add(node)
        genes_per_genome = [len(allele) for allele in alleles_of_genomes.values()]
        num_genomes_present = len(alleles_of_genomes)   # the genomes that carry the region
        weight = num_genomes_present / num_genomes

        # Expansion: the genes of the genome with the most genes in the region,
        # also on a log scale, log(1 + genes) / log(1 + log_size_cap), capped at 1.
        max_expansion = max(genes_per_genome)
        min_expansion = min(genes_per_genome) if num_genomes_present == num_genomes else 0
        expansion_log = min(np.log1p(max_expansion) / np.log1p(log_size_cap), 1.0)

        # Complexity (proxy for structural variation): the SynGCs needed to cover
        # all genomes present, taking their genome sets smallest-first and counting
        # each that adds a genome not yet covered (num_covering_syngcs), scaled to
        # 0 (one SynGC covers them) to 1 (one SynGC per genome of num_genomes).
        covered, num_covering_syngcs = set(), 0
        for genomes in sorted(map(set, genomes_sets), key=len):
            if not genomes <= covered:
                num_covering_syngcs += 1
                covered |= genomes
        complexity = (num_covering_syngcs - 1) / (num_genomes - 1) if num_genomes > 1 else 0

        # Diversity: the haplotype metric ``haplotype_type`` over the alleles of N
        # genomes, 1 - the mean similarity of their pairs. The similarity of two
        # alleles is
        #   'dissimilarity'  their Jaccard index: graded haplotype diversity, the
        #                    expected Jaccard distance between the alleles of two
        #                    random genomes. Equal to Nei's where the alleles share
        #                    nothing, lower where they overlap.
        #   'diversity'      1 if they are the same, else 0: Nei's haplotype
        #                    diversity, N/(N-1) * (1 - sum p^2), the probability
        #                    that two random genomes differ in the region at all.
        # Two genomes without the region have the same (empty) allele. N is
        # ``num_genomes``; with ``count_empty_allel`` False only the genomes
        # present are compared (N is num_genomes_present), so a single-variant VR
        # has 0. Both are 0 to 1; a backbone region has one allele, so 0.
        N = num_genomes if count_empty_allel else num_genomes_present
        diversity = 0.0
        if N > 1:
            genomes_per_allele = Counter(map(frozenset, alleles_of_genomes.values()))
            allele_counts = np.fromiter(genomes_per_allele.values(), dtype=float)
            # the pairs with the same allele, the empty one included
            similarity = ((allele_counts * (allele_counts - 1) / 2).sum()
                          + (N - num_genomes_present) * (N - num_genomes_present - 1) / 2)
            if haplotype_type == 'dissimilarity':
                # plus the Jaccard index of every two different alleles, from the SynGCs
                # they share (one sparse matrix product)
                alleles = list(genomes_per_allele)
                size = np.array([len(allele) for allele in alleles])
                node_index = {node: i for i, node in enumerate(nodes)}
                rows = np.repeat(np.arange(len(alleles)), size)
                cols = [node_index[node] for allele in alleles for node in allele]
                A = csr_matrix((np.ones(len(cols)), (rows, cols)), shape=(len(alleles), len(nodes)))
                shared = (A @ A.T).tocoo()
                pair = shared.row < shared.col
                a, b, n = shared.row[pair], shared.col[pair], shared.data[pair]
                similarity += (allele_counts[a] * allele_counts[b] * n / (size[a] + size[b] - n)).sum()
            diversity = float(1 - similarity / (N * (N - 1) / 2))

        # Composite variability score: the geometric and the arithmetic mean of the
        # four 0-to-1 metrics.
        metrics = np.array([complexity, diversity, expansion_log, weight])
        composite_variability_score_geometric = float(metrics.prod() ** (1 / len(metrics)))
        composite_variability_score_arithmetic = float(metrics.mean())

        return {
            'region_id': region_id,
            'region': 'BR' if is_backbone else 'VR',   # follows the region boundaries, not the metrics
            'x_min': min(region_x),
            'x_max': max(region_x),
            'num_synteny_gene_clusters': len(nodes),
            'num_gene_clusters': len({node.rsplit('_', 1)[0] for node in nodes}),
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
            'composite_variability_score_geometric': composite_variability_score_geometric,
            'composite_variability_score_arithmetic': composite_variability_score_arithmetic,
        }


    def _build_gene_calls_df(self, G):
        """Flatten ``{node -> {genome -> gene_caller_id}}`` into a long
        DataFrame indexed by (genome, gene_caller_id)."""
        rows = []
        for node, gene_calls in G.nodes(data='gene_calls'):
            if not gene_calls:
                continue
            for genome, gid in gene_calls.items():
                rows.append({'syn_cluster': node, 'genome': genome, 'gene_caller_id': gid})
        if not rows:
            return pd.DataFrame()
        return pd.DataFrame(rows).set_index(['genome', 'gene_caller_id'])


    def set_node_positions(self, node_positions):
        # Normal callers go through layout_all_components and produce a
        # position for every node. The membership check is defensive --
        # if a caller ever lays out only a subset, nodes outside that
        # subset keep their default (0, 0) from node_standard_attributes.
        for node in self.graph.nodes():
            if node in node_positions:
                self.graph.nodes()[node]['position'] = node_positions[node]


    def set_node_groups(self, node_groups):
        for node in self.graph.nodes():
            if node in node_groups:
                self.graph.nodes()[node]['group'] = node_groups[node]
            else:
                self.graph.nodes()[node]['group'] = ''


    def set_edge_positions(self, edge_positions):
        # Defensive: normal callers cover every edge via
        # layout_all_components. Edges that ever slip through stay at
        # the default route=[], length=0 from edge_standard_attributes.
        for edge_i, edge_j in self.graph.edges():
            if (edge_i, edge_j) not in edge_positions:
                continue
            route = edge_positions[(edge_i, edge_j)]
            self.graph[edge_i][edge_j]['route'] = route
            if route:
                length = route[-1][0] - route[0][0] + 1
            else:
                length = 0
            self.graph[edge_i][edge_j]['length'] = length


    def cut_edges(self, max_edge_length_filter):
        cutted_edges = []
        for edge_i, edge_j in self.graph.edges():
            length = self.graph[edge_i][edge_j]['length']
            if max_edge_length_filter != -1 and max_edge_length_filter < length:
                self.graph[edge_i][edge_j]['active'] = False
                cutted_edges += [(edge_i, edge_j)]
            else:
                self.graph[edge_i][edge_j]['active'] = True

        return(cutted_edges)


    def calculate_graph_distance(self):
        """Compute graph-based pairwise genome distances.

           Returns a (newick, distance_df, genome_names) tuple. `distance_df` is a square
           pandas DataFrame indexed and columned by genome names. `newick` may be an empty
           string if all pairwise distances are zero (i.e. genomes are graph-identical).
        """
        self.run.warning("Anvi'o now derives a pairwise distance between every pair of genomes directly from "
                         "the graph structure: genomes that traverse the same synteny gene clusters in the "
                         "same order are close, while genomes that diverge (through rearrangements, "
                         "insertions, or accessory content) are further apart. The full distance matrix is "
                         "written into the pan-graph-db and used to order genomes in downstream displays.",
                         header="COMPUTING GRAPH-BASED GENOME DISTANCES", lc="green")
        # sorted (not just list(set(...))) so the distance-matrix row/column order and the
        # newick leaf order below are deterministic across runs -- otherwise set iteration
        # order (process-dependent) would reshuffle the genome ordering every run.
        genome_names = sorted(set(it.chain(*[list(d.keys()) for node, d in self.graph.nodes(data='gene_calls')])))

        # The distance between two genomes is the Jaccard distance over the set of graph
        # elements, and the entire matrix is simply a single dot product, as in, there is
        # no reason to loop over pairs of genomes here (see commit logs for these lines
        # for details).
        genome_index = {genome_name: i for i, genome_name in enumerate(genome_names)}

        # one row per graph element (nodes first, then edges), one column per genome
        memberships = [d or {} for _, d in self.graph.nodes(data='gene_calls')] \
                    + [d.get('genomes') or [] for _, _, d in self.graph.edges(data=True)]

        M = np.zeros((len(memberships), len(genome_names)))
        for row, members in enumerate(memberships):
            columns = [genome_index[genome_name] for genome_name in members if genome_name in genome_index]
            if columns:
                M[row, columns] = 1.0

        num_shared = M.T @ M                                    # elements both genomes are in
        num_elements = M.sum(axis=0)                            # elements each genome is in
        num_total = num_elements[:, None] + num_elements[None, :]
        num_union = num_total - num_shared                      # similar + unsimilar

        # building a sparse distance matrix here
        X = np.divide(num_total - 2 * num_shared, num_union, out=np.zeros_like(num_union), where=num_union > 0)
        np.fill_diagonal(X, 0.0)

        if len(genome_names) > 1:
            upper_i, upper_j = np.triu_indices(len(genome_names), k=1)
            distances = X[upper_i, upper_j]
            argmin, argmax = distances.argmin(), distances.argmax()
            min_dist, max_dist = distances[argmin], distances[argmax]
            min_pair = (genome_names[upper_i[argmin]], genome_names[upper_j[argmin]])
            max_pair = (genome_names[upper_i[argmax]], genome_names[upper_j[argmax]])

            self.run.info('Smallest distance', f"d({min_pair[0]},{min_pair[1]}) = {round(min_dist, 3)}")
            self.run.info('Largest distance',  f"d({max_pair[0]},{max_pair[1]}) = {round(max_dist, 3)}")
            self.run.info_single("full matrix written to the pan-graph-db", cut_after=None)

        distance_matrix = pd.DataFrame(X, index=genome_names, columns=genome_names)

        # Check if all distances are zero (identical genomes)
        if np.all(X == 0):
            self.run.warning("All pairwise distances are zero (genomes are identical). No dendrogram will be generated.")
            return('', distance_matrix, genome_names)

        condensed_X = squareform(X)
        Z = linkage(condensed_X, 'ward')

        tree = to_tree(Z, False)
        newick = clustering.get_newick(tree, tree.dist, genome_names)

        return(newick, distance_matrix, genome_names)

    # TODO still empty..
    def generate_hybrid_genome(self, output_dir):
        self.run.info_single("The output hybrid genome function exists on paper but we haven't implemented it yet. sorry.")
