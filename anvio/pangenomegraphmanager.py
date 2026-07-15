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
import statistics as stat
from collections import Counter

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


def composite_variability_score(df):
    """Min-max-scale the four variability components of a region summary
    DataFrame and combine them into a composite score per region.

    ``df`` must contain ``complexity_normalized``, ``diversity_normalized``,
    ``max_expansion``, and ``weight_normalized``. Each is min-max scaled
    across the rows of ``df`` (constant columns collapse to 0). The
    composite score is the geometric mean of the four scaled components.

    Returns a DataFrame indexed like ``df`` with five new columns:
    ``composite_variability_score``, ``complexity_mm_scaled``,
    ``diversity_mm_scaled``, ``expansion_mm_scaled``, ``weight_mm_scaled``.
    """
    def _minmax(series):
        s_min, s_max = series.min(), series.max()
        if s_min == s_max:
            return pd.Series(0.0, index=series.index)
        return (series - s_min) / (s_max - s_min)

    c_mm = _minmax(df['complexity_normalized'])
    d_mm = _minmax(df['diversity_normalized'])
    e_mm = _minmax(df['max_expansion'])
    w_mm = _minmax(df['weight_normalized'])
    cvs = (c_mm * d_mm * e_mm * w_mm) ** (1 / 4)

    return pd.DataFrame({
        'composite_variability_score': cvs,
        'complexity_mm_scaled': c_mm,
        'diversity_mm_scaled': d_mm,
        'expansion_mm_scaled': e_mm,
        'weight_mm_scaled': w_mm,
    }, index=df.index)


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
            pos, epos, grps = TopologicalLayout().run_synteny_layout_algorithm(
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


    def summarize_all_components(self, scope='global'):
        """Run summarize() for every weakly-connected component and return
        the concatenated ``region_sides_df`` plus a ``{node -> backbone}``
        mapping (1 if the node sits in a BR region, 0 if VR).

        ``region_id`` values are plain ``"0"``, ``"1"``, ... per component
        (see ``_assign_region_ids``); the concatenated frame may carry
        duplicate ``region_id`` index values across components. Each row
        carries a ``component_id`` column so consumers (the pan_graph_regions
        DB table, ``self.regions``, ``self.region_sides_info``) can key by
        ``(component_id, region_id)``.

        ``scope`` is forwarded to ``summarize`` and controls whether the
        genome-count denominator (used for BR/VR classification and the
        normalized metrics) spans the whole graph (``'global'``) or only
        the current weakly-connected component (``'component'``).
        """
        n_components = sum(1 for _ in nx.weakly_connected_components(self.graph))

        region_dfs = []
        backbone_by_node = {}
        for i in range(n_components):
            cid = f"CP_{i + 1:04d}"
            region_sides_df, nodes_df, _ = self.summarize(component_id=cid, scope=scope)
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


    def summarize(self, component_id='CP_0001', scope='global'):
        """Compute region-level summaries for one weakly connected component.

        The component is selected by the ``component_id`` attribute set on
        every node (see Task 3 in HANDOFF.md): ``component_id="CP_0001"`` is the
        largest. ``region_id`` embeds its component and a 1-based index
        (``"CP_0001_1"``, ``"CP_0001_2"``, ...), so it is globally unique on its
        own; the ``component_id`` column is still attached to ``region_sides_df``
        and stored alongside on each region row in the DB.

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
            to the coverage rule."""
            backbone = (d.get('layer') or {}).get('backbone')
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
            cvs_cols = composite_variability_score(region_sides_df)
            region_sides_df = pd.concat([region_sides_df, cvs_cols], axis=1)
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


    def _summarize_one_region(self, region_id, values_list, num_genomes, is_backbone):
        """Compute one row of the region summary table (BR / VR + the
        complexity / diversity / expansion / weight metrics).

        ``is_backbone`` comes straight from the region-boundary partition
        (``_assign_region_ids``): a region built from core columns is BR, any
        other region is VR. The variability metrics below are still computed for
        display, but no longer decide the BR/VR label -- that follows the
        boundaries, so whatever fed them (the stored backbone flag when present,
        else coverage) is what classifies the region."""
        genomes_sets = [item[3] for item in values_list]
        genomes_involved = set(it.chain(*genomes_sets))
        nodes_sets = [item[2] for item in values_list]
        region_x = [item[0] for item in values_list]

        x_min = min(region_x)
        x_max = max(region_x)

        weight = len(genomes_involved)
        num_gene_clusters = len({n.rsplit('_', 1)[0] for n in nodes_sets})
        num_synteny_gene_clusters = len(nodes_sets)
        genome_occurences = list(it.chain(*genomes_sets))
        num_gene_calls = len(genome_occurences)

        counts = Counter(genome_occurences)
        values = [counts.get(g, 0) for g in genomes_involved]
        max_expansion = max(values)
        min_expansion = min(values) if len(genomes_involved) == num_genomes else 0

        # Complexity = number of genome-sets needed to cover all genomes_involved
        # when consumed smallest-first (proxy for structural variation).
        sum_of_genomes = set()
        complexity = 0
        for genome_set in sorted(genomes_sets, key=len):
            if set(genome_set).issubset(sum_of_genomes):
                continue
            complexity += 1
            sum_of_genomes |= set(genome_set)

        complexity_normalized = (complexity - 1) / (num_genomes - 1) if num_genomes > 1 else 0

        if len(genomes_sets) > 1 and not is_backbone:
            diversity = (stat.pvariance([1, num_genomes])
                         - stat.pvariance([len(s) for s in genomes_sets]))
            diversity_normalized = (stat.pvariance([1 / num_genomes, 1])
                                - stat.pvariance([len(s) / num_genomes for s in genomes_sets]))
        else:
            diversity = 0
            diversity_normalized = 0

        weight_normalized = weight / num_genomes if num_genomes > 0 else 0

        # BR/VR follows the region boundaries, not the variability metrics.
        region = 'BR' if is_backbone else 'VR'

        return {
            'region_id': region_id,
            'region': region,
            'x_min': x_min,
            'x_max': x_max,
            'num_synteny_gene_clusters': num_synteny_gene_clusters,
            'num_gene_clusters': num_gene_clusters,
            'num_gene_calls': num_gene_calls,
            'max_expansion': max_expansion,
            'min_expansion': min_expansion,
            'complexity': complexity,
            'complexity_normalized': complexity_normalized,
            'diversity': diversity,
            'diversity_normalized': diversity_normalized,
            'weight': weight,
            'weight_normalized': weight_normalized,
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

    # @staticmethod
    # def composite_variability_score(df, N):

    #     C_max = df['complexity'].max()
    #     # Option 2:
    #     # C_min = df['complexity'].min()

    #     E_max = df['expansion'].max()
    #     # Option 2:
    #     # E_min = df['max_expansion'].min()

    #     W_median = df['weight'].median()

    #     def log_min_max_normalize(X, X_min, X_max):
    #         if X_min == X_max:
    #             return 0.0
    #         X_norm = (math.log(1 + X) - math.log(1 + X_min)) / (math.log(1 + X_max) - math.log(1 + X_min))
    #         return(X_norm)

    #     def fixed_baseline_log_scaling(X, X_max):
    #         if X_max == 0.0:
    #             return(0.0)
    #         else:
    #             X_norm = (math.log(1 + X)) / (math.log(1 + X_max))
    #             return(X_norm)

    #     def func(row):

    #         C = row['complexity']
    #         E = row['expansion']
    #         W = row['weight']
    #         D = row['diversity']

    #         C_norm = fixed_baseline_log_scaling(C, C_max)
    #         # Option 2:
    #         # C_norm = log_min_max_normalize(C, C_min, C_max)

    #         E_norm = fixed_baseline_log_scaling(E, E_max)
    #         # Option 2:
    #         # E_norm = log_min_max_normalize(E, E_min, E_max)

    #         W_f = W/N
    #         # Option 2:
    #         # W_f = (1 - math.e**(-W/W_median))
    #         # Option 3:
    #         # W_f = (1 - math.e**(-W/W_median)) / (1 - math.e**(-N/W_median))

    #         CVS = (C_norm * D * E_norm)**(1/3) * W_f

    #         return([CVS, C_norm, E_norm, W_f])

    #     return(func)


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

        X = np.zeros([len(genome_names), len(genome_names)])
        min_pair = max_pair = None
        min_dist = float('inf')
        max_dist = float('-inf')
        for genome_i, genome_j in it.combinations(genome_names, 2):
            nodes_similar = 0
            edges_similar = 0
            nodes_unsimilar = 0
            edges_unsimilar = 0
            for _, data in self.graph.nodes(data=True):
                if genome_i in data['gene_calls'].keys() and genome_j in data['gene_calls'].keys():
                    nodes_similar += 1
                elif genome_i in data['gene_calls'].keys() or genome_j in data['gene_calls'].keys():
                    nodes_unsimilar += 1
            for _, _, data in self.graph.edges(data=True):
                edge_genomes = set(data.get('genomes', []))
                if genome_i in edge_genomes and genome_j in edge_genomes:
                    edges_similar += 1
                elif genome_i in edge_genomes or genome_j in edge_genomes:
                    edges_unsimilar += 1

            i = genome_names.index(genome_i)
            j = genome_names.index(genome_j)

            elements_similar = nodes_similar + edges_similar
            elements_unsimilar = nodes_unsimilar + edges_unsimilar

            d = elements_unsimilar / (elements_similar + elements_unsimilar)
            X[i][j] = d
            X[j][i] = d

            if d < min_dist:
                min_dist, min_pair = d, (genome_i, genome_j)
            if d > max_dist:
                max_dist, max_pair = d, (genome_i, genome_j)

        if min_pair is not None:
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
