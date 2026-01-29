# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import math
import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
import matplotlib.pyplot as plt
from collections import Counter
from statistics import mean, multimode

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree

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

    Definition of PangenomeGraphManager object edges:
    i: syn_cluster ----> 'GC1_1'
    j: syn_cluster ----> 'GC2_1'
    attributes:
        weight: --------> 2
        active: --------> True
        directions: ----> {G1: 'R', G2: 'R'}
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
            'alignment': ''
        }
        self.edge_standard_attributes = {
            'name': '',
            'weight': 0.0,
            'active': True,
            'directions': {},
            'route': [],
            'length': 0
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
            self.graph[syn_cluster_i][syn_cluster_j]['weight'] += attributes['weight']
            self.graph[syn_cluster_i][syn_cluster_j]['directions'].update(attributes['directions'])


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


    def summarize(self):
        """This is the main summary function for pangenome graphs. Input is the self

        Parameters
        ==========
        self.graph

        Returns
        =======
        pd.DataFrame
        """

        genome_names = set(it.chain(*[list(d.keys()) for node, d in self.graph.nodes(data='gene_calls')]))

        if not nx.is_directed_acyclic_graph(self.graph):
            raise ConfigError("Cyclic graphs, are not implemented.")

        all_positions = sorted(set([data['position'][0] for node, data in self.graph.nodes(data=True)]))
        all_positions_min = min(all_positions)
        all_positions_max = max(all_positions)

        regions_dict = {}
        edge_pos_dict = {}
        core_pos_num_y = {}
        for edge_i, edge_j, data in self.graph.edges(data=True):

            position_tuples = []
            edge_i_x, edge_i_y = self.graph.nodes[edge_i]['position']
            edge_j_x, edge_j_y = self.graph.nodes[edge_j]['position']

            position_tuples += [(edge_i_x, edge_i_y)]
            position_tuples += [(edge_j_x, edge_j_y)]

            for route_x, route_y in data['route']:
                position_tuples += [(route_x, route_y)]

            for pos_x, pos_y in position_tuples:
                if pos_x not in edge_pos_dict:
                    edge_pos_dict[pos_x] = set([pos_y])
                else:
                    edge_pos_dict[pos_x].add(pos_y)

        core_positions = sorted(set([data['position'][0] for node, data in self.graph.nodes(data=True) if len(data['gene_calls'].keys()) == len(genome_names)]))

        for core_position in core_positions:
            core_pos_num_y[core_position] = len(edge_pos_dict[core_position])

        core_pos_num_y_values = list(core_pos_num_y.values())
        mode_position = max(set(core_pos_num_y_values), key=core_pos_num_y_values.count)

        for core_position, num_y_positions in core_pos_num_y.items():
            if num_y_positions > mode_position:
                if core_position in core_positions:
                    core_positions.remove(core_position)
                    self.run.info_single(f"Position {core_position} is probably falsely annotated as core and might be a rearranged synteny gene cluster")

        if core_positions:
            core_position_min = min(core_positions)
            core_position_max = max(core_positions)

            overlap = [i for i in range(all_positions_min, core_position_min)] + [i for i in range(core_position_max+1, all_positions_max+1)]
            regions_dict |= {pos: 0 for pos in overlap}
            regions_id = 1

            # print(overlap)

            # regions_dict |= {pos:-1 for pos in core_positions}

            core_positions_pairs = map(tuple, zip(core_positions, core_positions[1:]))
        else:
            regions_id = 0
            core_positions_pairs = [(all_positions_min, all_positions_max)]

        core_chain = []
        for (i, j) in core_positions_pairs:
            if i != j-1:

                core_chain += [i]
                regions_dict |= {m: regions_id for m in core_chain}
                regions_id += 1

                regions_dict |= {k: regions_id for k in range(i+1, j)}
                regions_id += 1

                if j == core_positions[-1]:
                    core_chain = [j]
                else:
                    core_chain = []
            else:
                core_chain += [i]

                if j == core_positions[-1]:
                    core_chain += [j]

        regions_dict |= {n: regions_id for n in core_chain}

        regions_info_dict = {}
        node_regions_dict = {}
        i = 0
        for node, data in self.graph.nodes(data=True):
            node_x_position = data['position'][0]
            node_y_position = data['position'][1]
            genomes = list(data['gene_calls'].keys())
            if node_x_position in regions_dict.keys():
                region_id = regions_dict[node_x_position]
                node_regions_dict[i] = {
                    'syn_cluster': node,
                    'x': node_x_position,
                    'y': node_y_position,
                    'region_id': region_id
                }
                i += 1

                if region_id not in regions_info_dict:
                    regions_info_dict[region_id] = [(node_x_position, node_y_position, node, genomes)]
                else:
                    regions_info_dict[region_id] += [(node_x_position, node_y_position, node, genomes)]

        i = 0
        regions_summary_dict = {}
        num_genomes = len(genome_names)

        for region_id, values_list in dict(sorted(regions_info_dict.items())).items():

            genomes_sets = [item[3] for item in values_list]
            genomes_involved = set(it.chain(*genomes_sets))

            weight = len(genomes_involved)

            nodes_sets = [item[2] for item in values_list]

            num_gene_clusters = len(set([item.rsplit('_', 1)[0] for item in nodes_sets]))

            region_x_positions = [item[0] for item in values_list]

            region_x_positions_min = min(region_x_positions)
            region_x_positions_max = max(region_x_positions)

            if region_x_positions_min in set(core_positions) and region_x_positions_max in set(core_positions):
                P = 0
            elif region_x_positions_min not in set(core_positions) and region_x_positions_max not in set(core_positions):

                if region_x_positions_min != 0:
                    prior_core_region_id = regions_dict[region_x_positions_min - 1]
                    prior_core_region = regions_info_dict[prior_core_region_id]

                    predecessor = [item[2] for item in prior_core_region if item[0] == region_x_positions_min - 1][0]

                    P = sum([len(list(self.graph.successors(node))) - 1 for node in nodes_sets + [predecessor]])
                else:
                    P = 0
            else:
                print('Summary error.')

            genome_occurences = list(it.chain(*genomes_sets))
            counts = Counter(genome_occurences)

            K = len(nodes_sets) #Num synteny gene clusters
            T = len(genome_occurences)

            p = (1/K) * sum([len(genome_set)/weight for genome_set in genomes_sets])
            # Option 2:
            # p = (1/K) * sum([len(genome_set) for genome_set in genomes_sets])

            var = (1/K) * sum([(len(genome_set)/weight - p)**2 for genome_set in genomes_sets])
            # Option 2:
            # var = (1/K) * sum([(len(genome_set) - p)**2 for genome_set in genomes_sets])

            if p != 1:
                diversity = 1 - math.sqrt(var / (p * (1 - p)))
                # Option 2:
                # diversity = (weight - p) / (weight - 1) * (1 - math.sqrt(var / ((weight - 1)**2)))

            else:
                diversity = 0

            genome_counts = {item: counts.get(item, 0) for item in genomes_involved}
            values = list(genome_counts.values())

            max_expansion = max(values)
            min_expansion = 0 if len(genomes_involved) != len(genome_names) else min(values)
            mean_expansion = sum(values) / len(values)

            expansion = math.sqrt(max_expansion * mean_expansion)
            # Option 2:
            # expansion = max_expansion

            if P == 0:
                motif = 'BB'
            elif P == 1 and min_expansion == 0:
                motif = 'INDEL'
            else:
                motif = 'VR'

            complexity = P / weight
            # Option 2:
            # complexity = P / num_genomes

            regions_summary_dict[i] = {
                'region_id': region_id,
                'motif': motif,
                'x_min': region_x_positions_min,
                'x_max': region_x_positions_max,
                'num_synteny_gene_clusters': K,
                'num_gene_clusters': num_gene_clusters,
                'num_gene_calls': T,
                'pathways': P,
                'complexity': complexity,
                'expansion': expansion,
                'mean_expansion': mean_expansion,
                'max_expansion': max_expansion,
                'min_expansion': min_expansion,
                'prevalence': p,
                'prevalence_variance': var,
                'diversity': diversity,
                'weight': weight
            }
            i += 1

        i = 0
        gene_calls_dict = {}
        for node, data in self.graph.nodes(data='gene_calls'):
            for genome, gene_call in data.items():
                gene_calls_dict[i] = {'syn_cluster': node, 'genome': genome, 'gene_caller_id': gene_call}
                i += 1

        gene_calls_df = pd.DataFrame.from_dict(gene_calls_dict, orient='index').set_index(['genome', 'gene_caller_id'])
        nodes_df = pd.DataFrame.from_dict(node_regions_dict, orient='index').set_index('syn_cluster')
        region_sides_df = pd.DataFrame.from_dict(regions_summary_dict, orient='index').set_index('region_id')
        region_sides_df[['composite_variability_score', 'complexity_normalized', 'expansion_normalized', 'weight_factor']] = region_sides_df.apply(PangenomeGraphManager.composite_variability_score(region_sides_df, num_genomes), axis=1, result_type='expand')

        return(region_sides_df, nodes_df, gene_calls_df)


    @staticmethod
    def composite_variability_score(df, N):

        C_max = df['complexity'].max()
        # Option 2:
        # C_min = df['complexity'].min()

        E_max = df['expansion'].max()
        # Option 2:
        # E_min = df['max_expansion'].min()

        W_median = df['weight'].median()

        def log_min_max_normalize(X, X_min, X_max):
            if X_min == X_max:
                return 0.0
            X_norm = (math.log(1 + X) - math.log(1 + X_min)) / (math.log(1 + X_max) - math.log(1 + X_min))
            return(X_norm)

        def fixed_baseline_log_scaling(X, X_max):
            if X_max == 0.0:
                return(0.0)
            else:
                X_norm = (math.log(1 + X)) / (math.log(1 + X_max))
                return(X_norm)

        def func(row):

            C = row['complexity']
            E = row['expansion']
            W = row['weight']
            D = row['diversity']

            C_norm = fixed_baseline_log_scaling(C, C_max)
            # Option 2:
            # C_norm = log_min_max_normalize(C, C_min, C_max)

            E_norm = fixed_baseline_log_scaling(E, E_max)
            # Option 2:
            # E_norm = log_min_max_normalize(E, E_min, E_max)

            W_f = W/N
            # Option 2:
            # W_f = (1 - math.e**(-W/W_median))
            # Option 3:
            # W_f = (1 - math.e**(-W/W_median)) / (1 - math.e**(-N/W_median))

            CVS = (C_norm * E_norm * D)**(1/3) * W_f

            return([CVS, C_norm, E_norm, W_f])

        return(func)


    def reverse_edges(self, changed_edges):
        for (edge_i, edge_j) in changed_edges:
            directions = {genome:'L' for genome, direction in self.graph[edge_i][edge_j]['directions'].items()}
            # weight = self.graph[edge_i][edge_j]['weight']

            edge_attributes_ij = self.graph[edge_i][edge_j]
            edge_attributes_ij['directions'] = directions

            if self.graph.has_edge(edge_j, edge_i):
                edge_attributes_ji = self.graph[edge_j][edge_i]

                print('one double sided edge')
                edge_attributes_ji['weight'] += edge_attributes_ij['weight']
                edge_attributes_ji['directions'].update(edge_attributes_ij['directions'])
            else:
                self.add_edge_to_graph(edge_j, edge_i, edge_attributes_ij)
                self.graph.remove_edge(edge_i, edge_j)


    def set_node_positions(self, node_positions):
        for node in self.graph.nodes():
            self.graph.nodes()[node]['position'] = node_positions[node]


    def set_node_groups(self, node_groups):
        for node in self.graph.nodes():
            if node in node_groups:
                self.graph.nodes()[node]['group'] = node_groups[node]
            else:
                self.graph.nodes()[node]['group'] = ''


    def set_edge_positions(self, edge_positions):
        for edge_i, edge_j in self.graph.edges():
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


    def calculate_graph_distance(self, output_dir=''):
        self.run.warning(None, header="COMPUTING GRAPH-BASED GENOME DISTANCES", lc="green")
        genome_names = list(set(it.chain(*[list(d.keys()) for node, d in self.graph.nodes(data='gene_calls')])))

        X = np.zeros([len(genome_names), len(genome_names)])
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
                if genome_i in data['directions'].keys() and genome_j in data['directions'].keys():
                    edges_similar += 1
                elif genome_i in data['directions'].keys() or genome_j in data['directions'].keys():
                    edges_unsimilar += 1

            i = genome_names.index(genome_i)
            j = genome_names.index(genome_j)

            elements_similar = nodes_similar + edges_similar
            elements_unsimilar = nodes_unsimilar + edges_unsimilar

            X[i][j] = elements_unsimilar / (elements_similar + elements_unsimilar)
            X[j][i] = elements_unsimilar / (elements_similar + elements_unsimilar)

            self.run.info_single(f"d({genome_i},{genome_j}) = {round(X[i][j], 3)}", cut_after=None)

        # Check if all distances are zero (identical genomes)
        if np.all(X == 0):
            self.run.warning("All pairwise distances are zero (genomes are identical). No dendrogram will be generated.")
            self.run.info_single("Done.")
            return ''

        condensed_X = squareform(X)
        Z = linkage(condensed_X, 'ward')

        if output_dir:
            fig = plt.figure(figsize=(25, 10))
            ax = plt.axes()
            dendrogram(Z, ax=ax, labels=genome_names, orientation='right')
            plt.tight_layout()
            fig.savefig(os.path.join(output_dir, 'synteny_distance_dendrogram.svg'))
            plt.close(fig)

            distance_matrix = pd.DataFrame(X, index=genome_names, columns=genome_names)
            distance_matrix.to_csv(os.path.join(output_dir, 'synteny_distance_matrix.tsv'), sep='\t')

            self.run.info_single(f"Exported distance dendrogram to {os.path.join(output_dir, 'synteny_distance_dendrogram.svg')}.")
            self.run.info_single(f"Exported distance matrix to {os.path.join(output_dir, 'synteny_distance_matrix.tsv')}.")

        tree = to_tree(Z, False)
        newick = clustering.get_newick(tree, tree.dist, genome_names)

        with open(os.path.join(output_dir, 'synteny_distance_dendrogram.newick'), "w") as text_file:
            text_file.write(newick)
        self.run.info_single(f"Exported newick tree to {os.path.join(output_dir, 'synteny_distance_dendrogram.tree')}.")
        self.run.info_single("Done.")

        return(newick)

    # TODO still empty..
    def generate_hybrid_genome(self, output_dir):
        self.run.info_single("The output hybrid genome function exists on paper but we haven't implemented it yet. sorry.")
