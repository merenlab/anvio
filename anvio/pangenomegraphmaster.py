# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import numpy as np
import pandas as pd
import networkx as nx
import itertools as it
import matplotlib.pyplot as plt

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
            'type': '',
            'group': '',
            'layer': {}
        }
        self.edge_standard_attributes = {
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

        self.run.warning(None, header="Running connectivity check on pangenome graph.", lc="green")

        connectivity = nx.is_connected(self.graph.to_undirected())

        if connectivity == False:

            weakly_components = list(nx.weakly_connected_components(self.graph))
            self.run.info_single(f"The pangenome graph contains {len(weakly_components)} "
                                 f"independant components. The reason is likely to be the "
                                 f"fragmented nature of at least one of the genomes "
                                 f"present in the dataset. This happens quite often if the dataset "
                                 f"contains draft genomes and should not bother you to much.")

            subgraph = self.graph.subgraph(max(weakly_components, key=len))
            self.graph = nx.DiGraph(subgraph)
            self.run.info_single(f"Keeping the longest conntected subgraph with {len(self.graph.nodes())} nodes and {len(self.graph.edges())} edges.")

            connectivity = nx.is_connected(self.graph.to_undirected())
            if connectivity == True:
                self.run.info_single("The pangenome graph is now a connected cyclic graph.")
            else:
                raise ConfigError("Looks like the graph is still fragmented, please check the data"
                                  "for at least some level of consistency.")

        self.run.info_single("Done.")


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


    def summarize(self):
        """This is the main summary function for pangenome graphs. Input is the self

        Parameters
        ==========
        self.graph

        Returns
        =======
        pd.DataFrame
        """

        # self.run.warning(None, header="Generate pangenome graph summary tables", lc="green")

        genome_names = set(it.chain(*[list(d.keys()) for node, d in self.graph.nodes(data='gene_calls')]))

        if not nx.is_directed_acyclic_graph(self.graph):
            raise ConfigError("Cyclic graphs, are not implemented.")

        all_positions = sorted(set([data['position'][0] for node, data in self.graph.nodes(data=True)]))
        all_positions_min = min(all_positions)
        all_positions_max = max(all_positions)

        core_positions = sorted(set([data['position'][0] for node, data in self.graph.nodes(data=True) if len(data['gene_calls'].keys()) == len(genome_names)]))

        regions_dict = {}

        if core_positions:
            core_position_min = min(core_positions)
            core_position_max = max(core_positions)

            overlap = [i for i in range(all_positions_min, core_position_min)] + [i for i in range(core_position_max+1, all_positions_max+1)]
            regions_dict |= {pos:0 for pos in overlap}
            regions_id = 1

            # print(overlap)

            regions_dict |= {pos:-1 for pos in core_positions}

            core_positions_pairs = map(tuple, zip(core_positions, core_positions[1:]))
        else:
            regions_id = 0
            core_positions_pairs = [(all_positions_min, all_positions_max)]

        for (i,j) in core_positions_pairs:
            if i != j-1:
                regions_dict |= {i: regions_id for i in range(i+1,j)}
                regions_id += 1

        # print(regions_dict)

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

        # print(regions_info_dict)
        # print(node_regions_dict)
        # all_ghost_max_y_positions = {}

        # TODO height depends on the position of empty edges (should be fixed already)
        i = 0
        regions_summary_dict = {}
        for region_id, values_list in regions_info_dict.items():
            if region_id != -1:
                region_x_positions = [item[0] for item in values_list]
                region_y_positions = [item[1] for item in values_list]

                genomes_sets = [item[3] for item in values_list]
                genomes_involved = set(it.chain(*genomes_sets))
                weight = len(genomes_involved)

                region_x_positions_min = min(region_x_positions)
                region_x_positions_max = max(region_x_positions)

                length = region_x_positions_max - region_x_positions_min + 1
                quantity = len(values_list)
                height = len(set(region_y_positions))
                max_density = (height*length)

                if height <= 1:
                    if weight < len(genome_names)*0.5:
                        motif = 'INS'
                    elif weight == len(genome_names)*0.5:
                        motif = 'INDEL'
                    else:
                        motif = 'DEL'
                else:
                    motif = 'HVR'

                regions_summary_dict[i] = {
                    'region_id': region_id,
                    'height': height,
                    'length': length,
                    'weight': weight,
                    'motif': motif,
                    'x_min': region_x_positions_min,
                    'x_max': region_x_positions_max,
                    'quantity': len(values_list),
                    'density': quantity/max_density if max_density != 0 else 0
                }
                i += 1
            else:
                regions_summary_dict[i] = {
                    'region_id': -1,
                    'height': 1,
                    'length': -1,
                    'weight': len(genome_names),
                    'motif': 'BB',
                    'x_min': -1,
                    'x_max': -1,
                    'quantity': -1,
                    'density': 1.0
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

        return(region_sides_df, nodes_df, gene_calls_df)


    def reverse_edges(self, changed_edges):
        for (edge_i, edge_j) in changed_edges:
            directions = {genome:'L' for genome, direction in self.graph[edge_i][edge_j]['directions'].items()}
            weight = self.graph[edge_i][edge_j]['weight']

            edge_attributes = {
                'weight': weight,
                'directions': directions
            }

            self.add_edge_to_graph(edge_j, edge_i, edge_attributes)
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
        self.run.warning(None, header="Calculate synteny distance dendrogram", lc="green")
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

            self.run.info_single(f"d({genome_i},{genome_j}) = {round(X[i][j], 3)}")

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
