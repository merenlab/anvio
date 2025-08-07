# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import copy
import random
import networkx as nx

import anvio
import anvio.terminal as terminal

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

# ANCHOR - DirectedForce
class DirectedForce():
    """The first step is looking for open entrance points into the graph e.g. contigs beginning with singleton genes. Per definition
    of the maximum branching every node has to have a exactly single predecessor, but not necessarily a successor. We can exploit this
    definition by adding a artifical starting point connecting to all nodes in the graph. By adding weight to exactly one edge involving
    this starting point we force the algorithm in the direction of taking a specific node after the start as root, while also artificially
    closing all open entrance points into the graph. Since one edge from starting point is higher weighted the maximum branching has to take
    this route to created a MAXIMUM weighted tree representation. The tree to flow part uses this tree represenation starting from the last
    node (artificial stop node )of the longest branch (main branch) and tracking backwards connecting open leaves to the main branch.
    Connections that created loops in the original graph are reversed and marked to keep a steady flow in direction of the stop node. As soon
    as all nodes are fixed and connected, the starting and stop nodes are removed leaving a flow graph with a number of reversed edges. This
    calculation can be repeated by a number defined by the user to decrease the reversed edges optimizing the resulting flow graph. The simple
    example on the bottom shows green and red edges representing original (green) and artificially reversed (red) edges.

    While this procedure can hardly called a real graph algorithm as edges are reversed thereby heavily manipulating the original graph, it offers
    various advantages. Standard algorithms like topological sorting can be used to define x positions before reversing the changes edges back to
    the original direction, therefore functioning as a layout algorithm. Aside from that, I might be wrong but I think node to node distance can
    be estimated with reverse shortest path in this originally circular graph. When keeping the original direction of the edges in mind, while
    carefully tracing it should give a pretty good estimation of the distance (but it's not useful for me right now).
    """

    def __init__(self, seed=None, r=run, p=progress):

        self.seed = seed
        self.run = r
        self.progress = p


    def return_optimum_complexity(self, H, start_node=[], max_iterations=1):
        """One of the main functions it aims to decide which (in the best case)
        optimal set of edges need to be reversed to let the graph flow in one
        direction. It is a heuristic, therefore, suboptimal results are possible.
        The function can be used with parameters to e.g. increase the number of
        iterations to find a better solution.

        Parameters
        ==========
        nx.DiGraph()

        Returns
        =======
        list
        """

        max_weight = max([d for i, j, d in H.edges(data='weight')])
        changed_edges = []
        selfloops = list(nx.selfloop_edges(H))

        if selfloops:
            raise ConfigError("Looped graphs are not implemented in the algorithm please run remove edges from the networkx package first.")

        # TODO Give a warning if the node that is set or randomly picked is NOT a core synGC, this will likely mess up the graph a bit
        if max_iterations == 1 and not start_node:
            try:
                self.run.info_single("Low ressource mode: No start node was picked and the maximum number of iterations was set to 1.")
                # Suboptimal run, trying to find a sufficient starting point without a lot of ressources.
                G = nx.DiGraph(H)
                add_start = [node for node in G.nodes() if len(list(G.predecessors(node))) == 0]

                if not add_start:
                    G, M = self.find_maximum_branching(G)
                    add_start = [node for node in M.nodes() if len(list(M.predecessors(node))) == 0]
                    G = self.add_node_to_connector('START', add_start, G, max_weight)
                    M = self.add_node_to_connector('START', add_start, M, max_weight)
                else:
                    G = self.add_node_to_connector('START', add_start, G, max_weight)
                    # G, M, removed_nodes, removed_edges = self.find_maximum_branching(G)
                    G, M = self.find_maximum_branching(G)

                self.run.info_single("Solving complex graph.")
                changed_edges = self.run_tree_to_flow_network_algorithm(G, M, max_weight)
                # return(changed_edges, removed_nodes, removed_edges)
            except Exception as e:
                self.run.info_single(f"An error occured somewhere, we will try to continue anyway, but here is the error for posterity: '{e}'")

            return(changed_edges)

        else:
            # This area is reserved for multiple execution to find the optimum
            # I plan to add multithreading here
            # Will also be executed if you give a starting node to use in the graph
            all_starts = list(H.nodes())

            if set(start_node).issubset(all_starts):
                self.run.info_single("Low ressource mode: A start node was picked.")
                starting_list = start_node
            elif max_iterations >= len(all_starts):
                self.run.info_single("Medium ressource mode: No start node was picked but the maximum number of iterations is set less or equal then the number of nodes.")
                starting_list = all_starts
            else:
                self.run.info_single("High ressource mode: No start node was picked and the maximum number of iterations is set higher then number of nodes.")
                random.seed(self.seed)
                starting_list = random.sample(all_starts, max_iterations)

            min_complexity = len(H.nodes())

            iteration = 0
            self.run.info_single("Solving complex graph.")
            for start in starting_list:
                try:
                    self.run.info_single(f"Iteration {iteration}.")

                    G = nx.DiGraph(H)

                    # I Changed this for testing!
                    ebunch = [(node, start) for node in list(G.predecessors(start))]
                    G.remove_edges_from(ebunch)

                    add_start = [node for node in G.nodes() if len(list(G.predecessors(node))) == 0] + [start]
                    G = self.add_node_to_connector('START', add_start, G, max_weight)
                    G['START'][start]['weight'] = max_weight+10

                    # G, M, removed_nodes_remp, removed_edges_temp = self.find_maximum_branching(G)
                    G, M = self.find_maximum_branching(G)
                    changed_edges_temp = self.run_tree_to_flow_network_algorithm(G, M, max_weight)

                    complexity = len(changed_edges_temp)

                    if complexity < min_complexity:
                        min_complexity = complexity
                        changed_edges = changed_edges_temp + ebunch
                        # removed_nodes = removed_nodes_remp
                        # removed_edges = removed_edges_temp

                except Exception as e:
                    self.run.info_single(f"An error occured somewhere, we will try to continue anyway, but here is the error for posterity: '{e}'")

                iteration += 1

            # return(changed_edges, removed_nodes, removed_edges)
            return(changed_edges)


    def mean_M_path_weight(self, M, source, target):
        path = nx.shortest_path(G=M, source=source, target=target, weight='weight')
        path_weight = nx.path_weight(G=M, path=path, weight='weight')
        return(path_weight)


    def get_leaves(self, current_branch_successor, edmonds_graph_successors):

        leaves = set()
        while True:
            successors = edmonds_graph_successors[current_branch_successor]

            if successors:
                if len(successors) > 1:
                    for successor in successors:
                        leaves.update(self.get_leaves(successor, edmonds_graph_successors))
                    break
                else:
                    current_branch_successor = successors[0]
            else:
                leaves.update(set([current_branch_successor]))
                break

        return(leaves)


    def find_maximum_branching(self, G):
        """Slightly modified version of the find maximum aborescence algorithm in networkx
        by accessing the Edmonds class and therefore being able to access smaller aborescence
        in case no maximum aborescence is findable by sacrificing some information.

        Parameters
        ==========
        nx.DiGraph()

        Returns
        =======
        nx.DiGraph()
        """

        nx.set_edge_attributes(G, {(i, j): {'weight': d} for i, j, d in G.edges(data='weight')})
        nx.set_node_attributes(G, {k: {} for k in G.nodes()})

        edmonds = nx.algorithms.tree.branchings.Edmonds(G, seed=self.seed)
        M = edmonds.find_optimum(
            attr="weight",
            default=1,
            kind="max",
            style="arborescence",
            preserve_attrs=False,
            partition=None,
        )

        if not nx.algorithms.tree.recognition.is_arborescence(M):
            self.run.info_single('No maximum aborescence. Entering fallback mode.')
            M_sub_components = max(nx.weakly_connected_components(M), key=len)
            M_sub_graph = nx.DiGraph(M.subgraph(M_sub_components))
            if nx.algorithms.tree.recognition.is_arborescence(M_sub_graph):
                self.run.info_single('Found aborescence.')
                M = nx.DiGraph(M_sub_graph)

                removed_nodes = set(G.nodes()) - set(M.nodes())
                removed_edges = set(G.edges()) - set(M.edges())

                self.run.info_single(f'{len(removed_nodes)} nodes and {len(removed_edges)} edges removed to capture strongest signal.')
                G = nx.DiGraph(G.subgraph(M_sub_components))
            else:
                raise ConfigError("I'm very sorry to inform you that your data is not solvable by the current version of"
                                  "maximum flow. The fallback mode tried to solve your dataset by sacrificing some"
                                  "of the included information, but at this scale it will not lead to a acceptable result :(")
        # else:
        #     removed_nodes = []
        #     removed_edges = []

        nx.set_edge_attributes(M, {(i, j): {'weight': d} for i, j, d in G.edges(data='weight') if (i, j) in M.edges()})
        nx.set_node_attributes(M, {k: {} for k in G.nodes() if k in M.nodes()})

        # return(G, M, removed_nodes, removed_edges)
        return(G, M)


    def add_node_to_connector(self, new, connectors, G, weight, forward=True):
        for u in connectors:
            if forward == True:
                G.add_edge(
                    *(new, u),
                    weight=weight
                )
            else:
                G.add_edge(
                    *(u, new),
                    weight=weight
                )

        return(G)


    # TODO There is a bug somewhere, that breaks the graph with a node coming from reversing that is not connected to the start anymore. Is it time for a reimplementation already?
    def run_tree_to_flow_network_algorithm(self, G, M, max_weight):
        """Main algorithm function. It used the Maxmimum aborescence graph M and the original
        graph G to find the optimal edges to reverse. More documentation on this algoritm will
        be available as a picture.

        Parameters
        ==========
        nex.DiGraph(), nx.DiGraph(), float

        Returns
        =======
        list
        """

        M_edges = set(M.edges())
        M_nodes = set(M.nodes())

        G_edges = set(G.edges())
        G_nodes = set(G.nodes())

        M_removed_edges = G_edges - M_edges
        M_end = max([(self.mean_M_path_weight(M, 'START', node), node) for node in M_nodes if len(list(M.successors(node))) == 0])[1]

        G = self.add_node_to_connector('STOP', [M_end], G, max_weight, forward=False)
        M = self.add_node_to_connector('STOP', [M_end], M, max_weight, forward=False)

        M_nodes.add('STOP')
        M_predecessors = {M_node: list(M.predecessors(M_node))[0] for M_node in M_nodes if M_node != 'START'}
        M_successors = {M_node: list(M.successors(M_node)) for M_node in M_nodes}
        G_successors = {G_node: list(G.successors(G_node)) for G_node in G_nodes}

        # TODO multithreading?
        self.progress.new("Calculating initial distances")
        M_distances = {}
        for z, M_node in enumerate(M_nodes):

            M_distances[M_node] = self.mean_M_path_weight(M, 'START', M_node)
            self.progress.update(f"{str(z+1).rjust(len(str(len(M_nodes))), ' ')} / {len(M_nodes)}")

        self.progress.end()

        i = 0
        resolved_nodes = set(['STOP'])

        next_node = ''
        x = 0

        self.progress.new("Solving complex graph")
        while len(resolved_nodes) != len(G_nodes) + 1 or x < 1:
            if len(resolved_nodes) == len(G_nodes) + 1:
                x += 1

            i += 1

            visited_nodes = set(['STOP'])
            current_node = 'STOP'
            while len(visited_nodes) != len(G_nodes) + 1:

                self.progress.update(f"{str(len(resolved_nodes)).rjust(len(str(len(G_nodes) + 1)), ' ')} / {len(G_nodes) + 1}")

                # if next_node:
                #     current_branch_root = next_node
                #     next_node = ''
                # else:
                current_branch_root = M_predecessors[current_node]

                current_forward_connected = []
                current_backward_connected = []
                successor_branch_leaves = set()

                for current_branch_successor in M_successors[current_branch_root]:
                    if current_branch_successor not in visited_nodes and current_branch_successor != current_node:
                        successor_branch_leaves.update(self.get_leaves(current_branch_successor, M_successors))

                if not successor_branch_leaves:
                    current_node = current_branch_root
                else:
                    current_node = max([(M_distances[successor_branch_leaf], successor_branch_leaf) for successor_branch_leaf in successor_branch_leaves])[1]

                if current_node in resolved_nodes:
                    connected = True
                else:
                    connected = False

                if connected != True or x > 0:
                    for current_node_successor in G_successors[current_node]:

                        if current_node_successor in resolved_nodes:

                            if current_node_successor in nx.ancestors(M, current_node) or (current_node_successor not in visited_nodes and current_node_successor not in resolved_nodes):
                                if (current_node, current_node_successor) in M_removed_edges:
                                    current_backward_connected.append(current_node_successor)

                            else:
                                if (current_node, current_node_successor) in M_removed_edges:
                                    current_forward_connected.append(current_node_successor)
                                    connected = True

                                else:
                                    connected = True

                    if connected == False:
                        if len(list(G.successors(current_node))) == 0:
                            data = {
                                'weight': max_weight
                            }

                            if M.has_edge(current_node, 'STOP'):
                                M[current_node]['STOP']['weight'] += data['weight']
                            else:
                                M.add_edge(current_node, 'STOP', **data)
                            connected = True

                    if connected == True:
                        for current_forward in current_forward_connected:
                            data = G.get_edge_data(current_node, current_forward)

                            if M.has_edge(current_node, current_forward):
                                M[current_node][current_forward]['weight'] += data['weight']
                            else:
                                M.add_edge(current_node, current_forward, **data)

                            M_removed_edges.remove((current_node, current_forward))

                        for current_backward in current_backward_connected:
                            data = G.get_edge_data(current_node, current_backward)

                            if M.has_edge(current_backward, current_node):
                                M[current_backward][current_node]['weight'] += data['weight']
                            else:
                                M.add_edge(current_backward, current_node, **data)

                            M_removed_edges.remove((current_node, current_backward))

                        resolved_nodes.add(current_node)

                    else:
                        if current_backward_connected:

                            number = max([(G.get_edge_data(current_node, backward)['weight'], i) for (i, backward) in enumerate(current_backward_connected)])[1]
                            current_connector = current_backward_connected[number]
                            next_node = current_node

                            # while current_connector:
                            # while next_node not in resolved_nodes:
                            while True:

                                # print(next_node, current_connector)

                                if G.has_edge(M_predecessors[current_node], current_node):

                                    current_node = next_node

                                    data = G.get_edge_data(current_node, current_connector)

                                    M.remove_edge(M_predecessors[current_node], current_node)

                                    if M.has_edge(current_connector, current_node):
                                        M[current_connector][current_node]['weight'] += data['weight']
                                    else:
                                        M.add_edge(current_connector, current_node, **data)

                                    M_removed_edges.remove((current_node, current_connector))
                                    M_removed_edges.add((M_predecessors[current_node], current_node))

                                    M_successors[M_predecessors[current_node]].remove(current_node)
                                    M_successors[current_connector] += [current_node]

                                    next_node = M_predecessors[current_node]

                                    M_predecessors.pop(current_node, None)
                                    M_predecessors[current_node] = current_connector

                                    M_distances[current_node] = self.mean_M_path_weight(M, 'START', current_node)

                                    if current_node in resolved_nodes:
                                        break
                                    else:
                                        resolved_nodes.add(current_node)
                                        current_connector = current_node

                                else:
                                    break

                                # if M.out_degree(current_node) == 0:
                                # current_connector = current_node
                                # else:
                                # current_connector = ''
                                # print(next_node)

                            break

                visited_nodes.add(current_node)

                if not nx.is_directed_acyclic_graph(M):
                    raise ConfigError("Oh no. It looks like your graph is so complex or includes a motif I haven't seen before"
                                      "therefore the reattachement algorithm itself included a loop to the graph. We had multiple"
                                      "sanity checks to prevent this but unfortunatly nobody is perfect. We will include more"
                                      "checks in the next version. Sorry :/")

        self.progress.end()

        remaining_stops = [node for node in M.nodes() if M.out_degree(node) == 0 and node != 'STOP']

        for stop in remaining_stops:

            data = {
                'weight': max_weight
            }

            if M.has_edge(stop, 'STOP'):
                M[stop]['STOP']['weight'] += data['weight']
            else:
                M.add_edge(stop, 'STOP', **data)

            M_successors[stop] += ['STOP']

        if not nx.is_directed_acyclic_graph(M):
            raise ConfigError("Oh no. It looks like your graph is so complex or includes a motif I haven't seen before"
                              "therefore the reattachement algorithm itself included a loop to the graph. We had multiple"
                              "sanity checks to prevent this but unfortunatly nobody is perfect. We will include more"
                              "checks in the next version. Sorry :/")
        else:
            self.run.info_single("No loops. Roger roger and ready to go.")

        v = 0
        w = 0
        x = 0
        y = 0
        changed_edges_new = []
        G_all_edges = set(G.edges())
        for i, j in M.edges():
            if G.has_edge(i,j) and G.has_edge(j,i):
                # print('double sided.')
                G_all_edges.remove((i,j))
                G_all_edges.remove((j,i))
                changed_edges_new += [(j,i)]
                v += 1
            elif G.has_edge(j,i):
                # print('reverse sided.')
                G_all_edges.remove((j,i))
                changed_edges_new += [(j,i)]
                w += 1
            elif G.has_edge(i,j):
                # print('no change.')
                G_all_edges.remove((i,j))
                x += 1
            else:
                # print('problem.')
                y += 1

        self.run.info_single(f"{x} edges can be kept in the original direction.")
        self.run.info_single(f"{w} edges have to reversed to capture maximum force.")
        self.run.info_single(f"{v} edges seem to be double sided. Double sided edges are a sign for potential inversions.")
        self.run.info_single(f"{y} problems.")

        return(changed_edges_new)