# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Class for the pangenome graph topological layout algorithm.
"""

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

# ANCHOR - TopologicalLayout
class TopologicalLayout():
    """A class to calculate x,y positions on the nodes of a graph as well as
    group those nodes together in case they follow a one to one connection
    pattern.
    """

    def __init__(self, r=run, p=progress):

        self.run = r
        self.progress = p


    # def run_synteny_layout_algorithm(self, F, gene_cluster_grouping_threshold=-1, groupcompress=1.0, ungroup_open=[], ungroup_close=[]):
    def run_synteny_layout_algorithm(self, F, gene_cluster_grouping_threshold=-1, groupcompress=1.0):
        """One of the main functions. Describe!

        Parameters
        ==========

        Returns
        =======

        """

        L = nx.DiGraph(F)

        nx.set_edge_attributes(L, {(i, j): {'weight': -1} for i, j in L.edges()})
        nx.set_node_attributes(L, {k: {'genomes': list(d.keys())} for k, d in L.nodes(data='gene_calls')})

        if not nx.is_directed_acyclic_graph(L):
            raise ConfigError("Cyclic graphs, are not implemented, sorry.")

        x_list = {}
        positions = {}
        removed = set()
        edges = {}
        grouping = {}
        offset = {}

        add_start = set()
        add_stop = set()
        for node in L.nodes():
            if len(list(L.successors(node))) == 0:
                add_stop.add(node)
            if len(list(L.predecessors(node))) == 0:
                add_start.add(node)

        for start in add_start:
            L.add_edge(*('START', start), weight=-1)

        for stop in add_stop:
            L.add_edge(*(stop, 'STOP'), weight=-1)

        for x, generation in enumerate(nx.topological_generations(L)):
            x_list[x] = generation
            for node in generation:
                positions[node] = (x, -1)

        global_x = x
        global_y = 0

        layout_graph_nodes = list(L.nodes())
        layout_graph_successors = {layout_graph_node: list(L.successors(layout_graph_node)) for layout_graph_node in layout_graph_nodes}

        ghost = 0
        for x in range(global_x-1, 0, -1):
            for node in x_list[x]:
                node_x_position = positions[node][0]

                change = []
                for successor in layout_graph_successors[node]:
                    if successor != 'STOP':
                        successor_x_position = positions[successor][0]

                        if successor_x_position <= node_x_position:
                            raise ConfigError(f"The node {node} is succeded by the node {successor} which is in front of the node {node}"
                                              f"that does not make a lot of sense and should not happen. We are sorry for the inconvenience :/")
                        else:
                            change.append((successor_x_position, successor))

                if change:
                    if min(change)[0] > 1:
                        new_x_position = min([(x, n) for (x, n) in change if x > 1])[0] - 1
                        positions[node] = (new_x_position, -1)

                    for (position, extend_successor) in change:
                        x_difference = position - new_x_position

                        path_list = [node]

                        for i in range(1, x_difference):
                            path_list += ['GHOST_' + str(ghost)]
                            positions['GHOST_' + str(ghost)] = (new_x_position + i, -1)
                            ghost += 1

                        path_list += [extend_successor]

                        edges[(node, extend_successor)] = (path_list, L.edges[node, extend_successor])
                        L.remove_edge(node, extend_successor)

                        if len(path_list) == 2:
                            L.add_edges_from(map(tuple, zip(path_list, path_list[1:])), weight=-1)
                        else:
                            L.add_edges_from(map(tuple, zip(path_list, path_list[1:])), weight=-0.5)

        for i, j in L.edges():
            if positions[j][0] - positions[i][0] != 1 and i != 'START' and j != 'STOP' and (i,j) not in removed:
                raise ConfigError(f"Hmmm. This situation would create a very weird looking connection."
                                  f"The ede {(i, j)} is longer than it should be. I don't know what created"
                                  f"this but we will work on a solution on the next release. Sorry :(")

        longest_path = nx.bellman_ford_path(G=L, source='START', target='STOP', weight='weight')
        m = set(longest_path)

        dfs_list = list(nx.dfs_edges(L, source='START'))

        group = 0
        groups = {}
        groups_rev = {}
        # keep = False

        # TODO Currently no seperation between unequal genome context, but is it needed? I guess this is fine.
        if gene_cluster_grouping_threshold == -1:
            self.run.info_single("Setting algorithm to 'no grouping'")
        else:
            self.run.info_single(f"Setting algorithm to 'Grouping single connected chains size > {gene_cluster_grouping_threshold}'")

        for node_v, node_w in dfs_list:

            # if node_v in ungroup_open:
            #     keep = True
            #     if not keep:

            if node_v != 'START' and node_w != 'STOP' and L.in_degree(node_v) == 1 and L.out_degree(node_v) == 1 and L.in_degree(node_w) == 1 and L.out_degree(node_w) == 1:
                if node_v not in groups_rev.keys():
                    group_name = group
                    groups[group_name] = [node_v, node_w]
                    groups_rev[node_v] = group_name
                    groups_rev[node_w] = group_name
                    group += 1
                else:
                    group_name = groups_rev[node_v]
                    groups[group_name] += [node_w]
                    groups_rev[node_w] = group_name

            # if node_w in ungroup_close:
            #     keep = False

        for group_name, group_nodes in groups.items():

            group_nodes = [node for node in group_nodes if not node.startswith('GHOST_')]
            former_genome = set()
            former_type = ''
            label = 'GCG_' + str(group_name).zfill(8)
            condense_nodes = []
            for node in group_nodes:

                current_genome = set(L.nodes()[node]['gene_calls'].keys())
                current_type = L.nodes()[node]['type']

                if not former_genome and not former_type:
                    condense_nodes += [node]
                elif former_genome == current_genome and former_type == current_type:
                    condense_nodes += [node]
                else:
                    if len(condense_nodes) >= gene_cluster_grouping_threshold and gene_cluster_grouping_threshold != -1:
                        grouping[label] = condense_nodes

                    condense_nodes = [node]
                    label = 'GCG_' + str(group).zfill(8)
                    group += 1

                former_genome = current_genome
                former_type = current_type

            if len(condense_nodes) >= gene_cluster_grouping_threshold and gene_cluster_grouping_threshold != -1:
                grouping[label] = condense_nodes

            condense_nodes = [node]
            label = 'GCG_' + str(group).zfill(8)
            group += 1

        self.run.info_single(f"Grouped {len(sum(grouping.values(), []))} nodes in {len(grouping.keys())} groups")

        L.remove_nodes_from(['START', 'STOP'])

        branches = {}
        sortable = []
        for g in groups.keys():
            branch = groups[g]

            if not set(branch).issubset(m):
                start = positions[branch[0]][0]
                length = len(branch)
                impact = - sum([len(L.nodes()[br]['gene_calls'].keys()) if 'gene_calls' in L.nodes()[br].keys() else 0 for br in branch])

                if start in branches.keys():
                    if length in branches[start].keys():
                        if impact in branches[start][length].keys():
                            if branches[start][length][impact].keys():
                                num = max(branches[start][length][impact].keys()) + 1
                            else:
                                num = 1

                            branches[start][length][impact][num] = branch
                        else:
                            num = 1
                            branches[start][length][impact] = {num: branch}
                    else:
                        num = 1
                        branches[start][length] = {impact: {num: branch}}
                else:
                    num = 1
                    branches[start] = {length: {impact: {num: branch}}}

                sortable += [(start, length, impact, num)]

        left_nodes = set(L.nodes()) - set(groups_rev.keys())
        for n in left_nodes:

            if not set([n]).issubset(m):
                start = positions[n][0]
                length = 1
                impact = - len(L.nodes()[n]['gene_calls'].keys()) if 'gene_calls' in L.nodes()[n].keys() else 0

                if start in branches.keys():
                    if length in branches[start].keys():
                        if impact in branches[start][length].keys():
                            if branches[start][length].keys():
                                num = max(branches[start][length][impact].keys()) + 1
                            else:
                                num = 1

                            branches[start][length][impact][num] = [n]
                        else:
                            num = 1
                            branches[start][length][impact] = {num: [n]}
                    else:
                        num = 1
                        branches[start][length] = {impact: {num: [n]}}
                else:
                    num = 1
                    branches[start] = {length: {impact: {num: [n]}}}

                sortable += [(start, length, impact, num)]

        used = set()
        finished = set()

        y_new = 0
        for node in longest_path:
            x_pos = positions[node][0]
            positions[node] = (x_pos, y_new)
            used.add((x_pos, y_new))
            finished.add(node)

        stack = [longest_path]
        while stack:

            current = stack[0]

            remove = True
            for i,j,k,l in sorted(sortable, key=lambda x: (x[1], x[0], x[2]), reverse = False):
                branch = branches[i][j][k][l]

                # print(branch)

                branch_pred = set(L.predecessors(branch[0]))
                branch_succ = set(L.successors(branch[-1]))
                if (not branch_pred.isdisjoint(set(current))) or (not branch_succ.isdisjoint(set(current))) or (not branch_pred.isdisjoint(set(current)) and not branch_succ.isdisjoint(set(current))):

                    remove = False
                    sortable.remove((i,j,k,l))
                    y_new = max(sum([[positions[ypred][1] for ypred in branch_pred], [positions[ysucc][1] for ysucc in branch_succ]], []))

                    stack = [branch] + stack
                    while True:
                        repeat = False
                        for xnew in range(i, i+j):
                            if (xnew, y_new) in used:
                                repeat = True

                        if repeat == False:
                            break
                        else:
                            y_new += 1

                    for node in branch:
                        x_pos = positions[node][0]
                        positions[node] = (x_pos, y_new)
                        used.add((x_pos, y_new))
                        finished.add(node)

                        global_y = y_new if y_new > global_y else global_y

                    break

            if remove == True:
                stack.remove(current)

        if len(set(positions.values())) != len(positions.values()):
            raise ConfigError("No no no no. Something went very wrong here. Some nodes overlap in the UI."
                              "We don't want this, we definitely don't want this...")

        edge_positions = {}
        for edge_i, edge_j in edges.keys():
            path_list, infos = edges[(edge_i, edge_j)]

            L.add_edge(edge_i, edge_j, **infos)
            L[edge_i][edge_j]['weight'] = F[edge_i][edge_j]['weight']
            edge_positions[(edge_i, edge_j)] = [positions[p] for p in path_list[1:-1]]
            L.remove_nodes_from(path_list[1:-1])

        for node in L.nodes():
            # L.nodes[node]['pos'] = positions[node]
            offset[node] = positions[node][0]

        x_positions_list = []

        for group_name in grouping.keys():

            group = grouping[group_name]

            group_size = len(group)
            group_size_compressed = round(group_size * groupcompress)
            compressed_factor = group_size_compressed if group_size_compressed == 0 else group_size_compressed - 1

            node_distance_factor = compressed_factor / (group_size - 1)

            start_x = positions[group[0]][0]

            for i, node in enumerate(group):
                offset[node] = round(start_x + i * node_distance_factor)

            compressed_length = int(offset[group[-1]] - offset[group[0]])

            x_positions_list += [i for i in range(start_x, start_x + compressed_length + 1)]

        x_positions_list += list(offset.values())
        x_positions = set(x_positions_list)
        empty_spots = []

        for x in range(global_x+1):
            if x not in x_positions:
                empty_spots.append(x)

        node_positions = {}
        for node in L.nodes():
            decrease = 0
            x = offset[node]
            for e in empty_spots:
                if x > e:
                    decrease += 1
                    if e == empty_spots[-1]:
                        offset[node] = x - decrease
                        break
                else:
                    offset[node] = x - decrease
                    break
            node_positions[node] = (offset[node], positions[node][1])

        # edge_positions = {}
        for edge_i, edge_j, data in L.edges(data=True):
            if (edge_i, edge_j) in edge_positions:
                for i, (x, y) in enumerate(edge_positions[(edge_i, edge_j)]):
                    decrease = 0
                    for e in empty_spots:
                        if x > e-1:
                            decrease += 1
                        else:
                            edge_positions[(edge_i, edge_j)][i] = (x - decrease, y)
                            break
                    else:
                        edge_positions[(edge_i, edge_j)][i] = (x - decrease, y)
            else:
                if offset[edge_j] - offset[edge_i] != 1:
                    y = positions[edge_i][1]
                    x = offset[edge_i] + 1
                    while x < offset[edge_j]:
                        edge_positions[(edge_i, edge_j)].append((x, y))
                        x += 1
            # edge_positions[(edge_i, edge_j)] = data['route']

        node_groups = {}
        for label, nodes in grouping.items():
            for node in nodes:
                node_groups[node] = label

        return(node_positions, edge_positions, node_groups)
