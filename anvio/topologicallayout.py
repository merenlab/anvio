# pylint: disable=line-too-long
"""
    Topological layout for the pangenome graph.

    Ports the Nemesis notebook layout (Restart.ipynb cell 9) into
    anvi'o, plus re-adds the gene_cluster_grouping_threshold /
    groupcompress UI features from the previous implementation.
    See HANDOFF.md task 6.
"""

from collections import defaultdict
from dataclasses import dataclass
from itertools import chain

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


GHOST_PREFIX = "GHOST_"
START_NODE = "START"
STOP_NODE = "STOP"


class LayoutInvariantError(ConfigError):
    """Raised when the layout violates the no-duplicate-cell invariant."""


@dataclass(frozen=True)
class Branch:
    """A non-spine subtree to be placed."""
    nodes: tuple
    x_start: int
    length: int
    impact: int
    seq: int


def _impact_of(nodes, gc_size):
    return -sum(gc_size.get(n, 0) for n in nodes if not n.startswith(GHOST_PREFIX))


def _build_layout_graph(G):
    """Working copy of G with virtual START/STOP brackets added.
    Every original edge gets weight=-1 so Bellman-Ford yields the
    longest path."""
    L = nx.DiGraph()
    L.add_nodes_from(G.nodes(data=True))
    for u, v, data in G.edges(data=True):
        L.add_edge(u, v, **data)
        L[u][v]['weight'] = -1

    sources = [n for n in L.nodes() if L.in_degree(n) == 0]
    sinks = [n for n in L.nodes() if L.out_degree(n) == 0]
    for s in sources:
        L.add_edge(START_NODE, s, weight=-1)
    for s in sinks:
        L.add_edge(s, STOP_NODE, weight=-1)
    return L


def _inject_bracket_edges(L, spine, gc_size, total_genomes):
    """Bracket non-spine components by their flanking core GCs (or
    START/STOP).  Virtual weight=0 edges constrain the upcoming
    topological-generations pass; they never materialise as drawn
    edges."""
    spine_pos = {n: i for i, n in enumerate(spine)}

    def is_core(n):
        return (n not in (START_NODE, STOP_NODE)
                and gc_size.get(n, 0) == total_genomes)

    preceding = [None] * len(spine)
    last_core = None
    for i, n in enumerate(spine):
        if is_core(n):
            last_core = i
        preceding[i] = last_core
    following = [None] * len(spine)
    next_core = None
    for i in range(len(spine) - 1, -1, -1):
        if is_core(spine[i]):
            next_core = i
        following[i] = next_core

    spine_set = set(spine)
    # spine includes START/STOP, but for the "is this neighbour on the spine?"
    # test we want only the real spine. Without this, dangling-start nodes (whose
    # only L-predecessor is START) would have pred_on_spine=[0], flipping
    # is_second_source to False and skipping the bracket they were meant to get.
    # Symmetric problem for dangling-end nodes via STOP.
    spine_inner = spine_set - {START_NODE, STOP_NODE}
    non_spine = {n for n in L.nodes()
                 if n not in spine_set and n not in (START_NODE, STOP_NODE)}
    sub = nx.Graph()
    sub.add_nodes_from(non_spine)
    for u, v in L.edges():
        if u in non_spine and v in non_spine:
            sub.add_edge(u, v)

    n_added = 0
    for comp in nx.connected_components(sub):
        comp_set = set(comp)
        local_sources = []
        local_sinks = []
        pred_on_spine = []
        succ_on_spine = []
        for n in comp:
            preds = list(L.predecessors(n))
            succs = list(L.successors(n))
            if not any(p in comp_set for p in preds):
                local_sources.append(n)
            if not any(s in comp_set for s in succs):
                local_sinks.append(n)
            for p in preds:
                if p in spine_inner:
                    pred_on_spine.append(spine_pos[p])
            for s in succs:
                if s in spine_inner:
                    succ_on_spine.append(spine_pos[s])

        is_second_source = not pred_on_spine
        is_dangling_sink = not succ_on_spine
        if not (is_second_source or is_dangling_sink):
            continue

        if is_second_source:
            if succ_on_spine:
                merge_pos = min(succ_on_spine)
                if preceding[merge_pos] == merge_pos:
                    idx = preceding[merge_pos - 1] if merge_pos > 0 else None
                else:
                    idx = preceding[merge_pos]
                left_bracket = spine[idx] if idx is not None else START_NODE
            else:
                left_bracket = START_NODE
            for src in local_sources:
                if left_bracket != src and not L.has_edge(left_bracket, src):
                    L.add_edge(left_bracket, src, weight=0)
                    n_added += 1

        if is_dangling_sink:
            if pred_on_spine:
                div_pos = max(pred_on_spine)
                if following[div_pos] == div_pos:
                    idx = following[div_pos + 1] if div_pos + 1 < len(spine) else None
                else:
                    idx = following[div_pos]
                right_bracket = spine[idx] if idx is not None else STOP_NODE
            else:
                right_bracket = STOP_NODE
            for sink in local_sinks:
                if sink != right_bracket and not L.has_edge(sink, right_bracket):
                    L.add_edge(sink, right_bracket, weight=0)
                    n_added += 1

    return n_added


def _initial_x_by_generations(L):
    positions = {}
    global_x = 0
    for x, generation in enumerate(nx.topological_generations(L)):
        global_x = x
        for n in generation:
            positions[n] = x
    return positions, global_x


def _contract_x_with_ghosts(L, positions, global_x):
    """Walk x right-to-left; push every node as far right as possible
    and fill remaining gaps with GHOST_* placeholders so every edge in
    L spans exactly one column."""
    edges_saved = {}
    by_x = defaultdict(set)
    for n, x in positions.items():
        by_x[x].add(n)

    succ_cache = {n: list(L.successors(n)) for n in L.nodes()}

    ghost_id = 0
    for x in range(global_x - 1, 0, -1):
        for node in list(by_x[x]):
            node_x = positions[node]
            # Brackets (weight=0) need asymmetric treatment in this pass:
            #   * They must *cap* the rightward pull -- otherwise the anchor
            #     core of a dangling-start would overshoot past its dangling
            #     sibling and collide at the same x.
            #   * They must NOT themselves *induce* a pull -- otherwise a
            #     dangling-end whose only successor is its closing bracket
            #     would be dragged to next-core.x - 1 instead of staying one
            #     column after its real predecessor.
            # And ghost chains only ever materialise from real edges.
            change_real = []
            change_all = []
            for nb in succ_cache[node]:
                if nb == STOP_NODE:
                    continue
                nb_x = positions[nb]
                if nb_x <= node_x:
                    raise ConfigError(
                        f"Topological violation: {node} (x={node_x}) -> "
                        f"{nb} (x={nb_x}). This shouldn't happen on a DAG.")
                change_all.append(nb_x)
                if L.edges[node, nb].get("weight") != 0:
                    change_real.append((nb_x, nb))
            if not change_real:
                continue

            min_next_x = min(change_all)
            if min_next_x > 1:
                new_x = min_next_x - 1
                if new_x != node_x:
                    positions[node] = new_x
                    by_x[node_x].discard(node)
                    by_x[new_x].add(node)
                    node_x = new_x

            for nb_x, succ in change_real:
                gap = nb_x - node_x
                if gap == 1:
                    continue
                path = [node]
                for i in range(1, gap):
                    g = f"{GHOST_PREFIX}{ghost_id}"
                    ghost_id += 1
                    positions[g] = node_x + i
                    by_x[node_x + i].add(g)
                    path.append(g)
                path.append(succ)

                edge_data = dict(L.edges[node, succ])
                edges_saved[(node, succ)] = (path, edge_data)
                L.remove_edge(node, succ)
                w = -0.5 if len(path) > 2 else -1
                L.add_edges_from(
                    ((u, v) for u, v in zip(path, path[1:])),
                    weight=w,
                )

                cache = succ_cache[node]
                cache[cache.index(succ)] = path[1]
                for i in range(1, len(path) - 1):
                    succ_cache[path[i]] = [path[i + 1]]

    return positions, edges_saved


def _find_chain_groups(L):
    """DFS from START; collapse consecutive degree-1-1 nodes into
    chain groups."""
    groups = {}
    groups_rev = {}
    group_id = 0
    in_deg = L.in_degree
    out_deg = L.out_degree
    for u, v in nx.dfs_edges(L, source=START_NODE):
        if u == START_NODE or v == STOP_NODE:
            continue
        if in_deg(u) != 1 or out_deg(u) != 1:
            continue
        if in_deg(v) != 1 or out_deg(v) != 1:
            continue
        if u not in groups_rev:
            groups[group_id] = [u, v]
            groups_rev[u] = group_id
            groups_rev[v] = group_id
            group_id += 1
        else:
            gid = groups_rev[u]
            groups[gid].append(v)
            groups_rev[v] = gid
    return groups, groups_rev


def _make_branch_list(L, groups, groups_rev, spine_set, positions, gc_size):
    """Every non-spine chain plus every non-spine singleton, sorted
    by (length, x_start, impact)."""
    out = []
    seq = 0
    for chain_nodes in groups.values():
        if set(chain_nodes).issubset(spine_set):
            continue
        out.append(Branch(
            nodes=tuple(chain_nodes),
            x_start=positions[chain_nodes[0]],
            length=len(chain_nodes),
            impact=_impact_of(chain_nodes, gc_size),
            seq=seq,
        ))
        seq += 1

    in_groups = set(groups_rev)
    for n in L.nodes():
        if n in in_groups or n in spine_set:
            continue
        if n in (START_NODE, STOP_NODE):
            continue
        out.append(Branch(
            nodes=(n,),
            x_start=positions[n],
            length=1,
            impact=_impact_of((n,), gc_size),
            seq=seq,
        ))
        seq += 1
    out.sort(key=lambda b: (b.length, b.x_start, b.impact))
    return out


def _assign_y(L, spine, branches, positions_x):
    """Greedy y-assignment.  Spine at y=0; each branch placed at the
    lowest free row >= max(neighbour_y)."""
    y = {}
    used_ys = defaultdict(set)
    for n in spine:
        y[n] = 0
        used_ys[positions_x[n]].add(0)

    pending = list(branches)
    # Drop weight=0 brackets from the per-branch adjacency. Brackets only
    # exist to constrain topological generations; if we let them through
    # here, a dangling-end whose bracket successor lands on the spine would
    # be eligible immediately with adj_ys=[spine_y], landing one row above
    # the spine instead of staying with its real predecessor chain.
    def _real_preds(n):
        return {p for p in L.predecessors(n) if L.edges[p, n].get("weight") != 0}
    def _real_succs(n):
        return {s for s in L.successors(n) if L.edges[n, s].get("weight") != 0}
    branch_pred = [_real_preds(b.nodes[0]) for b in branches]
    branch_succ = [_real_succs(b.nodes[-1]) for b in branches]

    stack = [spine]
    global_y = 0

    while stack:
        current_nodes = stack[-1]
        current_set = set(current_nodes)
        placed_one = False

        for i, b in enumerate(pending):
            if b is None:
                continue
            bp = branch_pred[i]
            bs = branch_succ[i]
            if bp.isdisjoint(current_set) and bs.isdisjoint(current_set):
                continue

            adj_ys = [y[n] for n in chain(bp, bs) if n in y]
            y_new = max(adj_ys) if adj_ys else 0

            x_lo = b.x_start
            x_hi = x_lo + b.length - 1
            while any(y_new in used_ys[x] for x in range(x_lo, x_hi + 1)):
                y_new += 1

            for n in b.nodes:
                y[n] = y_new
                used_ys[positions_x[n]].add(y_new)
            if y_new > global_y:
                global_y = y_new

            pending[i] = None
            stack.append(list(b.nodes))
            placed_one = True
            break

        if not placed_one:
            stack.pop()

    return y, global_y


def _squeeze_empty_columns(positions_x, edge_paths, global_x):
    """Pre-compute a cumulative decrement[x] so each cell shifts left
    in O(1)."""
    used_x = set(positions_x.values())
    for path in edge_paths.values():
        for x, _ in path:
            used_x.add(x)
    decrement = [0] * (global_x + 2)
    running = 0
    for x in range(global_x + 2):
        if x not in used_x:
            running += 1
        decrement[x] = running
    new_x = {n: positions_x[n] - decrement[positions_x[n]] for n in positions_x}
    new_paths = {e: [(x - decrement[x], y) for (x, y) in path]
                 for e, path in edge_paths.items()}
    return new_x, new_paths


def _compress_y_rows(positions_y, edge_paths):
    """Dense remap of used y values (placement may leave intermediate
    rows empty)."""
    used_y = set(positions_y.values())
    for path in edge_paths.values():
        for _, y in path:
            used_y.add(y)
    if not used_y:
        return positions_y, edge_paths
    max_y = max(used_y)
    shift = [0] * (max_y + 1)
    running = 0
    for y in range(max_y + 1):
        if y not in used_y:
            running += 1
        shift[y] = running
    new_pos = {n: positions_y[n] - shift[positions_y[n]] for n in positions_y}
    new_paths = {e: [(x, y - shift[y]) for (x, y) in path]
                 for e, path in edge_paths.items()}
    return new_pos, new_paths


def _assert_no_cell_duplicates(node_positions, edge_paths):
    node_cells = {}
    for n, xy in node_positions.items():
        if xy in node_cells:
            raise LayoutInvariantError(
                f"Two nodes at the same cell {xy}: "
                f"{node_cells[xy]!r} and {n!r}")
        node_cells[xy] = n

    edge_cell_owner = {}
    for edge, path in edge_paths.items():
        for cell in path:
            if cell in node_cells:
                raise LayoutInvariantError(
                    f"Edge {edge!r} passes through node "
                    f"{node_cells[cell]!r} at cell {cell}")
            if cell in edge_cell_owner:
                raise LayoutInvariantError(
                    f"Two edges share an interior cell {cell}: "
                    f"{edge_cell_owner[cell]!r} and {edge!r}")
            edge_cell_owner[cell] = edge


def _build_node_groups(L, groups, gene_cluster_grouping_threshold):
    """DFS-chain → sub-split by (genome-set, type) → threshold filter
    → label.  Ghost nodes are filtered out before sub-splitting.
    Returns (grouping, node_groups)."""
    if gene_cluster_grouping_threshold == -1:
        return {}, {}

    grouping = {}
    next_id = max(groups.keys(), default=-1) + 1

    for group_name, group_nodes in groups.items():
        real_nodes = [n for n in group_nodes if not n.startswith(GHOST_PREFIX)]
        if not real_nodes:
            continue

        former_genome = None
        former_type = None
        label = f"GCG_{group_name:08d}"
        condense = []

        for node in real_nodes:
            current_genome = set(L.nodes[node].get('gene_calls', {}).keys())
            current_type = L.nodes[node].get('type', '')

            if former_genome is None and former_type is None:
                condense.append(node)
            elif former_genome == current_genome and former_type == current_type:
                condense.append(node)
            else:
                if len(condense) >= gene_cluster_grouping_threshold:
                    grouping[label] = condense
                condense = [node]
                label = f"GCG_{next_id:08d}"
                next_id += 1
            former_genome = current_genome
            former_type = current_type

        if len(condense) >= gene_cluster_grouping_threshold:
            grouping[label] = condense

    node_groups = {}
    for label, nodes in grouping.items():
        for node in nodes:
            node_groups[node] = label
    return grouping, node_groups


def _apply_groupcompress(positions_x, grouping, groupcompress):
    """Rewrite x-coords inside each grouped chain so it spans
    round(size * groupcompress) columns.  Multiple group nodes may
    end up sharing a column — that is intentional (the UI renders
    grouped runs as visually collapsed); the squeeze pass then
    collapses the now-empty trailing columns."""
    if not grouping or groupcompress == 1.0:
        return positions_x

    for nodes in grouping.values():
        size = len(nodes)
        if size < 2:
            continue
        size_compressed = round(size * groupcompress)
        compressed_factor = (size_compressed
                             if size_compressed == 0
                             else size_compressed - 1)
        step = compressed_factor / (size - 1)
        start_x = positions_x[nodes[0]]
        for i, node in enumerate(nodes):
            positions_x[node] = round(start_x + i * step)
    return positions_x


class TopologicalLayout():
    """Topological layout for the pangenome graph (Nemesis port).

    The layout algorithm itself lives in module-level helpers; this
    class wraps them with anvi'o logging, component selection, and
    the UI grouping / groupcompress pipeline.
    """

    def __init__(self, r=run, p=progress):
        self.run = r
        self.progress = p

    def run_synteny_layout_algorithm(self, F,
                                     gene_cluster_grouping_threshold=-1,
                                     groupcompress=1.0,
                                     component=0):
        """Lay out the ``component``-th largest weakly-connected
        component of ``F``.  Returns ``(node_positions, edge_positions,
        node_groups)``."""

        components = sorted(nx.weakly_connected_components(F),
                            key=len, reverse=True)
        if not components:
            raise ConfigError("The pangenome graph has no nodes — nothing to lay out.")
        if component < 0 or component >= len(components):
            raise ConfigError(
                f"Component index {component} is out of range — the graph "
                f"has {len(components)} component(s) (valid indices: "
                f"0..{len(components) - 1}).")

        sub = F.subgraph(components[component]).copy()

        if not nx.is_directed_acyclic_graph(sub):
            raise ConfigError(
                "The pangenome graph contains a cycle. The topological "
                "layout requires a DAG; the Nemesis engine should never "
                "produce a cycle — please report this.")

        gc_size = {n: len(d.get('gene_calls', {}))
                   for n, d in sub.nodes(data=True)}
        original_edges = set(sub.edges())

        self.run.info("Layout component", component)
        self.run.info("Component nodes", sub.number_of_nodes())
        self.run.info("Component edges", sub.number_of_edges())

        L = _build_layout_graph(sub)

        spine = nx.bellman_ford_path(L, START_NODE, STOP_NODE, weight="weight")
        spine_set = set(spine) - {START_NODE, STOP_NODE}
        self.run.info("Spine length", len(spine_set))

        total_genomes = max(gc_size.values(), default=0)
        n_brackets = _inject_bracket_edges(L, spine, gc_size, total_genomes)
        self.run.info("Bracket edges injected", n_brackets)

        positions_x, global_x = _initial_x_by_generations(L)
        positions_x, edges_saved = _contract_x_with_ghosts(L, positions_x, global_x)
        n_ghost = sum(1 for n in positions_x if n.startswith(GHOST_PREFIX))
        self.run.info("Ghost columns", n_ghost)

        groups, groups_rev = _find_chain_groups(L)
        L.remove_nodes_from([START_NODE, STOP_NODE])

        branches = _make_branch_list(L, groups, groups_rev, spine_set,
                                     positions_x, gc_size)
        y_map, _global_y = _assign_y(L, spine, branches, positions_x)

        edge_paths = {}
        for (i, j), (path, _attrs) in edges_saved.items():
            edge_paths[(i, j)] = [
                (positions_x[g], y_map[g]) for g in path[1:-1]
            ]
            for g in path[1:-1]:
                positions_x.pop(g, None)
                y_map.pop(g, None)

        for (s, t) in original_edges:
            if (s, t) not in edge_paths:
                edge_paths[(s, t)] = []

        positions_x.pop(START_NODE, None)
        positions_x.pop(STOP_NODE, None)
        y_map.pop(START_NODE, None)
        y_map.pop(STOP_NODE, None)

        # Verify the Nemesis layout's no-duplicate invariant before
        # the UI grouping pass mutates x's (groupcompress intentionally
        # collapses cells, which would trip this check).
        node_positions_pre = {n: (positions_x[n], y_map[n]) for n in positions_x}
        _assert_no_cell_duplicates(node_positions_pre, edge_paths)

        if gene_cluster_grouping_threshold == -1:
            self.run.info_single("Setting algorithm to 'no grouping'")
        else:
            self.run.info_single(
                f"Setting algorithm to 'Grouping single connected chains "
                f"size > {gene_cluster_grouping_threshold}'")

        grouping, node_groups = _build_node_groups(
            L, groups, gene_cluster_grouping_threshold)
        self.run.info_single(
            f"Grouped {sum(len(v) for v in grouping.values())} nodes "
            f"in {len(grouping)} groups")

        positions_x = _apply_groupcompress(positions_x, grouping, groupcompress)
        positions_x, edge_paths = _squeeze_empty_columns(
            positions_x, edge_paths, global_x)
        y_map, edge_paths = _compress_y_rows(y_map, edge_paths)

        node_positions = {n: (positions_x[n], y_map[n]) for n in positions_x}

        return node_positions, edge_paths, node_groups
