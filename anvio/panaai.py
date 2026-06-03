# pylint: disable=line-too-long
"""PangenomeAAI engine.

Builds a pangenome graph from the artifacts produced *upstream* of a pan-db:
the genomes-storage db (genome hash -> name), each genome's CONTIGS.db (gene
calls), and DIAMOND's raw cross-genome bitscore output.  Replaces the legacy
SyntenyGeneCluster + DirectedForce pipeline.

The engine fuses gene endpoints into super-nodes by AAI minbit, picking
fusions with a frontier-growth (Prim-style) walk over a stochastic top-K
bucket of the ranked edges.  Super-nodes are guaranteed to contain at most
one gene per genome.  The output graph is a DAG by construction; no
edge-reversal step is needed.

See HANDOFF.md (Task 4) for the migration story behind this module.
"""

import csv
import os
import random
from argparse import Namespace
from collections import Counter, defaultdict, deque

import networkx as nx

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.genomestorage as genomestorage

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
P = terminal.pluralize


# ---------------------------------------------------------------------------
# Pure helpers — no I/O, no logging.
# ---------------------------------------------------------------------------


def parse_external_genomes(path):
    """Read external-genomes.txt → dict ``{genome_name: contigs_db_path}``."""
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["name"]] = row["contigs_db_path"]
    return out


def split_anvio_id(seq_id):
    """``hash<8>_<gene_callers_id>`` → (hash, gene_callers_id_int).

    Returns ``(None, None)`` if the format doesn't match.
    """
    if "_" not in seq_id:
        return None, None
    hash_part, gid_part = seq_id.rsplit("_", 1)
    if not hash_part.startswith("hash") or not gid_part.isdigit():
        return None, None
    return hash_part, int(gid_part)


def build_gene_to_contig(genome_calls):
    """Build ``{(genome, gene_callers_id) -> contig}`` lookup."""
    out = {}
    for genome, calls in genome_calls.items():
        for gc in calls:
            out[(genome, gc["gene_callers_id"])] = gc["contig"]
    return out


def _cmp_window_topn(W_a, W_b, sing_minbit):
    """Top-n best-match window score.

    For each gene in ``W_a`` record its maximum minbit against any gene in
    ``W_b`` (zero if no candidate matches); sort descending; keep top
    ``n = min(|W_a|, |W_b|)``; return their mean.  Boundary truncation falls
    out for free near contig ends.  Returns ``None`` if either window is
    empty.
    """
    if not W_a or not W_b:
        return None
    scores = []
    for a in W_a:
        best = 0.0
        for b in W_b:
            key = (a, b) if a < b else (b, a)
            w = sing_minbit.get(key, 0.0)
            if w > best:
                best = w
        scores.append(best)
    n = min(len(W_a), len(W_b))
    scores.sort(reverse=True)
    return sum(scores[:n]) / n


def _aggregate_line_pairs(lines, line_names, edges, K):
    """Walk every AAI edge, accumulate (forward, reverse) sums per
    ordered line-pair ``(li_a, li_b)`` with ``li_a < li_b``, and record
    the per-edge ``(forward, reverse)`` signals.

    Returns ``(fwd_sum, fwd_cnt, rev_sum, rev_cnt, total, edge_signals)``.
    """
    line_idx_by_name = {nm: i for i, nm in enumerate(line_names)}
    line_lengths = [len(line) for line in lines]
    pos_offsets = [0]
    for ll in line_lengths:
        pos_offsets.append(pos_offsets[-1] + ll)

    token_to_id = {}
    for li, line in enumerate(lines):
        for pos, tok in enumerate(line):
            token_to_id.setdefault((li, tok), pos_offsets[li] + pos)

    sing_minbit = {}
    for (u_str, v_str, w) in edges:
        try:
            u_line, u_tok = u_str.rsplit(":", 1)
            v_line, v_tok = v_str.rsplit(":", 1)
        except ValueError:
            continue
        li_u = line_idx_by_name.get(u_line)
        li_v = line_idx_by_name.get(v_line)
        if li_u is None or li_v is None:
            continue
        u_id = token_to_id.get((li_u, u_tok))
        v_id = token_to_id.get((li_v, v_tok))
        if u_id is None or v_id is None:
            continue
        key = (u_id, v_id) if u_id < v_id else (v_id, u_id)
        prev = sing_minbit.get(key)
        if prev is None or w > prev:
            sing_minbit[key] = float(w)

    fwd_sum = defaultdict(float)
    fwd_cnt = defaultdict(int)
    rev_sum = defaultdict(float)
    rev_cnt = defaultdict(int)
    total = defaultdict(int)
    edge_signals = {}

    for (u_str, v_str, w) in edges:
        try:
            u_line, u_tok = u_str.rsplit(":", 1)
            v_line, v_tok = v_str.rsplit(":", 1)
        except ValueError:
            continue
        li_u = line_idx_by_name.get(u_line)
        li_v = line_idx_by_name.get(v_line)
        if li_u is None or li_v is None or li_u == li_v:
            continue
        u_id = token_to_id.get((li_u, u_tok))
        v_id = token_to_id.get((li_v, v_tok))
        if u_id is None or v_id is None:
            continue

        pos_u = u_id - pos_offsets[li_u]
        pos_v = v_id - pos_offsets[li_v]
        len_u = line_lengths[li_u]
        len_v = line_lengths[li_v]
        off_u = pos_offsets[li_u]
        off_v = pos_offsets[li_v]

        L_u = [off_u + (pos_u - d) for d in range(1, K + 1) if pos_u - d >= 0]
        R_u = [off_u + (pos_u + d) for d in range(1, K + 1) if pos_u + d < len_u]
        L_v = [off_v + (pos_v - d) for d in range(1, K + 1) if pos_v - d >= 0]
        R_v = [off_v + (pos_v + d) for d in range(1, K + 1) if pos_v + d < len_v]

        LL = _cmp_window_topn(L_u, L_v, sing_minbit)
        RR = _cmp_window_topn(R_u, R_v, sing_minbit)
        LR = _cmp_window_topn(L_u, R_v, sing_minbit)
        RL = _cmp_window_topn(R_u, L_v, sing_minbit)

        fwd_terms = [x for x in (LL, RR) if x is not None]
        rev_terms = [x for x in (LR, RL) if x is not None]
        forward = sum(fwd_terms) / len(fwd_terms) if fwd_terms else None
        reverse = sum(rev_terms) / len(rev_terms) if rev_terms else None

        edge_signals[(u_str, v_str)] = (forward, reverse)

        key = (li_u, li_v) if li_u < li_v else (li_v, li_u)
        total[key] += 1
        if forward is not None:
            fwd_sum[key] += forward
            fwd_cnt[key] += 1
        if reverse is not None:
            rev_sum[key] += reverse
            rev_cnt[key] += 1

    return fwd_sum, fwd_cnt, rev_sum, rev_cnt, total, edge_signals


def _tree_path_edges(a, b, parent):
    """Canonical tree edges along the spanning-tree path from ``a`` to ``b``."""
    ancestors_a = {a}
    cur = a
    while parent[cur] is not None:
        cur = parent[cur][0]
        ancestors_a.add(cur)
    path_b = []
    cur = b
    while cur not in ancestors_a:
        p = parent[cur][0]
        path_b.append((cur, p) if cur < p else (p, cur))
        cur = p
    lca = cur
    path_a = []
    cur = a
    while cur != lca:
        p = parent[cur][0]
        path_a.append((cur, p) if cur < p else (p, cur))
        cur = p
    return path_a + path_b


def compute_line_orientations(lines, line_names, edges,
                              K=10, min_hits=10,
                              tie_threshold=0.02,
                              min_score=0.0,
                              demotion_strategy="slimmest-margin",
                              precomputed=None):
    """Compute a per-line flip vector from the singleton-locality scan.

    Pipeline: (1) build per-pair labels gated by ``min_hits``,
    ``tie_threshold`` AND ``min_score``; (2) BFS spanning tree on the
    labeled graph; (3) detect odd fundamental cycles via non-tree edges
    and demote the weakest edge per ``demotion_strategy`` until balanced;
    (4) final BFS on the cleaned graph for flips + components.

    Returns ``(flips, components, inconsistencies, pair_label)``.
    """
    n_lines = len(lines)
    if precomputed is None:
        agg = _aggregate_line_pairs(lines, line_names, edges, K)
    else:
        agg = precomputed
    fwd_sum, fwd_cnt, rev_sum, rev_cnt, total = agg[:5]

    pair_label = {}
    pair_evidence = {}
    for key, n_total in total.items():
        if n_total < min_hits:
            continue
        avg_fwd = fwd_sum[key] / n_total
        avg_rev = rev_sum[key] / n_total
        diff = avg_fwd - avg_rev
        if abs(diff) < tie_threshold:
            continue
        if max(avg_fwd, avg_rev) < min_score:
            continue
        pair_label[key] = "same" if diff > 0 else "flip"
        pair_evidence[key] = (n_total, abs(diff))

    adj = defaultdict(list)
    for (a, b), label in pair_label.items():
        sign = +1 if label == "same" else -1
        adj[a].append((b, sign))
        adj[b].append((a, sign))

    color = [None] * n_lines
    parent = [None] * n_lines
    for start in range(n_lines):
        if color[start] is not None:
            continue
        color[start] = +1
        q = deque([start])
        while q:
            cur = q.popleft()
            for nb, sign in adj[cur]:
                if color[nb] is None:
                    color[nb] = color[cur] * sign
                    parent[nb] = (cur, sign)
                    q.append(nb)

    tree_edges = set()
    for node in range(n_lines):
        if parent[node] is None:
            continue
        p = parent[node][0]
        tree_edges.add((node, p) if node < p else (p, node))

    odd_cycles = []
    inconsistencies = []
    for (a, b), label in pair_label.items():
        if (a, b) in tree_edges:
            continue
        sign = +1 if label == "same" else -1
        if color[a] * sign == color[b]:
            continue
        cycle = _tree_path_edges(a, b, parent) + [(a, b)]
        odd_cycles.append(cycle)
        inconsistencies.append(("contradiction", a, b, cycle))

    if demotion_strategy == "slimmest-margin":
        weakness_key = lambda e: (pair_evidence[e][1], pair_evidence[e][0])
    elif demotion_strategy == "fewest-edges":
        weakness_key = lambda e: (pair_evidence[e][0], pair_evidence[e][1])
    else:
        raise ConfigError(f"Unknown demotion strategy {demotion_strategy!r}; "
                          f"expected 'slimmest-margin' or 'fewest-edges' :/")

    remaining = list(odd_cycles)
    while remaining:
        candidates = {edge for cycle in remaining for edge in cycle
                      if edge in pair_label}
        if not candidates:
            break
        weakest = min(candidates, key=weakness_key)
        del pair_label[weakest]
        remaining = [c for c in remaining if weakest not in c]

    adj_clean = defaultdict(list)
    for (a, b), label in pair_label.items():
        sign = +1 if label == "same" else -1
        adj_clean[a].append((b, sign))
        adj_clean[b].append((a, sign))

    color_final = [None] * n_lines
    components = []
    for start in range(n_lines):
        if color_final[start] is not None:
            continue
        color_final[start] = +1
        comp = {start}
        q = deque([start])
        while q:
            cur = q.popleft()
            for nb, sign in adj_clean[cur]:
                if color_final[nb] is None:
                    color_final[nb] = color_final[cur] * sign
                    comp.add(nb)
                    q.append(nb)
        components.append(frozenset(comp))

    flips = [c == -1 for c in color_final]
    return flips, components, inconsistencies, pair_label


def apply_line_flips(lines, flips):
    """Return a new lines list with every line marked ``True`` in
    ``flips`` reversed (token order)."""
    return [list(reversed(line)) if flip else list(line)
            for line, flip in zip(lines, flips)]


# ---------------------------------------------------------------------------
# Pangenome graph construction primitives.
# ---------------------------------------------------------------------------


def _uf_find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def _add_line_to_graph(G, line_idx, lines, line_names, flip,
                       super_of_gene, uf_parent,
                       in_g_lines, in_g_flip):
    """Add a line's gene-order chain to G in the chosen orientation."""
    tokens = lines[line_idx]
    name = line_names[line_idx]
    seq = list(reversed(tokens)) if flip else list(tokens)
    nodes = [f"{name}:{t}" for t in seq]
    for node in nodes:
        G.add_node(node, lines={line_idx}, genes={node})
        super_of_gene[node] = node
        uf_parent[node] = node
    for a, b in zip(nodes, nodes[1:]):
        G.add_edge(a, b, kind="gene_order")
    in_g_lines.add(line_idx)
    in_g_flip[line_idx] = flip


def _fuse_genes(G, u, v, w, score, mean, label,
                super_of_gene, uf_parent, rejected_edges,
                line_to_genome):
    """Try to contract the super-nodes containing ``u`` and ``v``.

    Records ``transitive-cycle`` in ``rejected_edges`` if the contraction
    would break the DAG, or ``same-genome-conflict`` if it would pull two
    genes from the same genome into one super-node.  No-op when the two
    genes are already in the same super-node.
    """
    A = _uf_find(uf_parent, super_of_gene[u])
    B = _uf_find(uf_parent, super_of_gene[v])
    if A == B:
        return
    if nx.has_path(G, A, B) or nx.has_path(G, B, A):
        rejected_edges.append((u, v, mean, "transitive-cycle"))
        return
    a_attrs = G.nodes[A]
    b_attrs = G.nodes[B]
    a_genomes = {line_to_genome[g.rsplit(":", 1)[0]] for g in a_attrs.get("genes", ())}
    b_genomes = {line_to_genome[g.rsplit(":", 1)[0]] for g in b_attrs.get("genes", ())}
    if a_genomes & b_genomes:
        rejected_edges.append((u, v, mean, "same-genome-conflict"))
        return
    merged_lines = a_attrs.get("lines", set()) | b_attrs.get("lines", set())
    merged_genes = a_attrs.get("genes", set()) | b_attrs.get("genes", set())
    fuse_log = (a_attrs.get("fuses", []) + b_attrs.get("fuses", [])
                + [{"u": u, "v": v, "minbit": w,
                    "decision_score": score, "mean": mean,
                    "label": label if label else "tie"}])
    nx.contracted_nodes(G, A, B, self_loops=False, copy=False)
    G.nodes[A]["lines"] = merged_lines
    G.nodes[A]["genes"] = merged_genes
    G.nodes[A]["fuses"] = fuse_log
    uf_parent[B] = A


def _try_commit_tie_line(li, lines, line_names, anchors, pending_merges, G,
                         super_of_gene, uf_parent,
                         in_g_lines, in_g_flip, rejected_edges,
                         line_to_genome):
    """Commit a tie-only line if it has >=2 anchors with distinct
    line-positions AND distinct graph topological positions."""
    li_anchors = anchors[li]
    if len(li_anchors) < 2:
        return False
    topo = {n: i for i, n in enumerate(nx.topological_sort(G))}
    for i in range(len(li_anchors)):
        for j in range(i + 1, len(li_anchors)):
            a1, a2 = li_anchors[i], li_anchors[j]
            if a1["pos"] == a2["pos"]:
                continue
            sup1 = _uf_find(uf_parent, super_of_gene[a1["in_gene"]])
            sup2 = _uf_find(uf_parent, super_of_gene[a2["in_gene"]])
            t1 = topo.get(sup1)
            t2 = topo.get(sup2)
            if t1 is None or t2 is None or t1 == t2:
                continue
            line_order = 1 if a1["pos"] > a2["pos"] else -1
            graph_order = 1 if t1 > t2 else -1
            flip = (line_order != graph_order)
            _add_line_to_graph(G, li, lines, line_names, flip,
                               super_of_gene, uf_parent,
                               in_g_lines, in_g_flip)
            for (u, v, w, score, mean, label) in pending_merges[li]:
                _fuse_genes(G, u, v, w, score, mean, label,
                            super_of_gene, uf_parent, rejected_edges,
                            line_to_genome)
            return True
    return False


def rename_pangenome_nodes(G, gene_clusters=None, line_to_genome=None):
    """Rename super-nodes deterministically in topological order.

    Panmode (``gene_clusters`` given): each super-node is named after its
    parent GC; multiple super-nodes from the same GC get ``_1``, ``_2``,
    ... suffixes in topological order.

    Non-panmode: synthesize sequential names ``GC_00000001_1``,
    ``GC_00000002_1``, ... in topological order.

    Node attributes (``lines``, ``genes``, ``fuses``) are preserved.
    Returns a new graph; the original is not mutated.
    """
    topo = list(nx.topological_sort(G))

    if gene_clusters is not None:
        seen = Counter()
        mapping = {}
        for n in topo:
            parent_gc = None
            for g in G.nodes[n].get("genes", ()):
                line_name, sep, gid_str = g.rpartition(":")
                if not sep or not gid_str.isdigit():
                    continue
                genome = line_to_genome.get(line_name) if line_to_genome else None
                if genome is None:
                    continue
                parent_gc = gene_clusters.get((genome, int(gid_str)))
                if parent_gc is not None:
                    break
            if parent_gc is None:
                parent_gc = "GC_UNKNOWN"
            seen[parent_gc] += 1
            mapping[n] = f"{parent_gc}_{seen[parent_gc]}"
    else:
        mapping = {n: f"GC_{i + 1:08d}_1" for i, n in enumerate(topo)}

    return nx.relabel_nodes(G, mapping, copy=True)


def compute_node_types(G, line_to_genome, gene_clusters=None):
    """Assign each super-node a ``type`` attribute.

    Per parent GC, compute ``splits`` (>=2 super-nodes came from this GC)
    and ``parent_multi_copy`` (the input gene-cluster table had >=2 genes
    from at least one genome in this GC).  Then per super-node:

      * splits AND parent_multi_copy        -> ``duplication``
      * splits AND NOT parent_multi_copy    -> ``rearrangement``
      * NOT splits, all genomes covered     -> ``core``
      * NOT splits, >1 genome covered       -> ``accessory``
      * NOT splits, exactly 1 genome        -> ``singleton``

    ``rna`` is a valid value in the schema but never assigned here.  In
    non-panmode (``gene_clusters=None``) every super-node is its own
    parent GC; only ``core``/``accessory``/``singleton`` are produced.
    """
    all_genomes = set(line_to_genome.values())

    nodes_by_parent_gc = defaultdict(list)
    for n in G.nodes():
        parent_gc = n.rsplit("_", 1)[0] if gene_clusters is not None else n
        nodes_by_parent_gc[parent_gc].append(n)

    parent_multi_copy = {}
    if gene_clusters is not None:
        gc_to_genome_counts = defaultdict(Counter)
        for (genome, _gid), gc_name in gene_clusters.items():
            gc_to_genome_counts[gc_name][genome] += 1
        for gc_name, gc_counts in gc_to_genome_counts.items():
            parent_multi_copy[gc_name] = any(v >= 2 for v in gc_counts.values())

    for parent_gc, nodes in nodes_by_parent_gc.items():
        splits = len(nodes) >= 2
        pmc = parent_multi_copy.get(parent_gc, False)

        if splits and pmc:
            for n in nodes:
                G.nodes[n]["type"] = "duplication"
        elif splits and not pmc:
            for n in nodes:
                G.nodes[n]["type"] = "rearrangement"
        else:
            for n in nodes:
                genomes = {line_to_genome[g.rsplit(":", 1)[0]]
                           for g in G.nodes[n].get("genes", ())}
                if genomes == all_genomes:
                    G.nodes[n]["type"] = "core"
                elif len(genomes) > 1:
                    G.nodes[n]["type"] = "accessory"
                else:
                    G.nodes[n]["type"] = "singleton"
    return G


def compute_component_ids(G):
    """Assign ``component_id`` per node.

    Components are weakly connected components sorted size-descending so
    ``component_id=0`` is the largest.  Mutates G in place and returns it.
    """
    components = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    for idx, comp in enumerate(components):
        for n in comp:
            G.nodes[n]["component_id"] = idx
    return G


# ---------------------------------------------------------------------------
# Debug table writers (gated on PangenomeAAIEngine.tables_dir).
# ---------------------------------------------------------------------------


def write_lines_file(genome_calls, savepath):
    """One row per (genome, contig): ``genome_contig: gid, gid, ...``."""
    with open(savepath, "w") as f:
        for genome in sorted(genome_calls):
            by_contig = defaultdict(list)
            for gc in genome_calls[genome]:
                by_contig[gc["contig"]].append(gc)
            for contig in sorted(by_contig):
                items = sorted(by_contig[contig], key=lambda x: x["start"])
                tokens = [str(gc["gene_callers_id"]) for gc in items]
                f.write(f"{genome}_{contig}: {', '.join(tokens)}\n")


def write_genome_map(genome_calls, savepath):
    """TSV ``line_name\\tgenome_name\\tcontig``."""
    with open(savepath, "w") as f:
        f.write("line_name\tgenome_name\tcontig\n")
        for genome in sorted(genome_calls):
            contigs = sorted({gc["contig"] for gc in genome_calls[genome]})
            for contig in contigs:
                line_name = f"{genome}_{contig}"
                f.write(f"{line_name}\t{genome}\t{contig}\n")


def write_edges_file(edges_dict, savepath):
    """TSV ``u\\tv\\tminbit`` sorted lexicographically by (u, v)."""
    rows = []
    for pair, w in edges_dict.items():
        u, v = sorted(pair)
        rows.append((u, v, w))
    rows.sort()
    with open(savepath, "w") as f:
        f.write("u\tv\tminbit\n")
        for u, v, w in rows:
            f.write(f"{u}\t{v}\t{w:.6f}\n")


def write_orientation_tsv(rows, savepath):
    """TSV of per-line-pair orientation diagnostics.  ``rows`` is a list of
    11-tuples produced by :py:meth:`PangenomeAAIEngine._compute_orientations`."""
    with open(savepath, "w") as f:
        f.write("line_a\tline_b\tgenome_a\tgenome_b\tn_edges\tn_fwd"
                "\tavg_forward\tn_rev\tavg_reverse\tdiff\torientation\n")
        for (line_a, line_b, genome_a, genome_b, n_edges, n_fwd,
             avg_fwd, n_rev, avg_rev, diff, orient) in rows:
            f.write(f"{line_a}\t{line_b}\t{genome_a}\t{genome_b}\t{n_edges}"
                    f"\t{n_fwd}\t{avg_fwd:.4f}\t{n_rev}\t{avg_rev:.4f}"
                    f"\t{diff:+.4f}\t{orient}\n")


def write_ranking_tsv(ranking, ranking_mean, savepath):
    """TSV of the final scored edge ranking."""
    with open(savepath, "w") as f:
        f.write(f"rank\tu\tv\tminbit\tdecision_score\tsupport"
                f"\tmean_{ranking_mean}\n")
        for rank, (u, v, w, score, support, mean) in enumerate(ranking, 1):
            f.write(f"{rank}\t{u}\t{v}\t{w:.6f}\t{score:.6f}"
                    f"\t{support:.6f}\t{mean:.6f}\n")


# ---------------------------------------------------------------------------
# PangenomeAAIEngine — orchestrates the full pipeline.
# ---------------------------------------------------------------------------


# ANCHOR - PangenomeAAIEngine
class PangenomeAAIEngine():
    """Build a pangenome DAG from raw DIAMOND results + CONTIGS.dbs via
    AAI minbit fusion.

    Replaces the legacy ``SyntenyGeneCluster`` + ``DirectedForce`` pipeline.
    Call :py:meth:`process` to run the full pipeline; it returns the
    constructed DiGraph plus the line bookkeeping that downstream consumers
    (``PangenomeGraph.add_layers``, ``calculate_graph_distance``, the
    layout, ...) need.
    """

    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        # Input artifacts.
        self.external_genomes = A('external_genomes')
        self.genomes_storage = A('genomes_storage')
        self.diamond_search_results = A('diamond_search_results')

        # AAI engine parameters.
        self.locality_window = A('locality_window')
        self.min_line_pair_hits = A('min_line_pair_hits')
        self.orientation_tie_threshold = A('orientation_tie_threshold')
        self.min_orientation_score = A('min_orientation_score')
        self.orientation_demotion_strategy = A('orientation_demotion_strategy')
        self.ranking_components = A('ranking_components')
        self.ranking_mean = A('ranking_mean')
        self.minbit_floor = A('minbit_floor')
        self.fusion_top_bucket_k = A('fusion_top_bucket_k')
        self.fusion_seed = A('fusion_seed')

        # Filtering.
        self.min_contig_chain = A('min_contig_chain')

        # Optional debug output.
        self.tables_dir = A('tables_dir')

        # Populated by process().
        self.hash_to_genome = None
        self.genome_calls = None
        self.lines = None
        self.line_names = None
        self.line_to_genome = None


    # -----------------------------------------------------------------------
    # Sanity check (called at the top of process()).
    # -----------------------------------------------------------------------

    def sanity_check(self):
        if not self.external_genomes:
            raise ConfigError("PangenomeAAIEngine needs an external-genomes file (`--external-genomes`) "
                              "to know which CONTIGS.dbs to read. None was given :/")
        if not self.genomes_storage:
            raise ConfigError("PangenomeAAIEngine needs a genomes-storage db (`--genomes-storage`) to "
                              "resolve DIAMOND hash prefixes back to genome names :/")
        if not self.diamond_search_results:
            raise ConfigError("PangenomeAAIEngine needs the DIAMOND search results file "
                              "(`--diamond-search-results`) to build cross-genome AAI edges :/")
        if not os.path.exists(self.external_genomes):
            raise ConfigError(f"external-genomes file not found: '{self.external_genomes}' :/")
        if not os.path.exists(self.genomes_storage):
            raise ConfigError(f"genomes-storage db not found: '{self.genomes_storage}' :/")
        if not os.path.exists(self.diamond_search_results):
            raise ConfigError(f"DIAMOND search results not found: '{self.diamond_search_results}' :/")

        if int(self.locality_window or 0) < 1:
            raise ConfigError("`--locality-window` must be >= 1 :/")

        if self.tables_dir is not None:
            os.makedirs(self.tables_dir, exist_ok=True)


    # -----------------------------------------------------------------------
    # Input loading.
    # -----------------------------------------------------------------------

    def _load_genome_hash_map(self):
        """Build ``{genome_hash: genome_name}`` from the genomes-storage db
        via :py:class:`anvio.genomestorage.GenomeStorage`."""
        gs = genomestorage.GenomeStorage(
            self.genomes_storage,
            run=terminal.Run(verbose=False),
            progress=terminal.Progress(verbose=False),
        )
        hash_to_genome = {info['genome_hash']: name
                          for name, info in gs.genomes_info.items()}
        self.run.info('Genomes in storage', len(hash_to_genome))
        return hash_to_genome


    def _load_genome_calls(self, genome_paths):
        """Read every CONTIGS.db listed in ``genome_paths`` via
        :py:class:`anvio.dbops.ContigsSuperclass`.

        Returns ``{genome_name: [{'gene_callers_id', 'contig', 'start',
        'stop', 'direction'}, ...]}``.
        """
        genome_calls = {}
        total_calls = 0

        self.progress.new("Reading gene calls", progress_total_items=len(genome_paths))
        for genome, db_path in sorted(genome_paths.items()):
            self.progress.update(f"{genome} ...")
            self.progress.increment()
            if not os.path.exists(db_path):
                self.progress.end()
                raise ConfigError(f"CONTIGS.db for genome '{genome}' not found at '{db_path}' :/")
            cs_args = Namespace(contigs_db=db_path)
            cs = dbops.ContigsSuperclass(
                cs_args,
                r=terminal.Run(verbose=False),
                p=terminal.Progress(verbose=False),
            )
            calls = []
            for gid, info in cs.genes_in_contigs_dict.items():
                calls.append({
                    "gene_callers_id": gid,
                    "contig": info["contig"],
                    "start": info["start"],
                    "stop": info["stop"],
                    "direction": info["direction"],
                })
            genome_calls[genome] = calls
            total_calls += len(calls)
        self.progress.end()

        self.run.info('Total gene calls loaded', pp(total_calls))
        return genome_calls


    # -----------------------------------------------------------------------
    # DIAMOND parsing.
    # -----------------------------------------------------------------------

    def _collect_self_bitscores(self):
        """First DIAMOND pass: ``{seqid: max self-bitscore}`` from rows where
        ``qseqid == sseqid``.  Used to normalize cross-genome bitscores."""
        self.progress.new("DIAMOND pass 1/2")
        self.progress.update("collecting self-bitscores ...")
        self_bs = {}
        n_rows = 0
        with open(self.diamond_search_results) as f:
            for line in f:
                n_rows += 1
                if n_rows % 1_000_000 == 0:
                    self.progress.update(f"collecting self-bitscores ({pp(n_rows)} rows) ...")
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                qseqid, sseqid = parts[0], parts[1]
                if qseqid != sseqid:
                    continue
                try:
                    bs = float(parts[11])
                except ValueError:
                    continue
                if bs <= 0:
                    continue
                prev = self_bs.get(qseqid)
                if prev is None or bs > prev:
                    self_bs[qseqid] = bs
        self.progress.end()
        self.run.info('Self-bitscores collected', pp(len(self_bs)))
        return self_bs


    def _parse_diamond_minbit(self, hash_to_genome, gene_to_contig,
                              self_bs, gene_clusters=None):
        """Second DIAMOND pass: minbit = bitscore / min(self[q], self[s]),
        reciprocal-averaged across both directions, dropped if below
        ``self.minbit_floor`` or (in panmode) if the endpoints are in
        different gene clusters.

        Returns ``{frozenset({u, v}): minbit}`` keyed by gene endpoint
        strings ``"<genome>_<contig>:<gene_callers_id>"``.
        """
        self.progress.new("DIAMOND pass 2/2")
        self.progress.update("computing minbit per cross-genome pair ...")

        dir_pair = {}
        n_rows = 0
        n_self_row = 0
        n_self_genome = 0
        n_no_self_bs = 0
        n_unknown_hash = 0
        n_unknown_contig = 0
        n_bad_format = 0
        n_kept_dir = 0
        n_cluster_mismatch = 0
        n_cluster_unknown = 0

        with open(self.diamond_search_results) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 12:
                    n_bad_format += 1
                    continue
                n_rows += 1
                if n_rows % 1_000_000 == 0:
                    self.progress.update(f"computing minbit ({pp(n_rows)} rows) ...")
                qseqid, sseqid = parts[0], parts[1]
                if qseqid == sseqid:
                    n_self_row += 1
                    continue
                try:
                    bs = float(parts[11])
                except ValueError:
                    n_bad_format += 1
                    continue
                if bs <= 0:
                    continue

                q_hash, q_gid = split_anvio_id(qseqid)
                s_hash, s_gid = split_anvio_id(sseqid)
                if q_hash is None or s_hash is None:
                    n_bad_format += 1
                    continue
                q_genome = hash_to_genome.get(q_hash)
                s_genome = hash_to_genome.get(s_hash)
                if q_genome is None or s_genome is None:
                    n_unknown_hash += 1
                    continue
                if q_genome == s_genome:
                    n_self_genome += 1
                    continue

                if gene_clusters is not None:
                    q_cluster = gene_clusters.get((q_genome, q_gid))
                    s_cluster = gene_clusters.get((s_genome, s_gid))
                    if q_cluster is None or s_cluster is None:
                        n_cluster_unknown += 1
                        continue
                    if q_cluster != s_cluster:
                        n_cluster_mismatch += 1
                        continue

                q_self = self_bs.get(qseqid)
                s_self = self_bs.get(sseqid)
                if q_self is None or s_self is None:
                    n_no_self_bs += 1
                    continue
                denom = min(q_self, s_self)
                if denom <= 0:
                    continue
                minbit = bs / denom

                q_contig = gene_to_contig.get((q_genome, q_gid))
                s_contig = gene_to_contig.get((s_genome, s_gid))
                if q_contig is None or s_contig is None:
                    n_unknown_contig += 1
                    continue

                u = f"{q_genome}_{q_contig}:{q_gid}"
                v = f"{s_genome}_{s_contig}:{s_gid}"
                prev = dir_pair.get((u, v))
                if prev is None or minbit > prev:
                    dir_pair[(u, v)] = minbit
                    if prev is None:
                        n_kept_dir += 1

        edges = {}
        seen = set()
        for (u, v), w in dir_pair.items():
            key = frozenset({u, v})
            if key in seen:
                continue
            rev = dir_pair.get((v, u))
            w_avg = (w + rev) / 2.0 if rev is not None else w
            if w_avg >= float(self.minbit_floor or 0.0):
                edges[key] = w_avg
            seen.add(key)

        self.progress.end()

        self.run.info('DIAMOND rows read', pp(n_rows))
        self.run.info('Self-row hits skipped', pp(n_self_row))
        self.run.info('Self-genome rows skipped', pp(n_self_genome))
        if gene_clusters is not None:
            self.run.info('Rows dropped (endpoint not in any GC)', pp(n_cluster_unknown))
            self.run.info('Rows dropped (endpoints in different GCs)', pp(n_cluster_mismatch))
        if n_no_self_bs:
            self.run.info('Rows with missing self-bitscore', pp(n_no_self_bs), mc='red')
        if n_unknown_hash:
            self.run.info('Rows with unknown genome hash', pp(n_unknown_hash), mc='red')
        if n_unknown_contig:
            self.run.info('Rows with unknown gene_callers_id', pp(n_unknown_contig), mc='red')
        if n_bad_format:
            self.run.info('Rows with bad format', pp(n_bad_format), mc='red')
        self.run.info('Unique directed pairs kept', pp(n_kept_dir))
        self.run.info('Unique undirected pairs (after reciprocal-avg)', pp(len(seen)))
        self.run.info('Passed minbit floor', pp(len(edges)), mc='green')

        return edges


    # -----------------------------------------------------------------------
    # Build lines, orientations, ranking, graph.
    # -----------------------------------------------------------------------

    def _build_lines(self, genome_calls):
        """Build ``lines``, ``line_names``, ``line_to_genome`` from
        ``genome_calls``, mirroring the on-disk layout produced by
        :py:func:`write_lines_file`."""
        lines = []
        line_names = []
        line_to_genome = {}
        for genome in sorted(genome_calls):
            by_contig = defaultdict(list)
            for gc in genome_calls[genome]:
                by_contig[gc["contig"]].append(gc)
            for contig in sorted(by_contig):
                items = sorted(by_contig[contig], key=lambda x: x["start"])
                line_name = f"{genome}_{contig}"
                line_names.append(line_name)
                lines.append([str(gc["gene_callers_id"]) for gc in items])
                line_to_genome[line_name] = genome
        return lines, line_names, line_to_genome


    def _filter_short_contigs(self, lines, line_names, line_to_genome):
        """Drop lines (contigs) shorter than ``self.min_contig_chain``."""
        threshold = int(self.min_contig_chain or 0)
        if threshold <= 1:
            return lines, line_names, line_to_genome

        dropped_per_genome = Counter()
        kept_lines = []
        kept_names = []
        kept_l2g = {}
        for line, name in zip(lines, line_names):
            if len(line) < threshold:
                dropped_per_genome[line_to_genome[name]] += 1
                continue
            kept_lines.append(line)
            kept_names.append(name)
            kept_l2g[name] = line_to_genome[name]

        n_dropped = sum(dropped_per_genome.values())
        if n_dropped:
            self.run.info(f"Contigs dropped (length < {threshold})", pp(n_dropped), mc='red')
            for genome, n in dropped_per_genome.most_common():
                self.run.info_single(f"  {genome}: {n}")
        return kept_lines, kept_names, kept_l2g


    def _compute_orientations(self, lines, line_names, line_to_genome, edges_list):
        """Run the locality scan + odd-cycle demotion, returning the data
        needed downstream plus a per-pair diagnostic ``rows`` list (used by
        :py:func:`write_orientation_tsv` when ``tables_dir`` is set)."""
        K = int(self.locality_window)
        self.run.info_single(f"Running locality scan (K={K}) on "
                             f"{P('line', len(lines))} / "
                             f"{P('edge', len(edges_list))} ...")

        fwd_sum, fwd_cnt, rev_sum, rev_cnt, total, edge_signals = (
            _aggregate_line_pairs(lines, line_names, edges_list, K))

        flips, components, inconsistencies, pair_label = compute_line_orientations(
            lines, line_names, edges_list,
            K=K,
            min_hits=int(self.min_line_pair_hits),
            tie_threshold=float(self.orientation_tie_threshold),
            min_score=float(self.min_orientation_score),
            demotion_strategy=self.orientation_demotion_strategy,
            precomputed=(fwd_sum, fwd_cnt, rev_sum, rev_cnt, total),
        )

        largest = max((len(c) for c in components), default=0)
        self.run.info('Orientation components',
                      f"{len(components)} (largest covers {largest}/{len(lines)} lines)")
        self.run.info('Odd-cycle contradictions detected',
                      f"{len(inconsistencies)} "
                      f"({'globally consistent' if not inconsistencies else 'resolved by demotion'})")

        # Build per-pair diagnostic rows for the orientation TSV.
        rows = []
        for key in sorted(total.keys()):
            li_a, li_b = key
            n_edges = total[key]
            n_fwd = fwd_cnt[key]
            n_rev = rev_cnt[key]
            avg_fwd = fwd_sum[key] / n_edges
            avg_rev = rev_sum[key] / n_edges
            diff = avg_fwd - avg_rev
            if n_edges < int(self.min_line_pair_hits):
                orient = "low-hits"
            elif abs(diff) < float(self.orientation_tie_threshold):
                orient = "tie"
            elif max(avg_fwd, avg_rev) < float(self.min_orientation_score):
                orient = "weak"
            elif key in pair_label:
                orient = pair_label[key]
            else:
                orient = "demoted"
            line_a = line_names[li_a]
            line_b = line_names[li_b]
            rows.append((line_a, line_b,
                         line_to_genome.get(line_a, line_a),
                         line_to_genome.get(line_b, line_b),
                         n_edges, n_fwd, avg_fwd, n_rev, avg_rev, diff, orient))

        return (fwd_sum, fwd_cnt, rev_sum, rev_cnt, total, edge_signals,
                pair_label, rows)


    def _compute_ranking(self, lines, line_names, edges_list,
                         total, edge_signals, pair_label):
        """Attach per-edge decision_score + per-pair support, combine the
        ranking components, and sort.

        Returns ``ranking`` (list of 6-tuples)."""
        components = [c.strip() for c in (self.ranking_components or "").split(",") if c.strip()]
        valid = ("minbit", "decision", "support")
        if not components:
            raise ConfigError(f"`--ranking-components` must list at least one of {valid} :/")
        bad = [c for c in components if c not in valid]
        if bad:
            raise ConfigError(f"`--ranking-components`: unknown component(s) {bad}; "
                              f"choices: {list(valid)} :/")

        support_by_pair = {}
        for (li_a, li_b), n in total.items():
            denom = min(len(lines[li_a]), len(lines[li_b]))
            support_by_pair[(li_a, li_b)] = (n / denom) if denom > 0 else 0.0

        line_idx_by_name = {nm: i for i, nm in enumerate(line_names)}
        scored = []
        for (u, v, w) in edges_list:
            u_line = u.rsplit(":", 1)[0] if ":" in u else u
            v_line = v.rsplit(":", 1)[0] if ":" in v else v
            li_u = line_idx_by_name.get(u_line)
            li_v = line_idx_by_name.get(v_line)
            if li_u is None or li_v is None:
                scored.append((u, v, w, 0.0, 0.0))
                continue
            pair_key = (li_u, li_v) if li_u < li_v else (li_v, li_u)
            label = pair_label.get(pair_key)
            fwd, rev = edge_signals.get((u, v), (None, None))
            if label == "same":
                score = fwd if fwd is not None else 0.0
            elif label == "flip":
                score = rev if rev is not None else 0.0
            else:
                score = 0.0
            support = support_by_pair.get(pair_key, 0.0)
            scored.append((u, v, w, score, support))

        pickers = {
            "minbit":   lambda w, sc, sp: w,
            "decision": lambda w, sc, sp: sc,
            "support":  lambda w, sc, sp: sp,
        }
        sel = [pickers[c] for c in components]
        if self.ranking_mean == "geometric":
            def mean_fn(w, sc, sp):
                prod = 1.0
                for p_fn in sel:
                    v = p_fn(w, sc, sp)
                    if v <= 0:
                        return 0.0
                    prod *= v
                return prod ** (1.0 / len(sel))
        elif self.ranking_mean == "arithmetic":
            def mean_fn(w, sc, sp):
                return sum(p_fn(w, sc, sp) for p_fn in sel) / len(sel)
        else:
            raise ConfigError(f"Unknown `--ranking-mean`: {self.ranking_mean!r}; "
                              f"expected 'geometric' or 'arithmetic' :/")

        ranking = sorted(
            ((u, v, w, score, support, mean_fn(w, score, support))
             for (u, v, w, score, support) in scored),
            key=lambda r: r[5],
            reverse=True,
        )

        self.run.info('Edges ranked',
                      f"{len(ranking)} (components={','.join(components)}, "
                      f"mean={self.ranking_mean})")
        return ranking


    def _build_pangenome_graph(self, lines, line_names, line_to_genome,
                               pair_label, ranking):
        """Build the pangenome DiGraph via frontier-growth fusion.

        Walks ``ranking`` top-down with a stochastic top-K bucket over
        edges whose endpoint-lines overlap the current frontier of
        ``in_g_lines`` (Prim-style growth).  For each picked edge:

        * Both lines already in G -> attempt a fuse.
        * One line off-graph, labeled pair -> add the off line in the
          implied orientation, then fuse.
        * One line off-graph, tie pair -> stash an anchor; commit the off
          line once it has >=2 anchors with distinct line- and graph-topo
          positions.
        * Both lines off-graph -> seed two new lines in the natural
          orientation (Line_X) and flipped if labeled flip (Line_Y).

        Returns ``(G, rejected_edges, uncommitted_lines, in_g_flip)``.
        """
        rng = random.Random(int(self.fusion_seed) if self.fusion_seed is not None else None)
        top_k = max(1, int(self.fusion_top_bucket_k))
        min_score = 0.0  # honor minbit_floor + ranking pipeline already

        G = nx.DiGraph()
        super_of_gene = {}
        uf_parent = {}
        in_g_lines = set()
        in_g_flip = {}
        anchors = defaultdict(list)
        pending_merges = defaultdict(list)
        rejected_edges = []

        line_idx_by_name = {nm: i for i, nm in enumerate(line_names)}

        pool = []
        for (u, v, w, score, support, mean) in ranking:
            if mean < min_score:
                continue
            u_line = u.rsplit(":", 1)[0] if ":" in u else u
            v_line = v.rsplit(":", 1)[0] if ":" in v else v
            li_u = line_idx_by_name.get(u_line)
            li_v = line_idx_by_name.get(v_line)
            if li_u is None or li_v is None or li_u == li_v:
                continue
            pool.append([u, v, w, score, support, mean, li_u, li_v])
        pool.sort(key=lambda r: r[5], reverse=True)

        def take_bucket(must_touch_g):
            bucket, indices = [], []
            for idx, row in enumerate(pool):
                if must_touch_g and (row[6] not in in_g_lines
                                     and row[7] not in in_g_lines):
                    continue
                bucket.append(row)
                indices.append(idx)
                if len(bucket) >= top_k:
                    break
            return bucket, indices

        total_pool = len(pool)
        self.progress.new("Pangenome fusion", progress_total_items=total_pool)
        processed = 0
        while pool:
            if not in_g_lines:
                bucket, indices = take_bucket(must_touch_g=False)
            else:
                bucket, indices = take_bucket(must_touch_g=True)
                if not bucket:
                    bucket, indices = take_bucket(must_touch_g=False)
                    if not bucket:
                        break
            pick = rng.randrange(len(bucket))
            u, v, w, score, support, mean, li_u, li_v = bucket[pick]
            pool.pop(indices[pick])
            processed += 1
            if processed & 0x3FF == 0:
                self.progress.update(f"lines {len(in_g_lines)}/{len(lines)}  "
                                     f"nodes {G.number_of_nodes()}  "
                                     f"rejects {len(rejected_edges)}")
            self.progress.increment()

            pair_key = (li_u, li_v) if li_u < li_v else (li_v, li_u)
            label = pair_label.get(pair_key)
            u_in_g = li_u in in_g_lines
            v_in_g = li_v in in_g_lines

            if u_in_g and v_in_g:
                _fuse_genes(G, u, v, w, score, mean, label,
                            super_of_gene, uf_parent, rejected_edges,
                            line_to_genome)
            elif u_in_g or v_in_g:
                if u_in_g:
                    off_li, off_gene, in_li, in_gene = li_v, v, li_u, u
                else:
                    off_li, off_gene, in_li, in_gene = li_u, u, li_v, v

                if label in ("same", "flip"):
                    off_flip = in_g_flip[in_li] ^ (label == "flip")
                    _add_line_to_graph(G, off_li, lines, line_names, off_flip,
                                       super_of_gene, uf_parent,
                                       in_g_lines, in_g_flip)
                    _fuse_genes(G, u, v, w, score, mean, label,
                                super_of_gene, uf_parent, rejected_edges,
                                line_to_genome)
                else:
                    off_tok = off_gene.rsplit(":", 1)[1]
                    try:
                        pos = lines[off_li].index(off_tok)
                    except ValueError:
                        rejected_edges.append((u, v, mean, "line-position-unknown"))
                        continue
                    anchors[off_li].append({
                        "pos": pos,
                        "in_gene": in_gene,
                        "off_gene": off_gene,
                    })
                    pending_merges[off_li].append((u, v, w, score, mean, label))
                    if _try_commit_tie_line(off_li, lines, line_names,
                                            anchors, pending_merges, G,
                                            super_of_gene, uf_parent,
                                            in_g_lines, in_g_flip,
                                            rejected_edges,
                                            line_to_genome):
                        anchors.pop(off_li, None)
                        pending_merges.pop(off_li, None)
            else:
                fu = False
                fv = True if label == "flip" else False
                _add_line_to_graph(G, li_u, lines, line_names, fu,
                                   super_of_gene, uf_parent,
                                   in_g_lines, in_g_flip)
                _add_line_to_graph(G, li_v, lines, line_names, fv,
                                   super_of_gene, uf_parent,
                                   in_g_lines, in_g_flip)
                _fuse_genes(G, u, v, w, score, mean, label,
                            super_of_gene, uf_parent, rejected_edges,
                            line_to_genome)

        self.progress.end()

        uncommitted_lines = sorted(li for li in range(len(lines))
                                   if li not in in_g_lines)
        for li in uncommitted_lines:
            for (u, v, w, score, mean, label) in pending_merges.get(li, []):
                rejected_edges.append((u, v, mean, "insufficient-anchors"))

        return G, rejected_edges, uncommitted_lines, in_g_flip


    # -----------------------------------------------------------------------
    # Orchestration.
    # -----------------------------------------------------------------------

    def process(self, gene_clusters=None):
        """Run the full pipeline.

        Parameters
        ==========
        gene_clusters : dict or None
            ``{(genome_name, gene_callers_id): gene_cluster_name}``.  In
            panmode this is loaded from the pan-db by the caller and
            passed in; in non-panmode pass ``None``.

        Returns
        =======
        G : networkx.DiGraph
            The pangenome DAG with named super-nodes.  Each node has
            ``genes``, ``lines``, ``type``, ``component_id`` attributes.
        lines : list[list[str]]
            Per-contig ordered gene_caller_id token sequences (post
            short-contig filter).
        line_names : list[str]
            ``{genome}_{contig}`` line identifiers, parallel to ``lines``.
        line_to_genome : dict[str, str]
            ``line_name -> genome_name`` lookup.
        in_g_flip : dict[int, bool]
            For each committed line index, whether the line was added in
            reversed orientation.
        """
        self.run.warning(None, header="BUILDING PANGENOME GRAPH (AAI ENGINE)", lc="green")

        self.sanity_check()

        self.run.info('External genomes', self.external_genomes)
        self.run.info('Genomes storage', self.genomes_storage)
        self.run.info('DIAMOND search results', self.diamond_search_results)
        if self.tables_dir:
            self.run.info('Intermediate tables dir', self.tables_dir, mc='green')
        if gene_clusters is None:
            self.run.info_single("No pan-db given — alignments will be unavailable downstream "
                                 "and node types collapse to core / accessory / singleton.",
                                 mc='yellow')
        else:
            self.run.info('Gene-cluster assignments', pp(len(gene_clusters)))

        # 1. Read GENOMES.db.
        hash_to_genome = self._load_genome_hash_map()

        # 2. Read external-genomes.txt + cross-check against GENOMES.db.
        genome_paths = parse_external_genomes(self.external_genomes)
        self.run.info('Genomes in external-genomes file', len(genome_paths))
        in_storage = set(hash_to_genome.values())
        missing = sorted(set(genome_paths) - in_storage)
        if missing:
            head = ', '.join(missing[:3]) + ('...' if len(missing) > 3 else '')
            self.run.warning(f"{len(missing)} genome(s) in external-genomes are not in the "
                             f"genomes-storage db and their DIAMOND hits will be ignored: {head}",
                             lc='yellow')

        # 3. Read CONTIGS.dbs.
        genome_calls = self._load_genome_calls(genome_paths)
        gene_to_contig = build_gene_to_contig(genome_calls)

        # 4. DIAMOND parsing.
        self_bs = self._collect_self_bitscores()
        edges_dict = self._parse_diamond_minbit(
            hash_to_genome, gene_to_contig, self_bs,
            gene_clusters=gene_clusters)

        # Dump tables if requested.
        if self.tables_dir:
            self._dump_tables_input_stage(genome_calls, edges_dict)

        # 5. Build lines, then drop short contigs.
        lines, line_names, line_to_genome = self._build_lines(genome_calls)
        lines, line_names, line_to_genome = self._filter_short_contigs(
            lines, line_names, line_to_genome)
        self.run.info('Lines retained', len(lines))

        # 6. Convert edges dict to the (u, v, w) list the orientation /
        #    ranking pipeline expects.
        edges_list = [(*sorted(pair), w) for pair, w in edges_dict.items()]

        # 7. Orientation scan + per-pair labeling.
        (fwd_sum, fwd_cnt, rev_sum, rev_cnt, total, edge_signals,
         pair_label, orientation_rows) = self._compute_orientations(
            lines, line_names, line_to_genome, edges_list)

        if self.tables_dir:
            write_orientation_tsv(
                orientation_rows,
                os.path.join(self.tables_dir, "orientation.tsv"))

        # 8. Edge ranking.
        ranking = self._compute_ranking(
            lines, line_names, edges_list, total, edge_signals, pair_label)

        if self.tables_dir:
            write_ranking_tsv(
                ranking, self.ranking_mean,
                os.path.join(self.tables_dir, "ranking.tsv"))

        # 9. Build the pangenome graph.
        self.run.info_single(
            f"Building pangenome graph (top-bucket-k={self.fusion_top_bucket_k}, "
            f"seed={self.fusion_seed}) ...")
        G, rejected_edges, uncommitted_lines, in_g_flip = self._build_pangenome_graph(
            lines, line_names, line_to_genome, pair_label, ranking)

        reasons = Counter(reason for (_u, _v, _m, reason) in rejected_edges)
        self.run.info('Graph nodes', pp(G.number_of_nodes()), mc='green')
        self.run.info('Graph edges', pp(G.number_of_edges()))
        self.run.info('Lines committed',
                      f"{len(in_g_flip)} / {len(lines)} "
                      f"({len(uncommitted_lines)} uncommitted)")
        self.run.info('Fuses rejected',
                      f"{len(rejected_edges)} ({dict(reasons) if reasons else '{}'})")

        # 10. Post-processing: rename, type, component id.
        G = rename_pangenome_nodes(
            G, gene_clusters=gene_clusters, line_to_genome=line_to_genome)
        G = compute_node_types(G, line_to_genome, gene_clusters=gene_clusters)
        G = compute_component_ids(G)

        type_counts = Counter(d.get("type", "") for _, d in G.nodes(data=True))
        self.run.info('Node types', dict(type_counts))
        component_counts = Counter(d["component_id"] for _, d in G.nodes(data=True))
        largest = component_counts.most_common(1)[0][1] if component_counts else 0
        self.run.info('Components',
                      f"{len(component_counts)} (largest = {largest} nodes)")

        # Cache for callers that may want to re-inspect.
        self.hash_to_genome = hash_to_genome
        self.genome_calls = genome_calls
        self.lines = lines
        self.line_names = line_names
        self.line_to_genome = line_to_genome

        return G, lines, line_names, line_to_genome, in_g_flip


    # -----------------------------------------------------------------------
    # Debug-table dumping helper (called when self.tables_dir is set).
    # -----------------------------------------------------------------------

    def _dump_tables_input_stage(self, genome_calls, edges_dict):
        """Write the post-parse, pre-orientation tables to ``self.tables_dir``.
        Orientation + ranking dumps happen later, once those stages run."""
        d = self.tables_dir
        write_lines_file(genome_calls, os.path.join(d, "lines.txt"))
        write_edges_file(edges_dict, os.path.join(d, "edges.tsv"))
        write_genome_map(genome_calls, os.path.join(d, "genome_map.tsv"))
        self.run.info_single(f"Intermediate tables written to {d}", mc='green')
