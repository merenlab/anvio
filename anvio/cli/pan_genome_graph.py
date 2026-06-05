#!/usr/bin/env python
"""Compute a graph representation of a pangenome"""

import argparse
import sys

import anvio
import anvio.panops as panops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ahenoch', 'meren']
__requires__ = ['pan-db']
__can_use__ = ['genomes-storage-db', 'external-genomes']
__provides__ = ['pan-graph-db']
__description__ = ("An anvi'o program to compute a graph representation of pangenomes. It will do its magic, and store it into your "
                   "pan-db, or report a JSON formatted graph file, for downstream visualization and analyses with `anvi-display-pan-graph`")
__resources__ = []


@terminal.time_program
def main():
    args = get_args()

    try:
        pangraph = panops.PangenomeGraph(args)
        pangraph.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', "Anvi'o artifacts for the pan graph to be computed.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db', {'required': False}))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': True}))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': True}))
    groupA.add_argument(*anvio.A('genomes-names'), **anvio.K('genomes-names', {'required': False}))
    groupA.add_argument(*anvio.A('diamond-search-results'), **anvio.K('diamond-search-results', {'required': True}))

    groupB = parser.add_argument_group('OUTPUT', "Where the resulting pan-graph-db should be written. If you don't specify a path, "
                                "anvi'o will write the database into the current working directory using the project name as the "
                                "filename (i.e., './<project_name>-PAN-GRAPH.db').")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file', {'required': False,
                                'help': "Output file path for the pan-graph-db this program will generate."}))
    groupB.add_argument(*anvio.A('tables-dir'), **anvio.K('tables-dir'))

    groupC = parser.add_argument_group('GRAPH BUILDING', "Controls how the gene-order graph is constructed from genome inputs.")

    groupC.add_argument('--min-contig-chain', default=5, type=int, help = "Skip contigs with fewer than this many genes "
                    "(filters very short/fragmented contigs).")
    groupC.add_argument('--no-remerge', default=False, action='store_true',
                    help = "Skip the post-engine remerge pass. By default (no flag), anvi'o walks every parent "
                    "gene cluster that produced multiple synteny clusters and merges pairs that look like "
                    "engine over-splits: same weakly connected component, not in-series in the DAG, and "
                    "disjoint genome coverage. Pass this flag to disable that pass and keep the engine's "
                    "original splits intact.")
    groupC.add_argument('--remerge-max-length', default=10, type=int,
                    help = "Asymmetry cap (in DAG edges) for remerge candidates. For each candidate pair, "
                    "the remerge pass measures shortest-path distances from BOTH the lowest common ancestor "
                    "(LCA) and the lowest common descendant (LCD), where they exist, to each of the two "
                    "nodes. If either |dist_a - dist_b| exceeds this value the merge is rejected. This "
                    "blocks topologically lopsided merges where one candidate sits much deeper/shallower "
                    "than the other relative to the shared ancestor or descendant. Default -1 disables the "
                    "check. When only one of LCA or LCD exists, only that side's check applies; pairs with "
                    "neither (extremely rare on a connected DAG) are not rejected by this guard.")
    groupC.add_argument('--no-include-non-coding-genes', default=False, action='store_true',
                    help = "Drop non-coding gene calls (rRNA, tRNA, etc.) from the graph. By default (no flag), "
                    "non-coding gene calls are kept as singleton super-nodes typed `rna`. Non-coding genes are "
                    "detected via each genome's CONTIGS.db `call_type` column (anything that is not CODING). "
                    "They have no DIAMOND hits (DIAMOND runs on translated proteins) and therefore never merge "
                    "across genomes -- each non-coding gene becomes its own node. Pass this flag to drop them "
                    "before component-id and layout computation; flanking gene-order edges are re-stitched so "
                    "line continuity is preserved.")

    groupD = parser.add_argument_group('AAI ENGINE PARAMETERS', "Parameters controlling the AAI-based gene-endpoint fusion engine. "
                    "These rarely need tuning; defaults work well for most datasets.")

    groupD.add_argument('--locality-window', default=10, type=int, help = "Flanking window size (in genes) on each side of a candidate "
                    "line pair, used to decide whether two contigs are co-oriented or flipped.")
    groupD.add_argument('--min-window-completeness', default=0, type=int, help = "Minimum number of genes that each of the four flanking "
                    "windows (left/right of each endpoint) must contain for an AAI edge to be kept. Must be 0..--locality-window. "
                    "0 disables the filter; --locality-window requires a fully populated window on both sides of both endpoints. "
                    "Edges whose any flanking window is shorter than this are trashed before orientation scoring and fusion.")
    groupD.add_argument('--min-line-pair-hits', default=100, type=int, help = "Minimum number of DIAMOND hits between two contig lines "
                    "required to consider them for orientation scoring.")
    groupD.add_argument('--orientation-tie-threshold', default=0.2, type=float, help = "Score margin under which a line-pair orientation "
                    "call is considered a tie and demoted (see --orientation-demotion-strategy).")
    groupD.add_argument('--min-orientation-score', default=0.7, type=float, help = "Minimum orientation score for a line pair to be "
                    "accepted; pairs below this are dropped.")
    groupD.add_argument('--orientation-demotion-strategy', default='slimmest-margin',
                    choices=['slimmest-margin', 'fewest-edges'], help = "How to break ties when committing line-pair "
                    "orientations during the spanning-tree walk.")
    groupD.add_argument('--ranking-components', default='minbit,decision,support', type=str, help = "Comma-separated list "
                    "of components used to rank AAI edges before fusion (any subset of minbit, decision, support).")
    groupD.add_argument('--ranking-mean', default='geometric', choices=['geometric', 'arithmetic'], help = "Mean used to "
                    "combine ranking components into a single score per edge.")
    groupD.add_argument('--minbit-floor', default=0.0, type=float, help = "Edges with AAI minbit below this value are dropped "
                    "before ranking.")
    groupD.add_argument('--min-ranking-score', default=0.05, type=float, help = "Edges whose combined ranking score (the mean of "
                    "--ranking-components) is below this value are skipped during fusion. Independent of --minbit-floor, which gates "
                    "on raw minbit; this gates on the aggregated score.")
    groupD.add_argument('--fusion-top-bucket-k', default=5, type=int, help = "Number of top-ranked edges sampled per fusion step "
                    "(Prim-style frontier growth with stochastic top-K bucketing).")
    groupD.add_argument('--fusion-seed', default=42, type=int, help = "Random seed for the stochastic top-K fusion sampler.")

    groupE = parser.add_argument_group('LAYOUT & SIMPLIFICATION', "Controls how the graph is compressed and long edges are filtered.")

    groupE.add_argument('--component', default=0, type=int, help = "Which weakly connected component to lay out and summarize. "
                    "Components are indexed largest-first; the default (0) selects the largest. All components are still persisted "
                    "in the pan-graph-db.")
    groupE.add_argument('--gene-cluster-grouping-threshold', default=-1, type=int, help = "Compress linear chains of nodes of "
                    "this length or longer into groups (-1 disables grouping; useful to simplify long conserved runs).")
    groupE.add_argument('--grouping-compression', default=1.0, type=float, help = "Compression factor for grouped chains "
                    "(1.0 = no compression; lower values squeeze grouped nodes closer in the layout).")
    groupE.add_argument('--max-edge-length-filter', default=-1, type=int, help = "Remove edges longer than this length threshold "
                    "after layout (-1 disables filtering; lower values prune long jumps/repeats).")

    groupF = parser.add_argument_group('METADATA & LAYERS', "Display and metadata options for the resulting pan-graph.")

    groupF.add_argument(*anvio.A('description'), **anvio.K('description'))
    groupF.add_argument(*anvio.A('project-name'), **anvio.K('project-name'))
    groupF.add_argument('--load-state', default='default', type=str, help="Initial display state name to store/use in the pan-graph-db.")
    groupF.add_argument('--import-values', default='start,stop,partial,call_type', type=str, help = "Comma-separated "
                    "numeric columns from the pangenome table to copy as node layers (e.g., start, stop, partial); "
                    "non-numeric/missing columns are ignored.")

    groupG = parser.add_argument_group('ROBUSTNESS & DEBUGGING', "Controls safety checks and noisy-data handling.")

    groupG.add_argument('--just-do-it', default=False, action="store_true", help = "Bypass safety checks (e.g., "
                    "over-splitting) and continue even when warnings would normally abort; use to diagnose tough datasets.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
