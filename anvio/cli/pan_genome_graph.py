#!/usr/bin/env python
# -*- coding: utf-8
"""Compute a graph representation of a pangenome"""

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
__authors__ = ['ahenoch']
__requires__ = ['pan-db', 'genomes-storage-db', 'external-genomes']
__provides__ = ['pan-graph-json']
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
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': False}))
    groupA.add_argument(*anvio.A('genomes-names'), **anvio.K('genomes-names', {'required': False}))
    groupA.add_argument(*anvio.A('pan-graph-db'), **anvio.K('pan-graph-db', {'required': False}))
    groupA.add_argument(*anvio.A('pan-graph-yaml'), **anvio.K('pan-graph-yaml', {'required': False}))

    groupB = parser.add_argument_group('OUTPUT', "By default, this program will store the resulting pangraph into the anvi'o pan-db "
                                "you have provided as input above, so it is accessible to downstream analyses seamlessly. Using the "
                                "the parameters below, you can ask anvi'o to store the resulting graph into a text output file (which "
                                "may be useful for developers for debugging purposes) or you can ask anvi'o to skip adding the results "
                                "to the pan-db.")
    groupB.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))
    groupB.add_argument('--output-synteny-gene-cluster-dendrogram', default=False, action="store_true", help="Write a dendrogram (.svg) "
                                "and distance matrix (.tsv) for each split SynGC to the output directory (for debugging "
                                "the paralog splitter).")

    groupC = parser.add_argument_group('GRAPH BUILDING & SPLITTING', "Controls how SynGC are generated and filtered.")

    groupC.add_argument('--circularize', default=False, action="store_true", help = "Connect contig ends back to starts "
                    "(only sensible for single-contig/circular genomes); can add cycles otherwise.")
    groupC.add_argument('--min-contig-chain', default=5, type=int, help = "Skip contigs with fewer than this many SynGCs "
                    "(filters very short/fragmented contigs).")
    groupC.add_argument('--min-k', default=1, type=int, help = "Minimum k-mer window size around each gene before splitting paralogs "
                    "(will auto-increase until each k-mer is genome-unique; raise to demand more context).")
    groupC.add_argument('--alpha', default=0.5, type=float, help = "Global context similarity cutoff (single-copy-core flank similarity "
                    "must be >= alpha to be considered same context; lower alpha merges more, higher splits more).")
    groupC.add_argument('--n', default=50, type=int, help = "Max number of single-copy-core neighbors to scan on each side when scoring "
                    "context similarity (0 disables global-context filtering).")
    groupC.add_argument('--beta', default=0.5, type=float, help = "Penalty weight for orientation mismatches in local k-mer comparisons "
                    "(higher beta penalizes strand flips more).")
    groupC.add_argument('--gamma', default=0.25, type=float, help = "Penalty weight for gaps/contig-edge markers in local k-mer comparisons "
                    "(higher gamma penalizes contig ends/missing neighbors more).")
    groupC.add_argument('--delta', default=0.75, type=float, help = "Maximum allowed distance for k-mer similarity "
                    "(values above delta are treated as mismatches; raise to be more permissive).")
    groupC.add_argument('--inversion-aware', default=False, action="store_true", help = "Also compare reversed k-mers, allowing inverted "
                    "contexts to cluster together (helps when inversions are common 🤞).")

    groupD = parser.add_argument_group('ANCHORING & PRIORITY', "Choose a reference genome or anchor gene for layout.")

    groupD.add_argument('--priority-genome', default='', type=str, help = "Genome name to prioritize when building edges "
                    "(its path gets extra weight so layout favors its ordering).")
    groupD.add_argument('--start-gene', default=None, type=str, help = "Regex/text to pick a SynGC as the "
                    "starting node for layout (looked up in --start-column; helpful to anchor the graph). "
                    "If not specified, anvi'o will use the synteny cluster with the smallest avg. gene caller ID "
                    "across all genomes to serve as the starting node automatically. If you have reoriented your "
                    "genomes, this should give you best results where the graph starts at the SynGC that describes "
                    "the first genes in your genomes.")
    groupD.add_argument('--start-column', default='COG24_FUNCTION_TEXT', type=str, help = "Annotation column to search for --start-gene.")

    groupE = parser.add_argument_group('LAYOUT & SIMPLIFICATION', "Controls how the graph is compressed and long edges are filtered.")

    groupE.add_argument('--gene-cluster-grouping-threshold', default=-1, type=int, help = "Compress linear chains of nodes of "
                    "this length or longer into groups (-1 disables grouping; useful to simplify long conserved runs).")
    groupE.add_argument('--grouping-compression', default=1.0, type=float, help = "Compression factor for grouped chains "
                    "(1.0 = no compression; lower values squeeze grouped nodes closer in the layout).")
    groupE.add_argument('--max-edge-length-filter', default=-1, type=int, help = "Remove edges longer than this length threshold "
                    "after layout (-1 disables filtering; lower values prune long jumps/repeats).")

    groupF = parser.add_argument_group('METADATA & LAYERS', "Display and metadata options for the resulting pan-graph.")

    groupF.add_argument('--project-name', default=None, help = "Optional name stored in the pan-graph-db metadata (for display/export).")
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
