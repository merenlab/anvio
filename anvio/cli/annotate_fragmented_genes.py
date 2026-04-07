#!/usr/bin/env python
# -*- coding: utf-8
"""Identify and annotate fragmented genes (pseudogenes) in pangenomes"""

import sys

import anvio
import anvio.panops as panops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__credits__ = ["Sean Crosson"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db', 'genomes-storage-db', 'external-genomes']
__provides__ = ['functions']
__description__ = ("Identify fragmented genes (pseudogenes) in a pangenome by finding adjacent genes from "
                   "the same genome within the same gene cluster, and annotate them under a 'PSEUDO_GENES' "
                   "function source in the relevant contigs databases")
__resources__ = []


@terminal.time_program
def main():
    args = get_args()

    try:
        annotator = panops.FragmentedGeneAnnotator(args)
        annotator.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', "Provide the pangenome database, genomes storage, and an external "
                    "genomes file. The program will scan gene clusters for fragmentation events and annotate "
                    "the relevant contigs databases.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))

    groupB = parser.add_argument_group('PARAMETERS', "Fine-tune the fragmentation detection logic.")
    groupB.add_argument('--min-full-length-ratio', default=0.50, type=float,
                        help="Minimum ratio of the longest fragment length to the full-length reference "
                             "gene length for that fragment to be labeled 'fragmented_gene' rather than "
                             "'gene_fragment'. If the longest fragment in a genome is shorter than this "
                             "fraction of the reference, ALL fragments in that genome will be labeled "
                             "'gene_fragment'. Default: %(default).2f")
    groupB.add_argument('--find-stray-fragments', default=False, action='store_true',
                        help="Also look for out-of-frame gene fragments that ended up in different gene "
                             "clusters. When a premature stop codon splits a gene and the downstream "
                             "fragment is in a different reading frame, MCL places the fragment in a "
                             "separate gene cluster. This flag enables a second scan that detects such "
                             "cases by finding truncated genes whose adjacent neighbor on the same contig "
                             "belongs to a different gene cluster and together they approximate the "
                             "full-length reference.")

    groupC = parser.add_argument_group('REPORTING', "Control what this program reports and whether it gets to update "
                        "anything at all.")
    groupC.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source', {'help': "If "
                             "provided, the consensus function from this annotation source will be "
                             "shown next to each gene cluster name in the terminal report. You can "
                             "always run `anvi-db-info` on your genome storage db to remember what "
                             "function annotation sources available and what they are called."}))
    groupC.add_argument('--skip-reporting', default=False, action='store_true',
                        help="Do not print per-gene-cluster visualizations to the terminal. Annotations "
                             "will still be written to contigs databases.")
    groupC.add_argument('--report-only', default=False, action='store_true',
                        help="Only report fragmentation events to the terminal. Do not write any "
                             "annotations to contigs databases.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
