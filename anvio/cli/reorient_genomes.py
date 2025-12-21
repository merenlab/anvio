#!/usr/bin/env python
# coding: utf-8

import sys

import anvio
import anvio.terminal as terminal

from anvio.argparse import ArgumentParser
from anvio.genomereorientation import GenomeReorienter
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2025, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['fasta-txt']
__provides__ = ['fasta']
__description__ = ("Reorient circular genomes in a fasta-txt so their coordinates match a chosen reference genome "
                   "using minimap2 + seqkit rotation.")


@terminal.time_program
def main():
    args = get_args()

    try:
        reorienter = GenomeReorienter(args)
        reorienter.process()

    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES')
    groupA.add_argument('--fasta-txt', required=True,
                        help="Two-column TAB-delimited file with genome name and FASTA file path.")
    groupA.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir', {'required': True}))

    groupRef = parser.add_argument_group('SELECTION OF REFERENCE',
                                         "If you do not set --reference, the program will auto-select a reference by "
                                         "choosing the genome with the fewest contigs; ties are broken by the largest "
                                         "total length. This mirrors how anvi'o tools often prioritize more complete assemblies.")
    groupRef.add_argument('--reference', required=False,
                          help="Genome name in fasta-txt to use as the reference orientation. If omitted, auto-selection applies.")
    groupRef.add_argument('--use-dnaa-for-reference-orientation', action='store_true',
                          help="Use DnaA gene location to orient the reference genome. The program will identify the DnaA "
                               "gene using an HMM profile (Bac_DnaA_C from Pfam), and rotate the reference to start near "
                               "the DnaA gene, which typically marks the origin of replication in bacterial genomes. This "
                               "option is useful for bacterial genomes but may not work well for plasmids or viral genomes "
                               "without DnaA.")

    groupViz = parser.add_argument_group('VISUALIZATION',
                                         "By default, the program will generate synteny ribbon plots showing "
                                         "alignment patterns between the reference and each query genome before "
                                         "and after reorientation/scaffolding. These visualizations help assess "
                                         "the quality of the reorientation process.")
    groupViz.add_argument('--skip-visualizing-alignments', action='store_true',
                          help="Skip generating alignment visualizations. This can speed up processing when "
                               "you only need the reoriented FASTA files without visual inspection.")
    groupViz.add_argument('--plot-width', type=int,
                          help="Width of alignment plots in characters. By default, plots automatically scale "
                               "to use the full terminal width. Use this to override with a specific width. "
                               "Note: widths below 100 characters may not display properly.")
    groupViz.add_argument('--plot-height', type=int, default=20,
                          help="Height of alignment plots in characters. Default: %(default)s")

    groupB = parser.add_argument_group('ADVANCED')
    groupB.add_argument('--threads', type=int, default=1,
                        help="Number of threads for minimap2.")
    groupB.add_argument('--log-file-path',
                        help="Write a detailed log to this file (otherwise only concise reporting is printed).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
