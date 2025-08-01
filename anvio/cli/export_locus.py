#!/usr/bin/env python
# -*- coding: utf-8
"""anvi-export-locus, first of its name, splitter of loci, and generator of mini contigs databases

  You want to play with it? This is how you could quickly test it:

  ```
  anvi-self-test --suite export-locus

  cd export_locus_test

  anvi-export-locus -c CONTIGS.db \
                    -O OUTPUT \
                    --use-hmm \
                    --hmm-sources Bacteria_71 \
                    --search-term "Exonuc_VII_L" \
                    -n 500,500 \
                    -O P_marinus_CCMP1375
  ```
"""

import sys
from anvio.argparse import ArgumentParser

import anvio

from anvio.splitter import LocusSplitter
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mschecht', 'ekiefl', 'ShaiberAlon']
__requires__ = ['contigs-db']
__provides__ = ['locus-fasta']
__description__ = ("This program helps you cut a 'locus' from a larger genetic context (e.g., contigs, "
                   "genomes). By default, anvi'o will locate a user-defined anchor gene, "
                   "extend its selection upstream and downstream based on the --num-genes "
                   "argument, then extract the locus to create a new contigs database. The "
                   "anchor gene must be provided as --search-term, --gene-caller-ids, or --hmm-sources. "
                   "If --flank-mode is designated, you MUST provide TWO flanking genes that "
                   "define the locus region (Please see --flank-mode help for more information). "
                   "If everything goes as plan, anvi'o will give you individual locus contigs "
                   "databases for every matching anchor gene found in the original contigs "
                   "database provided. Enjoy your mini contigs databases!")


def main():
    args = get_args()

    try:
        locus_splitter = LocusSplitter(args)
        locus_splitter.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    # Essential input
    groupA = parser.add_argument_group('INPUT DATABASE')
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    # Query options for anchor gene
    groupB =  parser.add_argument_group('ANCHOR GENE IDENTIFICATION OPTION #1', "Search based on gene function annotations or HMMs")
    groupB.add_argument('-s', '--search-term', help='Search term to be found in function annotation sources.')
    groupB.add_argument(*anvio.A('case-sensitive'), **anvio.K('case-sensitive'))
    groupB.add_argument(*anvio.A('exact-match'), **anvio.K('exact-match'))
    groupB.add_argument(*anvio.A('hmm-sources'), **anvio.K('hmm-sources'))
    groupB.add_argument(*anvio.A('annotation-sources'), **anvio.K('annotation-sources'))
    groupB.add_argument('--use-hmm', default = False, action='store_true', help='Use HMM hits instead of functional annotations. \
                            In other words, --search-term will be queried against HMM source annotations, NOT functional annotations. \
                            If you choose this option, you must also say which HMM source to use.')
    groupB.add_argument(*anvio.A('list-hmm-sources'), **anvio.K('list-hmm-sources'))

    groupC =  parser.add_argument_group('ANCHOR GENE IDENTIFICATION OPTION #2', "Search based on gene caller id(s)")
    groupC.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))

    # Output options
    groupD = parser.add_argument_group('THE OUTPUT', "Where should the output go. It will be one FASTA file with all matches \
                                       or one FASTA per match (see --separate-fasta)")
    groupD.add_argument(*anvio.A('output-dir'), **anvio.K('output-dir'))
    groupD.add_argument(*anvio.A('output-file-prefix'), **anvio.K('output-file-prefix', {'required': True}))

    # Additional
    groupE = parser.add_argument_group('ADDITIONAL STUFF', "Flags and parameters you can set according to your need")
    groupE.add_argument(*anvio.A('flank-mode'), **anvio.K('flank-mode'))
    groupE.add_argument(*anvio.A('num-genes'), **anvio.K('num-genes'))
    groupE.add_argument(*anvio.A('remove-partial-hits'), **anvio.K('remove-partial-hits'))
    groupE.add_argument(*anvio.A('delimiter'), **anvio.K('delimiter'))
    groupE.add_argument(*anvio.A('never-reverse-complement'), **anvio.K('never-reverse-complement'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
