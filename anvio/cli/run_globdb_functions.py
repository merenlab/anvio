#!/usr/bin/env python

import sys

import anvio

from anvio.globdb import GlobDBFunctions
from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'dspeth']
__requires__ = ['globdb-data', 'contigs-db']
__provides__ = ['functions']
__description__ = ("Annotate genes with gene family functions derived from GlobDB with "
                   "gene family-level cutoffs determined by empirical Local Alignment "
                   "Score Ratio (LASR) thresholds.")


@time_program
def main():
    args = get_args()

    try:
        globdb = GlobDBFunctions(args)
        aa_file = args.fasta_file or None
        globdb.process(aa_sequences_file_path=aa_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group("INPUT OPTION #1: AN ANVI'O CONTIGS DATABASE",
                                       "Use this option to annotate every gene in an anvi'o contigs "
                                       "database. Annotations are stored directly in the database.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group('INPUT OPTION #2: A BORING FASTA FILE',
                                       "Use this option to annotate amino acid sequences in a FASTA "
                                       "file. Results are written as a TAB-delimited output file.")
    groupB.add_argument(*anvio.A('fasta-file'), **anvio.K('fasta-file',
                        {'help': 'A FASTA-formatted file containing amino acid sequences to annotate.'}))
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    groupC = parser.add_argument_group('GENERAL OPTIONS')
    groupC.add_argument(*anvio.A('globdb-data-dir'), **anvio.K('globdb-data-dir'))
    groupC.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupC.add_argument(*anvio.A('temporary-dir-path'), **anvio.K('temporary-dir-path'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
