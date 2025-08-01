#!/usr/bin/env python
# -*- coding: utf-8
"""Return counts of AAs in bins, contigs, or gene caller ids"""

import sys

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.dbops import AA_counts


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["aa-frequencies-txt"]
__requires__ = ["splits-txt", "contigs-db", "profile-db", "collection"]
__description__ = ("Fetches the number of times each amino acid occurs from a contigs database in a "
                   "given bin, set of contigs, or set of genes")


@terminal.time_program
def main():
    try:
        args = get_args()
        aa_counts = AA_counts(args)
        aa_counts.report()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    group0 = parser.add_argument_group('MANDATORY STUFF', 'You have to set the following two parameters, then\
                                                           you will select one set of parameters from the\
                                                           following optional sections. If you select nothing\
                                                           from those sets, AA counts for everything in the\
                                                           contigs database will be reported.')
    group0.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': True}))
    group0.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    groupA = parser.add_argument_group('OPTIONAL PARAMS FOR BINS')
    groupA.add_argument(*anvio.A('profile-db'), **anvio.K('profile-db', {'required': False}))
    groupA.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupA.add_argument(*anvio.A('bin-ids-file'), **anvio.K('bin-ids-file'))

    groupB = parser.add_argument_group('OPTIONAL PARAMS FOR CONTIGS')
    groupB.add_argument(*anvio.A('contigs-of-interest'), **anvio.K('contigs-of-interest', {'help': "A file with\
                                                         contig names. There should be only one column in the file,\
                                                         and each line should correspond to a unique split name."}))

    groupC = parser.add_argument_group('OPTIONAL PARAMS FOR GENE CALLS')
    groupC.add_argument(*anvio.A('gene-caller-ids'), **anvio.K('gene-caller-ids'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
