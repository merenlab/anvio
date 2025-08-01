#!/usr/bin/env python
# -*- coding: utf-8
"""A program to generate a functions across groups stats output."""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomedescriptions import AggregateFunctions


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'adw96']
__requires__ = ['functions', 'genomes-storage-db', 'internal-genomes', 'external-genomes']
__provides__ = ['interactive']
__description__ = "Generate a TAB delimited file for the distribution of functions across groups of genomes/metagenomes"


@terminal.time_program
def main():
    args = get_args()

    try:
        facc = AggregateFunctions(args)
        facc.report_functions_per_group_stats(args.output_file)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('GENOMES', "Tell anvi'o where your genomes are.")
    groupA.add_argument(*anvio.A('internal-genomes'), **anvio.K('internal-genomes'))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage'))

    groupB = parser.add_argument_group('FUNCTIONS', "Tell anvi'o which functional annotation source you like above all, and other "
                                "important details you like about your analysis.")
    groupB.add_argument(*anvio.A('annotation-source'), **anvio.K('annotation-source', {'required': True}))
    groupB.add_argument(*anvio.A('aggregate-based-on-accession'), **anvio.K('aggregate-based-on-accession'))
    groupB.add_argument(*anvio.A('aggregate-using-all-hits'), **anvio.K('aggregate-using-all-hits'))
    groupB.add_argument('--min-occurrence', metavar="NUM GENOMES", default=1, help=("The minimum number of occurrence of any "
                                "given function accross genomes. If you set a value, those functions that occur in less number "
                                "of genomes will be excluded."), type=int)

    groupC = parser.add_argument_group('GROUPS', "How should anvi'o divide your genomes into groups?")
    groupC.add_argument(*anvio.A('groups-txt'), **anvio.K('groups-txt'))

    groupD = parser.add_argument_group('OUTPUT', "A.k.a., what you're really here for")
    groupD.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
