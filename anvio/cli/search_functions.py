#!/usr/bin/env python
"""Searches contigs database for a given function"""

import sys

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.functionops import FunctionsSearchResultReporter


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ivagljiva']
__requires__ = ['contigs-db']
__can_use__ = ['genomes-storage-db']
__provides__ = ['functions-txt']
__description__ = ("Search functions in an anvi'o contigs database or genomes storage. Basically, this program "
                   "searches for one or more search terms you define in functional annotations of "
                   "genes in an anvi'o contigs database, and generates multiple reports. The "
                   "default report simply tells you which contigs contain genes with functions "
                   "matching to serach terms you used, useful for viewing in the interface. You "
                   "can also request a much more comprehensive report, which gives you anything "
                   "you might need to know for each hit and serach term")


def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        FunctionsSearchResultReporter(args, run=run, progress=progress).process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('SEARCH IN', 'Relevant source databases')
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db', {'required': False}))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))
    groupA.add_argument(*anvio.A('external-genomes'), **anvio.K('external-genomes', {'required': False}))

    groupB = parser.add_argument_group('SEARCH FOR', 'Relevant terms')
    groupB.add_argument(*anvio.A('search-terms'), **anvio.K('search-terms', {'required': True}))
    groupB.add_argument(*anvio.A('case-sensitive'), **anvio.K('case-sensitive'))
    groupB.add_argument(*anvio.A('exact-match'), **anvio.K('exact-match'))
    groupB.add_argument(*anvio.A('delimiter'), **anvio.K('delimiter'))
    groupB.add_argument(*anvio.A('annotation-sources'), **anvio.K('annotation-sources'))
    groupB.add_argument(*anvio.A('list-annotation-sources'), **anvio.K('list-annotation-sources'))

    groupC = parser.add_argument_group('REPORT', "Anvi'o can report the hits in multiple ways. The output file will be a very simple 2-column\
                                                  TAB-delimited output that is compatible with anvi'o additional data format (so you can give\
                                                  it to the `anvi-interactive` to see which splits contained genes that were matching to your\
                                                  search terms). You can also ask anvi'o to generate a full-report, that contains much more and\
                                                  much helpful information about each hit. Optionally you can even ask the gene sequences to\
                                                  appear in this report.")
    groupC.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupC.add_argument(*anvio.A('full-report'), **anvio.K('full-report'))
    groupC.add_argument(*anvio.A('include-sequences'), **anvio.K('include-sequences'))
    groupC.add_argument(*anvio.A('verbose'), **anvio.K('verbose'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
