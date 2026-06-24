#!/usr/bin/env python
# -*- coding: utf-8
"""Compare two pan databases made from the same genomes storage"""

import sys

import anvio
import anvio.panops as panops
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError, HDF5Error


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['pan-db', 'genomes-storage-db']
__provides__ = ['misc-data-items']
__description__ = ("Compare two pan databases made from the same genomes storage to identify "
                   "fragmentation, combination, and compositional differences between gene clusters")
__resources__ = []


run = terminal.Run()
progress = terminal.Progress()


def main():
    args = get_args()

    try:
        comparer = panops.ComparePan(args, run, progress)
        comparer.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)
    except HDF5Error as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT', "The two pan databases to compare and the genomes storage they share.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('compared-pan-db'), **anvio.K('compared-pan-db', {'required': True}))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': True}))

    groupB = parser.add_argument_group('OUTPUT', "Where to store the comparison results.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupB.add_argument('--skip-output-files', default=False, action='store_true',
                        help="Skip writing comparison results to a text file.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
