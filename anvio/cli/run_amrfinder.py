#!/usr/bin/env python

import sys

import anvio

from anvio.amrfinder import AMRFinder
from anvio.terminal import time_program
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']
__requires__ = ['contigs-db', 'amrfinder-data']
__provides__ = ['functions']
__description__ = ("You can use this to annotate genes in a contigs-db with NCBI's AMRFinderPlus to find antimicrobial "
                   "resistance genes, stress response genes (such as metal and biocide resistance), and virulence genes.")


@time_program
def main():
    args = get_args()

    try:
        amrfinder = AMRFinder(args)
        amrfinder.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('amrfinder-data-dir'), **anvio.K('amrfinder-data-dir'))
    parser.add_argument(*anvio.A('full-amrfinder-report'), **anvio.K('full-amrfinder-report'))
    parser.add_argument(*anvio.A('skip-e-values'), **anvio.K('skip-e-values'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it', {'help': "Using this flag, anvi'o will overwrite "
                        "previous AMRFinderPlus annotations in your contigs db."}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
