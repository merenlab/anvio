#!/usr/bin/env python

import sys

import anvio

from anvio.amrfinder import AMRFinderSetup
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2026, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['jessika-fuessel']
__provides__ = ['amrfinder-data']
__description__ = "Downloads and sets up NCBI's AMRFinderPlus database for `anvi-run-amrfinder`"


def main():
    args = get_args()

    try:
        setup = AMRFinderSetup(args)
        setup.setup()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('amrfinder-data-dir'), **anvio.K('amrfinder-data-dir'))
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset', {'help': "Use this flag if you want to re-download the "
                        "AMRFinderPlus database if you already have a version on your computer."}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
