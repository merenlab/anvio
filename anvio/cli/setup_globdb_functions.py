#!/usr/bin/env python

import sys

import anvio

from anvio.globdb import GlobDBFunctionsSetup
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'dspeth']
__provides__ = ['globdb-data']
__description__ = "Download and set up the GlobDB gene family database for functional annotation"


def main():
    args = get_args()

    try:
        setup = GlobDBFunctionsSetup(args)
        setup.create()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    parser.add_argument(*anvio.A('globdb-data-dir'), **anvio.K('globdb-data-dir'))
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
