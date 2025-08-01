#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio

from anvio.cogs import COGsSetup
from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ['cogs-data']
__description__ = "Download and setup NCBI's Clusters of Orthologous Groups database"


def main():
    args = get_args()

    try:
        setup = COGsSetup(args)
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

    parser.add_argument(*anvio.A('cog-version'), **anvio.K('cog-version'))
    parser.add_argument('--cog-data-dir', default=None, type=str, help="The directory for COG data to be stored. If you leave it\
                        as is without specifying anything, the default destination for the data directory will be used to set things\
                        up. The advantage of it is that everyone will be using a single data directory, but then you may need\
                        superuser privileges to do it. Using this parameter you can choose the location of the data directory somewhere\
                        you like. However, when it is time to run COGs, you will need to remember that path and provide it to the program.")
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset'))
    parser.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
