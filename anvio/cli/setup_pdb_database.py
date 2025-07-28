#!/usr/bin/env python
# -*- coding: utf-8

import sys
import anvio
from anvio.argparse import ArgumentParser

import anvio.terminal as terminal
import anvio.structureops as structops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__provides__ = ['pdb-db']
__description__ = ("Setup or update an offline database of representative PDB structures clustered at 95%")


@terminal.time_program
def main():
    args = get_args()

    try:
        structops.PDBDatabase(args).process()
    except ConfigError as e:
        print(e)
        sys.exit(1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(2)


def get_args():
    parser = ArgumentParser(description=__description__)
    parser.add_argument(*anvio.A('pdb-database-path'), **anvio.K('pdb-database-path'))
    parser.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    parser.add_argument('--update', action='store_true',
                        help="Use this flag if you would like to update your current database.")
    parser.add_argument('--skip-modeller-update', action='store_true',
                        help="By default, MODELLER's search DB is updated when this program is ran "
                             "so that if MODELLER finds a protein, its structure is guaranteed "
                             "to be in this database. If you don't want to touch the MODELLER "
                             "database, use this flag.")
    parser.add_argument(*anvio.A('reset'), **anvio.K('reset'))

    return parser.get_args(parser)

if __name__ == '__main__':
    main()
