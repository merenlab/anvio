#!/usr/bin/env python
# -*- coding: utf-8
"""A script to upgrade the description in an anvi'o db"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.dbops as dbops

from anvio.errors import ConfigError, FilesNPathsError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ["pan-db", "profile-db", "contigs-db", "genomes-storage-db"]
__description__ = "Update the description in an anvi'o database"


def main():
    args = get_args()

    try:
        dbops.update_description_in_db_from_file(args.db_path, args.description)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('db_path', metavar='DB', help="An anvi'o database.")
    parser.add_argument(*anvio.A('description'), **anvio.K('description', {'required': True}))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
