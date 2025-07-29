#!/usr/bin/env python
# -*- coding: utf-8
"""Adds a DEAFULT collection with EVERYTHING in it into a pan or profile database."""

import sys

import anvio

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.collections import TablesForCollections


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__provides__ = ["collection", "bin"]
__requires__ = ["pan-db", "profile-db", "contigs-db"]
__description__ = ("A script to add a 'DEFAULT' collection in an anvi'o pan or profile database with "
                   "either (1) a single bin that describes all items available in the profile database, "
                   "or (2) as many bins as there are items in the profile database wher every item has its "
                   "own bin. The former is the default behavior that will be useful in most instances "
                   "where you need to use this script. The latter is most useful if you are Florian and/or "
                   "have something very specific in mind.")


def main():
    args = get_args()

    try:
        TablesForCollections(args.pan_or_profile_db).add_default_collection_to_db(contigs_db_path=args.contigs_db,
                                                                                  collection_name=args.collection_name,
                                                                                  bin_name=args.bin_id,
                                                                                  bin_each_item_separately=args.bin_each_item_separately)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('DATABASE INPUTS', "Provide a profile-, gene-, or pan-db, and if relevant, "
                                       "also a contigs-db here")
    groupA.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db'))
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    groupB = parser.add_argument_group('COLLECTION NAME', "Go with the flow or be a rebel")
    groupB.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name', {'required': False,
                                                               'help': "Name for the new collection. If you don't provide any then \
                                                                        it will be named \"DEFAULT\".",
                                                               'default': "DEFAULT"}))

    groupC = parser.add_argument_group('BINS', "You either provide a bin name for everything, or let every "
                                        "item be their own bins.")
    groupC.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id', {'required': False, 'help': "Name for the new bin. "
                "If you don't provide any then it will be named \"EVERYTHING\".", 'default': "EVERYTHING"}))
    groupC.add_argument('--bin-each-item-separately', default=False, action="store_true", help="If you declare this flag "
                "every item in your database will have its own bin (which is crazy, but may be necessary).")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
