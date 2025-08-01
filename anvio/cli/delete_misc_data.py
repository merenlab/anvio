#!/usr/bin/env python
# -*- coding: utf-8
"""Remove stuff from misc data tables"""

import sys
from anvio.argparse import ArgumentParser

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import MiscDataTableFactory


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren', 'ekiefl']
__requires__ = ["pan-db", "profile-db", "misc-data-items", "misc-data-layers", "misc-data-layer-orders", "misc-data-nucleotides", "misc-data-amino-acids"]
__description__ = ("Remove stuff from 'additional data' or 'order' tables for either items or layers in either "
                   "pan or profile databases. OR, remove stuff from the 'additional data' tables for nucleotides "
                   "or amino acids in contigs databases")
__resources__ = [("Working with anvi'o additional data tables", "http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology")]


@terminal.time_program
def main():
    args = get_args()

    A = lambda x: args.__dict__[x] if x in args.__dict__ else None

    try:
        if not A('contigs_db') and not A('pan_or_profile_db'):
            raise ConfigError("Please provide either a contigs database (--contigs-db) or a profile/pan "
                              "database (--pan-or-profile-db)")

        if A('contigs_db') and A('pan_or_profile_db'):
            raise ConfigError("You provided a contigs db (--contigs-db) and a profile/pan "
                              "db (--pan-or-profile-db). Please provide only one.")

        table_for_additional_data = MiscDataTableFactory(args)

        if args.list_available_keys:
            table_for_additional_data.list_data_keys()
            sys.exit(0)

        keys_to_remove = [k for k in args.keys_to_remove.split(',')] if args.keys_to_remove else []
        groups_to_remove = [k for k in args.groups_to_remove.split(',')] if args.groups_to_remove else []

        table_for_additional_data.remove(data_keys_list=keys_to_remove, data_groups_list=groups_to_remove)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('Database input', "Provide 1 of these")
    group2 = parser.add_argument_group('Details', "Everything else.")

    group1.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db', {'required': False}))
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    group2.add_argument(*anvio.A('target-data-table'), **anvio.K('target-data-table', {'required': True}))
    group2.add_argument('--keys-to-remove', help="A comma-separated list of data keys to remove from the database. If you\
                                                  do not use this parameter, anvi'o will simply remove everything from the\
                                                  target data table immediately. Please note that you should not use this\
                                                  parameter together with `--groups-to-remove` in a single command.")
    group2.add_argument('--groups-to-remove', help="A comma-separated list of data groups to remove from the database. If you\
                                                  do not use this parameter, anvi'o will simply remove everything from the\
                                                  target data table immediately. Please note that you should not use this\
                                                  parameter together with `--keys-to-remove` in a single command.")
    group2.add_argument('--list-available-keys', action='store_true', help="Using this flag will list available data keys in\
                                                  the target data table and quit without doing anything else.")
    group2.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
