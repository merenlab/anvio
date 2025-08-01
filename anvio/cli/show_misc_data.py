#!/usr/bin/env python
# -*- coding: utf-8
"""List misc data keys in pan or profile databases"""

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
__requires__ = ["pan-db", "profile-db", "contigs-db"]
__description__ = "Show all misc data keys in all misc data tables"


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def run_program():
    args = get_args()

    d = args.__dict__
    if not d['contigs_db'] and not d['pan_or_profile_db']:
        raise ConfigError("Please provide either a contigs database (--contigs-db) or a profile/pan "
                          "database (--pan-or-profile-db)")

    if d['contigs_db'] and d['pan_or_profile_db']:
        raise ConfigError("You provided a contigs db (--contigs-db) and a profile/pan "
                          "db (--pan-or-profile-db). Please provide only one.")

    if args.target_data_table:
        table_for_additional_data = MiscDataTableFactory(args)
        table_for_additional_data.list_data_keys()
    else:
        if args.pan_or_profile_db:
            table_names = ['items', 'layers', 'layer_orders']
        else:
            table_names = ['nucleotides', 'amino_acids']

        for table_name in table_names:
            args.target_data_table = table_name
            table_for_additional_data = MiscDataTableFactory(args)
            table_for_additional_data.list_data_keys()


def get_args():
    parser = ArgumentParser(description=__description__)

    group1 = parser.add_argument_group('Database input', "Provide 1 of these")
    group2 = parser.add_argument_group('Details', "Everything else.")

    group1.add_argument(*anvio.A('pan-or-profile-db'), **anvio.K('pan-or-profile-db', {'required': False}))
    group1.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db', {'required': False}))

    group2.add_argument(*anvio.A('target-data-table'), **anvio.K('target-data-table'))
    group2.add_argument(*anvio.A('target-data-group'), **anvio.K('target-data-group'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
