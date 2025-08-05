#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_kegg_modules_db


current_version, next_version = [x[1:] for x in __name__.split('_to_')]

brite_table_name      = "brite_hierarchies"
brite_table_structure = ['hierarchy_accession', 'hierarchy_name', 'ortholog_accession', 'ortholog_name', 'categorization']
brite_table_types     = [        'str'        ,       'str'     ,         'str'       ,      'str'     ,      'str'      ]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    is_kegg_modules_db(db_path)

    modules_db = db.DB(db_path, None, ignore_version = True)

    if str(modules_db.get_version()) != current_version:
        modules_db.disconnect()
        raise ConfigError(f"Version of this modules database is not {current_version} (hence, this script cannot really do anything).")

    progress.new("Migrating")
    progress.update("...")

    # just to be on the safe side.
    try:
        modules_db.drop_table(brite_table_name)
    except:
        pass

    try:
        modules_db.remove_meta_key_value_pair('num_brite_hierarchies')
        modules_db.remove_meta_key_value_pair('total_brite_entries')
        modules_db.remove_meta_key_value_pair('is_brite_setup')
    except:
        pass

    modules_db.set_meta_value('num_brite_hierarchies', None)
    modules_db.set_meta_value('total_brite_hierarchies', None)
    modules_db.set_meta_value('is_brite_setup', 0)

    progress.update("Creating a new table for KEGG BRITE categorizations")

    modules_db.create_table(brite_table_name, brite_table_structure, brite_table_types)

    progress.update("Updating version")
    modules_db.remove_meta_key_value_pair('version')
    modules_db.set_version(next_version)

    modules_db.disconnect()
    progress.end()

    run.info_single(f"The modules database is now {next_version}. An empty table of KEGG BRITE hierarchy "
                    "categorizations of all orthologs was created, and related self table attributes were added. "
                    "We suggest generating a new modules database from `anvi-setup-kegg-data` and re-running "
                    "`anvi-run-kegg-kofams` on your contigs databases to benefit from these useful annotations :)",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade KEGG modules database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('modules_db', metavar = 'MODULES_DB', help = 'Modules database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.modules_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
