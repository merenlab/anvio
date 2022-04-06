#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


current_version, next_version = [x[1:] for x in __name__.split('_to_')]

kegg_brite_info_table_name      = 'kegg_brite_info'
kegg_brite_info_table_structure = ['hierarchy_accession', 'max_depth', 'max_depth_excluding_subcategories']
kegg_brite_info_table_types     = [        'text'       ,  'numeric' ,             'numeric'              ]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)

    if str(contigs_db.get_version()) != current_version:
        contigs_db.disconnect()
        raise ConfigError(f"Version of this contigs database is not {current_version} (hence, this script cannot really do anything).")

    progress.new("Migrating")
    progress.update("Creating a new table to support manipulation of KEGG BRITE categorizations")

    # just to be on the safe side.
    try:
        contigs_db.drop_table(kegg_brite_info_table_name)
    except:
        pass

    contigs_db.create_table(kegg_brite_info_table_name, kegg_brite_info_table_structure, kegg_brite_info_table_types)

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    contigs_db.disconnect()
    progress.end()

    run.info_single(f"The contigs database is now {next_version}. This upgrade added an empty table of the "
                    "maximum depths of KEGG BRITE hierarchies used in KEGG gene annotation -- pretty obscure. "
                    "But this should draw your attention to the fact that `anvi-run-kegg-kofams` will now categorize "
                    "your genes in the useful BRITE system when run with an up-to-date modules database. "
                    "We suggest running `anvi-setup-kegg-kofams` (with the `-D` flag in anvi'o v7 or earlier), "
                    "and running (or rerunning) `anvi-run-kegg-kofams` on your contigs database to benefit "
                    "from BRITE annotations :)",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade contigs database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
