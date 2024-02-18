#!/usr/bin/env python
# -*- coding: utf-8

import sys
import sqlite3
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    contigs_db_info = dbinfo.ContigsDBInfo(db_path)
    if str(contigs_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided contigs database is {contigs_db_info.version}, "
            f"not the required version, {current_version}, so this script cannot upgrade the database.")

    progress.new("Migrating")
    progress.update("Removing old SCGs...")

    self_table = contigs_db_info.get_self_table()
    contigs_db = contigs_db_info.load_db()
    
    if self_table['scg_taxonomy_was_run']:
        with sqlite3.connect(db_path) as connection:
            cursor = connection.cursor()

            # Construct the parameterized query for SCG removal
            placeholders = ', '.join(['?'] * len(constants.default_scgs_for_taxonomy))
            where_clause = f"gene_name NOT IN ({placeholders})"

            # Execute the DELETE statement to remove old SCGs
            cursor.execute(f"DELETE FROM scg_taxonomy WHERE {where_clause}", constants.default_scgs_for_taxonomy)

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()

    message = (
        "Congratulations! Your contigs database is now version 23. This update removed the old SCGs "
        "that anvi-estimate-scg-taxonomy used to use, and updated the version of the database."
    )

    run.info_single(message, nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"A simple script to upgrade an anvi'o contigs database from version {current_version} to version {next_version}")
    parser.add_argument("contigs_db", metavar="CONTIGS_DB", help=f"An anvi'o contigs database of version {current_version}")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)