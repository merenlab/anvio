#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

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
    progress.update("Adding any missing database metavariables")

    self_table = contigs_db_info.get_self_table()
    contigs_db = contigs_db_info.load_db()
    for key in [
        'reaction_network_ko_annotations_hash',
        'reaction_network_kegg_database_release',
        'reaction_network_modelseed_database_sha'
    ]:
        if key not in self_table:
            contigs_db.set_meta_value(key, None)

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()

    message = (
        "Congratulations! Your contigs database is now version 22. This update fixes a bug that prevented "
        "construction and storage of a metabolic reaction network in version 21 databases that were created "
        "from scratch and not migrated from version 20. You can blame the author of this script for the "
        "oversight responsible for the bug and this annoying migration. If the version 21 database was "
        "created from scratch, then placeholders for three 'metavariables' should now have been added to the "
        "database; if the version 21 database migrated from version 20, then nothing at all should have "
        "changed in the database except the version number."
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
