#!/usr/bin/env python

"""Migration script to add the modifications table to profile databases."""

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal
import anvio.tables as t

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    profile_db_info = dbinfo.ProfileDBInfo(db_path)
    if str(profile_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided profile database is {profile_db_info.version}, not the "
            f"required version, {current_version}, so this script cannot upgrade the database."
        )

    progress.new("Migrating")
    profile_db = profile_db_info.load_db()

    existing_tables = profile_db.get_table_names()

    if t.modifications_table_name not in existing_tables:
        progress.update("Creating modifications table")
        profile_db.create_table(t.modifications_table_name, t.modifications_table_structure, t.modifications_table_types)

    progress.update("Updating version")
    profile_db._exec("DELETE FROM self WHERE key = 'version'")
    profile_db._exec("INSERT INTO self VALUES ('version', ?)", (next_version,))

    profile_db.disconnect()
    progress.end()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('db_path', help='Path to a profile database to migrate')
    args = parser.parse_args()

    try:
        migrate(args.db_path)
    except ConfigError as e:
        run.error(e)
        sys.exit(-1)
