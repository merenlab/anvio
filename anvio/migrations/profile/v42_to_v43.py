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
    else:
        progress.update("Updating modifications table")
        # If the existing table does not have the desired columns (count and strand),
        # aggregate the old per-read rows into the new aggregated table.
        columns = [row[1] for row in profile_db._exec('PRAGMA TABLE_INFO(%s)' % t.modifications_table_name)]
        if 'count' not in columns or 'strand' not in columns:
            progress.update("Aggregating existing modifications into new schema")
            profile_db._exec('ALTER TABLE modifications RENAME TO modifications_old')
            profile_db.create_table(t.modifications_table_name, t.modifications_table_structure, t.modifications_table_types)
            # If old table had a strand column, preserve it; otherwise default to '.'
            if 'strand' in columns:
                profile_db._exec(
                    'INSERT INTO modifications (sample_id, split_name, pos, pos_in_contig, modification, strand, count) '
                    'SELECT sample_id, split_name, pos, pos_in_contig, lower(modification), strand, COUNT(*) '
                    'FROM modifications_old GROUP BY sample_id, split_name, pos, pos_in_contig, lower(modification), strand')
            else:
                profile_db._exec(
                    "INSERT INTO modifications (sample_id, split_name, pos, pos_in_contig, modification, strand, count) "
                    "SELECT sample_id, split_name, pos, pos_in_contig, lower(modification), '.' as strand, COUNT(*) "
                    "FROM modifications_old GROUP BY sample_id, split_name, pos, pos_in_contig, lower(modification)")
            profile_db._exec('DROP TABLE modifications_old')

    # Add new meta keys for modifications support with safe defaults
    progress.update("Setting modifications metadata")
    
    # Check existing keys by querying the self table
    existing_keys_rows = profile_db._exec("SELECT key FROM self WHERE key IN ('modifications_profiled', 'min_coverage_for_modifications', 'modification_filters')").fetchall()
    existing_keys = set([row[0] for row in existing_keys_rows])
    
    if 'modifications_profiled' not in existing_keys:
        profile_db.set_meta_value('modifications_profiled', 0)
    
    if 'min_coverage_for_modifications' not in existing_keys:
        profile_db.set_meta_value('min_coverage_for_modifications', 0)
    
    if 'modification_filters' not in existing_keys:
        profile_db.set_meta_value('modification_filters', '{"filters": {}, "default": 0.0}')

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
        run.warning(str(e))
        sys.exit(-1)
