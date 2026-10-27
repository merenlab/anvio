#!/usr/bin/env python

"""Migration script to add the number of reads mapped tables, and to move remaining zero-value rows
from view tables into zero_coverage tables.

This migration:
  1. Creates the (empty) tables for the number of reads mapped per split and contig, and sets the
     `num_reads_mapped_stored` flag to False, since this information cannot be recovered from the
     existing data.
  2. Moves (item, layer) pairs that have value=0 across all view tables into the zero_coverage tables.
     Profile databases generated since v42 store explicit zero rows for contigs that were profiled
     but ended up with zero coverage (i.e., none of their reads survived filtering). Such contigs
     are now recorded in the zero_coverage tables, just like those skipped by the BAM index prefilter.

IMPORTANT: the second step is destructive -- it modifies view tables in place. That part is wrapped
in a single SQLite transaction so that a crash at any point leaves the view tables unchanged.
"""

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.tables as t
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

essential_data_fields = ['std_coverage',
                         'mean_coverage',
                         'mean_coverage_Q2Q3',
                         'detection',
                         'abundance',
                         'variability']

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    profile_db_info = dbinfo.ProfileDBInfo(db_path)
    if str(profile_db_info.version) != current_version:
        raise ConfigError(f"The version of the provided profile database is {profile_db_info.version}, not the "
                          f"required version, {current_version}, so this script cannot upgrade the database.")

    progress.new("Migrating")
    profile_db = profile_db_info.load_db()

    progress.update("Creating tables for the number of reads mapped per split and contig")
    existing_tables = profile_db.get_table_names()
    for table_name, structure, types in [(t.num_reads_assigned_splits_table_name, t.view_table_structure, t.view_table_types),
                                         (t.num_reads_mapped_contigs_table_name, t.view_table_structure, t.view_table_types),
                                         (t.zero_num_reads_assigned_splits_table_name, t.zero_num_reads_assigned_splits_table_structure, t.zero_num_reads_assigned_splits_table_types)]:
        if table_name not in existing_tables:
            profile_db.create_table(table_name, structure, types)

    # the number of reads cannot be recovered from the existing data, so these tables stay empty
    profile_db.remove_meta_key_value_pair('num_reads_mapped_stored')
    profile_db.set_meta_value('num_reads_mapped_stored', False)

    existing_tables = profile_db.get_table_names()

    # When SNVs were not profiled, the variability tables may be empty or contain None values, and
    # should be excluded from both validation and zero-row deletion.
    self_table = profile_db.get_table_as_dict('self')
    SNVs_profiled = int(self_table['SNVs_profiled']['value']) if 'SNVs_profiled' in self_table else False

    # We accumulate all SQL operations and execute them in a single transaction.
    queued_operations = []
    num_pairs_moved = {}

    for target in ['splits', 'contigs']:
        zero_cov_table = f'zero_coverage_{target}'
        progress.update(f"Processing {target} view tables")

        table_names = []
        for field in essential_data_fields:
            if field == 'variability' and not SNVs_profiled:
                continue
            table_name = f'{field}_{target}'
            if table_name in existing_tables:
                table_names.append(table_name)

        # Use detection as the reference field to identify zero-coverage items. If detection = 0, all
        # other fields must also be 0.
        detection_table = f'detection_{target}'
        if detection_table not in existing_tables or zero_cov_table not in existing_tables:
            continue

        progress.update(f"Reading {detection_table} for zero-detection pairs")
        zero_pairs = set()
        for entry_id, item, layer, value in profile_db.get_all_rows_from_table(detection_table):
            if value == 0 or value == 0.0:
                zero_pairs.add((item, layer))

        if not zero_pairs:
            continue

        # Validate consistency across all other view tables. Pairs that have nonzero values in any
        # other field, or that are missing rows in some view tables, stay in the view tables as-is.
        progress.update(f"Validating {len(zero_pairs)} zero-coverage {target} across all fields")
        inconsistent_pairs = set()
        for table_name in table_names:
            if table_name == detection_table:
                continue

            zero_pairs_found = set()
            for entry_id, item, layer, value in profile_db.get_all_rows_from_table(table_name):
                pair = (item, layer)
                if pair in zero_pairs:
                    zero_pairs_found.add(pair)
                    if value != 0 and value != 0.0:
                        inconsistent_pairs.add(pair)

            inconsistent_pairs.update(zero_pairs - inconsistent_pairs - zero_pairs_found)

        if inconsistent_pairs:
            run.warning(f"This profile database had {len(inconsistent_pairs)} zero-detection {target[:-1]} pair(s) "
                        f"with non-zero values in other fields, or with missing rows in some view tables. These pairs "
                        f"will be left in the view tables as-is and will NOT be moved into the zero_coverage table. "
                        f"The migration will proceed normally.")
            zero_pairs -= inconsistent_pairs

        # pairs that are somehow already listed in the zero_coverage table only need their rows deleted
        already_listed = set(profile_db.get_all_rows_from_table(zero_cov_table))
        pairs_to_insert = list(zero_pairs - already_listed)

        if not zero_pairs:
            continue

        if pairs_to_insert:
            queued_operations.append((f'INSERT INTO "{zero_cov_table}" (item, layer) VALUES (?, ?)', pairs_to_insert))

        # Queue deletions from all view tables for these pairs. We use a temp table with an index
        # to avoid O(n*m) full scans.
        temp_table = f'_zero_pairs_{target}_TEMP'
        queued_operations.append((f'CREATE TEMP TABLE "{temp_table}" (item text, layer text)', None))
        queued_operations.append((f'INSERT INTO "{temp_table}" (item, layer) VALUES (?, ?)', list(zero_pairs)))
        queued_operations.append((f'CREATE INDEX "_idx_{temp_table}" ON "{temp_table}" (item, layer)', None))

        for table_name in table_names:
            queued_operations.append((f'DELETE FROM "{table_name}" WHERE EXISTS '
                                      f'(SELECT 1 FROM "{temp_table}" z WHERE z.item = "{table_name}".item AND z.layer = "{table_name}".layer)',
                                      None))

        queued_operations.append((f'DROP TABLE IF EXISTS "{temp_table}"', None))
        num_pairs_moved[target] = len(zero_pairs)

    # Execute everything in a single transaction
    progress.update("Applying changes in a single transaction")

    conn = profile_db.conn
    conn.commit()
    cursor = conn.cursor()

    try:
        cursor.execute('BEGIN IMMEDIATE')

        for sql, params in queued_operations:
            if params is not None:
                cursor.executemany(sql, params)
            else:
                cursor.execute(sql)

        # Update version inside the same transaction
        cursor.execute("DELETE FROM self WHERE key = 'version'")
        cursor.execute("INSERT INTO self VALUES ('version', ?)", (next_version,))

        conn.commit()
    except Exception:
        conn.rollback()
        progress.end()
        raise

    if num_pairs_moved:
        progress.update("Reclaiming disk space")
        profile_db._exec('VACUUM')

    progress.update("Done")
    profile_db.disconnect()

    progress.end()

    moved_msg = ""
    if num_pairs_moved:
        moved_msg = (f" In addition, {', '.join([f'{n} zero-coverage {target}' for target, n in num_pairs_moved.items()])} "
                     f"(item/sample pairs) were moved from the view tables into the zero_coverage tables. These pairs were "
                     f"left behind during the previous migration step because they did have reads mapping to them, yet all "
                     f"those reads were filtered out so effectively their coverage (and other stats) are still zeros.")

    run.info_single(f"Congratulations! Your profile database is now version {next_version}. This update added new tables "
                    f"to store the number of reads that map to each split and contig, which is needed by some downstream "
                    f"analyses such as the calculation of distribution of coverage (discov) scores. Please note that these tables are EMPTY "
                    f"in your database, since this information cannot be recovered from the existing data. If you need it, "
                    f"you will have to re-profile your BAM files with `anvi-profile`.{moved_msg}", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"A simple script to upgrade an anvi'o profile database from version {current_version} to version {next_version}")
    parser.add_argument("profile_db", metavar="PROFILE_DB", help=f"An anvi'o profile database of version {current_version}")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
