#!/usr/bin/env python

"""Migration script to remove zero-value rows from view tables into zero_coverage tables.

Previously, the profile database stored explicit zero-value rows in all 12 view tables
(6 fields x splits/contigs) for every contig/split that had no coverage. This is extremely
wasteful when mapping to large reference collections where most items have no coverage.

This migration:
  1. Creates two new tables: zero_coverage_splits and zero_coverage_contigs
  2. Identifies all (item, layer) pairs that have value=0 across ALL view tables
  3. Validates consistency: if any field is 0 for a pair, all fields must be 0
  4. Moves those pairs into the zero_coverage tables
  5. Deletes the zero-value rows from all view tables
  6. Reclaims disk space via VACUUM

IMPORTANT: this migration is destructive — it modifies view tables in place.
The entire operation is wrapped in a single SQLite transaction so that a crash
at any point leaves the database unchanged.
"""

import sys
import argparse

import anvio.dbinfo as dbinfo
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
        raise ConfigError(
            f"The version of the provided profile database is {profile_db_info.version}, not the "
            f"required version, {current_version}, so this script cannot upgrade the database."
        )

    progress.new("Migrating")
    profile_db = profile_db_info.load_db()

    existing_tables = profile_db.get_table_names()

    # Check whether SNVs were profiled. When they were not, the variability tables
    # may be empty or contain None values, and should be excluded from both
    # validation and zero-row deletion.
    self_table = profile_db.get_table_as_dict('self')
    SNVs_profiled = int(self_table['SNVs_profiled']['value']) if 'SNVs_profiled' in self_table else False

    # We accumulate all SQL operations and execute them in a single transaction.
    queued_operations = []

    # The zero_coverage tables use the same column name 'layer' as the view tables,
    # matching the naming convention used throughout the profile database.
    queued_operations.append(('CREATE TABLE IF NOT EXISTS "zero_coverage_splits" (item text, layer text)', None))
    queued_operations.append(('CREATE TABLE IF NOT EXISTS "zero_coverage_contigs" (item text, layer text)', None))

    for target in ['splits', 'contigs']:
        zero_cov_table = f'zero_coverage_{target}'
        progress.update(f"Processing {target} view tables")

        # Collect all view table names for this target that exist, skipping
        # variability tables when SNVs were not profiled (they may be empty
        # or contain None values that should not be treated as zero-coverage).
        table_names = []
        for field in essential_data_fields:
            if field == 'variability' and not SNVs_profiled:
                continue
            table_name = f'{field}_{target}'
            if table_name in existing_tables:
                table_names.append(table_name)

        if not table_names:
            continue

        # Use detection as the reference field to identify zero-coverage items.
        # If detection = 0, all other fields must also be 0.
        detection_table = f'detection_{target}'
        if detection_table not in existing_tables:
            continue

        # Step 1: Read the detection table to find zero-detection (item, layer) pairs
        progress.update(f"Reading {detection_table} for zero-detection pairs")
        detection_rows = profile_db.get_all_rows_from_table(detection_table)

        # Build set of zero-detection (item, layer) pairs
        zero_pairs = set()
        for entry_id, item, layer, value in detection_rows:
            if value == 0 or value == 0.0:
                zero_pairs.add((item, layer))

        if not zero_pairs:
            continue

        progress.update(f"Validating {len(zero_pairs)} zero-coverage {target} across all fields")

        # Step 2: Validate consistency across all other view tables.
        # For each zero-detection (item, layer) pair, every other field should also be 0.
        # Extremely old databases may violate this invariant, but we keep them around for
        # historical reasons (also see migration script from 40 to 41). So here we collect
        # inconsistent pairs and exclude them from the zero_coverage movement rather than
        # aborting. They stay in the view tables as-is.
        inconsistent_pairs = set()

        for field in essential_data_fields:
            if field == 'detection':
                continue
            if field == 'variability' and not SNVs_profiled:
                continue

            table_name = f'{field}_{target}'
            if table_name not in existing_tables:
                continue

            progress.update(f"Validating {table_name}")
            rows = profile_db.get_all_rows_from_table(table_name)

            # Track which zero-detection pairs we find in this table
            zero_pairs_found = set()

            for entry_id, item, layer, value in rows:
                pair = (item, layer)
                if pair in zero_pairs:
                    zero_pairs_found.add(pair)
                    if value != 0 and value != 0.0:
                        inconsistent_pairs.add(pair)

            # Check that all zero-detection pairs have a row in this table.
            # Missing rows mean the data is incomplete; exclude those pairs too.
            truly_zero = zero_pairs - inconsistent_pairs
            missing_pairs = truly_zero - zero_pairs_found
            if missing_pairs:
                inconsistent_pairs.update(missing_pairs)

        if inconsistent_pairs:
            run.warning(f"This profile database had a historical data inconsistency: {len(inconsistent_pairs)} "
                        f"zero-detection {target[:-1]} pair(s) had non-zero values in other fields, or were "
                        f"missing rows in some view tables. This is a known bug from very early versions of "
                        f"anvi'o. These pairs will be left in the view tables as-is and will NOT be moved "
                        f"into the zero_coverage table. The migration will proceed normally.")
            zero_pairs -= inconsistent_pairs

        if not zero_pairs:
            continue

        # Step 3: Queue inserts into the zero_coverage table
        zero_cov_rows = list(zero_pairs)
        queued_operations.append((
            f'INSERT INTO "{zero_cov_table}" (item, layer) VALUES (?, ?)',
            zero_cov_rows
        ))

        # Step 4: Queue deletions from all view tables for these pairs.
        # We use a temp table with an index to avoid O(n*m) full scans.
        temp_table = f'_zero_pairs_{target}_TEMP'
        queued_operations.append((f'CREATE TEMP TABLE "{temp_table}" (item text, layer text)', None))
        queued_operations.append((
            f'INSERT INTO "{temp_table}" (item, layer) VALUES (?, ?)',
            zero_cov_rows
        ))
        queued_operations.append((f'CREATE INDEX "_idx_{temp_table}" ON "{temp_table}" (item, layer)', None))

        for table_name in table_names:
            queued_operations.append((
                f'DELETE FROM "{table_name}" WHERE EXISTS '
                f'(SELECT 1 FROM "{temp_table}" z WHERE z.item = "{table_name}".item AND z.layer = "{table_name}".layer)',
                None
            ))

        queued_operations.append((f'DROP TABLE IF EXISTS "{temp_table}"', None))

        progress.update(f"Queued removal of {len(zero_pairs)} zero-coverage {target}")

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

    progress.update("Reclaiming disk space")
    profile_db._exec('VACUUM')

    progress.update("Done")
    profile_db.disconnect()

    progress.end()

    run.info_single(
        f"Congratulations! Your profile database is now version {next_version}. Zero-value "
        f"rows have been moved from the view tables into the new zero_coverage_splits and "
        f"zero_coverage_contigs tables, significantly reducing database size for sparse "
        f"coverage profiles.",
        nl_after=1, nl_before=1, mc='green'
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=
        f"A simple script to upgrade an anvi'o profile database from version "
        f"{current_version} to version {next_version}"
    )
    parser.add_argument(
        "profile_db",
        metavar="PROFILE_DB",
        help=f"An anvi'o profile database of version {current_version}"
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
