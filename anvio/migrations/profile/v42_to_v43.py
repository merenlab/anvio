#!/usr/bin/env python

"""Migration script that adds an index on variable_nucleotides(split_name).

The variable_nucleotides table is registered in tables/__init__.py with the
'indexed' flag set to True, but no column-specific index was ever created.
For typical anvi'o profile databases (low millions of SNV rows) this has not
mattered. For large read-recruitment workflows (hundreds of millions of SNV
rows in a single merged profile db) any query of the form
'WHERE split_name IN (...)' degenerates into a full table scan, which makes
per-split or per-bin SNV access prohibitively slow.

This migration creates that index. Nothing else changes. The operation is
idempotent (IF NOT EXISTS) and wrapped in a single transaction so a crash
leaves the database at version 42.

Heads up: on a 700M-row variable_nucleotides table, the CREATE INDEX itself
takes ~30-60 minutes and grows the database file by ~5-15 GB. There is no
incremental progress reporting available from SQLite for this operation.
"""

import sys
import argparse

import anvio.tables as t
import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

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

    # Heads up to the user before we start the long part. We print this through
    # `run` rather than `progress.update` because progress messages get
    # overwritten and we want this one to stick around in the log.
    progress.reset()
    run.warning(
        f"Anvi'o is about to build an index on the 'split_name' column of the "
        f"'{t.variable_nts_table_name}' table. For typical profile databases this finishes in "
        f"seconds. For very large databases (hundreds of millions of SNV rows -- the kind you "
        f"get from massive read-recruitment workflows) it can take 30-60 minutes and grow the "
        f"profile db file by several gigabytes. The operation runs in a single transaction, so "
        f"if you interrupt it the database stays at version {current_version} and you can "
        f"safely rerun this migration later.",
        header="ABOUT TO BUILD AN INDEX -- THIS MAY TAKE A WHILE",
        lc='yellow'
    )

    conn = profile_db.conn
    conn.commit()
    cursor = conn.cursor()

    try:
        cursor.execute('BEGIN IMMEDIATE')

        progress.update(f"Creating index on {t.variable_nts_table_name}(split_name)...")
        cursor.execute(
            f'CREATE INDEX IF NOT EXISTS variable_nucleotides_split_name_idx '
            f'ON "{t.variable_nts_table_name}"(split_name)'
        )

        # Update version inside the same transaction so a partial state cannot
        # leave us at v43 without the index, or vice versa.
        cursor.execute("DELETE FROM self WHERE key = 'version'")
        cursor.execute("INSERT INTO self VALUES ('version', ?)", (next_version,))

        conn.commit()
    except Exception:
        conn.rollback()
        progress.end()
        raise

    profile_db.disconnect()
    progress.end()

    run.info_single(
        f"Your profile database is now version {next_version}. SNV-table queries that filter on "
        f"'split_name' will now use an index instead of a full table scan -- which is the "
        f"difference between subsecond and tens-of-minutes on very large profile databases.",
        nl_after=1, nl_before=1, mc='green'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            f"A simple script to upgrade an anvi'o profile database from version "
            f"{current_version} to version {next_version}."
        )
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
