#!/usr/bin/env python

"""Migration script to add the `clippings` table to profile databases.

The `clippings` table stores read-edge soft- and hard-clip events detected
during BAM profiling. It is structurally similar to `indels`, but with an
additional `side` column (L/R) to disambiguate which end of the read clipped
at a given breakpoint position.

This migration is purely additive: the `clippings` table is created empty.
Existing profile databases will simply have no clipping data until they are
re-profiled. Single, merged, and blank profile databases are all handled
identically.
"""

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

clippings_table_name      = 'clippings'
clippings_table_structure = ['sample_id', 'split_name', 'pos', 'pos_in_contig', 'corresponding_gene_call', 'in_noncoding_gene_call', 'in_coding_gene_call', 'base_pos_in_codon', 'codon_order_in_gene', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'reference', 'type', 'side', 'sequence', 'length', 'count', 'coverage']
clippings_table_types     = ['text'     , 'text'      , 'integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'text', 'text', 'text', 'text', 'integer', 'integer', 'integer']


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
    progress.update("Creating empty `clippings` table")

    profile_db = profile_db_info.load_db()

    existing_tables = profile_db.get_table_names()
    if clippings_table_name in existing_tables:
        # extremely unlikely, but be defensive: a previous half-applied migration
        # could leave the table behind. We refuse to clobber it.
        progress.end()
        raise ConfigError(
            f"A table named '{clippings_table_name}' already exists in this profile database, "
            f"which is unexpected at version {current_version}. Anvi'o refuses to overwrite it. "
            f"If you know what you are doing, drop the table manually and re-run the migration."
        )

    profile_db.create_table(clippings_table_name, clippings_table_structure, clippings_table_types)

    progress.update("Updating version")
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    profile_db.disconnect()

    progress.end()

    run.info_single(
        f"Your profile database is now version {next_version}. An empty `clippings` table has "
        f"been added. To populate it with clip events, you will need to re-profile the BAM file "
        f"with the current version of anvi-profile.",
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
