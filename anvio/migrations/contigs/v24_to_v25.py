#!/usr/bin/env python

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.tables.sequencefeatures import create_sequence_features_tables

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(
            f"The version of the provided contigs database is {contigs_db.get_version()}, not the "
            f"required version, {current_version}, so this script cannot upgrade the database."
        )

    progress.new("Adding sequence-features tables")
    progress.update("...")

    # delegate to the single source of truth for the v25 sequence-features
    # schema. `create_sequence_features_tables` creates the five new tables,
    # their indexes, and seeds the builtin `feature_types` registry. The
    # same function is invoked by `ContigsDatabase.touch` so fresh
    # databases and migrated databases end up in the identical end state.
    create_sequence_features_tables(contigs_db)

    progress.end()

    progress.new("Updating version")
    progress.update("...")
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    contigs_db.disconnect()
    progress.end()

    message = (
        f"Congratulations! Your contigs database is now version {next_version}. Five additive "
        f"tables for arbitrary sequence features (`contigs_sequence_features`, `feature_types`, "
        f"`feature_relationships`, `feature_qualifiers`, and `CDS_features`) have been created "
        f"and the builtin feature-type registry has been populated. No existing rows in any "
        f"other table were touched. You can now use `anvi-import-genbank-features` to populate "
        f"the new tables from a GenBank file."
    )
    run.info_single(message, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar='CONTIGS_DB', help='Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
