#!/usr/bin/env python

"""Migration script to add the (empty) table for read identity histograms to auxiliary databases.

This migration creates the table that stores histograms of the percent identity of reads to the
reference per split and sample, and sets the `identity_histograms_stored` flag to False, since
this information cannot be recovered from the existing data.
"""

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

# defined here rather than imported from anvio.tables, so this script keeps working if the table
# definition changes in the future
identity_histograms_table_name = 'split_identity_histograms'
identity_histograms_table_structure = ['split_name', 'sample_name', 'num_below_min', 'histogram']
identity_histograms_table_types = ['str', 'str', 'integer', 'blob']

run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    auxiliary_db_info = dbinfo.AuxiliaryDBInfo(db_path)
    if str(auxiliary_db_info.version) != current_version:
        raise ConfigError(f"The version of the provided auxiliary database is {auxiliary_db_info.version}, not the "
                          f"required version, {current_version}, so this script cannot upgrade the database.")

    progress.new("Migrating")
    auxiliary_db = auxiliary_db_info.load_db()

    progress.update("Creating the table for read identity histograms")
    if identity_histograms_table_name not in auxiliary_db.get_table_names():
        auxiliary_db.create_table(identity_histograms_table_name, identity_histograms_table_structure, identity_histograms_table_types)

    auxiliary_db._exec(f"""CREATE INDEX IF NOT EXISTS identity_histograms_index ON {identity_histograms_table_name}(split_name, sample_name)""")

    # set the self table flag to indicate that the new table is empty to downstream programs
    auxiliary_db.set_meta_value('identity_histograms_stored', False)

    progress.update("Updating the version")
    auxiliary_db.remove_meta_key_value_pair('version')
    auxiliary_db.set_version(next_version)

    auxiliary_db.disconnect()
    progress.end()

    run.info_single(f"Congratulations! Your auxiliary database is now version {next_version}. This update added a new table "
                    f"to store histograms of the percent identity of reads to the reference genomes for each split and sample. "
                    f"Please note that this table is EMPTY in your database, since this information cannot be recovered from "
                    f"the existing data. If you need it, you will have to re-profile your BAM files with `anvi-profile`.",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=f"A simple script to upgrade an anvi'o auxiliary database from version {current_version} to version {next_version}")
    parser.add_argument("auxiliary_db", metavar="AUXILIARY_DB", help=f"An anvi'o auxiliary database of version {current_version}")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.auxiliary_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
