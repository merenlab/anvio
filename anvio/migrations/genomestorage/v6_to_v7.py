#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

# tables that are specific to genome storage databases that require updatin'
tables = {
    "gene_function_calls": {
        "structure": [
            "genome_name",
            "gene_callers_id",
            "source",
            "accession",
            "function",
            "e_value",
        ],
        "types": ["str", "numeric", "text", "text", "text", "numeric"],
    },
}


def drop_entry_id_column_from_table(db_path, table_name, table_properties):
    progress.new("Modifying '%s'" % table_name)

    structure = table_properties["structure"]
    types = table_properties["types"]
    db_fields = ", ".join(["%s %s" % (t[0], t[1]) for t in zip(structure, types)])
    temp_table_name = table_name + "_TEMP"

    _db = db.DB(db_path, None, ignore_version=True)

    progress.update("Creating a temporary table")
    _db._exec("""CREATE TABLE %s (%s)""" % (temp_table_name, db_fields))

    progress.update("Copying data into the temporary table")
    _db._exec(
        """INSERT INTO %s SELECT %s FROM %s"""
        % (temp_table_name, ", ".join(structure), table_name)
    )

    progress.update("Dropping the original table")
    _db._exec("""DROP TABLE IF EXISTS %s""" % (table_name))

    progress.update("Renaming temporary table to become the original")
    _db._exec("""ALTER TABLE %s RENAME TO %s""" % (temp_table_name, table_name))

    progress.update("Committing changes")
    _db.disconnect()

    progress.end()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    genomes_db = db.DB(db_path, None, ignore_version=True)
    if str(genomes_db.get_version()) != current_version:
        genomes_db.disconnect()
        raise ConfigError(
            "Version of this genome storage is not %s (hence, this script cannot really do anything)."
            % current_version
        )
    genomes_db.disconnect()

    # drop entry ids one by one
    for table_name in tables:
        drop_entry_id_column_from_table(
            db_path, table_name, table_properties=tables[table_name]
        )

    # set the version
    genomes_db = db.DB(db_path, None, ignore_version=True)
    genomes_db.remove_meta_key_value_pair("version")
    genomes_db.set_version(next_version)
    genomes_db.disconnect()

    progress.end()

    run.info_single(
        "Your genomes storage is now version %s. All good." % (next_version),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade a genome storage from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "genomes_storage",
        metavar="GENOMES_STORAGE",
        help="An anvi'o genomes storage of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genomes_storage)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
