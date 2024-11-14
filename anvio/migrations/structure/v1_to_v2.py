#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

# tables that are specific to structure databases that require updatin'
tables = {
    "templates": {
        "structure": ["corresponding_gene_call", "pdb_id", "chain_id", "ppi"],
        "types": ["integer", "text", "text", "real"],
    },
    "models": {
        "structure": [
            "corresponding_gene_call",
            "molpdf",
            "GA341_score",
            "DOPE_score",
            "picked_as_best",
        ],
        "types": ["integer", "real", "real", "real", "integer"],
    },
    "residue_info": {
        "structure": [
            "corresponding_gene_call",
            "codon_order_in_gene",
            "contact_numbers",
            "codon",
            "amino_acid",
            "codon_number",
            "sec_struct",
            "rel_solvent_acc",
            "phi",
            "psi",
            "NH_O_1_index",
            "NH_O_1_energy",
            "O_NH_1_index",
            "O_NH_1_energy",
            "NH_O_2_index",
            "NH_O_2_energy",
            "O_NH_2_index",
            "O_NH_2_energy",
        ],
        "types": [
            "integer",
            "integer",
            "text",
            "text",
            "text",
            "integer",
            "text",
            "real",
            "real",
            "real",
            "integer",
            "real",
            "integer",
            "real",
            "integer",
            "real",
            "integer",
            "real",
        ],
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

    # make sure someone is not being funny
    utils.is_structure_db(db_path)

    structure_db = db.DB(db_path, None, ignore_version=True)
    if str(structure_db.get_version()) != current_version:
        structure_db.disconnect()
        raise ConfigError(
            "Version of this structure database is not %s (hence, this script cannot really do anything)."
            % current_version
        )
    structure_db.disconnect()

    # drop entry ids one by one
    for table_name in tables:
        drop_entry_id_column_from_table(
            db_path, table_name, table_properties=tables[table_name]
        )

    # set the version
    structure_db = db.DB(db_path, None, ignore_version=True)
    structure_db.remove_meta_key_value_pair("version")
    structure_db.set_version(next_version)
    structure_db.disconnect()

    progress.end()

    run.info_single(
        "Your structure db is now version %s. %d of its tables were cleaned from a historical "
        "design artifact." % (next_version, len(tables)),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade structure database from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "structure_db",
        metavar="STRUCTURE_DB",
        help="An anvi'o structure database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.structure_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
