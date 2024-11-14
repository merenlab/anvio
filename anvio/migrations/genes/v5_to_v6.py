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

# tables that are specific to genes databases that require updatin'
tables = {
    "item_additional_data": {
        "structure": ["item_name", "data_key", "data_value", "data_type", "data_group"],
        "types": ["text", "text", "text", "text", "text"],
    },
    "layer_additional_data": {
        "structure": ["item_name", "data_key", "data_value", "data_type", "data_group"],
        "types": ["text", "text", "text", "text", "text"],
    },
    "gene_level_coverage_stats": {
        "structure": [
            "gene_callers_id",
            "sample_name",
            "mean_coverage",
            "detection",
            "non_outlier_mean_coverage",
            "non_outlier_coverage_std",
            "gene_coverage_values_per_nt",
            "non_outlier_positions",
        ],
        "types": [
            "numeric",
            "text",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "blob",
            "blob",
        ],
    },
    "gene_level_inseq_stats": {
        "structure": [
            "gene_callers_id",
            "sample_name",
            "mean_coverage",
            "insertions",
            "insertions_normalized",
            "mean_disruption",
            "below_disruption",
            "gene_coverage_values_per_nt",
        ],
        "types": [
            "numeric",
            "text",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "numeric",
            "blob",
        ],
    },
    "collections_bins_info": {
        "structure": ["collection_name", "bin_name", "source", "html_color"],
        "types": ["text", "text", "text", "text"],
    },
    "collections_of_contigs": {
        "structure": ["collection_name", "contig", "bin_name"],
        "types": ["text", "text", "text"],
    },
    "collections_of_splits": {
        "structure": ["collection_name", "split", "bin_name"],
        "types": ["text", "text", "text"],
    },
}


def drop_entry_id_column_from_table(db_path, table_name, table_properties):
    progress.new("Modifying '%s'" % table_name)

    structure = table_properties["structure"]
    types = table_properties["types"]
    db_fields = ", ".join(["%s %s" % (t[0], t[1]) for t in zip(structure, types)])
    temp_table_name = table_name + "_TEMP"

    _db = db.DB(db_path, None, ignore_version=True)

    # if table doesn't exist (due to some historical glitch), generate it, and return
    if table_name not in _db.get_table_names():
        _db._exec("""CREATE TABLE %s (%s)""" % (table_name, db_fields))
        _db.disconnect()
        progress.end()

        return

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
    utils.is_genes_db(db_path)

    genes_db = db.DB(db_path, None, ignore_version=True)
    if str(genes_db.get_version()) != current_version:
        genes_db.disconnect()
        raise ConfigError(
            "Version of this genes database is not %s (hence, this script cannot really do anything)."
            % current_version
        )
    genes_db.disconnect()

    # drop entry ids one by one
    for table_name in tables:
        drop_entry_id_column_from_table(
            db_path, table_name, table_properties=tables[table_name]
        )

    # set the version
    genes_db = db.DB(db_path, None, ignore_version=True)
    genes_db.remove_meta_key_value_pair("version")
    genes_db.set_version(next_version)
    genes_db.disconnect()

    progress.end()

    run.info_single(
        "Your genes db is now version %s. %d of its tables were cleaned from a historical "
        "design artifact." % (next_version, len(tables)),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade genes database from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "genes_db",
        metavar="GENES_DB",
        help="An anvi'o genes database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genes_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
