#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version = "10"
next_version = "11"

run = terminal.Run()
progress = terminal.Progress()

genes_in_splits_table_name = "genes_in_splits"
genes_in_splits_table_structure = [
    "entry_id",
    "split",
    "gene_callers_id",
    "start_in_split",
    "stop_in_split",
    "percentage_in_split",
]
genes_in_splits_table_types = [
    "numeric",
    "text",
    "numeric",
    "numeric",
    "numeric",
    "numeric",
]

genes_in_splits_summary_table_name = "genes_in_splits_summary"
genes_in_splits_summary_table_structure = [
    "split",
    "num_genes",
    "avg_gene_length",
    "ratio_coding",
]
genes_in_splits_summary_table_types = ["text", "numeric", "numeric", "numeric"]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this contigs database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    progress.new("Removing '" + genes_in_splits_summary_table_name + "'")
    contigs_db._exec("DROP TABLE %s;" % genes_in_splits_summary_table_name)
    progress.end()

    progress.new("Upgrading '" + genes_in_splits_table_name + "'")

    progress.update("Creating temporary table")
    contigs_db.create_table(
        genes_in_splits_table_name + "_temp",
        genes_in_splits_table_structure,
        genes_in_splits_table_types,
    )

    progress.update("Moving unique records")
    contigs_db._exec(
        "INSERT INTO %s SELECT * FROM %s GROUP BY %s;"
        % (
            genes_in_splits_table_name + "_temp",
            genes_in_splits_table_name,
            ", ".join(genes_in_splits_table_structure[1:]),
        )
    )
    progress.update("Updating entry_id")
    contigs_db._exec(
        "UPDATE %s SET entry_id = rowid - 1;" % (genes_in_splits_table_name + "_temp")
    )

    progress.update("Swapping temporary table with actual table")
    contigs_db._exec(
        "ALTER TABLE %s RENAME TO %s;"
        % (genes_in_splits_table_name, genes_in_splits_table_name + "_old")
    )
    contigs_db._exec(
        "ALTER TABLE %s RENAME TO %s;"
        % (genes_in_splits_table_name + "_temp", genes_in_splits_table_name)
    )

    progress.update("Removing old table")
    contigs_db._exec("DROP TABLE %s;" % (genes_in_splits_table_name + "_old"))

    try:
        progress.update("Optimizing the database")
        contigs_db.conn.isolation_level = None
        contigs_db._exec("VACUUM;")
    except:
        pass

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair("version")
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()
    run.info_single(
        "The contigs database is now %s." % (next_version),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade CONTIGS.db from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "contigs_db",
        metavar="CONTIGS_DB",
        help="Contigs database at version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
