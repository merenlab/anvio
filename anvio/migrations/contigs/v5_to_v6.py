#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_contigs_db(db_path)

    # make sure the version is 2
    contigs_db = db.DB(db_path, None, ignore_version=True)
    if str(contigs_db.get_version()) != "5":
        raise ConfigError(
            "Version of this contigs database is not 5 (hence, this script cannot really do anything)."
        )

    progress.new("Trying to upgrade the contigs database")
    progress.update("...")

    # drop the old tables:
    try:
        contigs_db._exec("""DROP TABLE %s""" % (t.splits_taxonomy_table_name))
        contigs_db._exec("""DROP TABLE %s""" % (t.taxon_names_table_name))
        contigs_db._exec("""DROP TABLE %s""" % (t.genes_taxonomy_table_name))
    except:
        pass
    contigs_db.commit()

    # create new empty ones
    contigs_db.create_table(
        t.splits_taxonomy_table_name,
        t.splits_taxonomy_table_structure,
        t.splits_taxonomy_table_types,
    )
    contigs_db.create_table(
        t.taxon_names_table_name,
        t.taxon_names_table_structure,
        t.taxon_names_table_types,
    )
    contigs_db.create_table(
        t.genes_taxonomy_table_name,
        t.genes_taxonomy_table_structure,
        t.genes_taxonomy_table_types,
    )

    # set the version
    contigs_db.remove_meta_key_value_pair("version")
    contigs_db.set_version("6")

    # bye
    contigs_db.disconnect()

    # bye
    progress.end()
    run.info_single("The contigs database successfully upgraded from version 5 to 6!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade contigs database from version 5 to version 6"
    )
    parser.add_argument("contigs_db", metavar="CONTIGS_DB", help="Contigs database")
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
