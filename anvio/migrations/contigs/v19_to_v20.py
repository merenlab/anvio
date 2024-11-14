#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils


import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split("_to_")]

run = terminal.Run()
progress = terminal.Progress()


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

    tables_updated = False

    progress.new("Migrating things")

    gene_function_sources = contigs_db.get_meta_value("gene_function_sources")

    # there is no need to do anything if COG functions are not in the database
    if gene_function_sources and "COG_FUNCTION" in gene_function_sources:
        progress.update("Dealing with the functions table")

        functions_updated = []
        functions = contigs_db.get_some_rows_from_table("gene_functions", "1")

        for entry in functions:
            entry = list(entry)
            # (46520, 8114, 'COG_CATEGORY', 'J', 'Translation, ribosomal structure and biogenesis', 0)
            if entry[2] == "COG_CATEGORY":
                entry[2] = "COG14_CATEGORY"
            elif entry[2] == "COG_FUNCTION":
                entry[2] = "COG14_FUNCTION"
            else:
                pass

            functions_updated.append(entry[1:])

        contigs_db._exec("""DELETE FROM gene_functions WHERE 1""")
        contigs_db._exec_many(
            """INSERT INTO gene_functions VALUES (?,?,?,?,?)""", functions_updated
        )

        progress.update("Dealing with the self table")
        gene_function_sources = contigs_db.get_meta_value("gene_function_sources")

        if gene_function_sources:
            gene_function_sources = gene_function_sources.replace(
                "COG_FUNCTION", "COG14_FUNCTION"
            ).replace("COG_CATEGORY", "COG14_CATEGORY")
            contigs_db.remove_meta_key_value_pair("gene_function_sources")
            contigs_db.set_meta_value("gene_function_sources", gene_function_sources)

        tables_updated = True

    progress.update("Updating version")
    contigs_db.remove_meta_key_value_pair("version")
    contigs_db.set_version(next_version)

    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()

    if tables_updated:
        message = (
            "Done! Your contigs db is now version 20. Here is a little info on why this update was necessary (in case you have nothing "
            "else to read): Recently the NCBI released a new version of the COGs. The last time they did that, it was 2014, so it is "
            "kind of a big deal. We've made the necessary changes in anvi'o to make sure by default you are using the 2020 release. "
            "But to make sure COG functional annotations in existing contigs dbs are not lost in the process, we decided to keep the "
            "version name in the functional annotation source. This migration replaced `COG_FUNCTION` with `COG14_FUNCTION`, & `COG_CATEGORY` "
            "with `COG14_CATEGORY`. So the old annotations stay. Please note that if you would like to run `anvi-run-ncbi-cogs` with "
            "the old COGs, you will need to include in your command line `--version COG14`. Otherwise anvi'o will run the default ones, "
            "which is going to be COGs of the 2020 relase."
        )
    else:
        message = (
            "You didn't have any COG annotations, so you get a free pass to version 20."
        )

    run.info_single(message, nl_after=1, nl_before=1, mc="green")


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
