#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split("_to_")]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_db_info = dbinfo.PanDBInfo(db_path)
    if str(pan_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided pan database is {pan_db_info.version}, not the required version, "
            f"{current_version}, so this script cannot upgrade the database."
        )

    pan_db = pan_db_info.load_db()

    if str(pan_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this contigs database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    progress.new("Migrating")
    progress.update("Updating the self table with a new variable")

    try:
        pan_db.remove_meta_key_value_pair("user_provided_gene_clusters_txt")
    except:
        pass

    pan_db.set_meta_value("user_provided_gene_clusters_txt", False)

    progress.update("Updating version")
    pan_db.remove_meta_key_value_pair("version")
    pan_db.set_version(next_version)

    progress.update("Committing changes")
    pan_db.disconnect()

    progress.end()

    message = (
        "Yes! Your pan database is now version 18. This wasn't a biggie and anvi'o just added a new "
        "variable into the self table to track where are the gene clusters coming from in a given "
        "pan-db."
    )
    run.info_single(message, nl_after=1, nl_before=1, mc="green")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade the pan database from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "pan_db",
        metavar="PAN_DB",
        help="An anvi'o pan database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
