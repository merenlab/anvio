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


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_genes_db(db_path)

    # make sure the version is accurate
    genes_db = db.DB(db_path, None, ignore_version=True)
    if str(genes_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this genes database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    try:
        genes_db.remove_meta_key_value_pair("mode")
    except:
        pass

    genes_db.update_meta_value("mode", "STANDARD")

    # set the version
    genes_db.remove_meta_key_value_pair("version")
    genes_db.set_version(next_version)

    # bye
    genes_db.disconnect()
    progress.end()

    run.info_single(
        "Your genes db is now %s (anvi'o set its mode to 'standard', so if you had been working with TN-Seq databases,\
                     you are in trouble (applies only to two people in the world at the time this message was being written, \
                     so PROBABLY it is OK if you ignore this message))."
        % next_version,
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
