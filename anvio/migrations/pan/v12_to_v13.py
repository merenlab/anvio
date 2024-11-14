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
    utils.is_pan_db(db_path)

    # make sure the version is accurate
    pan_db = db.DB(db_path, None, ignore_version=True)
    if str(pan_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this pan database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    # gene_clusters_ordered -> items_ordered
    pan_db.set_meta_value(
        "items_ordered", pan_db.get_meta_value("gene_clusters_ordered")
    )
    pan_db.remove_meta_key_value_pair("gene_clusters_ordered")
    pan_db.remove_meta_key_value_pair("gene_clusters_clustered")

    # set the version
    pan_db.remove_meta_key_value_pair("version")
    pan_db.set_version(next_version)

    # now bye for real!
    pan_db.disconnect()

    progress.end()

    run.info_single(
        "Your pan db is now %s (lucky you)." % next_version,
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade pan database from version %s to version %s"
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
