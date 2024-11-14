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

current_version = "20"
next_version = "21"

item_orders_table_name = "item_orders"
item_orders_table_structure = ["name", "type", "data"]
item_orders_table_types = ["text", "text", "text"]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version=True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError(
            "Version of this profile database is not %s (hence, this script cannot really do anything)."
            % current_version
        )

    progress.new("Trying to upgrade the profile database")
    progress.update("...")

    try:
        profile_db.create_table(
            item_orders_table_name, item_orders_table_structure, item_orders_table_types
        )
    except:
        pass

    clusterings = profile_db.get_table_as_dict("clusterings")

    # move clustering data into the new table
    for clustering in clusterings:
        newick = clusterings[clustering]["newick"]
        profile_db._exec(
            """INSERT INTO %s VALUES (?,?,?)""" % item_orders_table_name,
            tuple([clustering, "newick", newick]),
        )

    # update keys
    for old_key, new_key in [
        ("available_clusterings", "available_item_orders"),
        ("contigs_clustered", "contigs_ordered"),
        ("default_clustering", "default_item_order"),
    ]:
        try:
            profile_db.set_meta_value(new_key, profile_db.get_meta_value(old_key))
        except:
            pass

    # remove stuff that are not irrelevant
    try:
        profile_db._exec("DROP TABLE clusterings;")
        profile_db.remove_meta_key_value_pair("available_clusterings")
        profile_db.remove_meta_key_value_pair("contigs_clustered")
        profile_db.remove_meta_key_value_pair("default_clustering")
    except:
        pass

    # commit
    try:
        profile_db._exec("COMMIT")
    except:
        pass

    # cleanup
    try:
        profile_db._exec("vacuum")
    except:
        pass

    # set the version
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single(
        "Your profile db is now %s." % next_version, nl_after=1, nl_before=1, mc="green"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade profile database from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "profile_db",
        metavar="PROFILE_DB",
        help="An anvi'o profile database of version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
