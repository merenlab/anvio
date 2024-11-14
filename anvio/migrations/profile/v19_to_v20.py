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

current_version = "19"
next_version = "20"


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

    # drop the table
    try:
        profile_db._exec("DROP TABLE gene_coverages;")
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

    # remove irrelevant self table entry
    try:
        profile_db.remove_meta_key_value_pair("gene_coverages_computed")
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
