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

    utils.is_profile_db(db_path)

    progress.new("Durr Durr")
    progress.update("...")

    profile_db = db.DB(db_path, None, ignore_version=True)

    is_merged = profile_db.get_meta_value("merged")

    if is_merged:
        # merged profile
        try:
            profile_db._exec("DROP TABLE relative_abundance_contigs")
            profile_db._exec("DROP TABLE relative_abundance_splits")
            profile_db._exec("DROP TABLE max_normalized_ratio_contigs")
            profile_db._exec("DROP TABLE max_normalized_ratio_splits")
        except:
            pass

    else:
        # other profile
        pass

    # set the version
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)

    # 안녕
    profile_db.disconnect()

    progress.end()

    if is_merged:
        run.info_single(
            f"The profile database is now {next_version}. There were some unnecessary "
            f"tables in it, but they are no more.",
            nl_after=1,
            nl_before=1,
            mc="green",
        )
    else:
        run.info_single(
            f"The profile database is now {next_version}.",
            nl_after=1,
            nl_before=1,
            mc="green",
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade PROFILE.db from version %s to version %s"
        % (current_version, next_version)
    )
    parser.add_argument(
        "profile_db",
        metavar="PROFILE_DB",
        help="Profile database at version %s" % current_version,
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
