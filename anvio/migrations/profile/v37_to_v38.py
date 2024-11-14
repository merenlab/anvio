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

    profile_db = db.DB(db_path, None, ignore_version=True)

    is_blank = profile_db.get_meta_value("blank")
    is_merged = profile_db.get_meta_value("merged")

    progress.new("Durr Durr")
    progress.update("...")

    msg = ""

    if is_blank:
        # nothing to be done since we don't have
        # fetch filters for blank profiles
        msg = "But this was a blank profile database, so anvi'o did nothing."
        pass
    elif is_merged:
        # we need to update dis.
        try:
            profile_db.remove_meta_key_value_pair("fetch_filter")
        except:
            pass

        samples = profile_db.get_meta_value("samples").split(",")
        profile_db.set_meta_value("fetch_filter", ", ".join(["None"] * len(samples)))

        msg = (
            "This was a merged profile database, so anvi'o assumed none of the single profiles "
            "had any fetch filters, and marked them as such."
        )

        # BEWARE OF THIS CHEATING.
        # yes we are not here for this, but we will squeeze it in anyway. so far we have not been
        # tracking specific 'min percent identity' paramters set for individual single profiles
        # to filter short reads that are taken into consideration during profiling. but at this
        # stage of the codebase, the merger class does store that informaiton in merged profile
        # self tables, so we are also reprsenting that information in previous versions of
        # merged profile databases:
        profile_db.set_meta_value(
            "min_percent_identity", ", ".join(["0.0"] * len(samples))
        )
    else:
        try:
            profile_db.remove_meta_key_value_pair("fetch_filter")
        except:
            pass

        profile_db.set_meta_value("fetch_filter", "None")

        msg = (
            "This was a single profile database, so anvi'o marked it with a blank fetch filter (which "
            "really is the case since fetch filters are just being introduced in anvi'o, and any single "
            "profile database that was generated in previous versions do not have any fetch filters (unless "
            "you are Florian -- because if you are, you need to re-profile all your databases you had profiled "
            "profiled with a fetch filter)."
        )

    profile_db.set_version(next_version)

    #              خدا حافظ
    profile_db.disconnect()

    progress.end()

    run.info_single(
        f"The profile database is now {next_version}. We just added a very fancy feature in anvi'o, "
        f"'fetch filter', that enables you to define very specific filters regarding what to work with "
        f"from BAM files during profiling, and this update is all about that. {msg}",
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
