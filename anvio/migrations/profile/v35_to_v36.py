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

views_table_name = "views"
views_table_structure = ["view_id", "target_table"]
views_table_types = ["str", "str"]

view_table_structure = ["contig", "sample", "value"]
view_table_types = ["text", "text", "numeric"]

essential_data_fields_for_anvio_profiles = [
    "std_coverage",
    "mean_coverage",
    "mean_coverage_Q2Q3",
    "detection",
    "abundance",
    "variability",
]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version=True)
    is_merged = profile_db.get_meta_value("merged")
    is_blank = profile_db.get_meta_value("blank")
    sample_name = profile_db.get_meta_value("sample_id")
    tables_in_db = profile_db.get_table_names()
    profile_db.disconnect()

    is_full_profile = (
        "mean_coverage_splits" in tables_in_db or "atomic_data_splits" in tables_in_db
    )

    run.info("Profile db type", "Merged" if is_merged else "Single")
    run.info("Full profile", is_full_profile)
    run.info("Is blank", is_blank)

    progress.new("Durr Durr")
    progress.update("...")

    if is_blank:
        ########################
        #     BLANK PROFILE    #
        ########################

        pass

    elif is_full_profile and not is_merged:
        #########################
        #     SINGLE PROFILE    #
        #########################
        profile_db = db.DB(db_path, None, ignore_version=True)

        # remove the default view variable in self, and add it back with 'mean_coverage'
        profile_db.remove_meta_key_value_pair("default_view")
        profile_db.set_meta_value("default_view", "mean_coverage")

        for target in ["splits", "contigs"]:
            # get rid of the hideous view called 'single'.
            profile_db._exec('''DELETE FROM views WHERE view_id = "single"''')

            atomic_data = profile_db.get_table_as_dict(f"atomic_data_{target}")

            for view in essential_data_fields_for_anvio_profiles:
                table_name = f"{view}_{target}"

                # le creationeaux au de neuvo tabl
                profile_db._exec(
                    f"""CREATE TABLE {table_name} (item text, layer text, value numeric)"""
                )
                view_data = []
                for split_name in atomic_data:
                    view_data.append(
                        (split_name, sample_name, atomic_data[split_name][view]),
                    )

                # populate the new view table
                profile_db._exec_many(
                    """INSERT INTO %s VALUES (?,?,?)""" % (table_name), view_data
                )

                # update the views table
                if target == "splits":
                    profile_db._exec(
                        """INSERT INTO views VALUES (?,?)""", (view, table_name)
                    )

        # баяртай
        profile_db.disconnect()

    elif is_full_profile and is_merged:
        #########################
        #     MERGED PROFILE    #
        #########################

        # open the profile database without rowid prepend.
        profile_db = db.DB(db_path, None, ignore_version=True, skip_rowid_prepend=True)

        # learn your samples
        sample_names = [
            s.strip() for s in profile_db.get_meta_value("samples").split(",")
        ]

        # drop the contents of the view table.
        profile_db._exec("DELETE FROM views")

        for target in ["splits", "contigs"]:
            for view in essential_data_fields_for_anvio_profiles:
                table_name = f"{view}_{target}"

                progress.update(f"Working on table '{table_name} ...'")

                table_data = profile_db.get_table_as_dict(table_name)

                # drop the old view table
                profile_db._exec(f"DROP TABLE {table_name}")

                # create a new view table!
                profile_db._exec(
                    f"""CREATE TABLE {table_name} (item text, layer text, value numeric)"""
                )

                # fill in the new view data from the old format
                view_data = []
                for split_name in table_data:
                    for sample_name in sample_names:
                        view_data.append(
                            (
                                split_name,
                                sample_name,
                                table_data[split_name][sample_name],
                            ),
                        )

                # populate new view table
                profile_db._exec_many(
                    """INSERT INTO %s VALUES (?,?,?)""" % (table_name), view_data
                )

                # if splits, I sits
                if target == "splits":
                    profile_db._exec(
                        """INSERT INTO views VALUES (?,?)""", (view, table_name)
                    )

        # さようなら
        profile_db.disconnect()

    else:
        ###########################
        #     SURPRISE PROFILE    #
        ###########################

        raise ConfigError(
            "Anvi'o is confuse. Your profile database does not fit into anything we have "
            "anticipated to run into here. For full disclosure, [the rest of the sentence "
            "was left blank intentionally just to drive you mad as you drive anvi'o mad -- "
            "eye for an eye]."
        )

    # set the version
    profile_db = db.DB(db_path, None, ignore_version=True)
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)
    profile_db.disconnect()

    progress.end()
    run.info_single(
        f"The profile database is now {next_version}. This upgrade fixed one of the most annoying "
        f"early design decisions we have made (and when we say 'we', we actually mean 'Meren', and "
        f"the rest of us accept no blame for it). This design shortcoming prevented anvi'o to merge "
        f"more than 2,000 samples. The current update reflects a significant change in the structure "
        f"of the 'view' tables of anvi'o and not only removes this limitation, but also results in "
        f"significant speed and memory gains during `anvi-merge`. But this operation is similar to "
        f"changing the entire flooring of an apartment while having to make sure each piece of the "
        f"furniture put back to their place properly once the flooring is redone. The aim of this "
        f"migration script was to put the furniture back. If you are reading this message, you are "
        f"most likely ⭐",
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
