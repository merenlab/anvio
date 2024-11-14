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

current_version = "16"
next_version = "17"


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

    is_merged = profile_db.get_meta_value("merged")
    tables_in_db = profile_db.get_table_names()
    is_full_profile = (
        "portion_covered_splits" in tables_in_db or "atomic_data_splits" in tables_in_db
    )

    run.info("Profile db type", "Merged" if is_merged else "Single")
    run.info("Full profile", is_full_profile)

    progress.new("Trying to upgrade the profile database")
    progress.update("...")

    if is_full_profile and is_merged:
        profile_db._exec(
            "ALTER TABLE portion_covered_splits RENAME TO detection_splits;"
        )
        profile_db._exec(
            "ALTER TABLE portion_covered_contigs RENAME TO detection_contigs;"
        )
        profile_db._exec(
            "ALTER TABLE mean_coverage_Q1Q3_splits RENAME TO mean_coverage_Q2Q3_splits;"
        )
        profile_db._exec(
            "ALTER TABLE mean_coverage_Q1Q3_contigs RENAME TO mean_coverage_Q2Q3_contigs;"
        )

        profile_db._exec(
            'DELETE FROM %s WHERE view_id = "portion_covered"' % (t.views_table_name)
        )
        profile_db._exec(
            'INSERT INTO %s VALUES ("detection", "detection_splits")'
            % t.views_table_name
        )
        profile_db._exec(
            'DELETE FROM %s WHERE view_id = "mean_coverage_Q1Q3"' % (t.views_table_name)
        )
        profile_db._exec(
            'INSERT INTO %s VALUES ("mean_coverage_Q2Q3", "mean_coverage_Q2Q3_splits")'
            % t.views_table_name
        )

    elif is_full_profile and not is_merged:
        profile_db.cursor.execute(
            "ALTER TABLE atomic_data_contigs RENAME TO atomic_data_contigs_TEMP;"
        )
        profile_db.cursor.execute(
            "CREATE TABLE atomic_data_contigs (contig text, std_coverage numeric, mean_coverage numeric, mean_coverage_Q2Q3 numeric, max_normalized_ratio numeric, relative_abundance numeric, detection numeric, abundance numeric, variability numeric, __parent__ text);"
        )
        profile_db.cursor.execute(
            "INSERT INTO atomic_data_contigs(contig, std_coverage, mean_coverage, mean_coverage_Q2Q3, max_normalized_ratio, relative_abundance, detection, abundance, variability, __parent__) SELECT contig, std_coverage, mean_coverage, mean_coverage_Q1Q3, max_normalized_ratio, relative_abundance, portion_covered, abundance, variability, __parent__ FROM atomic_data_contigs_TEMP;"
        )
        profile_db.cursor.execute("DROP TABLE atomic_data_contigs_TEMP;")

        profile_db.cursor.execute(
            "ALTER TABLE atomic_data_splits RENAME TO atomic_data_splits_TEMP;"
        )
        profile_db.cursor.execute(
            "CREATE TABLE atomic_data_splits (contig text, std_coverage numeric, mean_coverage numeric, mean_coverage_Q2Q3 numeric, max_normalized_ratio numeric, relative_abundance numeric, detection numeric, abundance numeric, variability numeric, __parent__ text);"
        )
        profile_db.cursor.execute(
            "INSERT INTO atomic_data_splits(contig, std_coverage, mean_coverage, mean_coverage_Q2Q3, max_normalized_ratio, relative_abundance, detection, abundance, variability, __parent__) SELECT contig, std_coverage, mean_coverage, mean_coverage_Q1Q3, max_normalized_ratio, relative_abundance, portion_covered, abundance, variability, __parent__ FROM atomic_data_splits_TEMP;"
        )
        profile_db.cursor.execute("DROP TABLE atomic_data_splits_TEMP;")

    # update states
    states = profile_db.get_table_as_dict(t.states_table_name)
    for state in states:
        profile_db._exec(
            'DELETE FROM %s WHERE name = "%s"' % (t.states_table_name, state)
        )
        profile_db._exec(
            "INSERT INTO %s VALUES (?,?,?)" % (t.states_table_name),
            (
                state,
                states[state]["content"]
                .replace("portion_covered", "detection")
                .replace("mean_coverage_Q1Q3", "mean_coverage_Q2Q3"),
                states[state]["last_modified"],
            ),
        )

    # set the version
    profile_db.remove_meta_key_value_pair("version")
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single(
        "Database successfully upgraded to version 17!",
        nl_after=1,
        nl_before=1,
        mc="green",
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
