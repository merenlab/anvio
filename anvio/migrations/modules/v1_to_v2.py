#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse
import re

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

    utils.is_kegg_modules_db(db_path)

    # make sure the current version is 1
    modules_db = db.DB(db_path, None, ignore_version=True)
    if str(modules_db.get_version()) != current_version:
        modules_db.disconnect()
        raise ConfigError(
            "Version of this modules database is not %s (hence, "
            "this script cannot really do anything)." % current_version
        )

    # split any orthology entries with multiple KOs into separate lines
    where_clause = (
        "data_name = 'ORTHOLOGY' AND (data_value LIKE '%+%' OR data_value LIKE '%-%')"
    )
    bad_rows = modules_db.get_some_rows_from_table("modules", where_clause)

    modules_db.remove_some_rows_from_table("modules", where_clause)

    good_rows = []
    for mod, orth, kos, definition, line in bad_rows:
        if orth != "ORTHOLOGY":
            raise ConfigError(
                "Serious SQL error here. We asked for ORTHOLOGY but we got %s" % orth
            )
        split_kos = re.split("\+|\-", kos)
        for ko in split_kos:
            if len(ko) != 6 or ko[0] != "K":
                raise ConfigError(
                    "Uh oh. We split a KO that doesn't seem to have the right format. It looks like this: %s, "
                    "and the string that it came from is %s." % (ko, kos)
                )
            good_rows.append((mod, orth, ko, definition, line))

    modules_db.insert_many("modules", entries=good_rows)

    prev_total_entries = modules_db.get_meta_value("total_entries")
    modules_db.remove_meta_key_value_pair("total_entries")
    modules_db.set_meta_value(
        "total_entries", prev_total_entries - len(bad_rows) + len(good_rows)
    )
    modules_db.remove_meta_key_value_pair("version")
    modules_db.set_version(next_version)
    modules_db.disconnect()

    run.info_single(
        "Huzzah! Your modules db is now %s. This update fixed KOs that did not have individual "
        "entries in the modules table." % (next_version),
        nl_after=1,
        nl_before=1,
        mc="green",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="A simple script to upgrade KEGG Modules database from version 1 to version 2"
    )
    parser.add_argument(
        "modules_db", metavar="MODULES_DB", help="KEGG Modules database"
    )
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.modules_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
