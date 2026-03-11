#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.dbinfo as dbinfo
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    pan_db_info = dbinfo.PanDBInfo(db_path)
    if str(pan_db_info.version) != current_version:
        raise ConfigError(
            f"The version of the provided pan database is {pan_db_info.version}, not the required "
            f"version, {current_version}, so this script cannot upgrade the database."
        )

    pan_db = pan_db_info.load_db()

    progress.new("Migrating")
    progress.update("Updating the self table with one new variable")

    if 'db_variant' not in pan_db_info.get_self_table():
        pan_db.set_meta_value('db_variant', 'sequence-based')

    progress.update("Updating version")
    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)

    progress.update("Committing changes")
    pan_db.disconnect()

    progress.end()

    message = (f"Done! Your pan database is now version {current_version}. The purpose of this chahge was to make space "
               f"for different pan-db variants. So far we only had sequence-based pangeomes. Now anvi'o can generate "
               f"structure-informed pangenomes to better identify distant homologs, and with this change the code will be "
               f"able to recognize which variant a given pan-db file is.")
    run.info_single(message, nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade the pan database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('pan_db', metavar = 'PAN_DB', help = "An anvi'o pan database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
