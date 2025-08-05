#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import argparse

import anvio.db as db
import anvio.dictio as dictio
import anvio.terminal as terminal 

from anvio.errors import ConfigError
from anvio.dbinfo import is_profile_db


run = terminal.Run()
progress = terminal.Progress()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No profile database is given.")

    # make sure someone is not being funny
    is_profile_db(db_path)

    # make sure the version is 5
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != '14':
        raise ConfigError("Version of this profile database is not 14 (hence, this script cannot really do anything).")

    is_merged = profile_db.get_meta_value('merged')

    progress.new("Trying to upgrade the %s profile database" % 'merged' if is_merged else 'single')

    # update the runinfo.cp
    input_dir = os.path.dirname(os.path.abspath(db_path))
    P = lambda x: os.path.join(input_dir, x)
    E = lambda x: os.path.exists(x)

    runinfo_path = P('RUNINFO.cp') if E(P('RUNINFO.cp')) else None
    runinfo_path = P('RUNINFO.mcp') if E(P('RUNINFO.mcp')) else None

    if runinfo_path:
        runinfo = dictio.read_serialized_object(runinfo_path)
        if 'blank' not in runinfo:
            runinfo['blank'] = False

            dictio.write_serialized_object(runinfo, runinfo_path)

    # add the new value
    profile_db.set_meta_value('blank', False)

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version('15')

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single("Database successfully upgraded to version 15!", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database to from version 14 version 15')
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = 'Profile database (of version 14)')
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
