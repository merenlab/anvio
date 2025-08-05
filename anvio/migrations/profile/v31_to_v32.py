#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_profile_db

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

indels_table_name      = 'indels'
indels_table_structure = ['sample_id', 'split_name', 'type', 'sequence', 'start_in_contig', 'start_in_split', 'length' , 'coverage']
indels_table_types     = ['text'     , 'text'      , 'text', 'text'    , 'numeric'        , 'numeric'       , 'numeric', 'numeric' ]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    try:
        profile_db.create_table(indels_table_name, indels_table_structure, indels_table_types)
    except:
        pass

    profile_db.set_meta_value('INDELs_profiled', 0)
    profile_db.set_meta_value('min_percent_identity', 0)

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # full as opposed to a "blank profile"
    tables_in_db = profile_db.get_table_names()
    is_full_profile = 'mean_coverage_Q2Q3_splits' in tables_in_db or 'atomic_data_splits' in tables_in_db

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    if is_full_profile:
        run.info_single("Your profile db is now %s. We just added a bunch of new variables to the `self` table "
                        "of your database. All good now." % next_version, nl_after=1, nl_before=1, mc='green')
    else:
        run.info_single("Your profile db is now version %s. But essentially nothing really happened to your "
                        "database since it was a blank profile (which is OK, move along)." \
                                                            % next_version, nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = "An anvi'o profile database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
