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

current_version = '17'
next_version    = '18'


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    is_merged = profile_db.get_meta_value('merged')
    is_blank = profile_db.get_meta_value('blank')
    is_full_profile = 'portion_covered_splits' in  profile_db.get_table_names()

    run.info('Profile db type', 'Merged' if is_merged else 'Single')
    run.info('Full profile', is_full_profile)

    progress.new("Trying to upgrade the profile database")
    progress.update('...')

    if is_merged:
        num_samples = len(profile_db.get_meta_value('samples').split(','))
        profile_db.remove_meta_key_value_pair('total_reads_mapped')
        profile_db.set_meta_value('total_reads_mapped', ', '.join(['0'] * num_samples))
    elif is_blank:
        profile_db.remove_meta_key_value_pair('total_reads_mapped')
        profile_db.set_meta_value('total_reads_mapped', '0')
    else:
        pass

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single('Your profile db is now %s.' % next_version, nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = "An anvi'o profile database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
