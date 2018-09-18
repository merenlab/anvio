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

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    try:
        profile_db.remove_meta_key_value_pair('max_contig_length')
    except:
        pass

    profile_db.update_meta_value('max_contig_length', sys.maxsize)

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = "An anvi'o profile database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
