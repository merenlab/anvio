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


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No profile database is given.")

    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    # make sure the version is 15
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != '15':
        raise ConfigError("Version of this profile database is not 15 (hence, this script cannot really do anything).")

    is_merged = profile_db.get_meta_value('merged')

    progress.new("Trying to upgrade the %s profile database" % 'merged' if is_merged else 'single')

    available_clusterings = []
    clusterings = profile_db.get_table_as_dict('clusterings')

    if clusterings:
        profile_db._exec('''DELETE FROM clusterings''')

        for entry in clusterings:
            clustering_id = ':'.join([entry, 'euclidean', 'ward'])
            clustering_newick = clusterings[entry]['newick']
            profile_db._exec('''INSERT INTO clusterings VALUES (?,?)''', tuple([clustering_id, clustering_newick]))
            available_clusterings.append(clustering_id)

        profile_db.remove_meta_key_value_pair('available_clusterings')
        profile_db.set_meta_value('available_clusterings', ','.join(available_clusterings))

        default_clustering = profile_db.get_meta_value('default_clustering')
        profile_db.remove_meta_key_value_pair('default_clustering')
        profile_db.set_meta_value('default_clustering', ':'.join([default_clustering, 'euclidean', 'ward']))

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version('16')

    # bye
    profile_db.disconnect()
    progress.end()

    run.info_single("Database successfully upgraded to version 16!", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database to from version 15 version 16')
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = 'Profile database (of version 15)')
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
