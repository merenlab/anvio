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

taxonomy_estimation_bin_name             = 'taxonomy_estimation_bin'
taxonomy_estimation_bin_structure        = ['entry_id', 'collection_name', 'bin_name', 'source'  , 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
taxonomy_estimation_bin_types            = [ 'numeric',   'text'   ,        'text'  ,  'text',      'text',   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")


    # make sure someone is not being funny
    utils.is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version = True)


    profile_db.create_table(taxonomy_estimation_bin_name, taxonomy_estimation_bin_structure,taxonomy_estimation_bin_types)

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
