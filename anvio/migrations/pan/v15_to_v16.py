#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

def migrate(db_path):
    pan_db = db.DB(db_path, None, ignore_version = True, skip_rowid_prepend=True)

    progress.new("Bleep bloop")
    progress.update('...')

    try:
        is_sensitive = int(pan_db.get_meta_value('diamond_sensitive'))
    except:
        is_sensitive = 0

    try:
        pan_db.remove_meta_key_value_pair('diamond_sensitive')
        pan_db.remove_meta_key_value_pair('additional_params_for_seq_search')
    except:
        pass

    if is_sensitive:
        pan_db.set_meta_value('additional_params_for_seq_search', "--sensitive")
    else:
        pan_db.set_meta_value('additional_params_for_seq_search', "")

    pan_db.remove_meta_key_value_pair('version')
    pan_db.set_version(next_version)
    pan_db.disconnect()

    progress.end()

    run.info_single(f"The pan database is now {next_version}. This upgrade sync'd the database design so there is "
                    f"more room for activities. And by activities, we mean `--additional-params-for-seq-search`. "
                    f"See the help menu for the program `anvi-pan-genome` for details :)",
                    nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade pan database from version %s to version %s' % (current_version, next_version))
    parser.add_argument('pan_db', metavar = 'PAN_DB', help = "An anvi'o pan database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.pan_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
