#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.utils as utils

import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version = True)

    is_blank = profile_db.get_meta_value('blank')
    is_merged = profile_db.get_meta_value('merged')

    progress.new("Harr Harr")
    progress.update('...')

    msg = ""

    if is_blank:
        msg = "But this was a blank profile database, so anvi'o did nothing."
        pass
    else:
        try:
            profile_db.remove_meta_key_value_pair('skip_edges_for_variant_profiling')
        except:
            pass

        profile_db.set_meta_value('skip_edges_for_variant_profiling', 0)

    #              خدا حافظ
    profile_db.disconnect()

    progress.end()

    run.info_single(f"The profile database is now {next_version}. We recently added a parameter for "
                    f"ancient DNA-friendly profiling of SNVs, SCVs, and SAAVs. If you are interested "
                    f"to learn more, see https://github.com/merenlab/anvio/pull/2081. Version bump was "
                    f"necessary to make sure all anvi'o profile-dbs have a means to track whether "
                    f"an additional parameter was used by the user to skip the edges of short reads "
                    f"during profiling.", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade PROFILE.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = 'Profile database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
