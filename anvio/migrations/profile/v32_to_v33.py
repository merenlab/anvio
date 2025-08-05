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

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    # make sure someone is not being funny
    is_profile_db(db_path)

    # make sure the version is accurate
    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # full as opposed to a "blank profile"
    tables_in_db = profile_db.get_table_names()
    is_full_profile = 'mean_coverage_Q2Q3_splits' in tables_in_db or 'atomic_data_splits' in tables_in_db

    # ---------------------------------------------------------------------------------

    progress.new('Updating DB')

    variable_nts_table_name = 'variable_nucleotides'
    if variable_nts_table_name in tables_in_db:
        progress.update('Renaming columns in variability table...')
        profile_db._exec("""ALTER TABLE %s RENAME COLUMN in_partial_gene_call TO in_noncoding_gene_call;""" % variable_nts_table_name)
        profile_db._exec("""ALTER TABLE %s RENAME COLUMN in_complete_gene_call TO in_coding_gene_call;""" % variable_nts_table_name)

    # ---------------------------------------------------------------------------------

    # bye
    profile_db.disconnect()
    progress.end()

    if is_full_profile:
        run.info_single("Your profile db is now %s. This update renamed two column names in the `variabile_nucleotides` table "
                        "of your profile database (`in_partial_gene_call` has become `in_noncoding_gene_call`, and "
                        "`in_complete_gene_call` has become `in_complete_gene_call`" % next_version, nl_after=1, nl_before=1, mc='green')
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
