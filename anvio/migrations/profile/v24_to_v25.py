#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db
import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.constants import codons
from anvio.dbinfo import is_profile_db

run = terminal.Run()
progress = terminal.Progress()

current_version, next_version = [x[1:] for x in __name__.split('_to_')]


variable_codons_table_name           = 'variable_codons'
variable_codons_table_structure      = ['entry_id', 'sample_id', 'corresponding_gene_call', 'codon_order_in_gene', 'reference', 'departure_from_reference', 'coverage'] + codons
variable_codons_table_types          = [ 'numeric',    'text'  ,        'numeric'         ,       'numeric'      ,    'text'  ,          'numeric'        , 'numeric' ] + ['numeric'] * len(codons)


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
    tables_in_db = profile_db.get_table_names()
    is_full_profile = 'mean_coverage_Q2Q3_splits' in tables_in_db or 'atomic_data_splits' in tables_in_db

    run.info('Profile db type', 'Merged' if is_merged else 'Single')
    run.info('Full profile', is_full_profile)

    if is_full_profile:
        # add our new table
        profile_db.create_table(variable_codons_table_name, variable_codons_table_structure, variable_codons_table_types)

        # drop the old table
        profile_db._exec('DROP TABLE variable_amino_acid_frequencies;')

        # rename the table old
        profile_db._exec('ALTER TABLE variable_nucleotide_positions RENAME TO variable_nucleotides;')

        # remove stuff no longer necessary
        profile_db.remove_meta_key_value_pair('AA_frequencies_profiled')

        # add the sad fact
        profile_db.set_meta_value('SCVs_profiled', False)

        try:
            # clean after yourself
            profile_db.conn.isolation_level = None
            profile_db._exec('vacuum')
        except:
            pass

        full_upgrade = True
    else:
        full_upgrade = False

    # set the version
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)

    # bye
    profile_db.disconnect()
    progress.end()

    if full_upgrade:
        run.info_single("Your profile db is now version %s. If you had amino acids profiled for this database,\
                         you just lost all of that content :( The only option is to re-profile all your databases\
                         and merge them again. We are very sorry about the inconvenience. If you don't know what\
                         this message is talking about, then you have nothing to worry about." \
                                                            % next_version, nl_after=1, nl_before=1, mc='green')
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
