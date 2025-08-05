#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db

import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_profile_db

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

indels_table_name      = 'indels'
indels_table_structure = ['sample_id', 'split_name', 'pos'    , 'pos_in_contig', 'corresponding_gene_call', 'in_noncoding_gene_call', 'in_coding_gene_call' , 'base_pos_in_codon', 'codon_order_in_gene', 'cov_outlier_in_split', 'cov_outlier_in_contig', 'reference', 'type', 'sequence', 'length' , 'count'  , 'coverage']
indels_table_types     = ['text'     , 'text'      , 'integer', 'integer'      , 'integer'                , 'integer'               , 'integer'             , 'integer'          , 'integer'            , 'integer'             , 'integer'              , 'text'     , 'text', 'text'    , 'integer', 'integer', 'integer']

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    progress.new("Redefining indels table...")
    progress.update("...")

    # delete if exists
    try:
        profile_db.drop_table(indels_table_name)
    except:
        pass

    # create with new structure
    profile_db.create_table(indels_table_name, indels_table_structure, indels_table_types)

    progress.update("Updating self table")
    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)
    profile_db.set_meta_value('min_indel_fraction', '0.0')
    profile_db.set_meta_value('INDELs_profiled', '0')

    progress.update("Committing changes")
    profile_db.disconnect()

    progress.end()
    run.info_single("The profile database is now %s. This upgrade redefined the stored format of INDELS to "
                    "provide a more robust working framework. If you are upgrading this database from `v6`, "
                    "you don't have anything to worry about. But if you were using the active branch of anvi'o, "
                    "then you lost your INDELs now and you would need to re-profile your BAM files if you want "
                    "them back :)" % (next_version), nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade PROFILE.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = 'Profile database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
