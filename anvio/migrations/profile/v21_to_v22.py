#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import gzip
import h5py
import time
import argparse

import anvio.db as db
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version = "21"
next_version = "22"

run = terminal.Run()
progress = terminal.Progress()

split_coverages_table_name       = 'split_coverages'
split_coverages_table_structure  = ['split_name', 'sample_name', 'coverages']
split_coverages_table_types      = [    'str'   ,     'str'    ,   'blob'   ]

item_additional_data_table_name      = 'item_additional_data'
item_additional_data_table_structure = ['entry_id', 'sample_id', 'split_name', 'key', 'value', 'type']
item_additional_data_table_types     = [ 'numeric',    'text'  ,    'text'   , 'str',  'str' ,  'str']


def convert_numpy_array_to_binary_blob(array, compress=True):
    if compress:
        return gzip.compress(memoryview(array), compresslevel=1)
    else:
        return memoryview(array)


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    dbops.is_profile_db(db_path)

    profile_db = db.DB(db_path, None, ignore_version = True)
    if str(profile_db.get_version()) != current_version:
        raise ConfigError("Version of this profile database is not %s (hence, this script cannot really do anything)." % current_version)

    auxiliary_path = os.path.join(os.path.dirname(db_path), 'AUXILIARY-DATA.h5')
    new_auxiliary_path = os.path.join(os.path.dirname(db_path), 'AUXILIARY-DATA.db')

    if os.path.exists(auxiliary_path):
        fp = h5py.File(auxiliary_path, 'r')
        G = lambda x: fp.attrs[x].decode('utf-8') if isinstance(fp.attrs[x], bytes) else fp.attrs[x]
        auxiliary_db = db.DB(new_auxiliary_path, '2', new_database=True)

        auxiliary_db.set_meta_value('db_type', 'auxiliary data for coverages')
        auxiliary_db.set_meta_value('contigs_db_hash', G('hash'))
        auxiliary_db.set_meta_value('creation_date', time.time())
        auxiliary_db.create_table(split_coverages_table_name, split_coverages_table_structure, split_coverages_table_types)
        auxiliary_db._exec("""CREATE INDEX IF NOT EXISTS covering_index ON %s(split_name, sample_name)""" % (split_coverages_table_name))

        sample_names_in_db = set(list(list(fp['/data/coverages'].values())[0].keys()))
        split_names_in_db = list(fp['/data/coverages'].keys())

        run.info("Auxiliary data file found", auxiliary_path)
        run.info("Splits found", len(split_names_in_db))
        run.info("Samples found", len(sample_names_in_db))
        run.info("New auxiliary data path", new_auxiliary_path)

        progress.new('Processing auxiliary')
        counter, total = 0, len(sample_names_in_db)

        entries = []
        for sample_name in sample_names_in_db:
            for split_name in split_names_in_db:
                entries.append((split_name, sample_name, convert_numpy_array_to_binary_blob(fp['/data/coverages/%s/%s' % (split_name, sample_name)].value),))

            counter += 1
            progress.update('sample %d of %d ...' % (counter, total))

            if counter % 10 == 0:
                progress.update("Writing buffer into a new database file ...")
                auxiliary_db.insert_many(split_coverages_table_name, entries=entries)
                entries = []

        auxiliary_db.insert_many(split_coverages_table_name, entries=entries)

        progress.end()
        auxiliary_db.disconnect()
        fp.close()

        os.remove(auxiliary_path)

        run.info_single("(anvi'o just created a new, up-to-date auxiliary data file (which ends with extension .db), and deleted \
                         the old one (the one that ended with the extension .h5)).", nl_before=2, nl_after=1)
    else:
         run.warning("Althought he actual purpose of this script is to upgrade your AUXILIARY-DATA.h5 file, it doesn't seem to be where \
                      anvi'o expects it to be, Anvi'o will upgrade your Profile.db alone but as a consequence you will not be able to \
                      use its auxiliary data with this profile database.")

    # we also added a totally new table to this version:
    profile_db.create_table(item_additional_data_table_name, item_additional_data_table_structure, item_additional_data_table_types)

    profile_db.remove_meta_key_value_pair('version')
    profile_db.set_version(next_version)
    profile_db.disconnect()

    run.info_single('Done! Your profile db is now version %s.' % next_version, nl_after=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade profile database and AUXILIARY-DATA.h5 from version %s to version %s' % (current_version, next_version))
    parser.add_argument('profile_db', metavar = 'PROFILE_DB', help = "An anvi'o profile database of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.profile_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
