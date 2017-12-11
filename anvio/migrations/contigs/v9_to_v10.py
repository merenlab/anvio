#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import gzip
import h5py
import argparse

import anvio.db as db
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version = "9"
next_version = "10"

run = terminal.Run()
progress = terminal.Progress()

nt_position_info_table_name       = 'nt_position_info'
nt_position_info_table_structure  = ['contig_name', 'position_info']
nt_position_info_table_types      = [    'str'    ,      'blob'    ]

def convert_numpy_array_to_binary_blob(array, compress=True):
    if compress:
        return gzip.compress(memoryview(array), compresslevel=1)
    else:
        return memoryview(array)


def update_contigs_db(contigs_db_path, just_do_it=False):
    if contigs_db_path is None:
        raise ConfigError("No database path is given.")

    dbops.is_contigs_db(contigs_db_path)

    contigs_db = db.DB(contigs_db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    if not just_do_it:
        try:
            run.warning("This script will try to upgrade your profile database from v%s to v%s. Here is the idea behind this upgrade: \
                         We have been maintaining an axuiliary data file for anvi'o contigs databases. But we decided that this data \
                         could be stored in contigs databases, reducing the clutter on our file systems. If you think you are ready, \
                         just press ENTER to continue. If you want to cancel the upgrade and think more about it, press CTRL+C now. \
                         If you want to avoid this message the next time, use '--just-do-it'. PLEASE NOTE, running this script will\
                         remove the obsolete '.h5' file from your directory once it is successfully finished upgrading your contigs \
                         database. If you would like to save it for some reason, please consider backing it up first." % (current_version, next_version))
            input("Press ENTER to continue...\n")
        except:
            print()
            sys.exit()

    auxiliary_path = ''.join(contigs_db_path[:-3]) + '.h5'

    if not os.path.exists(auxiliary_path):
        raise ConfigError("%s, the target of this script does not seem to be where it should have been :/" % auxiliary_path)

    fp = h5py.File(auxiliary_path, 'r')

    contigs_db.create_table(nt_position_info_table_name, nt_position_info_table_structure, nt_position_info_table_types)

    contig_names_in_db = list(fp['/data/nt_position_info'].keys())

    run.info("Auxiliary data file found", auxiliary_path)
    run.info("Contigs found", len(contig_names_in_db))

    progress.new('Processing the auxiliary data file')
    counter, total = 0, len(contig_names_in_db)

    entries = []
    for contig_name in contig_names_in_db:
        entries.append((contig_name, convert_numpy_array_to_binary_blob(fp['/data/nt_position_info/%s' % (contig_name)].value),))

        counter += 1
        progress.update('contig %d of %d ...' % (counter, total))

        if counter % 10 == 0:
            progress.update("Writing buffer to the new database ...")
            contigs_db.insert_many(nt_position_info_table_name, entries=entries)
            entries = []

    contigs_db.insert_many(nt_position_info_table_name, entries=entries)

    progress.end()
    fp.close()

    # we also want to upgrade this table name, which was renamed within #654 re:
    # merenlab/pc_to_gene_cluster PR
    contigs_db._exec('ALTER TABLE gene_protein_sequences RENAME TO gene_amino_acid_sequences;')

    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)
    contigs_db.disconnect()

    os.remove(auxiliary_path)

    run.info_single('Done! Your contigs db is now %s, and the now-obsolete ".h5" file is gone!' % next_version, nl_after=1, nl_before=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db and CONTIGS.h5 from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    parser.add_argument('--just-do-it', default=False, action="store_true", help = "Do not bother me with warnings")
    args, unknown = parser.parse_known_args()

    try:
        update_contigs_db(args.contigs_db, just_do_it = args.just_do_it)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
