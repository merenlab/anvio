#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse

import anvio.db as db

import anvio.terminal as terminal

from anvio.errors import ConfigError
from anvio.dbinfo import is_contigs_db

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

run = terminal.Run()
progress = terminal.Progress()

# tables that are specific to contigs databases that require updatin'
tables = {
    'genes_in_splits':            {'structure': ['split', 'gene_callers_id', 'start_in_split', 'stop_in_split', 'percentage_in_split'],
                                       'types': ['text', 'numeric', 'numeric', 'numeric', 'numeric']},
    'gene_functions':             {'structure': ['gene_callers_id', 'source', 'accession', 'function', 'e_value'],
                                       'types': ['numeric', 'text', 'text', 'text', 'numeric']},
    'hmm_hits_in_splits':         {'structure': ['hmm_hit_entry_id', 'split', 'percentage_in_split', 'source'],
                                       'types': ['numeric', 'text', 'numeric', 'text']},
    'scg_taxonomy':               {'structure': ['gene_callers_id', 'gene_name', 'accession', 'percent_identity', 't_domain', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"],
                                       'types': ['numeric', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text']},
    'nucleotide_additional_data': {'structure': ['item_name', 'data_key', 'data_value', 'data_type', 'data_group'],
                                       'types': ['text', 'text', 'text', 'text', 'text']},
    'amino_acid_additional_data': {'structure': ['item_name', 'data_key', 'data_value', 'data_type', 'data_group'],
                                       'types': ['text', 'text', 'text', 'text', 'text']},
    'collections_bins_info':      {'structure': ['collection_name', 'bin_name', 'source', 'html_color'],
                                       'types': ['text', 'text', 'text', 'text']},
    'collections_of_contigs':     {'structure': ['collection_name', 'contig', 'bin_name'],
                                       'types': ['text', 'text', 'text']},
    'collections_of_splits':      {'structure': ['collection_name', 'split', 'bin_name'],
                                       'types': ['text', 'text', 'text']},
}


def drop_entry_id_column_from_table(db_path, table_name, table_properties):
    progress.new("Modifying '%s'" % table_name)

    structure = table_properties['structure']
    types = table_properties['types']
    db_fields = ', '.join(['%s %s' % (t[0], t[1]) for t in zip(structure, types)])
    temp_table_name = table_name + '_TEMP'

    _db = db.DB(db_path, None, ignore_version = True)

    progress.update("Creating a temporary table")
    _db._exec('''CREATE TABLE %s (%s)''' % (temp_table_name, db_fields))

    progress.update("Copying data into the temporary table")
    _db._exec('''INSERT INTO %s SELECT %s FROM %s''' % (temp_table_name, ', '.join(structure), table_name))

    progress.update("Dropping the original table")
    _db._exec('''DROP TABLE IF EXISTS %s''' % (table_name))

    progress.update("Renaming temporary table to become the original")
    _db._exec('''ALTER TABLE %s RENAME TO %s''' % (temp_table_name, table_name))

    progress.update("Committing changes")
    _db.disconnect()

    progress.end()


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        contigs_db.disconnect()
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)
    contigs_db.disconnect()

    # drop entry ids one by one
    for table_name in tables:
        drop_entry_id_column_from_table(db_path, table_name, table_properties=tables[table_name])

    contigs_db = db.DB(db_path, None, ignore_version = True)
    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)
    contigs_db.disconnect()

    run.info_single("Your contigs db is now %s. This update carried one more issue into the graveyard "
                    "of bad design decisions we've made years ago by altering %d tables in your database." \
                            % (next_version, len(tables)), nl_after=1, nl_before=1, mc='green')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
