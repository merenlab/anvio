#!/usr/bin/env python
# -*- coding: utf-8

import sys
import argparse
import numpy as np

import anvio.db as db
import anvio.utils as utils

import anvio.terminal as terminal

from anvio.errors import ConfigError

current_version, next_version = [x[1:] for x in __name__.split('_to_')]

genes_in_contigs_table_structure = ['gene_callers_id', 'contig', 'start' , 'stop'  , 'direction', 'partial', 'call_type', 'source', 'version']
genes_in_contigs_table_types     = [    'numeric'    ,  'text' ,'numeric','numeric',   'text'   , 'numeric',  'numeric' ,  'text' ,   'text' ]

run = terminal.Run()
progress = terminal.Progress()

def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    utils.is_contigs_db(db_path)

    contigs_db = db.DB(db_path, None, ignore_version = True)
    if str(contigs_db.get_version()) != current_version:
        raise ConfigError("Version of this contigs database is not %s (hence, this script cannot really do anything)." % current_version)

    aa_seqs = contigs_db.get_table_as_dataframe('gene_amino_acid_sequences')
    aa_seqs = dict(zip(aa_seqs['gene_callers_id'], aa_seqs['sequence']))

    gene_calls = contigs_db.get_table_as_dataframe('genes_in_contigs')

    progress.new("Adding coding type to each gene call", progress_total_items=gene_calls.shape[0])

    coding_num = 1
    noncoding_num = 2

    call_types = coding_num * np.ones(gene_calls.shape[0]).astype(int)
    for i, gene_call in gene_calls.iterrows():
        if i % 5000 == 0:
            progress.update("done %d of %d" % (i, gene_calls.shape[0]), increment=True)

        if not len(aa_seqs.get(gene_call['gene_callers_id'], '')):
            call_types[i] = noncoding_num

    gene_calls['call_type'] = call_types

    progress.end()

    contigs_db.drop_table('genes_in_contigs')
    contigs_db.create_table('genes_in_contigs', genes_in_contigs_table_structure, genes_in_contigs_table_types)
    contigs_db.insert_rows_from_dataframe('genes_in_contigs', gene_calls)

    # -------------------------------------

    progress.new("Adding self table values")
    progress.update("...")

    contigs_db.remove_meta_key_value_pair('version')
    contigs_db.set_version(next_version)

    contigs_db.set_meta_value('external_gene_calls', 0)
    contigs_db.set_meta_value('external_gene_amino_acid_seqs', 0)
    contigs_db.set_meta_value('skip_predict_frame', 0)


    progress.update("Committing changes")
    contigs_db.disconnect()

    progress.end()
    run.info_single("genes_in_contigs now has 'coding_type' as a column", nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade CONTIGS.db from version %s to version %s' % (current_version, next_version))
    parser.add_argument('contigs_db', metavar = 'CONTIGS_DB', help = 'Contigs database at version %s' % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.contigs_db)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
