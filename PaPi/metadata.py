# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.


'''Storing and retrieving metadata regarding contigs and splits'''


metadata_table_structure = ['contig', 'length', 'GC_content', 'std_coverage', 'mean_coverage', 'normalized_coverage', 'relative_abundance', 'portion_covered', 'abundance', 'variability', '__parent__']
metadata_table_types     = [ 'text' ,'numeric',   'numeric' ,   'numeric'   ,    'numeric'   ,       'numeric'      ,       'numeric'     ,     'numeric'    ,  'numeric' ,   'numeric'  ,    'text'   ]


import copy
import numpy

from PaPi.terminal import Progress
from PaPi.terminal import pretty_print as pp

# Mock progress object that will not report anything, for general clarity.
progress = Progress()
progress.verbose = False


class Metadata:
    def __init__(self, p=progress):
        self.metadata_contigs = {}
        self.metadata_splits = {}
        self.tnf_contigs = {}
        self.tnf_splits = {}
        self.progress = p


    def store_metadata_for_contigs_and_splits(self, sample_id, contigs, db):
        self.progress.new('Storing metadata')

        num_contigs = pp(len(contigs))
        cur_contig = 1

        kmers = sorted(contigs.values()[0].tnf.keys())

        # this loop will get both metadata and tnf information form Contig instanes and store them into the db
        # at once. this was broken down into about 10 functions, but this structure seems to be the most efficient
        # although it looks crappy:
        for contig_name in contigs:
            self.progress.update("Processing contig %s of %s" % (pp(cur_contig), num_contigs))
            contig = contigs[contig_name]
            contig_metadata = contig.get_metadata_dict()
            contig_tnf = contig.tnf

            self.metadata_contigs[contig.name] = {'contig': contig.name}
            for metadata_field in metadata_table_structure[1:]:
                self.metadata_contigs[contig.name][metadata_field] = contig_metadata[metadata_field]

            self.tnf_contigs[contig.name] = contig.tnf
            self.tnf_contigs[contig.name]['contig'] = contig.name

            # contig is done, deal with splits in it:
            for split in contig.splits:
                split_metadata = split.get_metadata_dict()
                self.metadata_splits[split.name] = {'contig': split.name}
                for metadata_field in metadata_table_structure[1:]:
                    self.metadata_splits[split.name][metadata_field] = split_metadata[metadata_field]

                self.tnf_splits[split.name] = split.tnf
                self.tnf_splits[split.name]['name'] = split.name
                self.tnf_splits[split.name]['__parent__'] = contig.name


        self.progress.update("Generating tables ...")
        gen_metadata_tables(self.metadata_splits, self.metadata_contigs, db)
        gen_tnf_tables(kmers, self.tnf_splits, self.tnf_contigs, db)
        self.progress.end()



def gen_metadata_tables(metadata_splits, metadata_contigs, db):
    # all objects are ready, creating tables next.
    db.create_table('metadata_splits', metadata_table_structure, metadata_table_types)
    db_entries = [tuple([split] + [metadata_splits[split][h] for h in metadata_table_structure[1:]]) for split in metadata_splits]
    db._exec_many('''INSERT INTO metadata_splits VALUES (?,?,?,?,?,?,?,?,?,?,?)''', db_entries)

    db.create_table('metadata_contigs', metadata_table_structure, metadata_table_types)
    db_entries = [tuple([split] + [metadata_contigs[metadata_splits[split]['__parent__']][h] for h in metadata_table_structure[1:]]) for split in metadata_splits]
    db._exec_many('''INSERT INTO metadata_contigs VALUES (?,?,?,?,?,?,?,?,?,?,?)''', db_entries)

    db.commit()


def gen_tnf_tables(kmers, tnf_splits, tnf_contigs, db):
    tnf_table_structure = ['contig'] + kmers
    tnf_table_types = ['text'] + ['numeric'] * len(kmers)

    db.create_table('tnf_splits', tnf_table_structure, tnf_table_types)
    db_entries = [tuple([split] + [tnf_splits[split][h] for h in tnf_table_structure[1:]]) for split in tnf_splits]
    db._exec_many('''INSERT INTO %s VALUES (%s)''' % ('tnf_splits', (','.join(['?'] * len(tnf_table_structure)))), db_entries)

    db.create_table('tnf_contigs', tnf_table_structure, tnf_table_types)
    db_entries = [tuple([split] + [tnf_contigs[tnf_splits[split]['__parent__']][h] for h in tnf_table_structure[1:]]) for split in tnf_splits]
    db._exec_many('''INSERT INTO %s VALUES (%s)''' % ('tnf_contigs', (','.join(['?'] * len(tnf_table_structure)))), db_entries)

    db.commit()

