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


metadata_table_structure = ['entry_id', 'sample_id', 'name', 'length', 'GC_content', 'std_coverage', 'mean_coverage', 'normalized_coverage', 'portion_covered', 'abundance', 'variability', '__parent__']
metadata_table_mapping   = [    int   ,     str    ,  str  ,    int  ,     float   ,     float     ,      float     ,         float        ,       float      ,    float   ,     float    ,      str    ]
metadata_table_types     = [ 'numeric',   'text'   ,'text' ,'numeric',   'numeric' ,   'numeric'   ,    'numeric'   ,       'numeric'      ,     'numeric'    ,  'numeric' ,   'numeric'  ,     'text'  ]


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
        self.entry_id_c = 0
        self.entry_id_s = 0


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

            self.metadata_contigs[self.entry_id_c] = {'name': contig.name, 'sample_id': sample_id}
            for metadata_field in metadata_table_structure[3:]:
                self.metadata_contigs[self.entry_id_c][metadata_field] = contig_metadata[metadata_field]

            self.tnf_contigs[self.entry_id_c] = contig.tnf
            self.tnf_contigs[self.entry_id_c]['name'] = contig.name
            self.tnf_contigs[self.entry_id_c]['sample_id'] = sample_id

            self.entry_id_c += 1

            # contig is done, deal with splits in it:
            for split in contig.splits:
                split_metadata = split.get_metadata_dict()
                self.metadata_splits[self.entry_id_s] = {'name': split.name, 'sample_id': sample_id}
                for metadata_field in metadata_table_structure[3:]:
                    self.metadata_splits[self.entry_id_s][metadata_field] = split_metadata[metadata_field]

                self.tnf_splits[self.entry_id_s] = split.tnf
                self.tnf_splits[self.entry_id_s]['name'] = split.name
                self.tnf_splits[self.entry_id_s]['sample_id'] = sample_id

                self.entry_id_s += 1


        # all objects are ready, creating tables next.
        self.progress.update("Writing into the database ...")
        db.create_table('metadata_splits', metadata_table_structure, metadata_table_types)
        db_entries = [tuple([entry_id] + [self.metadata_splits[entry_id][h] for h in metadata_table_structure[1:]]) for entry_id in self.metadata_splits]
        db._exec_many('''INSERT INTO metadata_splits VALUES (?,?,?,?,?,?,?,?,?,?,?,?)''', db_entries)

        db.create_table('metadata_contigs', metadata_table_structure, metadata_table_types)
        db_entries = [tuple([entry_id] + [self.metadata_contigs[entry_id][h] for h in metadata_table_structure[1:]]) for entry_id in self.metadata_contigs]
        db._exec_many('''INSERT INTO metadata_contigs VALUES (?,?,?,?,?,?,?,?,?,?,?,?)''', db_entries)

        # FIXME: This a pretty shitty way to do this; table structure is dynamic, but not accessible from other
        # modules:
        tnf_table_structure = ['entry_id', 'sample_id', 'name'] + kmers
        tnf_table_types = ['text'] * len(tnf_table_structure)

        for t, v in [('tnf_contigs', self.tnf_contigs), ('tnf_splits', self.tnf_splits)]:
            db.create_table(t, tnf_table_structure, tnf_table_types)
            db_entries = [tuple([entry_id] + [v[entry_id][h] for h in tnf_table_structure[1:]]) for entry_id in v]
            db._exec_many('''INSERT INTO %s VALUES (%s)''' % (t, (','.join(['?'] * len(tnf_table_structure)))), db_entries)

        db.commit()
        self.progress.end()
