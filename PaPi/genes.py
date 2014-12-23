# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

'''The purpose of this file is to provide a class to keep coverage values for each gene in contigs for a sample.
   Simply, you create an instance from it, keep sending contig instances from contig.py::Contig class along with
   a list of inferred start/stop locations for each reading frame. Once you are done, you export that information.'''


genes_table_structure = ['entry_id', 'prot', 'sample_id', 'mean_coverage']
genes_table_mapping   = [    int   ,  str  ,    str     ,      float     ]
genes_table_types     = [ 'numeric', 'text',  'text'    ,    'numeric'   ]


import copy
import numpy


from PaPi.terminal import Progress

# Mock progress object that will not report anything, for general clarity.
progress = Progress()
progress.verbose = False


class Genes:
    def __init__(self, p=progress):
        self.genes = {}
        self.progress = p
        self.entry_id = 0

        # we keep coverage values in contig.py/Contig instances only for splits, during the profiling,
        # coverage for contigs are temporarily calculated, and then discarded. probably that behavior
        # should change for good. but for now I will generate a dict to keep contig coverages to avoid
        # even more redundant computations:
        self.contig_coverages = {}


    def analyze_contig(self, contig, sample_id, start_stop_pos_list):
        if contig.name not in self.contig_coverages:
            contig_coverage = []
            for split in contig.splits:
                contig_coverage.extend(split.coverage.c)
            self.contig_coverages[contig.name] = contig_coverage

        for prot, start, stop in start_stop_pos_list:
            gene_coverage = numpy.mean(self.contig_coverages[contig.name][start:stop])
            self.add_gene_entry(prot, sample_id, gene_coverage)


    def add_gene_entry(self, prot, sample_id, coverage):
            self.genes[self.entry_id] = {'prot': prot,
                                         'sample_id': sample_id,
                                         'mean_coverage': coverage}
            self.entry_id += 1


    def create_genes_table(self, db):
        db.create_table('genes', genes_table_structure, genes_table_types)
        db_entries = [tuple([entry_id] + [self.genes[entry_id][h] for h in genes_table_structure[1:]]) for entry_id in self.genes]
        db._exec_many('''INSERT INTO genes VALUES (?,?,?,?)''', db_entries)
        db.commit()