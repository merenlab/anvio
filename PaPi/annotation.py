# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

"""
    Here is described the annotation class. Any parser implemented in parsers.py must generate headers
    matching 'header' variable.
"""

annotation_table_structure = ['prot', 'contig', 'start', 'stop'   , 'direction', 'figfam', 'function', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
annotation_table_mapping   = [ str  ,   str   ,  int   ,   int    ,     str    ,   str   ,    str    ,    str    ,   str    ,    str   ,    str    ,    str   ,     str    ]
annotation_table_types     = ['text',  'text' ,'numeric','numeric',   'text'   ,  'text' ,   'text'  ,   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

splits_table_structure = ['split', 'taxonomy', 'num_genes', 'num_tax_calls', 'num_function_calls', 'tax_accuracy']
splits_table_mapping   = [  str  ,     str   ,    int     ,      int       ,          int        ,      float    ]
splits_table_types     = [ 'text',   'text'  ,  'numeric' ,   'numeric'    ,       'numeric'     ,    'numeric'  ]

__version__ = "0.0.1"

import os
import sys
import operator
from collections import Counter

import PaPi.db as db
import PaPi.fastalib as u
import PaPi.utils as utils
import PaPi.dictio as dictio
import PaPi.filesnpaths as filesnpaths
from PaPi.utils import ConfigError
from PaPi.contig import Split


class Annotation:
    def __init__(self, db_path):
        self.db_path = db_path
        self.db = None

        # a dictionary to keep contigs
        self.contigs = {}
        self.split_length = None


    def create_new_database(self, contigs_fasta, source, split_length, parser="unknown"):
        if type(source) == type(dict()):
            self.matrix_dict = source
            self.check_keys(['prot'] + self.matrix_dict.values()[0].keys())
        if type(source) == type(str()):
            self.matrix_dict = utils.get_TAB_delimited_file_as_dictionary(source,
                                                                          column_names = annotation_table_structure,
                                                                          column_maping = annotation_table_mapping)

        # populate contigs dict with contig lengths
        fasta = u.SequenceSource(contigs_fasta)
        while fasta.next():
            self.contigs[fasta.id] = {'length': len(fasta.seq)}

        self.sanity_check()

        # init a new db
        self.db = db.DB(self.db_path, __version__, new_database = True)

        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)
        self.db.set_meta_value('annotation_source', parser)
        self.split_length = split_length

        # create annotation main table using fields in 'annotation_table_structure' variable
        self.db.create_table('annotation', annotation_table_structure, annotation_table_types)

        # push raw entries.
        db_entries = [tuple([prot] + [self.matrix_dict[prot][h] for h in annotation_table_structure[1:]]) for prot in self.matrix_dict]
        self.db._exec_many('''INSERT INTO annotation VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', db_entries)

        # compute and push split taxonomy information.
        self.init_splits_table()

        # bye.
        self.db.disconnect()


    def sanity_check(self):
        contig_names_in_matrix = set([v['contig'] for v in self.matrix_dict.values()])
        contig_names_in_fasta  = set(self.contigs.keys())

        for contig in contig_names_in_matrix:
            if contig not in contig_names_in_fasta:
                raise ConfigError, "We have a problem... Every contig name there is in your input files for annotation\
                                    must be found in the FASTA file. But it seems it is not the case. I did not check\
                                    all, but there there is at least one contig name ('%s') that appears in the\
                                    annotation files, but missing in the FASTA file. You may need to format the\
                                    names in your FASTA file to match contig names in your annotation files. Keep in\
                                    mind that contig names must match the ones in your BAM files later on. Even when\
                                    you use one software for assembly and mapping, disagreements between contig names\
                                    may arise. We know that it is the case with CLC, for instance. OK. Going back to the\
                                    issue. Here is one contig name from your FASTA file: '%s', and here is one from your\
                                    input files: '%s'. You should make them identical (and make sure whatever solution\
                                    you come up with will not make them incompatible with names in your BAM files\
                                    later on." % (contig, contig_names_in_fasta.pop(), contig_names_in_matrix.pop())


    def init_database(self):
        self.db = db.DB(self.db_path, __version__)


    def check_keys(self, keys):
        missing_keys = [key for key in annotation_table_structure if key not in keys]
        if len(missing_keys):
            raise ConfigError, "Your input lacks one or more header fields to generate a PaPi annotation db. Here is\
                                what you are missing: %s. The complete list (and order) of headers in your TAB\
                                delimited matrix file (or dictionary) must follow this: %s." % (', '.join(missing_keys),
                                                                                                ', '.join(annotation_table_structure))


    def init_splits_table(self):
        # build a dictionary for fast access to all proteins identified within a contig
        prots_in_contig = {}
        for prot in self.matrix_dict:
            contig = self.matrix_dict[prot]['contig']
            if prots_in_contig.has_key(contig):
                prots_in_contig[contig].add(prot)
            else:
                prots_in_contig[contig] = set([prot])

        splits_dict = {}
        for contig in self.contigs:
            chunks = utils.get_chunks(self.contigs[contig]['length'], self.split_length)
            for i in range(0, len(chunks)):
                split = Split(contig, i).name
                start = chunks[i][0]
                stop = chunks[i][1]

                taxa = []
                functions = []
                for prot in prots_in_contig[contig]:
                    if self.matrix_dict[prot]['stop'] > start and self.matrix_dict[prot]['start'] < stop:
                        taxa.append(self.matrix_dict[prot]['t_species'])
                        functions.append(self.matrix_dict[prot]['function'])

                taxonomy_strings = [t for t in taxa if t]
                function_strings = [f for f in functions if f]

                splits_dict[split] = {'taxonomy': None,
                                      'num_genes': len(taxa),
                                      'num_tax_calls': len(taxonomy_strings),
                                      'num_function_calls': len(function_strings),
                                      'tax_accuracy': 0.0}

                distinct_taxa = set(taxonomy_strings)

                if not len(distinct_taxa):
                    continue

                if len(distinct_taxa) == 1:
                    splits_dict[split]['taxonomy'] = distinct_taxa.pop()
                    splits_dict[split]['tax_accuracy'] = 1.0
                else:
                    d = Counter()
                    for taxon in taxonomy_strings:
                        d[taxon] += 1
                    consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
                    splits_dict[split]['taxonomy'] = distinct_taxa.pop()
                    splits_dict[split]['tax_accuracy'] = 1.0

        # create splits table using fields in 'splits_table_structure' info
        self.db.create_table('splits', splits_table_structure, splits_table_types)

        # push raw entries.
        db_entries = [tuple([split] + [splits_dict[split][h] for h in splits_table_structure[1:]]) for split in splits_dict]
        self.db._exec_many('''INSERT INTO splits VALUES (?,?,?,?,?,?)''', db_entries)


    def get_consensus_taxonomy_for_split(self, contig, t_level = 't_species', start = 0, stop = sys.maxint):
        """Returns (c, n, t, o) where,
            c: consensus taxonomy (the most common taxonomic call for each gene found in the contig),
            n: total number of genes found in the contig,
            t: total number of genes with known taxonomy,
            o: number of taxonomic calls that matches the consensus among t
        """

        response = self.db.cursor.execute("""SELECT %s FROM annotation WHERE contig='%s' and stop > %d and start < %d""" % (t_level, contig, start, stop))
        rows = response.fetchall()

        num_genes = len(rows)
        tax_str_list = [t[0] for t in rows if t[0]]
        distinct_taxa = set(tax_str_list)

        if not len(distinct_taxa):
            return None, num_genes, 0, 0

        if len(distinct_taxa) == 1:
            return distinct_taxa.pop(), num_genes, len(tax_str_list), len(tax_str_list)
        else:
            d = Counter()
            for t in tax_str_list:
                d[t] += 1
            consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
            return consensus, num_genes, len(tax_str_list), occurrence
