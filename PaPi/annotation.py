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

levels_of_taxonomy = ["t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
mapping_taxonomy   = [   str    ,   str    ,    str   ,    str    ,    str   ,     str    ]
db_types_taxonomy  = [  'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

header    = ['prot', 'contig', 'start',  'end'   , 'direction', 'figfam', 'function']
mapping   = [ str  ,   str   ,  int   ,   int    ,     str    ,    str  ,    str    ]
db_types  = ['text',  'text' ,'numeric','numeric',    'text'  ,  'text' ,   'text'  ]

header.extend(levels_of_taxonomy)
mapping.extend(mapping_taxonomy)
db_types.extend(db_types_taxonomy)

__version__ = "0.0.2"

import os
import sys
import operator
from collections import Counter

import PaPi.db as db
import PaPi.utils as utils
import PaPi.dictio as dictio
import PaPi.filesnpaths as filesnpaths
from PaPi.utils import ConfigError


class Annotation:
    def __init__(self, db_path):
        self.db_path = db_path
        self.version = None
        self.db = None


    def init_database_from_matrix(self, source):
        if type(source) == type(dict()):
            self.matrix_dict = source
            self.check_keys(['prot'] + self.matrix_dict.values()[0].keys())
        if type(source) == type(str()):
            self.matrix_dict = utils.get_TAB_delimited_file_as_dictionary(source,
                                                                          column_names = header,
                                                                          column_maping = mapping)

        self.db = db.DB(self.db_path, __version__, new_database = True)

        db_fields = ', '.join(['%s %s' % (t[0], t[1]) for t in zip(header, db_types)])
        self.db._exec('''CREATE TABLE annotation (%s)''' % db_fields)

        db_entries = [tuple([prot] + [self.matrix_dict[prot][h] for h in header[1:]]) for prot in self.matrix_dict]
        self.db._exec_many('''INSERT INTO annotation VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', db_entries)

        self.db.disconnect()


    def init_database(self):
        self.db = db.DB(self.db_path, __version__)


    def check_keys(self, keys):
        missing_keys = [key for key in header if key not in keys]
        if len(missing_keys):
            raise ConfigError, "Your input lacks one or more header fields to generate a PaPi annotation db. Here is\
                                what you are missing: %s. The complete list (and order) of headers in your TAB\
                                delimited matrix file (or dictionary) must follow this: %s." % (', '.join(missing_keys),
                                                                                                ', '.join(header))


    def get_consensus_taxonomy_for_contig(self, contig, t_level = 't_species', start = 0, stop = sys.maxint):
        """Returns (c, n, t, o) where,
            c: consensus taxonomy (the most common taxonomic call for each gene found in the contig),
            n: total number of genes found in the contig,
            t: total number of genes with known taxonomy,
            o: number of taxonomic calls that matches the consensus among t
        """

        response = self.db._exec("""SELECT %s FROM annotation WHERE contig='%s'""" % (t_level, contig, ))

        rows = response.fetchall()
        num_genes = len(rows)
        tax_str_list = [t[0] for t in rows if t[0]]
        distinct_taxa = set(tax_str_list)

        if not len(distinct_taxa):
            return None, num_genes, None

        if len(distinct_taxa) == 1:
            return distinct_taxa.pop(), num_genes, len(tax_str_list), len(tax_str_list)
        else:
            d = Counter()
            for t in tax_str_list:
                d[t] += 1
            consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
            return consensus, num_genes, len(tax_str_list), occurrence
