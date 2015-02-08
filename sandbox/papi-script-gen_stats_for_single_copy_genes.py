#!/usr/bin/env python
# -*- coding: utf-8

"""
Copyright (C) 2015, PaPi Authors

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.

Please read the COPYING file.
"""


import os
import sys
import argparse
from collections import Counter

import PaPi.fastalib as u
import PaPi.utils as utils
import PaPi.terminal as terminal 
import PaPi.annotation as annotation

from PaPi.utils import ConfigError
from PaPi.filesnpaths import FilesNPathsError

run = terminal.Run()
progress = terminal.Progress()

parser = argparse.ArgumentParser(description='A simple script to generate info from search tables')
parser.add_argument('contigs_fasta', metavar = 'FASTA',
                    help = 'A FASTA file with one or more contigs.')
parser.add_argument('annotation_db', metavar = 'ANNOTATION_DB',
                    help = 'Annotation database to read from.')
parser.add_argument('--list-sources', action='store_true', default=False,
                    help = 'Show available single-copy gene search results and exit.')
parser.add_argument('--source', default=None,
                    help = 'Source to focus on. If none declared, all single-copy gene sources\
                            are going to be listed.')

args = parser.parse_args()

contigs = set([])
contig_lengths = {}
contig_genes = {}
genes = {}

db = annotation.AnnotationDB(args.annotation_db)
search_contigs_dict = db.db.get_table_as_dict('search_contigs')
search_info_dict = db.db.get_table_as_dict('search_info')

sources = {}
for source in search_info_dict:
    if search_info_dict[source]['search_type'] == "singlecopy":
        sources[source] = [g.strip() for g in search_info_dict[source]['genes'].split(',')]

if args.list_sources:
    print sources.keys()
    sys.exit()

if args.source:
    if args.source not in sources:
        print 'bad --source. here are the available ones: %s' % ', '.join(sources)
    else:
        sources = {args.source: sources[source]}


fasta = u.SequenceSource(args.contigs_fasta)
while fasta.next():
    contigs.add(fasta.id)
    contig_lengths[fasta.id] = len(fasta.seq) - fasta.seq.count('N')

hits_output = open(args.annotation_db + '.hits', 'w')
hits_output.write('source\tcontig\tgene\te_value\n')
for entry in search_contigs_dict.values():
    hits_output.write('%s\t%s\t%s\t%.4g\n' % (entry['source'], entry['contig'], entry['gene_name'], entry['e_value'] ))
hits_output.close()

genes_output = open(args.annotation_db + '.genes', 'w')
genes_output.write('source\tgene\n')
for source in sources:
    for gene in sources[source]:
        genes_output.write('%s\t%s\n' % (source, gene))
genes_output.close()

#for e in es:
#    print
#    print '##################################################'
#    print e
#    print '##################################################'
#    print
#    for contig in contigs:
#        if contig_genes[e].has_key(contig) and len(contig_genes[e][contig]):
#            print contig, contig_lengths[contig], len(contig_genes[e][contig])
#    print
