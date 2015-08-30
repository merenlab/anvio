#!/usr/bin/env python
# -*- coding: utf-8

"""
Get sequences for a list of proteins.
"""

import textwrap
import argparse

import anvio.tables as t
import anvio.terminal as terminal 
import anvio.dbops as dbops

progress = terminal.Progress()

parser = argparse.ArgumentParser(description='A simple script to print out sequences for a given list of contigs database protein ids')
parser.add_argument('contigs_db', metavar = 'CONTIGS_DB',
                    help = 'Contigs database to read from.')
parser.add_argument('genes_list', metavar = 'PROT_IDs',
                    help = 'Protein IDs.')

args = parser.parse_args()

gene_ids = set([p.strip() for p in open(args.genes_list).readlines()])

# {'function': 'Cobalamin biosynthesis protein CbiG''direction': 'f', 't_phylum': '', 'figfam': '', 'stop': 5202, 't_order': '', 'start': 4162, 't_species': 'Ruminococcus sp.', 't_class': '', 'contig': '204_10M_MERGED.PERFECT.gz.keep_contig_6515', 't_family': ''}

db = dbops.ContigsDatabase(args.contigs_db, quiet=False)
contig_sequences = db.db.get_table_as_dict(t.contig_sequences_table_name)
genes_in_contigs = db.db.get_table_as_dict(t.genes_contigs_table_name)
db.disconnect()

for gene_id in gene_ids:
    g = genes_in_contigs[gene_id]
    g_sequence = contig_sequences[g['contig']]['sequence'][g['start']:g['stop']]
    print '>%s|%s' % (gene_id, '|'.join('%s:%s' % (k, g[k]) for k in ['contig', 'start', 'stop', 'direction', 'function', 't_species']))
    print textwrap.fill(g_sequence, 80)
