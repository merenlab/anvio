#!/usr/bin/env python
# -*- coding: utf-8

"""
Goes through the collections_* tables and creates a TAB delimited matrix based on
all available collections in a given contigs database.
"""

import argparse

import anvio.terminal as terminal 
import anvio.contigs as contigs


run = terminal.Run()
progress = terminal.Progress()

parser = argparse.ArgumentParser(description='A simple script to generate info from search tables')
parser.add_argument('contigs_db', metavar = 'CONTIGS_DB',
                    help = 'Contigs database to read from.')
parser.add_argument('-o', '--output', metavar = "OUTPUT_FILE.txt", default="COLLECTIONS.txt",
                    help = 'Output file name.')

args = parser.parse_args()

contigs = set([])
contig_lengths = {}

db = contigs.ContigsDatabase(args.contigs_db, quiet=False)
collections_splits_table = db.db.get_table_as_dict(contigs.collections_splits_table_name)
collections_info_table = db.db.get_table_as_dict(contigs.collections_info_table_name)
contig_lengths_table = db.db.get_table_as_dict(contigs.contig_lengths_table_name)
contig_lengths = dict([(c, contig_lengths_table[c]['length']) for c in contig_lengths_table])
db.disconnect()

sources = collections_info_table.keys()

splits = {}
for entry in collections_splits_table.values():
    split = entry['split']
    source = entry['source']
    cluster_id = entry['cluster_id']

    if splits.has_key(split):
        splits[split][source] = cluster_id
    else:
        splits[split] = {source: cluster_id}


output = open(args.output, 'w')
output.write('split\t%s\n' % '\t'.join(sources))
for split in splits.keys():
    line = [split]
    for source in sources:
        if splits[split].has_key(source):
            line.append(splits[split][source])
        else:
            line.append(None)
    output.write('\t'.join(line) + '\n')
output.close()

run.info('Collections matrix', args.output)
