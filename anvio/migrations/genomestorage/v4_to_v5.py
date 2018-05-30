#!/usr/bin/env python
# -*- coding: utf-8

import os
import sys
import time
import argparse

import anvio.db as db
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

with terminal.SuppressAllOutput():
    import h5py

current_version = '4'
next_version    = '5'

run = terminal.Run()
progress = terminal.Progress()

essential_genome_info = ['gc_content', 'num_contigs', 'num_splits', 'total_length', 'num_genes', 'percent_completion', 'percent_redundancy',
                         'genes_are_called', 'avg_gene_length', 'num_genes_per_kb', ]

gene_function_calls_table_name         = 'gene_functions'
gene_function_calls_table_structure    = ['entry_id', 'gene_callers_id', 'source', 'accession', 'function', 'e_value']
gene_function_calls_table_types        = [ 'numeric',     'numeric'    ,  'text' ,    'text'  ,   'text'  , 'numeric']

genome_info_table_name       = 'genome_info'
genome_info_table_structure  = ['genome_name', 'genome_hash', 'external_genome'] + essential_genome_info
genome_info_table_types      = [    'str'    ,      'str'   ,     'numeric'    ] + ['numeric'] * len(essential_genome_info)

gene_info_table_name       = 'gene_info'
gene_info_table_structure  = ['genome_name', 'gene_caller_id', 'aa_sequence', 'dna_sequence', 'partial', 'length' ]
gene_info_table_types      = [    'str'    ,     'numeric'   ,    'text'    ,     'text'    , 'numeric', 'numeric']

genome_gene_function_calls_table_name      = 'gene_function_calls'
genome_gene_function_calls_table_structure = ['genome_name', ] + gene_function_calls_table_structure[:]
genome_gene_function_calls_table_types     = [    'str'    , ] + gene_function_calls_table_types[:]


def migrate(db_path):
    if db_path is None:
        raise ConfigError("No database path is given.")

    fp = h5py.File(db_path, 'r')

    if int(fp.attrs['version']) != int(current_version):
        fp.close()
        raise ConfigError("Genome storage version is not %s." % current_version)

    old_storage_hash = str(fp.attrs['hash'])
    functions_are_available = fp.attrs['functions_are_available']

    run.info("Outdated genomes storage found (%s)" % old_storage_hash, db_path)

    genome_storage_db_path = db_path[:-3] + '.db'
    filesnpaths.is_output_file_writable(genome_storage_db_path, ok_if_exists=False)

    genomes_db = db.DB(genome_storage_db_path, next_version, new_database=True)
    genomes_db.create_table(genome_info_table_name, genome_info_table_structure, genome_info_table_types)
    genomes_db.create_table(gene_info_table_name, gene_info_table_structure, gene_info_table_types)
    genomes_db.create_table(genome_gene_function_calls_table_name, genome_gene_function_calls_table_structure, genome_gene_function_calls_table_types)
    genomes_db._exec("CREATE INDEX covering_index ON %s (gene_callers_id, genome_name);" % genome_gene_function_calls_table_name)

    genomes_db.set_meta_value('db_type', 'genomestorage')
    genomes_db.set_meta_value('creation_date', time.time())
    genomes_db.set_meta_value('hash', old_storage_hash)
    genomes_db.set_meta_value('functions_are_available', functions_are_available)

    I = lambda genome_name, key: fp['/info/genomes/%s/%s' % (genome_name, key)]

    genome_names = [d for d in fp['/info/genomes']]

    progress.new('Bleep bloop')
    progress.update('Adding genomes')
    genome_info_entries = []
    for genome_name in genome_names:
        values = (genome_name, )

        for column_name in genome_info_table_structure[1:]:
            # dirty workaround for backwards compatibility,
            # "percent_completion" may be "percent_complete" in some old genome storages, 
            # because ozcan forgot to add that into upgrade script :(
            if column_name == 'percent_completion' and '/info/genomes/%s/percent_completion' % genome_name not in fp:
                column_name = 'percent_complete'

            attr = I(genome_name, column_name)

            if attr.dtype == 'int64':
                values += (int(attr.value), )
            elif attr.dtype == 'float64':
                values += (float(attr.value), )
            else:
                values += ((attr.value), )

        genome_info_entries.append(values)
    genomes_db.insert_many(genome_info_table_name, entries=genome_info_entries)
    del genome_info_entries

    progress.update('Adding genes')
    gene_entries = []
    for genome_name in genome_names:
        for gene_callers_id in fp['/data/genomes/%s' % genome_name]:
            G = lambda key: fp['/data/genomes/%s/%s/%s' % (genome_name, gene_callers_id, key)].value
            gene_entries.append((genome_name, gene_callers_id, G('aa_sequence'), G('dna_sequence'), int(G('partial')), int(G('length')), ))
    genomes_db.insert_many(gene_info_table_name, entries=gene_entries)
    del gene_entries

    progress.update('Adding functions')
    functions_entries = []
    entry_id_counter = 0
    for genome_name in genome_names:
        for gene_callers_id in fp['/data/genomes/%s' % genome_name]:
            functions_path = '/data/genomes/%s/%s/functions' % (genome_name, gene_callers_id)
            if functions_path in fp:
                for source in fp[functions_path]:
                    annotation_list = str(fp['/data/genomes/%s/%s/functions/%s' % (genome_name, gene_callers_id, source)].value).split('|||')

                    functions_entries.append((genome_name, entry_id_counter, gene_callers_id, source, annotation_list[0], annotation_list[1], 0, ))
                    entry_id_counter += 1
    genomes_db.insert_many(genome_gene_function_calls_table_name, entries=functions_entries)

    genomes_db.disconnect()

    progress.end()

    os.remove(db_path)

    run.info_single("Your genomes storage is now at version %s. The new on is at %s, and anvi'o just removed\
                     the old one, which was at %s from your disk." % (next_version, genome_storage_db_path, db_path), nl_after=1, nl_before=1, mc='green')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A simple script to upgrade genomes storage from version %s to version %s' % (current_version, next_version))
    parser.add_argument('genomes_storage', metavar = 'GENOMES_STORAGE', help = "An anvi'o genomes storage of version %s" % current_version)
    args, unknown = parser.parse_known_args()

    try:
        migrate(args.genomes_storage)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
