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
    Classes to create and access the annotation database.
"""

contig_sequences_table_name          = 'contig_sequences'
contig_sequences_table_structure     = ['contig', 'sequence']
contig_sequences_table_types         = [  'str' ,   'str'   ]

contig_lengths_table_name            = 'contig_lengths'
contig_lengths_table_structure       = ['contig', 'length' ]
contig_lengths_table_types           = [  'str' , 'numeric']

splits_info_table_name               = 'splits_info'
splits_info_table_structure          = ['split', 'order_in_parent' , 'start' ,  'end'  , 'parent' ]
splits_info_table_types              = ['text' ,     'numeric     ','numeric','numeric',  'text'  ]

genes_contigs_table_name             = 'genes_in_contigs'
genes_contigs_table_structure        = ['prot', 'contig', 'start', 'stop'   , 'direction', 'figfam', 'function', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
genes_contigs_table_types            = ['text',  'text' ,'numeric','numeric',   'text'   ,  'text' ,   'text'  ,   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

genes_splits_summary_table_name      = 'genes_in_splits_summary'
genes_splits_summary_table_structure = ['split', 'taxonomy', 'num_genes', 'avg_gene_length', 'ratio_coding', 'ratio_hypothetical', 'ratio_with_tax', 'tax_accuracy']
genes_splits_summary_table_types     = [ 'text',   'text'  ,  'numeric' ,     'numeric'    ,   'numeric'   ,      'numeric'      ,     'numeric'   ,   'numeric'   ]

genes_splits_table_name              = 'genes_in_splits'
genes_splits_table_structure         = ['entry_id', 'split', 'prot', 'start_in_split', 'stop_in_split', 'percentage_in_split']
genes_splits_table_types             = [ 'numeric',  'text', 'text',    'numeric'    ,    'numeric'   ,       'numeric'      ]

hmm_hits_info_table_name             = 'hmm_hits_info'
hmm_hits_info_table_structure        = ['source', 'ref' , 'search_type', 'genes']
hmm_hits_info_table_types            = [ 'text' , 'text',    'text'    , 'text' ]

hmm_hits_contigs_table_name          = 'hmm_hits_in_contigs'
hmm_hits_contigs_table_structure     = ['entry_id', 'source', 'contig', 'start' , 'stop'  , 'gene_name', 'gene_id', 'e_value']
hmm_hits_contigs_table_types         = [ 'numeric',  'text' ,  'text' ,'numeric','numeric',   'text'   ,  'text'  , 'numeric']

hmm_hits_splits_table_name           = 'hmm_hits_in_splits'
hmm_hits_splits_table_structure      = ['entry_id', 'source', 'gene_unique_identifier', 'gene_name', 'split', 'percentage_in_split', 'e_value']
hmm_hits_splits_table_types          = [ 'numeric',  'text' ,          'text'         ,   'text'   ,  'text',       'numeric'      , 'numeric']

collections_info_table_name          = 'collections_info'
collections_info_table_structure     = ['source',  'ref']
collections_info_table_types         = [ 'text' , 'text']

collections_contigs_table_name       = 'collections_of_contigs'
collections_contigs_table_structure  = ['entry_id', 'source', 'contig', 'cluster_id']
collections_contigs_table_types      = [ 'numeric',  'text' ,  'text' ,    'text'   ]

collections_splits_table_name        = 'collections_of_splits'
collections_splits_table_structure   = ['entry_id', 'source', 'split', 'cluster_id']
collections_splits_table_types       = [ 'numeric',  'text' , 'text' ,    'text'   ]


__version__ = "0.4.2"


import os
import sys
import numpy
import random
import hashlib
import operator
from collections import Counter

import PaPi.db as db
import PaPi.fastalib as u
import PaPi.utils as utils
import PaPi.contig as contig
import PaPi.dictio as dictio
import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths

from PaPi.utils import ConfigError
from PaPi.commandline import HMMSearch
from PaPi.parsers import parser_modules

run = terminal.Run()
progress = terminal.Progress()


class AnnotationDatabase:
    """To create an empty annotation database and/or access one."""
    def __init__(self, db_path, run=run, progress=progress, quiet = True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if os.path.exists(self.db_path):
            self.db = db.DB(self.db_path, __version__)

            self.run.info('Annotation database', 'An existing database, %s, has been initiated.' % self.db_path, quiet = self.quiet)
            self.run.info('Number of contigs', self.db.get_meta_value('num_contigs'), quiet = self.quiet)
            self.run.info('Total number of nucleotides', self.db.get_meta_value('total_length'), quiet = self.quiet)
            self.run.info('Split length', self.db.get_meta_value('split_length'), quiet = self.quiet)
        else:
            self.db = None


    def create(self, contigs_fasta, split_length):
        if os.path.exists(self.db_path):
            raise ConfigError, "PaPi will not overwrite an existing annotation database. Please choose a different name\
                                or remove the existing database ('%s') first." % (self.db_path)

        if not split_length:
            raise ConfigError, "Creating a new annotation database requires split length information to be\
                                provided. But the AnnotationDatabase class was called to create one without this\
                                bit of information. Not cool."

        if not os.path.exists(contigs_fasta):
            raise ConfigError, "Creating a new annotation database requires a FASTA file with contigs to be provided."


        if not self.db_path.lower().endswith('.db'):
            raise ConfigError, "Please make sure your output file name has a '.db' extension. PaPi developers apologize\
                                for imposing their views on how local databases should be named, and are humbled by your\
                                cooperation."

        try:
            split_length = int(split_length)
        except:
            raise ConfigError, "Split size must be an integer."

        self.db = db.DB(self.db_path, __version__, new_database = True)

        # this will be the unique information that will be passed downstream whenever this db is used:
        self.db.set_meta_value('annotation_hash', '%08x' % random.randrange(16**8))
        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)

        self.db.create_table(contig_sequences_table_name, contig_sequences_table_structure, contig_sequences_table_types)
        self.db.create_table(contig_lengths_table_name, contig_lengths_table_structure, contig_lengths_table_types)
        self.db.create_table(splits_info_table_name, splits_info_table_structure, splits_info_table_types)

        # lets process and store the FASTA file.
        fasta = u.SequenceSource(contigs_fasta)
        num_contigs, total_length = 0, 0
        db_entries_contig_sequences = []
        db_entries_contig_lengths = []
        db_entries_splits_info = []

        while fasta.next():
            num_contigs += 1
            contig_length = len(fasta.seq)
            chunks = utils.get_chunks(contig_length, split_length)

            for order in range(0, len(chunks)):
                start, end = chunks[order]
                db_entries_splits_info.append((contig.gen_split_name(fasta.id, order), order, start, end, fasta.id), )

            db_entries_contig_sequences.append((fasta.id, fasta.seq), )
            db_entries_contig_lengths.append((fasta.id, contig_length), )
            total_length += contig_length

        self.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % contig_sequences_table_name, db_entries_contig_sequences)
        self.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % contig_lengths_table_name, db_entries_contig_lengths)
        self.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % splits_info_table_name, db_entries_splits_info)

        # set some useful meta values:
        self.db.set_meta_value('num_contigs', num_contigs)
        self.db.set_meta_value('total_length', total_length)
        self.db.set_meta_value('num_splits', len(db_entries_splits_info))
        self.db.set_meta_value('genes_annotation_source', None)

        # creating empty default tables
        self.db.create_table(hmm_hits_info_table_name, hmm_hits_info_table_structure, hmm_hits_info_table_types)
        self.db.create_table(hmm_hits_splits_table_name, hmm_hits_splits_table_structure, hmm_hits_splits_table_types)
        self.db.create_table(hmm_hits_contigs_table_name, hmm_hits_contigs_table_structure, hmm_hits_contigs_table_types)
        self.db.create_table(genes_contigs_table_name, genes_contigs_table_structure, genes_contigs_table_types)
        self.db.create_table(genes_splits_summary_table_name, genes_splits_summary_table_structure, genes_splits_summary_table_types)
        self.db.create_table(genes_splits_table_name, genes_splits_table_structure, genes_splits_table_types)
        self.db.create_table(collections_info_table_name, collections_info_table_structure, collections_info_table_types)
        self.db.create_table(collections_contigs_table_name, collections_contigs_table_structure, collections_contigs_table_types)
        self.db.create_table(collections_splits_table_name, collections_splits_table_structure, collections_splits_table_types)

        self.disconnect()

        self.run.info('Annotation database', 'A new database, %s, has been created.' % (self.db_path), quiet = self.quiet)
        self.run.info('Number of contigs', num_contigs, quiet = self.quiet)
        self.run.info('Total number of nucleotides', total_length, quiet = self.quiet)
        self.run.info('Split length', split_length, quiet = self.quiet)


    def disconnect(self):
        self.db.disconnect()


class Table(object):
    """Superclass for rudimentary table needs and operations"""
    def __init__(self, db_path, run=run, progress=progress):
        if not os.path.exists(db_path):
            raise ConfigError, "Annotation database ('%s') does not exist. You must create one first." % db_path

        self.db_path = db_path
        self.next_available_id = {}

        self.run = run
        self.progress = progress

        annotation_db = AnnotationDatabase(self.db_path, quiet = False)
        self.split_length = annotation_db.db.get_meta_value('split_length')
        contig_lengths_table = annotation_db.db.get_table_as_dict(contig_lengths_table_name)
        self.splits = annotation_db.db.get_table_as_dict(splits_info_table_name)
        annotation_db.disconnect()

        self.contig_name_to_splits = {}
        for split_name in self.splits:
            parent = self.splits[split_name]['parent']
            if self.contig_name_to_splits.has_key(parent):
                self.contig_name_to_splits[parent].append(split_name)
            else:
                self.contig_name_to_splits[parent] = [split_name]

        self.contig_lengths = dict([(c, contig_lengths_table[c]['length']) for c in contig_lengths_table])


    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError, "If you need unique ids, you must call 'set_next_available_id' first"

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1


    def set_next_available_id(self, table):
        annotation_db = AnnotationDatabase(self.db_path)
        table_content = annotation_db.db.get_table_as_dict(table)
        if table_content:
            self.next_available_id[table] = max(table_content.keys()) + 1
        else:
            self.next_available_id[table] = 0

        annotation_db.disconnect()


    def export_contigs_in_db_into_FASTA_file(self):
        annotation_db = AnnotationDatabase(self.db_path, quiet = True)
        contig_sequences_table = annotation_db.db.get_table_as_dict(contig_sequences_table_name)
        annotation_db.disconnect()

        self.progress.new('Exporting contigs into a FASTA file')
        self.progress.update('...')
        contigs_fasta_path = os.path.join(filesnpaths.get_temp_directory_path(), 'contigs.fa')
        contigs_fasta = u.FastaOutput(contigs_fasta_path)
        for contig in contig_sequences_table:
            contigs_fasta.write_id(contig)
            contigs_fasta.write_seq(contig_sequences_table[contig]['sequence'], split=False)

        self.progress.end()
        self.run.info('FASTA for contigs', contigs_fasta_path)

        return contigs_fasta_path


class TablesForCollections(Table):
    """Populates the collections_* tables, where clusters of contigs and splits are kept"""
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, run, progress)

        # set these dudes so we have access to unique IDs:
        self.set_next_available_id(collections_contigs_table_name)
        self.set_next_available_id(collections_splits_table_name)


    def append(self, source, clusters_dict):
        annotation_db = AnnotationDatabase(self.db_path)

        # FIXME: this check is being done on multiple places, merge them:
        collections_info_table = annotation_db.db.get_table_as_dict(collections_info_table_name)
        if source in collections_info_table:
            self.run.info('WARNING', 'Clustering data for "%s" will be replaced with the incoming data' % source, header = True, display_only = True)
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (collections_info_table_name, source))
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (collections_contigs_table_name, source))
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (collections_splits_table_name, source))

        # push information about this search result into serach_info table.
        db_entries = tuple([source, ''])
        annotation_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % collections_info_table_name, db_entries)
        # then populate serach_data table for each contig.
        db_entries = [tuple([self.next_id(collections_contigs_table_name), source] + [v[h] for h in collections_contigs_table_structure[2:]]) for v in clusters_dict.values()]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % collections_contigs_table_name, db_entries)

        db_entries = self.process_splits(source, clusters_dict)
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % collections_splits_table_name, db_entries)

        annotation_db.disconnect()


    def process_splits(self, source, clusters_dict):
        db_entries_for_splits = []

        contig_to_cluster_id = dict([(d['contig'], d['cluster_id']) for d in clusters_dict.values()])

        for contig in contig_to_cluster_id:
            for split_name in self.contig_name_to_splits[contig]:
                db_entry = tuple([self.next_id(collections_splits_table_name), source, split_name, contig_to_cluster_id[contig]])
                db_entries_for_splits.append(db_entry)

        return db_entries_for_splits


class TablesForSearches(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        self.debug = False

        Table.__init__(self, self.db_path, run, progress)

        self.set_next_available_id(hmm_hits_contigs_table_name)
        self.set_next_available_id(hmm_hits_splits_table_name)


    def populate_search_tables(self, sources = {}):
        if not len(sources):
            import PaPi.data.hmm
            sources = PaPi.data.hmm.sources

        if not sources:
            return

        commander = HMMSearch()
        contigs_fasta = self.export_contigs_in_db_into_FASTA_file()
        proteins_in_contigs_fasta = commander.run_prodigal(contigs_fasta)
        if not self.debug:
            os.remove(contigs_fasta)

        for source in sources:
            kind_of_search = sources[source]['kind']
            all_genes_searched_against = sources[source]['genes']
            hmm_model = sources[source]['model']
            reference = sources[source]['ref']
            hmm_scan_hits_txt = commander.run_hmmscan(source,
                                                      all_genes_searched_against,
                                                      hmm_model,
                                                      reference)

            if not hmm_scan_hits_txt:
                search_results_dict = {}
            else:
                parser = parser_modules['search']['hmmscan'](proteins_in_contigs_fasta, hmm_scan_hits_txt)
                search_results_dict = parser.get_search_results()

            self.append(source, reference, kind_of_search, all_genes_searched_against, search_results_dict)

        if not self.debug:
            commander.clean_tmp_dirs()


    def append(self, source, reference, kind_of_search, all_genes, search_results_dict):
        annotation_db = AnnotationDatabase(self.db_path)

        hmm_hits_info_table = annotation_db.db.get_table_as_dict(hmm_hits_info_table_name)
        if source in hmm_hits_info_table:
            self.run.info('WARNING', 'Data for "%s" will be replaced with the incoming data' % source, header = True, display_only = True)
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (hmm_hits_info_table_name, source))
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (hmm_hits_contigs_table_name, source))
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (hmm_hits_splits_table_name, source))

        # push information about this search result into serach_info table.
        db_entries = [source, reference, kind_of_search, ', '.join(all_genes)]
        annotation_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % hmm_hits_info_table_name, db_entries)
        # then populate serach_data table for each contig.
        db_entries = [tuple([self.next_id(hmm_hits_contigs_table_name), source] + [v[h] for h in hmm_hits_contigs_table_structure[2:]]) for v in search_results_dict.values()]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % hmm_hits_contigs_table_name, db_entries)

        db_entries = self.process_splits(source, search_results_dict)
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % hmm_hits_splits_table_name, db_entries)

        annotation_db.disconnect()


    def process_splits(self, source, search_results_dict):
        hits_per_contig = {}
        for hit in search_results_dict.values():
            if hits_per_contig.has_key(hit['contig']):
                hits_per_contig[hit['contig']].append(hit)
            else:
                hits_per_contig[hit['contig']] = [hit]

        db_entries_for_splits = []

        for contig in self.contig_lengths:
            if not hits_per_contig.has_key(contig):
                # no hits for this contig. pity!
                continue

            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits[split_name]['start']
                stop = self.splits[split_name]['end']

                # FIXME: this really needs some explanation.
                for hit in hits_per_contig[contig]:
                    if hit['stop'] > start and hit['start'] < stop:
                        gene_length = hit['stop'] - hit['start']
                        # if only a part of the gene is in the split:
                        start_in_split = (start if hit['start'] < start else hit['start']) - start
                        stop_in_split = (stop if hit['stop'] > stop else hit['stop']) - start
                        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length
                        
                        gene_unique_identifier = hashlib.sha224('_'.join([contig, hit['gene_name'], str(hit['start']), str(hit['stop'])])).hexdigest()
                        db_entry = tuple([self.next_id(hmm_hits_splits_table_name), source, gene_unique_identifier, hit['gene_name'], split_name, percentage_in_split, hit['e_value']])
                        db_entries_for_splits.append(db_entry)

        return db_entries_for_splits


class TablesForGenes(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, run, progress)

        # this class keeps track of genes that occur in splits, and responsible
        # for generating the necessary table in the annotation database
        self.genes_in_splits = GenesInSplits()


    def create(self, genes_dict, parser):
        self.genes_dict = genes_dict

        self.sanity_check()

        # oepn connection
        annotation_db = AnnotationDatabase(self.db_path)

        # test whether there are already genes tables populated
        genes_annotation_source = annotation_db.db.get_meta_value('genes_annotation_source')
        if genes_annotation_source:
            self.run.info('WARNING', 'Previous genes annotation data from "%s" will be replaced with the incoming data' % parser, header = True, display_only = True)
            annotation_db.db._exec('''DELETE FROM %s''' % (genes_contigs_table_name))
            annotation_db.db._exec('''DELETE FROM %s''' % (genes_splits_table_name))
            annotation_db.db._exec('''DELETE FROM %s''' % (genes_splits_summary_table_name))

        # set the parser
        annotation_db.db.remove_meta_key_value_pair('genes_annotation_source')
        annotation_db.db.set_meta_value('genes_annotation_source', parser)
        # push raw entries
        db_entries = [tuple([prot] + [self.genes_dict[prot][h] for h in genes_contigs_table_structure[1:]]) for prot in self.genes_dict]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''' % genes_contigs_table_name, db_entries)
        # disconnect like a pro.
        annotation_db.disconnect()


        # compute and push split taxonomy information.
        self.init_genes_splits_summary_table()


    def sanity_check(self):
        # check whether input matrix dict 
        keys_found = ['prot'] + self.genes_dict.values()[0].keys()
        missing_keys = [key for key in genes_contigs_table_structure if key not in keys_found]
        if len(missing_keys):
            raise ConfigError, "Your input lacks one or more header fields to generate a PaPi annotation db. Here is\
                                what you are missing: %s. The complete list (and order) of headers in your TAB\
                                delimited matrix file (or dictionary) must follow this: %s." % (', '.join(missing_keys),
                                                                                                ', '.join(genes_contigs_table_structure))


        contig_names_in_matrix = set([v['contig'] for v in self.genes_dict.values()])
        contig_names_in_db  = set(self.contig_lengths.keys())

        for contig in contig_names_in_matrix:
            if contig not in contig_names_in_db:
                raise ConfigError, "We have a problem... Every contig name found in the input file you provide\
                                    must be found in the annotation database. But it seems it is not the case. I did not check\
                                    all, but there there is at least one contig name ('%s') that appears in your\
                                    matrices, but missing in the database. You may need to format the contig\
                                    names in your FASTA file and regenerate the annotation database to match contig\
                                    names appear in your matrices. Keep in mind that contig names must also match the\
                                    ones in your BAM files later on. Even when you use one software for assembly and\
                                    mapping, disagreements between contig names may arise. We know that it is the case\
                                    with CLC for instance. OK. Going back to the issue. Here is one contig name from\
                                    the annotation database (which was originally in your contigs FASTA): '%s', and\
                                    here is one from your input files you just provided: '%s'. You should make them\
                                    identical (and make sure whatever solution you come up with will not make them\
                                    incompatible with names in your BAM files later on. Sorry about this mess, but\
                                    there is nothing much PaPi can do about this issue." %\
                                                    (contig, contig_names_in_db.pop(), contig_names_in_matrix.pop())


    def init_genes_splits_summary_table(self):
        # build a dictionary for fast access to all proteins identified within a contig
        prots_in_contig = {}
        for prot in self.genes_dict:
            contig = self.genes_dict[prot]['contig']
            if prots_in_contig.has_key(contig):
                prots_in_contig[contig].add(prot)
            else:
                prots_in_contig[contig] = set([prot])

        contigs_without_annotation = list(set(self.contig_lengths.keys()) - set(prots_in_contig.keys()))
        run.info('Percent of contigs annotated', '%.1f%%' % (len(prots_in_contig) * 100.0 / len(self.contig_lengths)))

        for contig in contigs_without_annotation:
            prots_in_contig[contig] = set([])

        splits_dict = {}
        split_to_prot = {}
        for contig in self.contig_lengths:
            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits[split_name]['start']
                stop = self.splits[split_name]['end']

                taxa = []
                functions = []
                gene_start_stops = []
                # here we go through all genes in the contig and identify the all the ones that happen to be in
                # this particular split to generate summarized info for each split. BUT one important that is done
                # in the following loop is self.genes_in_splits.add call, which populates GenesInSplits class.
                for prot in prots_in_contig[contig]:
                    if self.genes_dict[prot]['stop'] > start and self.genes_dict[prot]['start'] < stop:
                        taxa.append(self.genes_dict[prot]['t_species'])
                        functions.append(self.genes_dict[prot]['function'])
                        gene_start_stops.append((self.genes_dict[prot]['start'], self.genes_dict[prot]['stop']), )
                        self.genes_in_splits.add(split_name, start, stop, prot, self.genes_dict[prot]['start'], self.genes_dict[prot]['stop'])


                taxonomy_strings = [t for t in taxa if t]
                function_strings = [f for f in functions if f]

                # here we identify genes that are associated with a split even if one base of the gene spills into 
                # the defined start or stop of a split, which means, split N, will include genes A, B and C in this
                # scenario:
                #
                # contig: (...)------[ gene A ]--------[     gene B    ]----[gene C]---------[    gene D    ]-----(...)
                #         (...)----------x---------------------------------------x--------------------------------(...)
                #                        ^ (split N start)                       ^ (split N stop)
                #                        |                                       |
                #                        |<-              split N              ->|
                #
                # however, when looking at the coding versus non-coding nucleotide ratios in a split, we have to make
                # sure that only the relevant portion of gene A and gene C is counted:
                total_coding_nts = 0
                for gene_start, gene_stop in gene_start_stops:
                    total_coding_nts += (gene_stop if gene_stop < stop else stop) - (gene_start if gene_start > start else start)

                splits_dict[split_name] = {'taxonomy': None,
                                           'num_genes': len(taxa),
                                           'avg_gene_length': numpy.mean([(l[1] - l[0]) for l in gene_start_stops]) if len(gene_start_stops) else 0.0,
                                           'ratio_coding': total_coding_nts * 1.0 / (stop - start),
                                           'ratio_hypothetical': (len(functions) - len(function_strings)) * 1.0 / len(functions) if len(functions) else 0.0,
                                           'ratio_with_tax': len(taxonomy_strings) * 1.0 / len(taxa) if len(taxa) else 0.0,
                                           'tax_accuracy': 0.0}
                distinct_taxa = set(taxonomy_strings)

                if not len(distinct_taxa):
                    continue

                if len(distinct_taxa) == 1:
                    splits_dict[split_name]['taxonomy'] = distinct_taxa.pop()
                    splits_dict[split_name]['tax_accuracy'] = 1.0
                else:
                    d = Counter()
                    for taxon in taxonomy_strings:
                        d[taxon] += 1
                    consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
                    splits_dict[split_name]['taxonomy'] = consensus
                    splits_dict[split_name]['tax_accuracy'] = occurrence * 1.0 / len(taxonomy_strings)

        # open connection
        annotation_db = AnnotationDatabase(self.db_path)
        # push raw entries for splits table
        db_entries = [tuple([split] + [splits_dict[split][h] for h in genes_splits_summary_table_structure[1:]]) for split in splits_dict]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % genes_splits_summary_table_name, db_entries)
        # push entries for genes in splits table
        db_entries = [tuple([entry_id] + [self.genes_in_splits.splits_to_prots[entry_id][h] for h in genes_splits_table_structure[1:]]) for entry_id in self.genes_in_splits.splits_to_prots]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % genes_splits_table_name, db_entries)
        # disconnect
        annotation_db.disconnect()


    def get_consensus_taxonomy_for_split(self, contig, t_level = 't_species', start = 0, stop = sys.maxint):
        """Returns (c, n, t, o) where,
            c: consensus taxonomy (the most common taxonomic call for each gene found in the contig),
            n: total number of genes found in the contig,
            t: total number of genes with known taxonomy,
            o: number of taxonomic calls that matches the consensus among t
        """

        response = self.db.cursor.execute("""SELECT %s FROM %s WHERE contig='%s' and stop > %d and start < %d""" % (t_level, genes_contigs_table_name, contig, start, stop))
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


class GenesInSplits:
    def __init__(self):
        self.entry_id = 0
        self.splits_to_prots = {}

    def add(self, split_name, split_start, split_end, prot_id, prot_start, prot_end):

        gene_length = prot_end - prot_start

        if gene_length <= 0:
            raise ConfigError, "annotation.py/GeneInSplits: OK. There is something wrong. We have this gene, '%s',\
                                which starts at position %d and ends at position %d. Well, it doesn't look right,\
                                does it?" % (prot_id, prot_start, prot_end)

        # if only a part of the gene is in the split:
        start_in_split = (split_start if prot_start < split_start else prot_start) - split_start
        stop_in_split = (split_end if prot_end > split_end else prot_end) - split_start
        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length

        self.splits_to_prots[self.entry_id] = {'split': split_name,
                                               'prot': prot_id,
                                               'start_in_split': start_in_split,
                                               'stop_in_split': stop_in_split,
                                               'percentage_in_split': percentage_in_split}
        self.entry_id += 1

