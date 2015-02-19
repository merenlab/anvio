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

annotation_table_name      = 'annotation'
annotation_table_structure = ['prot', 'contig', 'start', 'stop'   , 'direction', 'figfam', 'function', "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]
annotation_table_mapping   = [ str  ,   str   ,  int   ,   int    ,     str    ,   str   ,    str    ,    str    ,   str    ,    str   ,    str    ,    str   ,     str    ]
annotation_table_types     = ['text',  'text' ,'numeric','numeric',   'text'   ,  'text' ,   'text'  ,   'text'  ,  'text'  ,  'text'  ,  'text'   ,  'text'  ,   'text'   ]

splits_table_name      = 'splits'
splits_table_structure = ['split', 'taxonomy', 'num_genes', 'avg_gene_length', 'ratio_coding', 'ratio_hypothetical', 'ratio_with_tax', 'tax_accuracy']
splits_table_mapping   = [  str  ,     str   ,    int     ,       float      ,    'float'    ,        float        ,       float     ,     float     ]
splits_table_types     = [ 'text',   'text'  ,  'numeric' ,     'numeric'    ,   'numeric'   ,      'numeric'      ,     'numeric'   ,   'numeric'   ]

splits_to_prots_table_name      = 'genes_in_splits'
splits_to_prots_table_structure = ['entry_id', 'split', 'prot', 'start_in_split', 'stop_in_split', 'percentage_in_split']
splits_to_prots_table_mapping   = [    int   ,   str  ,  str  ,       int       ,       int      ,         float        ]
splits_to_prots_table_types     = [ 'numeric',  'text', 'text',    'numeric'    ,    'numeric'   ,       'numeric'      ]

search_info_table_name      = 'search_info'
search_info_table_structure = ['source', 'ref' , 'search_type', 'genes']
search_info_table_mapping   = [  str   ,  str  ,      str     ,   str  ]
search_info_table_types     = [ 'text' , 'text',    'text'    , 'text' ]

search_contigs_table_name      = 'search_contigs'
search_contigs_table_structure = ['entry_id', 'source', 'contig', 'start' , 'stop'  , 'gene_name', 'gene_id', 'e_value']
search_contigs_table_mapping   = [    int   ,   str   ,    str  ,   int   ,   int   ,     str    ,    str   ,   float  ]   
search_contigs_table_types     = [ 'numeric',  'text' ,  'text' ,'numeric','numeric',   'text'   ,  'text'  , 'numeric']

search_splits_table_name      = 'search_splits'
search_splits_table_structure = ['entry_id', 'source', 'gene_unique_identifier', 'gene_name', 'split', 'percentage_in_split', 'e_value']
search_splits_table_mapping   = [    int   ,   str   ,            str          ,     str    ,   str  ,         float        ,   float  ]   
search_splits_table_types     = [ 'numeric',  'text' ,          'text'         ,   'text'   ,  'text',       'numeric'      , 'numeric']

collections_info_table_name      = 'collections_info'
collections_info_table_structure = ['source',  'ref']
collections_info_table_types     = [ 'text' , 'text']

collections_contigs_table_name      = 'collections_contigs'
collections_contigs_table_structure = ['entry_id', 'source', 'contig', 'cluster_id']
collections_contigs_table_types     = [ 'numeric',  'text' ,  'text' ,    'text'   ]

collections_splits_table_name      = 'collections_splits'
collections_splits_table_structure = ['entry_id', 'source', 'split', 'cluster_id']
collections_splits_table_types     = [ 'numeric',  'text' , 'text' ,    'text'   ]


__version__ = "0.4.0"


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
import PaPi.dictio as dictio
import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths

from PaPi.utils import ConfigError
from PaPi.contig import Split
from PaPi.commandline import HMMSearch
from PaPi.parsers import parser_modules

run = terminal.Run()
progress = terminal.Progress()


def magic_touch(annotation_db_path, split_length = None, run=run, progress=progress):
    """you call this function, and it creates a new database for you if the path you sent does not exists, or gets
       the split size if the annotation_db_path does exist. It is like 'touch' command in UNIX environments, hence
       my silly naming."""

    # FIXME: this is a risky convenience (maybe I should be more explicit): 
    if os.path.isdir(annotation_db_path):
        annotation_db_path = os.path.join(annotation_db_path, 'ANNOTATION.db')

    if os.path.exists(annotation_db_path):
        if split_length:
            raise ConfigError, "You can not declare a split length (-L/--split-length) when you are\
                                working with an already existing annotation database."

        annotation_db = AnnotationDB(annotation_db_path)
        split_length = annotation_db.db.get_meta_value('split_length')
        run.info('Database touch', 'An existing database, %s, has been initiated.' % annotation_db_path)
        run.info('Split length', split_length)
        annotation_db.disconnect()
    else:
        split_length = split_length
        if not split_length:
            raise ConfigError, "A new annotation database is about to be created, but the split length information\
                                is missing from the request. If you just run a PaPi program and got this error,\
                                you should try to take a look at the --help menu for relevant parameter."

        annotation_db = AnnotationDB(annotation_db_path, split_length, create_new = True)
        run.info('Database touch', 'A new database, %s, has been created.' % annotation_db_path)
        run.info('Split length', split_length)
        annotation_db.disconnect()

    return (annotation_db_path, split_length)


class PopulateAnnotationDB:
    def __init__(self, args = None):
        self.contigs_fasta = None
        self.parser = None
        self.split_length = 20000
        self.input_files = []
        self.skip_search_tables = False
        self.annotation_db_path = None
        self.debug = False

        if args:
            self.contigs_fasta = args.contigs_fasta
            self.parser = args.parser
            self.split_length = args.split_length
            self.input_files = args.input_files
            self.skip_search_tables = args.skip_search_tables
            self.annotation_db_path = args.output
            self.debug = args.debug

        self.sanity_checked = False
        self.contig_lengths = {}


    def sanity_check(self):
        if not self.contigs_fasta:
            raise ConfigError, "This is not going to work without a FASTA file of contigs. Please see the help menu :/"

        if not os.path.exists(self.contigs_fasta):
            raise ConfigError, "PaPi can't find a FASTA file here: '%s' :/" % self.contigs_fasta

        if not self.annotation_db_path:
            raise ConfigError, "Annotation database path cannot be None :/"

        # well maybe .. since we are here ...
        self.contig_lengths = utils.get_read_lengths_from_fasta(self.contigs_fasta)

        self.sanity_checked = True


    def _run(self):
        self.sanity_check()

        # get a new, empty annotation database is created with all tables and stuff:
        annotation_db = AnnotationDB(self.annotation_db_path, self.split_length, create_new = True)
        annotation_db.db.disconnect()

        self.populate_annotation_tables()

        self.populate_search_tables()


    def populate_annotation_tables(self):
        if not self.sanity_check:
            raise ConfigError, "You must first call sanity_check()"

        if type(self.parser) == type(None):
            raise ConfigError, "You must specify a parser. Please see --help menu or the documentation for a list."
        if self.parser not in parser_modules['annotation']:
            raise ConfigError, "I don't know what to do with '%s'. Please enter a valid parser. Here is a list of\
                                parsers available for annotation data: %s" % (self.parser, ', '.join(parser_modules['annotation']))

        parser = parser_modules['annotation'][self.parser](self.input_files, annotation_table_structure)
        annotations_dict = parser.get_annotations_dict()
        annotation_tables = AnnotationTables(self.annotation_db_path, self.contig_lengths)
        annotation_tables.create(annotations_dict, self.parser)


    def populate_search_tables(self, sources = {}):
        if self.skip_search_tables:
            return

        if not self.sanity_check:
            raise ConfigError, "You must first call sanity_check()"

        if not len(sources):
            import PaPi.data.hmm
            sources = PaPi.data.hmm.sources

        if not sources:
            return

        commander = HMMSearch()
        search_tables = SearchTables(self.annotation_db_path, self.contig_lengths)

        proteins_in_contigs_fasta = commander.run_prodigal(self.contigs_fasta)

        search_results_dict = {}
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
                continue

            parser = parser_modules['search']['hmmscan'](proteins_in_contigs_fasta, hmm_scan_hits_txt)
            search_results_dict = parser.get_search_results()

            search_tables.append(source, reference, kind_of_search, all_genes_searched_against, search_results_dict)

        if not self.debug:
            commander.clean_tmp_dirs()


class AnnotationDB:
    """This class will connect to an annotation database at db_path. It will create one if requested."""
    def __init__(self, db_path, split_length = None, create_new = False, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress
        self.split_length = None

        self.db = None

        if create_new:
            self.new_annotation_database(split_length)
        else:
            self.init_database()
            self.split_length = self.db.get_meta_value('split_length')


    def init_database(self):
        self.db = db.DB(self.db_path, __version__)


    def new_annotation_database(self, split_length):
        if not split_length:
            raise ConfigError, "Creating a new annotation database requires split length information to be\
                                provided. But the AnnotationDB class was called to create one, without the\
                                split length information. Not cool."

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

        self.db.create_table(search_info_table_name, search_info_table_structure, search_info_table_types)
        self.db.create_table(search_splits_table_name, search_splits_table_structure, search_splits_table_types)
        self.db.create_table(search_contigs_table_name, search_contigs_table_structure, search_contigs_table_types)
        self.db.create_table(annotation_table_name, annotation_table_structure, annotation_table_types)
        self.db.create_table(splits_table_name, splits_table_structure, splits_table_types)
        self.db.create_table(splits_to_prots_table_name, splits_to_prots_table_structure, splits_to_prots_table_types)
        self.db.create_table(collections_info_table_name, collections_info_table_structure, collections_info_table_types)
        self.db.create_table(collections_contigs_table_name, collections_contigs_table_structure, collections_contigs_table_types)
        self.db.create_table(collections_splits_table_name, collections_splits_table_structure, collections_splits_table_types)


    def disconnect(self):
        self.db.disconnect()


class Table(object):
    def __init__(self, db_path):
        self.db_path = db_path
        self.next_available_id = {}

        annotation_db = AnnotationDB(self.db_path)
        self.split_length = annotation_db.split_length
        annotation_db.disconnect()


    def next_id(self, table):
        if table not in self.next_available_id:
            raise ConfigError, "If you need unique ids, you must call 'set_next_available_id' first"

        self.next_available_id[table] += 1
        return self.next_available_id[table] - 1


    def set_next_available_id(self, table):
        annotation_db = AnnotationDB(self.db_path)
        table_content = annotation_db.db.get_table_as_dict(table)
        if table_content:
            self.next_available_id[table] = max(table_content.keys()) + 1
        else:
            self.next_available_id[table] = 0

        annotation_db.disconnect()



class CollectionsTables(Table):
    """Populates the collections_* tables, where clusters of contigs and splits are kept"""
    def __init__(self, db_path, contig_lengths):
        self.db_path = db_path
        self.contig_lengths = contig_lengths

        Table.__init__(self, self.db_path)

        # set these dudes so we have access to unique IDs:
        self.set_next_available_id(collections_contigs_table_name)
        self.set_next_available_id(collections_splits_table_name)


    def append(self, source, clusters_dict):
        annotation_db = AnnotationDB(self.db_path)

        # FIXME: this check is being done on multiple places, merge them:
        collections_info_table = annotation_db.db.get_table_as_dict(collections_info_table_name)
        if source in collections_info_table:
            run.info('WARNING', 'Clustering data for "%s" will be replaced with the incoming data' % source, header = True, display_only = True)
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
            chunks = utils.get_chunks(self.contig_lengths[contig], self.split_length)
            for i in range(0, len(chunks)):
                split = Split(contig, i).name

                db_entry = tuple([self.next_id(collections_splits_table_name), source, split, contig_to_cluster_id[contig]])
                db_entries_for_splits.append(db_entry)

        return db_entries_for_splits


class SearchTables(Table):
    def __init__(self, db_path, contig_lengths):
        self.db_path = db_path
        self.contig_lengths = contig_lengths

        Table.__init__(self, self.db_path)

        self.set_next_available_id(search_contigs_table_name)
        self.set_next_available_id(search_splits_table_name)


    def append(self, source, reference, kind_of_search, all_genes, search_results_dict):
        annotation_db = AnnotationDB(self.db_path)

        search_info_table = annotation_db.db.get_table_as_dict(search_info_table_name)
        if source in search_info_table:
            run.info('WARNING', 'Data for "%s" will be replaced with the incoming data' % source, header = True, display_only = True)
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (search_info_table_name, source))
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (search_contigs_table_name, source))
            annotation_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (search_splits_table_name, source))

        # push information about this search result into serach_info table.
        db_entries = [source, reference, kind_of_search, ', '.join(all_genes)]
        annotation_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % search_info_table_name, db_entries)
        # then populate serach_data table for each contig.
        db_entries = [tuple([self.next_id(search_contigs_table_name), source] + [v[h] for h in search_contigs_table_structure[2:]]) for v in search_results_dict.values()]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % search_contigs_table_name, db_entries)

        db_entries = self.process_splits(source, search_results_dict)
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % search_splits_table_name, db_entries)

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

            chunks = utils.get_chunks(self.contig_lengths[contig], self.split_length)
            for i in range(0, len(chunks)):
                split = Split(contig, i).name
                start = chunks[i][0]
                stop = chunks[i][1]

                # FIXME: this really needs some explanation.
                for hit in hits_per_contig[contig]:
                    if hit['stop'] > start and hit['start'] < stop:
                        gene_length = hit['stop'] - hit['start']
                        # if only a part of the gene is in the split:
                        start_in_split = (start if hit['start'] < start else hit['start']) - start
                        stop_in_split = (stop if hit['stop'] > stop else hit['stop']) - start
                        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length
                        
                        gene_unique_identifier = hashlib.sha224('_'.join([contig, hit['gene_name'], str(hit['start']), str(hit['stop'])])).hexdigest()
                        db_entry = tuple([self.next_id(search_splits_table_name), source, gene_unique_identifier, hit['gene_name'], split, percentage_in_split, hit['e_value']])
                        db_entries_for_splits.append(db_entry)

        return db_entries_for_splits



class AnnotationTables:
    def __init__(self, annotation_db_path, contig_lengths):
        self.annotation_db_path = annotation_db_path
        self.contig_lengths = contig_lengths

        # read split length
        annotation_db = AnnotationDB(self.annotation_db_path)
        self.split_length = annotation_db.db.get_meta_value('split_length')
        annotation_db.disconnect()

        # this class keeps track of genes that occur in splits, and responsible
        # for generating the necessary table in the annotation database
        self.genes_in_splits = GenesInSplits()


    def create(self, data_source, parser):
        if type(data_source) == type(dict()):
            self.matrix_dict = data_source
            if len(self.matrix_dict):
                self.check_keys(['prot'] + self.matrix_dict.values()[0].keys())
        if type(data_source) == type(str()):
            self.matrix_dict = utils.get_TAB_delimited_file_as_dictionary(data_source,
                                                                          column_names = annotation_table_structure,
                                                                          column_maping = annotation_table_mapping)

        self.sanity_check()

        # oepn connection
        annotation_db = AnnotationDB(self.annotation_db_path)
        # set the parser
        annotation_db.db.set_meta_value('annotation_source', parser)
        # push raw entries
        db_entries = [tuple([prot] + [self.matrix_dict[prot][h] for h in annotation_table_structure[1:]]) for prot in self.matrix_dict]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''' % annotation_table_name, db_entries)
        # disconnect like a pro.
        annotation_db.disconnect()


        # compute and push split taxonomy information.
        self.init_splits_table()


    def sanity_check(self):
        contig_names_in_matrix = set([v['contig'] for v in self.matrix_dict.values()])
        contig_names_in_fasta  = set(self.contig_lengths.keys())

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

        contigs_without_annotation = list(set(self.contig_lengths.keys()) - set(prots_in_contig.keys()))
        run.info('Num contigs in FASTA', len(self.contig_lengths))
        run.info('Num contigs w annotation', len(prots_in_contig))
        run.info('Num contigs w/o annotation', len(contigs_without_annotation))

        for contig in contigs_without_annotation:
            prots_in_contig[contig] = set([])

        splits_dict = {}
        split_to_prot = {}
        for contig in self.contig_lengths:
            chunks = utils.get_chunks(self.contig_lengths[contig], self.split_length)
            for i in range(0, len(chunks)):
                split = Split(contig, i).name
                start = chunks[i][0]
                stop = chunks[i][1]

                taxa = []
                functions = []
                gene_start_stops = []
                # here we go through all genes in the contig and identify the all the ones that happen to be in
                # this particular split to generate summarized info for each split. BUT one important that is done
                # in the following loop is self.genes_in_splits.add call, which populates GenesInSplits class.
                for prot in prots_in_contig[contig]:
                    if self.matrix_dict[prot]['stop'] > start and self.matrix_dict[prot]['start'] < stop:
                        taxa.append(self.matrix_dict[prot]['t_species'])
                        functions.append(self.matrix_dict[prot]['function'])
                        gene_start_stops.append((self.matrix_dict[prot]['start'], self.matrix_dict[prot]['stop']), )
                        self.genes_in_splits.add(split, start, stop, prot, self.matrix_dict[prot]['start'], self.matrix_dict[prot]['stop'])


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

                splits_dict[split] = {'taxonomy': None,
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
                    splits_dict[split]['taxonomy'] = distinct_taxa.pop()
                    splits_dict[split]['tax_accuracy'] = 1.0
                else:
                    d = Counter()
                    for taxon in taxonomy_strings:
                        d[taxon] += 1
                    consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
                    splits_dict[split]['taxonomy'] = consensus
                    splits_dict[split]['tax_accuracy'] = occurrence * 1.0 / len(taxonomy_strings)

        # open connection
        annotation_db = AnnotationDB(self.annotation_db_path)
        # push raw entries for splits table
        db_entries = [tuple([split] + [splits_dict[split][h] for h in splits_table_structure[1:]]) for split in splits_dict]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % splits_table_name, db_entries)
        # push entries for genes in splits table
        db_entries = [tuple([entry_id] + [self.genes_in_splits.splits_to_prots[entry_id][h] for h in splits_to_prots_table_structure[1:]]) for entry_id in self.genes_in_splits.splits_to_prots]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % splits_to_prots_table_name, db_entries)
        # disconnect
        annotation_db.disconnect()


    def get_consensus_taxonomy_for_split(self, contig, t_level = 't_species', start = 0, stop = sys.maxint):
        """Returns (c, n, t, o) where,
            c: consensus taxonomy (the most common taxonomic call for each gene found in the contig),
            n: total number of genes found in the contig,
            t: total number of genes with known taxonomy,
            o: number of taxonomic calls that matches the consensus among t
        """

        response = self.db.cursor.execute("""SELECT %s FROM %s WHERE contig='%s' and stop > %d and start < %d""" % (t_level, annotation_table_name, contig, start, stop))
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

