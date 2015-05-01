# -*- coding: utf-8
"""
    Classes to create, access, and/or populate annotation and profile databases.
"""

import os
import sys
import time
import copy
import numpy
import random
import hashlib
import operator
from collections import Counter

import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.kmers as kmers
import anvio.contig as contig
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections

from anvio.errors import ConfigError
from anvio.commandline import HMMSearch
from anvio.parsers import parser_modules
from anvio.tableops import Table


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class AnnotationSuperclass(object):
    def __init__(self, args, r = run, p = progress):
        self.run = r
        self.progress = p

        self.a_meta = {}
        self.genes_in_contigs_dict = {}
        self.genes_in_splits = {}
        self.genes_in_splits_summary_dict = {}
        self.genes_in_splits_summary_headers = []
        self.split_to_genes_in_splits_ids = {} # for fast access to all self.genes_in_splits entries for a given split
        self.contigs_basic_info = {}
        self.split_sequences = {}
        self.hmm_sources_info = {}
        self.singlecopy_gene_hmm_sources = set([])
        self.non_singlecopy_gene_hmm_sources = set([])

        self.hmm_searches_dict = {}   # <--- upon initiation, this dict only keeps hmm hits for non-singlecopy
        self.hmm_searches_header = [] #      gene searches... single-copy gene info is accessed through completeness.py

        try:
            self.annotation_db_path = args.annotation_db
        except:
            self.run.warning('AnnotationSuperclass class called with args without annotation_db member')
            return

        if not self.annotation_db_path:
            return

        filesnpaths.is_file_exists(self.annotation_db_path)

        self.progress.new('Loading the annotation DB')
        annotation_db = AnnotationDatabase(self.annotation_db_path)

        self.progress.update('Setting annotation metadata dict')
        self.a_meta = annotation_db.meta

        self.a_meta['creation_date'] = utils.get_time_to_date(self.a_meta['creation_date']) if self.a_meta.has_key('creation_date') else 'unknown'

        self.progress.update('Reading contigs basic info')
        self.contigs_basic_info = annotation_db.db.get_table_as_dict(t.contigs_info_table_name, string_the_key = True)

        self.progress.update('Reading splits basic info')
        self.splits_basic_info = annotation_db.db.get_table_as_dict(t.splits_info_table_name)

        self.progress.update('Reading genes in contigs table')
        self.genes_in_contigs_dict = annotation_db.db.get_table_as_dict(t.genes_contigs_table_name)

        self.progress.update('Reading genes in splits table')
        self.genes_in_splits = annotation_db.db.get_table_as_dict(t.genes_splits_table_name)

        self.progress.update('Reading genes in splits summary table')
        self.genes_in_splits_summary_dict = annotation_db.db.get_table_as_dict(t.genes_splits_summary_table_name)
        self.genes_in_splits_summary_headers = annotation_db.db.get_table_structure(t.genes_splits_summary_table_name)

        self.progress.update('Identifying HMM searches for single-copy genes and others')
        self.hmm_sources_info = annotation_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
        self.singlecopy_gene_hmm_sources = set([s for s in self.hmm_sources_info.keys() if self.hmm_sources_info[s]['search_type'] == 'singlecopy'])
        self.non_singlecopy_gene_hmm_sources = set([s for s in self.hmm_sources_info.keys() if self.hmm_sources_info[s]['search_type'] != 'singlecopy'])

        self.progress.update('Generating split to genes in splits mapping dict')
        for entry_id in self.genes_in_splits:
            split_name = self.genes_in_splits[entry_id]['split']
            if split_name in self.split_to_genes_in_splits_ids:
                self.split_to_genes_in_splits_ids[split_name].add(entry_id)
            else:
                self.split_to_genes_in_splits_ids[split_name] = set([entry_id])

        self.progress.end()

        annotation_db.disconnect()
        run.info('Annotation DB', 'Initialized: %s (v. %s)' % (self.annotation_db_path, annotation_db.db.version))


    def init_split_sequences(self, min_contig_length = 0):
        self.progress.new('Loading split sequences')

        annotation_db = AnnotationDatabase(self.annotation_db_path)

        self.progress.update('Identifying contigs shorter than M')
        contigs_shorter_than_M = set([c for c in self.contigs_basic_info if self.contigs_basic_info[c]['length'] < min_contig_length])

        self.progress.update('Reading contig sequences')
        contigs_sequences = annotation_db.db.get_table_as_dict(t.contig_sequences_table_name, string_the_key = True)

        self.progress.update('Filtering out shorter contigs')
        for contig_name in contigs_shorter_than_M:
            contigs_sequences.pop(contig_name)

        self.progress.update('Discarding split names coming from short contigs')
        split_names_to_discard = set([])
        for split_name in self.splits_basic_info:
            if self.splits_basic_info[split_name]['parent'] in contigs_shorter_than_M:
                split_names_to_discard.add(split_name)

        for split_name in split_names_to_discard:
            self.splits_basic_info.pop(split_name)

        # FIXME: THIS VARIABLE SHOULD BE REMOVED ONCE WE TAG version 1.0
        db_has_num_splits = 'num_splits' in self.contigs_basic_info.keys()

        self.progress.update('Generating split sequences dict')
        for split_name in self.splits_basic_info:
            split = self.splits_basic_info[split_name]

            if split['parent'] in contigs_shorter_than_M:
                contigs_shorter_than_M.remove(split['parent'])
                continue

            if db_has_num_splits and self.contigs_basic_info[split['parent']]['num_splits'] == 1:
                self.split_sequences[split_name] = contigs_sequences[split['parent']]['sequence']
            else:
                self.split_sequences[split_name] = contigs_sequences[split['parent']]['sequence'][split['start']:split['end']]

        self.progress.end()

        annotation_db.disconnect()


    def init_non_singlecopy_gene_hmm_sources(self, split_names_of_interest = None, return_each_gene_as_a_layer = False):
        if not self.annotation_db_path:
            return

        self.progress.new('Loading split sequences')
        self.progress.update('...')

        annotation_db = AnnotationDatabase(self.annotation_db_path)

        if len(self.non_singlecopy_gene_hmm_sources):
            non_singlecopy_gene_hmm_results_dict = utils.get_filtered_dict(annotation_db.db.get_table_as_dict(t.hmm_hits_splits_table_name), 'source', self.non_singlecopy_gene_hmm_sources)
            non_singlecopy_gene_hmm_info_dict = annotation_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
            for source in self.singlecopy_gene_hmm_sources:
                non_singlecopy_gene_hmm_info_dict.pop(source)
        else:
            return 

        if split_names_of_interest:
            non_singlecopy_gene_hmm_results_dict = utils.get_filtered_dict(non_singlecopy_gene_hmm_results_dict, 'split', set(split_names_of_interest))

        sources_tmpl = {}

        # the following conditional is pretty critical. here is more info about the difference:
        # https://github.com/meren/anvio/issues/123
        if return_each_gene_as_a_layer:
            for source in self.non_singlecopy_gene_hmm_sources:
                search_type = self.hmm_sources_info[source]['search_type']
                for gene_name in [g.strip() for g in non_singlecopy_gene_hmm_info_dict[source]['genes'].split(',')]:
                    search_term = 'hmmx_%s_%s' % (search_type, gene_name)
                    sources_tmpl[search_term] = 0
                    self.hmm_searches_header.append(search_term)

            # fill all splits with 0s, so this is treated as a numeric column:
            for split_name in split_names_of_interest if split_names_of_interest else self.splits_basic_info:
                self.hmm_searches_dict[split_name] = copy.deepcopy(sources_tmpl)

            for e in non_singlecopy_gene_hmm_results_dict.values():
                search_term = 'hmmx_%s_%s' % (self.hmm_sources_info[e['source']]['search_type'], e['gene_name'])
                self.hmm_searches_dict[e['split']][search_term] = 1
        else:
            for source in self.non_singlecopy_gene_hmm_sources:
                search_type = self.hmm_sources_info[source]['search_type']
                sources_tmpl[search_type] = {}
                self.hmm_searches_header.append(search_type)

            for e in non_singlecopy_gene_hmm_results_dict.values():
                if not e['split'] in self.hmm_searches_dict:
                    self.hmm_searches_dict[e['split']] = copy.deepcopy(sources_tmpl)

                search_type = self.hmm_sources_info[e['source']]['search_type']
                self.hmm_searches_dict[e['split']][search_type] = e['gene_name']

        self.progress.end()


class ProfileSuperclass(object):
    def __init__(self, args, r = run, p = progress):
        self.args = args
        self.run = r
        self.progress = p

        self.clusterings = {}
        self.views = {}

        try:
            self.profile_db_path = args.profile_db
        except:
            self.run.warning('ProfileSuperclass class called with args without profile_db member')
            return

        if not self.profile_db_path:
            return

        filesnpaths.is_file_exists(self.profile_db_path)

        self.progress.new('Loading the annotation DB')
        profile_db = ProfileDatabase(self.profile_db_path)

        self.progress.update('Setting profile metadata dict')
        self.p_meta = profile_db.meta

        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if self.p_meta.has_key('creation_date') else 'unknown'
        self.p_meta['samples'] = sorted([s.strip() for s in self.p_meta['samples'].split(',')])

        self.progress.update('Reading clusterings dict')
        self.clusterings = profile_db.db.get_table_as_dict(t.clusterings_table_name)

        self.progress.end()

        profile_db.disconnect()


    def load_views(self):
        profile_db = ProfileDatabase(self.profile_db_path)

        views_table = profile_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            self.views[view] = {'table_name': table_name,
                                'header': profile_db.db.get_table_structure(table_name)[1:],
                                'dict': profile_db.db.get_table_as_dict(table_name)}

        profile_db.disconnect()


class DatabasesMetaclass(ProfileSuperclass, AnnotationSuperclass, object):
    """Essential data to load for a given run"""
    def __init__(self, args, r = run, p = progress):
        self.args = args
        self.run = r
        self.progress = p

        filesnpaths.is_file_exists(args.annotation_db)
        filesnpaths.is_file_exists(args.profile_db)

        is_annotation_and_profile_dbs_compatible(args.annotation_db, args.profile_db)

        AnnotationSuperclass.__init__(self, self.args, self.run, self.progress)
        ProfileSuperclass.__init__(self, self.args, self.run, self.progress)

        self.init_split_sequences()



####################################################################################################
#
#     DATABASES
#
####################################################################################################


class ProfileDatabase:
    """To create an empty profile database and/or access one."""
    def __init__(self, db_path, run=run, progress=progress, quiet = True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if os.path.exists(self.db_path):
            is_profile_db(self.db_path)
            self.db = db.DB(self.db_path, t.profile_db_version)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])
            self.samples = set([s.strip() for s in self.meta['samples'].split(',')])

            self.run.info('Profile database', 'An existing database, %s, has been initiated.' % self.db_path, quiet = self.quiet)
            self.run.info('Samples', self.meta['samples'], quiet = self.quiet)
        else:
            self.db = None


    def create(self, meta_values = {}):
        if os.path.exists(self.db_path):
            raise ConfigError, "anvio will not overwrite an existing profile database. Please choose a different name\
                                or remove the existing database ('%s') first." % (self.db_path)

        if not self.db_path.lower().endswith('.db'):
            raise ConfigError, "Please make sure your output file name has a '.db' extension. anvio developers apologize\
                                for imposing their views on how local databases should be named, and are humbled by your\
                                cooperation."

        self.db = db.DB(self.db_path, t.profile_db_version, new_database = True)

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

        # creating empty default tables
        self.db.create_table(t.clusterings_table_name, t.clusterings_table_structure, t.clusterings_table_types)
        self.db.create_table(t.gene_coverages_table_name, t.gene_coverages_table_structure, t.gene_coverages_table_types)
        self.db.create_table(t.variable_positions_table_name, t.variable_positions_table_structure, t.variable_positions_table_types)
        self.db.create_table(t.views_table_name, t.views_table_structure, t.views_table_types)
        ccollections.create_blank_collections_tables(self.db)

        self.disconnect()

        self.run.info('Annotation database', 'A new database, %s, has been created.' % (self.db_path), quiet = self.quiet)


    def disconnect(self):
        self.db.disconnect()


class AnnotationDatabase:
    """To create an empty annotation database and/or access one."""
    def __init__(self, db_path, run=run, progress=progress, quiet = True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.meta = {}
        self.init()


    def init(self):
        if os.path.exists(self.db_path):
            is_annotation_db(self.db_path)
            self.db = db.DB(self.db_path, t.annotation_db_version)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            self.run.info('Annotation database', 'An existing database, %s, has been initiated.' % self.db_path, quiet = self.quiet)
            self.run.info('Number of contigs', self.meta['num_contigs'], quiet = self.quiet)
            self.run.info('Number of splits', self.meta['num_splits'], quiet = self.quiet)
            self.run.info('Total number of nucleotides', self.meta['total_length'], quiet = self.quiet)
            self.run.info('Split length', self.meta['split_length'], quiet = self.quiet)
        else:
            self.db = None


    def create(self, contigs_fasta, split_length, kmer_size = 4):
        if os.path.exists(self.db_path):
            raise ConfigError, "anvio will not overwrite an existing annotation database. Please choose a different name\
                                or remove the existing database ('%s') first." % (self.db_path)

        if not split_length:
            raise ConfigError, "Creating a new annotation database requires split length information to be\
                                provided. But the AnnotationDatabase class was called to create one without this\
                                bit of information. Not cool."

        if not os.path.exists(contigs_fasta):
            raise ConfigError, "Creating a new annotation database requires a FASTA file with contigs to be provided."


        if not self.db_path.lower().endswith('.db'):
            raise ConfigError, "Please make sure your output file name has a '.db' extension. anvio developers apologize\
                                for imposing their views on how local databases should be named, and are humbled by your\
                                cooperation."

        try:
            split_length = int(split_length)
        except:
            raise ConfigError, "Split size must be an integer."

        try:
            kmer_size = int(kmer_size)
        except:
            raise ConfigError, "K-mer size must be an integer."
        if kmer_size < 2 or kmer_size > 8:
            raise ConfigError, "We like our k-mer sizes between 2 and 8, sorry! (but then you can always change the\
                                source code if you are not happy to be told what you can't do, let us know how it goes!)."

        self.db = db.DB(self.db_path, t.annotation_db_version, new_database = True)

        # know thyself
        self.db.set_meta_value('db_type', 'annotation')
        # this will be the unique information that will be passed downstream whenever this db is used:
        self.db.set_meta_value('annotation_hash', '%08x' % random.randrange(16**8))
        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)

        self.db.create_table(t.contig_sequences_table_name, t.contig_sequences_table_structure, t.contig_sequences_table_types)

        # lets process and store the FASTA file.
        fasta = u.SequenceSource(contigs_fasta)
        db_entries_contig_sequences = []

        contigs_kmer_table = KMerTablesForContigsAndSplits('kmer_contigs', k=kmer_size)
        splits_kmer_table = KMerTablesForContigsAndSplits('kmer_splits', k=kmer_size)

        contigs_info_table = InfoTableForContigs(split_length)
        splits_info_table = InfoTableForSplits()

        while fasta.next():
            contig_length, split_start_stops, contig_gc_content = contigs_info_table.append(fasta.id, fasta.seq)

            contig_kmer_freq = contigs_kmer_table.get_kmer_freq(fasta.seq)

            for order in range(0, len(split_start_stops)):
                start, end = split_start_stops[order]
                split_name = contig.gen_split_name(fasta.id, order)

                # this is very confusing, because both contigs_kmer_table and splits_kmer_able in fact
                # holds kmer values for splits only. in one table, each split has a kmer value of their
                # contigs (to not lose the genomic context while clustering based on kmers), in the other
                # one each split holds its own kmer value.
                contigs_kmer_table.append(split_name, fasta.seq[start:end], kmer_freq = contig_kmer_freq)
                splits_kmer_table.append(split_name, fasta.seq[start:end])

                splits_info_table.append(split_name, fasta.seq[start:end], order, start, end, contig_gc_content, fasta.id)

            db_entries_contig_sequences.append((fasta.id, fasta.seq), )

        self.db.set_meta_value('kmer_size', kmer_size)
        contigs_kmer_table.store(self.db)
        splits_kmer_table.store(self.db)
        contigs_info_table.store(self.db)
        splits_info_table.store(self.db)

        self.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.contig_sequences_table_name, db_entries_contig_sequences)

        # set some useful meta values:
        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value('num_contigs', contigs_info_table.total_contigs)
        self.db.set_meta_value('total_length', contigs_info_table.total_nts)
        self.db.set_meta_value('num_splits', splits_info_table.total_splits)
        self.db.set_meta_value('genes_annotation_source', None)

        # creating empty default tables
        self.db.create_table(t.hmm_hits_info_table_name, t.hmm_hits_info_table_structure, t.hmm_hits_info_table_types)
        self.db.create_table(t.hmm_hits_splits_table_name, t.hmm_hits_splits_table_structure, t.hmm_hits_splits_table_types)
        self.db.create_table(t.hmm_hits_contigs_table_name, t.hmm_hits_contigs_table_structure, t.hmm_hits_contigs_table_types)
        self.db.create_table(t.genes_contigs_table_name, t.genes_contigs_table_structure, t.genes_contigs_table_types)
        self.db.create_table(t.genes_splits_summary_table_name, t.genes_splits_summary_table_structure, t.genes_splits_summary_table_types)
        self.db.create_table(t.genes_splits_table_name, t.genes_splits_table_structure, t.genes_splits_table_types)
        ccollections.create_blank_collections_tables(self.db)

        self.disconnect()

        self.run.info('Annotation database', 'A new database, %s, has been created.' % (self.db_path), quiet = self.quiet)
        self.run.info('Number of contigs', contigs_info_table.total_contigs, quiet = self.quiet)
        self.run.info('Number of contigs', splits_info_table.total_splits, quiet = self.quiet)
        self.run.info('Total number of nucleotides', contigs_info_table.total_nts, quiet = self.quiet)
        self.run.info('Split length', split_length, quiet = self.quiet)


    def disconnect(self):
        self.db.disconnect()


####################################################################################################
#
#     TABLES
#
####################################################################################################


class InfoTableForContigs:
    def __init__(self, split_length):
        self.db_entries = []
        self.total_nts = 0
        self.total_contigs = 0
        self.split_length = split_length


    def append(self, seq_id, sequence):
        sequence_length = len(sequence)
        gc_content = utils.get_GC_content_for_sequence(sequence)

        # how many splits will there be?
        split_start_stops = utils.get_split_start_stops(sequence_length, self.split_length)

        self.total_nts += sequence_length
        self.total_contigs += 1
        db_entry = tuple([seq_id, sequence_length, gc_content, len(split_start_stops)])
        self.db_entries.append(db_entry)

        return (sequence_length, split_start_stops, gc_content)


    def store(self, db):
        db.create_table(t.contigs_info_table_name, t.contigs_info_table_structure, t.contigs_info_table_types)
        if len(self.db_entries):
            db._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.contigs_info_table_name, (','.join(['?'] * len(self.db_entries[0])))), self.db_entries)


class InfoTableForSplits:
    def __init__(self):
        self.db_entries = []
        self.total_splits = 0


    def append(self, seq_id, sequence, order, start, end, parent_gc_content, parent):
        self.total_splits += 1
        sequence_length = len(sequence)
        db_entry = tuple([seq_id, order, start, end, sequence_length, utils.get_GC_content_for_sequence(sequence), parent_gc_content, parent])
        self.db_entries.append(db_entry)


    def store(self, db):
        db.create_table(t.splits_info_table_name, t.splits_info_table_structure, t.splits_info_table_types)
        if len(self.db_entries):
            db._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.splits_info_table_name, (','.join(['?'] * len(self.db_entries[0])))), self.db_entries)


class KMerTablesForContigsAndSplits:
    def __init__(self, table_name, k = 4):
        self.table_name = table_name
        self.kmers_class = kmers.KMers(k)
        self.kmers = sorted(list(self.kmers_class.kmers[k]))

        self.kmer_dict = {}
        self.db_entries = []

        self.kmers_table_structure = ['contig'] + self.kmers
        self.kmers_table_types = ['text'] + ['numeric'] * len(self.kmers)


    def get_kmer_freq(self, sequence):
        return self.kmers_class.get_kmer_frequency(sequence)


    def append(self, seq_id, sequence, kmer_freq = None):
        if not kmer_freq:
            kmer_freq = self.kmers_class.get_kmer_frequency(sequence)

        db_entry = tuple([seq_id] + [kmer_freq[kmer] for kmer in self.kmers])
        self.db_entries.append(db_entry)


    def store(self, db):
        db.create_table(self.table_name, self.kmers_table_structure, self.kmers_table_types)
        db._exec_many('''INSERT INTO %s VALUES (%s)''' % (self.table_name, (','.join(['?'] * len(self.kmers_table_structure)))), self.db_entries)


class TableForViews(Table):
    def __init__(self, db_path, version, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, version, run, progress)

        self.db_entries = []


    def append(self, view_id, target_table):
        self.db_entries.append((view_id, target_table),)


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        profile_db.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.views_table_name, self.db_entries)
        profile_db.disconnect()


class TableForVariability(Table):
    def __init__(self, db_path, version, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, version, run, progress)

        self.num_entries = 0
        self.db_entries = []
        self.set_next_available_id(t.variable_positions_table_name)


    def append(self, profile):
        db_entry = tuple([self.next_id(t.variable_positions_table_name)] + [profile[h] for h in t.variable_positions_table_structure[1:]])
        self.db_entries.append(db_entry)
        self.num_entries += 1
        if self.num_entries % 100 == 0:
            progress.update('Information for %d SNP sites have been added ...' % self.num_entries)


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        profile_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''' % t.variable_positions_table_name, self.db_entries)
        profile_db.disconnect()


class TableForGeneCoverages(Table):
    '''The purpose of this class is to keep coverage values for each gene in contigs for found in a sample.
       Simply, you create an instance from it, keep sending contig instances from contig.py::Contig class along with
       a list of inferred start/stop locations for each reading frame. Once you are done, you call create_gene_coverages_table.'''
    def __init__(self, db_path, version, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, version, run, progress)

        self.genes = []
        self.set_next_available_id(t.gene_coverages_table_name)

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
        self.genes.append({'prot': prot, 'sample_id': sample_id, 'mean_coverage': coverage})


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        db_entries = [tuple([self.next_id(t.gene_coverages_table_name)] + [gene[h] for h in t.gene_coverages_table_structure[1:]]) for gene in self.genes]
        profile_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.gene_coverages_table_name, db_entries)
        profile_db.disconnect()


class TablesForSearches(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        self.debug = False

        Table.__init__(self, self.db_path, t.annotation_db_version, run, progress)

        self.set_next_available_id(t.hmm_hits_contigs_table_name)
        self.set_next_available_id(t.hmm_hits_splits_table_name)


    def populate_search_tables(self, sources = {}):
        if not len(sources):
            import anvio.data.hmm
            sources = anvio.data.hmm.sources

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
        self.delete_entries_for_key('source', source, [t.hmm_hits_info_table_name, t.hmm_hits_contigs_table_name, t.hmm_hits_splits_table_name])

        annotation_db = AnnotationDatabase(self.db_path)

        # push information about this search result into serach_info table.
        db_entries = [source, reference, kind_of_search, ', '.join(all_genes)]
        annotation_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.hmm_hits_info_table_name, db_entries)
        # then populate serach_data table for each contig.
        db_entries = [tuple([self.next_id(t.hmm_hits_contigs_table_name), source] + [v[h] for h in t.hmm_hits_contigs_table_structure[2:]]) for v in search_results_dict.values()]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % t.hmm_hits_contigs_table_name, db_entries)

        db_entries = self.process_splits(source, search_results_dict)
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.hmm_hits_splits_table_name, db_entries)

        annotation_db.disconnect()


    def process_splits(self, source, search_results_dict):
        hits_per_contig = {}
        for hit in search_results_dict.values():
            if hits_per_contig.has_key(hit['contig']):
                hits_per_contig[hit['contig']].append(hit)
            else:
                hits_per_contig[hit['contig']] = [hit]

        db_entries_for_splits = []

        for contig in self.contigs_info:
            if not hits_per_contig.has_key(contig):
                # no hits for this contig. pity!
                continue

            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                # FIXME: this really needs some explanation.
                for hit in hits_per_contig[contig]:
                    if hit['stop'] > start and hit['start'] < stop:
                        gene_length = hit['stop'] - hit['start']
                        # if only a part of the gene is in the split:
                        start_in_split = (start if hit['start'] < start else hit['start']) - start
                        stop_in_split = (stop if hit['stop'] > stop else hit['stop']) - start
                        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length
                        
                        gene_unique_identifier = hashlib.sha224('_'.join([contig, hit['gene_name'], str(hit['start']), str(hit['stop'])])).hexdigest()
                        db_entry = tuple([self.next_id(t.hmm_hits_splits_table_name), source, gene_unique_identifier, hit['gene_name'], split_name, percentage_in_split, hit['e_value']])
                        db_entries_for_splits.append(db_entry)

        return db_entries_for_splits


class TablesForGenes(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, t.annotation_db_version, run, progress)

        # this class keeps track of genes that occur in splits, and responsible
        # for generating the necessary table in the annotation database
        self.genes_in_splits = GenesInSplits()


    def create(self, genes_dict, parser):
        self.genes_dict = genes_dict

        self.sanity_check()

        # oepn connection
        annotation_db = AnnotationDatabase(self.db_path)

        self.splits_info = annotation_db.db.get_table_as_dict(t.splits_info_table_name)

        # test whether there are already genes tables populated
        genes_annotation_source = annotation_db.meta['genes_annotation_source']
        if genes_annotation_source:
            self.run.warning('Previous genes annotation data from "%s" will be replaced with the incoming data' % parser)
            annotation_db.db._exec('''DELETE FROM %s''' % (t.genes_contigs_table_name))
            annotation_db.db._exec('''DELETE FROM %s''' % (t.genes_splits_table_name))
            annotation_db.db._exec('''DELETE FROM %s''' % (t.genes_splits_summary_table_name))

        # set the parser
        annotation_db.db.remove_meta_key_value_pair('genes_annotation_source')
        annotation_db.db.set_meta_value('genes_annotation_source', parser)
        # push raw entries
        db_entries = [tuple([prot] + [self.genes_dict[prot][h] for h in t.genes_contigs_table_structure[1:]]) for prot in self.genes_dict]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''' % t.genes_contigs_table_name, db_entries)
        # disconnect like a pro.
        annotation_db.disconnect()


        # compute and push split taxonomy information.
        self.init_genes_splits_summary_table()


    def sanity_check(self):
        # check whether input matrix dict 
        keys_found = ['prot'] + self.genes_dict.values()[0].keys()
        missing_keys = [key for key in t.genes_contigs_table_structure if key not in keys_found]
        if len(missing_keys):
            raise ConfigError, "Your input lacks one or more header fields to generate a anvio annotation db. Here is\
                                what you are missing: %s. The complete list (and order) of headers in your TAB\
                                delimited matrix file (or dictionary) must follow this: %s." % (', '.join(missing_keys),
                                                                                                ', '.join(t.genes_contigs_table_structure))


        contig_names_in_matrix = set([v['contig'] for v in self.genes_dict.values()])
        contig_names_in_db  = set(self.contigs_info.keys())

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
                                    there is nothing much anvio can do about this issue." %\
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

        contigs_without_annotation = list(set(self.contigs_info.keys()) - set(prots_in_contig.keys()))
        run.info('Percent of contigs annotated', '%.1f%%' % (len(prots_in_contig) * 100.0 / len(self.contigs_info)))

        for contig in contigs_without_annotation:
            prots_in_contig[contig] = set([])

        splits_dict = {}
        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

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


                taxonomy_strings = [tt for tt in taxa if tt]
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
        db_entries = [tuple([split] + [splits_dict[split][h] for h in t.genes_splits_summary_table_structure[1:]]) for split in splits_dict]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % t.genes_splits_summary_table_name, db_entries)
        # push entries for genes in splits table
        db_entries = [tuple([entry_id] + [self.genes_in_splits.splits_to_prots[entry_id][h] for h in t.genes_splits_table_structure[1:]]) for entry_id in self.genes_in_splits.splits_to_prots]
        annotation_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.genes_splits_table_name, db_entries)
        # disconnect
        annotation_db.disconnect()


    def get_consensus_taxonomy_for_split(self, contig, t_level = 't_species', start = 0, stop = sys.maxint):
        """Returns (c, n, t, o) where,
            c: consensus taxonomy (the most common taxonomic call for each gene found in the contig),
            n: total number of genes found in the contig,
            tt: total number of genes with known taxonomy,
            o: number of taxonomic calls that matches the consensus among tt
        """

        response = self.db.cursor.execute("""SELECT %s FROM %s WHERE contig='%s' and stop > %d and start < %d""" % (t_level, t.genes_contigs_table_name, contig, start, stop))
        rows = response.fetchall()

        num_genes = len(rows)
        tax_str_list = [tt[0] for tt in rows if tt[0]]
        distinct_taxa = set(tax_str_list)

        if not len(distinct_taxa):
            return None, num_genes, 0, 0

        if len(distinct_taxa) == 1:
            return distinct_taxa.pop(), num_genes, len(tax_str_list), len(tax_str_list)
        else:
            d = Counter()
            for tt in tax_str_list:
                d[tt] += 1
            consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
            return consensus, num_genes, len(tax_str_list), occurrence


class GenesInSplits:
    def __init__(self):
        self.entry_id = 0
        self.splits_to_prots = {}

    def add(self, split_name, split_start, split_end, prot_id, prot_start, prot_end):

        gene_length = prot_end - prot_start

        if gene_length <= 0:
            raise ConfigError, "dbops.py/GeneInSplits: OK. There is something wrong. We have this gene, '%s',\
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


####################################################################################################
#
#     HELPER FUNCTIONS
#
####################################################################################################


def is_annotation_db(db_path):
    if get_db_type(db_path) != 'annotation':
        raise ConfigError, '"%s" is not a anvio annotation database.' % db_path


def is_profile_db(db_path):
    if get_db_type(db_path) != 'profile':
        raise ConfigError, '"%s" is not a anvio profile database.' % db_path


def get_db_type(db_path):
    try:
        database = db.DB(db_path, None, ignore_version = True)
    except:
        raise ConfigError, 'Are you sure "%s" is a database file? Because, you know, probably\
                            it is not at all..' % db_path

    tables = database.get_table_names()
    if 'self' not in tables:
        database.disconnect()
        raise ConfigError, '"%s" does not seem to be a anvio database...' % db_path

    db_type = database.get_meta_value('db_type')
    database.disconnect()

    return db_type


def is_annotation_and_profile_dbs_compatible(annotation_db_path, profile_db_path):
    is_annotation_db(annotation_db_path)
    is_profile_db(profile_db_path)

    annotation_db = AnnotationDatabase(annotation_db_path)
    profile_db = ProfileDatabase(profile_db_path)

    a_hash = annotation_db.meta['annotation_hash']
    p_hash = profile_db.meta['annotation_hash']
    merged = profile_db.meta['merged']

    annotation_db.disconnect()
    profile_db.disconnect()

    if a_hash != p_hash:
        raise ConfigError, 'The annotation database and the profile database does not\
                            seem to be compatible. More specifically, this annotation\
                            database is not the one that was used when %s generated\
                            this profile database.'\
                                % 'anvi-merge' if merged else 'anvi-profile'

    return True
