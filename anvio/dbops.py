# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to create, access, and/or populate contigs and profile databases.
"""

import os
import sys
import time
import copy
import random
import hashlib
import datetime
import textwrap
from itertools import chain
from collections import Counter

import numpy

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.kmers as kmers
import anvio.terminal as terminal
import anvio.contigops as contigops
import anvio.samplesops as samplesops
import anvio.filesnpaths as filesnpaths
import anvio.genecalling as genecalling
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.parsers import parser_modules
from anvio.tableops import Table
from anvio.constants import codon_to_AA

from anvio.drivers.hmmer import HMMer


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class DBClassFactory:
    """Factory pattern to get the appropriate class for a given anvi'o db type"""
    def __init__(self):
        self.DB_CLASSES = {'profile': ProfileDatabase,
                           'contigs': ContigsDatabase,
                           'pan': PanDatabase}

    def get_db_class(self, db_path):
        db_type = get_db_type(db_path)

        if db_type not in self.DB_CLASSES:
            raise ConfigError, "DBClassFactory speaking. I do not know a class for database type\
                                %s :/ I can deal with these though: '%s'" % (', '.join(self.DB_CLASSES))

        return self.DB_CLASSES[db_type]


    def get_db_object(self, db_path):
        anvio_db_class = self.get_db_class(db_path)
        return anvio_db_class(db_path)


class ContigsSuperclass(object):
    def __init__(self, args, r=run, p=progress):
        self.run = r
        self.progress = p

        self.a_meta = {}

        self.splits_basic_info = {}
        self.splits_taxonomy_dict = {}
        self.split_sequences = {}
        self.contigs_basic_info = {}
        self.contig_sequences = {}

        self.genes_in_contigs_dict = {}
        self.contig_name_to_genes = {}
        self.genes_in_splits = {}
        self.genes_in_splits_summary_dict = {}
        self.genes_in_splits_summary_headers = []
        self.split_name_to_gene_caller_ids_dict = {} # for fast access to all self.genes_in_splits entries for a given split
        self.gene_callers_id_to_split_name_dict = {} # for fast access to a split name that contains a given gene callers id

        self.auxiliary_contigs_data_available = False
        self.nt_positions_info = None

        self.gene_function_call_sources = []
        self.gene_function_calls_dict = {}
        self.gene_function_calls_initiated = False

        self.hmm_sources_info = {}
        self.hmm_searches_dict = {}   # <--- upon initiation, this dict only keeps hmm hits for non-singlecopy
        self.hmm_searches_header = [] #      gene searches... single-copy gene info is accessed through completeness.py

        self.singlecopy_gene_hmm_sources = set([])
        self.non_singlecopy_gene_hmm_sources = set([])

        try:
            self.contigs_db_path = args.contigs_db
        except:
            # ContigsSuperclass class called with args without contigs_db member..
            return

        if not self.contigs_db_path:
            return

        filesnpaths.is_file_exists(self.contigs_db_path)

        self.progress.new('Loading the contigs DB')
        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)

        self.progress.update('Setting contigs self data dict')
        self.a_meta = contigs_db.meta

        self.a_meta['creation_date'] = utils.get_time_to_date(self.a_meta['creation_date']) if 'creation_date' in self.a_meta else 'unknown'

        self.progress.update('Reading contigs basic info')
        self.contigs_basic_info = contigs_db.db.get_table_as_dict(t.contigs_info_table_name, string_the_key=True)

        self.progress.update('Reading splits basic info')
        self.splits_basic_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)

        self.progress.update('Reading genes in contigs table')
        self.genes_in_contigs_dict = contigs_db.db.get_table_as_dict(t.genes_in_contigs_table_name)

        self.progress.update('Populating contig name to gene IDs dict')
        for contig_name in self.contigs_basic_info:
            self.contig_name_to_genes[contig_name] = set([])
        for gene_unique_id in self.genes_in_contigs_dict:
            e = self.genes_in_contigs_dict[gene_unique_id]
            self.contig_name_to_genes[e['contig']].add((gene_unique_id, e['start'], e['stop']), )

        self.progress.update('Reading genes in splits table')
        self.genes_in_splits = contigs_db.db.get_table_as_dict(t.genes_in_splits_table_name)

        self.progress.update('Reading genes in splits summary table')
        self.genes_in_splits_summary_dict = contigs_db.db.get_table_as_dict(t.genes_in_splits_summary_table_name)
        self.genes_in_splits_summary_headers = contigs_db.db.get_table_structure(t.genes_in_splits_summary_table_name)

        self.progress.update('Identifying HMM searches for single-copy genes and others')
        self.hmm_sources_info = contigs_db.db.get_table_as_dict(t.hmm_hits_info_table_name)
        for hmm_source in self.hmm_sources_info:
            self.hmm_sources_info[hmm_source]['genes'] = sorted([g.strip() for g in self.hmm_sources_info[hmm_source]['genes'].split(',')])

        self.singlecopy_gene_hmm_sources = set([s for s in self.hmm_sources_info.keys() if self.hmm_sources_info[s]['search_type'] == 'singlecopy'])
        self.non_singlecopy_gene_hmm_sources = set([s for s in self.hmm_sources_info.keys() if self.hmm_sources_info[s]['search_type'] != 'singlecopy'])

        self.progress.update('Generating "split name" to "gene entry ids" mapping dict')
        for entry_id in self.genes_in_splits:
            split_name = self.genes_in_splits[entry_id]['split']
            if split_name in self.split_name_to_gene_caller_ids_dict:
                self.split_name_to_gene_caller_ids_dict[split_name].add(entry_id)
            else:
                self.split_name_to_gene_caller_ids_dict[split_name] = set([entry_id])

        for split_name in self.splits_basic_info:
            if split_name not in self.split_name_to_gene_caller_ids_dict:
                self.split_name_to_gene_caller_ids_dict[split_name] = set([])

        self.progress.update('Generating "gene caller id" to "split name" mapping dict')
        for entry in self.genes_in_splits.values():
            self.gene_callers_id_to_split_name_dict[entry['gene_callers_id']] = entry['split']

        contigs_db.disconnect()

        self.progress.update('Accessing the auxiliary data file')
        auxiliary_contigs_data_path = ''.join(self.contigs_db_path[:-3]) + '.h5'
        if os.path.exists(auxiliary_contigs_data_path):
            self.auxiliary_contigs_data_available = True
            self.nt_positions_info = auxiliarydataops.AuxiliaryDataForNtPositions(auxiliary_contigs_data_path, self.a_meta['contigs_db_hash'])
            self.progress.end()
        else:
            self.progress.end()
            self.run.warning("Auxiliary contigs data ('%s') is not available. Some operations related to\
                              variability analyses will not be available." % auxiliary_contigs_data_path)

        if self.auxiliary_contigs_data_available:
            self.run.info('Auxiliary Data', 'Found: %s (v. %s)' % (auxiliary_contigs_data_path, anvio.__hdf5__version__))

        self.run.info('Contigs DB', 'Initialized: %s (v. %s)' % (self.contigs_db_path, anvio.__contigs__version__))


    def init_splits_taxonomy(self, t_level = 't_genus'):
        if not self.contigs_db_path:
            return

        if t_level not in t.taxon_names_table_structure[1:]:
            raise ConfigError, "Pretty close. But the taxonomic level '%s' is not known to anvi'o. How about\
                                one of these: %s." % (t_level, ','.join(t.taxon_names_table_structure[1:]))

        self.progress.new('Initializing splits taxonomy')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        splits_taxonomy_table = contigs_db.db.get_table_as_dict(t.splits_taxonomy_table_name)
        taxon_names_table = contigs_db.db.get_table_as_dict(t.taxon_names_table_name)

        for split_name in splits_taxonomy_table:
            taxon_id = splits_taxonomy_table[split_name]['taxon_id']
            if taxon_id:
                if t_level in taxon_names_table[taxon_id] and taxon_names_table[taxon_id][t_level]:
                    self.splits_taxonomy_dict[split_name] = taxon_names_table[taxon_id][t_level]

        contigs_db.disconnect()
        self.progress.end()

        if len(splits_taxonomy_table):
            self.run.info('Taxonomy', 'Initiated for taxonomic level for "%s"' % t_level)


    def init_contig_sequences(self, min_contig_length=0):
        self.progress.new('Loading contig sequences')

        self.progress.update('Identifying contigs shorter than M')
        contigs_shorter_than_M = set([c for c in self.contigs_basic_info if self.contigs_basic_info[c]['length'] < min_contig_length])

        self.progress.update('Reading contig sequences')
        contigs_db = ContigsDatabase(self.contigs_db_path)
        self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name, string_the_key=True)
        contigs_db.disconnect()

        self.progress.update('Filtering out shorter contigs')
        for contig_name in contigs_shorter_than_M:
            self.contig_sequences.pop(contig_name)

        self.progress.end()

        return contigs_shorter_than_M


    def init_split_sequences(self, min_contig_length=0):
        contigs_shorter_than_M = self.init_contig_sequences(min_contig_length)

        self.progress.new('Computing split sequences from contigs')

        self.progress.update('Discarding split names coming from short contigs')
        split_names_to_discard = set([])
        for split_name in self.splits_basic_info:
            if self.splits_basic_info[split_name]['parent'] in contigs_shorter_than_M:
                split_names_to_discard.add(split_name)

        for split_name in split_names_to_discard:
            self.splits_basic_info.pop(split_name)

        self.progress.update('Generating split sequences dict')
        for split_name in self.splits_basic_info:
            split = self.splits_basic_info[split_name]

            if split['parent'] in contigs_shorter_than_M:
                contigs_shorter_than_M.remove(split['parent'])
                continue

            if self.contigs_basic_info[split['parent']]['num_splits'] == 1:
                self.split_sequences[split_name] = self.contig_sequences[split['parent']]['sequence']
            else:
                self.split_sequences[split_name] = self.contig_sequences[split['parent']]['sequence'][split['start']:split['end']]

        self.progress.end()


    def init_non_singlecopy_gene_hmm_sources(self, split_names_of_interest=None, return_each_gene_as_a_layer=False):
        if not self.contigs_db_path or not len(self.non_singlecopy_gene_hmm_sources):
            return

        self.progress.new('Initializing non-single-copy HMM sources')
        self.progress.update('...')

        non_singlecopy_gene_hmm_info_dict = {}
        for source in self.non_singlecopy_gene_hmm_sources:
            non_singlecopy_gene_hmm_info_dict[source] = self.hmm_sources_info[source]

        contigs_db = ContigsDatabase(self.contigs_db_path)
        non_singlecopy_gene_hmm_results_dict = utils.get_filtered_dict(contigs_db.db.get_table_as_dict(t.hmm_hits_splits_table_name), 'source', self.non_singlecopy_gene_hmm_sources)
        hmm_hits_table = utils.get_filtered_dict(contigs_db.db.get_table_as_dict(t.hmm_hits_table_name), 'source', self.non_singlecopy_gene_hmm_sources)

        if split_names_of_interest:
            non_singlecopy_gene_hmm_results_dict = utils.get_filtered_dict(non_singlecopy_gene_hmm_results_dict, 'split', set(split_names_of_interest))

        sources_tmpl = {}

        # the following conditional is pretty critical. here is more info about the difference:
        # https://github.com/meren/anvio/issues/123
        if return_each_gene_as_a_layer:
            for source in self.non_singlecopy_gene_hmm_sources:
                search_type = self.hmm_sources_info[source]['search_type']
                for gene_name in non_singlecopy_gene_hmm_info_dict[source]['genes']:
                    search_term = 'hmmx_%s_%s' % (search_type, gene_name)
                    sources_tmpl[search_term] = 0
                    self.hmm_searches_header.append((search_term, source),)

            # fill all splits with 0s, so this is treated as a numeric column:
            for split_name in split_names_of_interest if split_names_of_interest else self.splits_basic_info:
                self.hmm_searches_dict[split_name] = copy.deepcopy(sources_tmpl)

            for e in non_singlecopy_gene_hmm_results_dict.values():
                hmm_hit = hmm_hits_table[e['hmm_hit_entry_id']]
                search_term = 'hmmx_%s_%s' % (self.hmm_sources_info[e['source']]['search_type'], hmm_hit['gene_name'])
                self.hmm_searches_dict[e['split']][search_term] = 1
        else:
            for source in self.non_singlecopy_gene_hmm_sources:
                search_type = 'hmms_%s' % self.hmm_sources_info[source]['search_type']
                sources_tmpl[source] = []
                self.hmm_searches_header.append((search_type, source),)

            for e in non_singlecopy_gene_hmm_results_dict.values():
                hmm_hit = hmm_hits_table[e['hmm_hit_entry_id']]
                if not e['split'] in self.hmm_searches_dict:
                    self.hmm_searches_dict[e['split']] = copy.deepcopy(sources_tmpl)

                search_type = 'hmms_%s' % self.hmm_sources_info[e['source']]['search_type']

                # populate hmm_searches_dict with hmm_hit and unique identifier (see #180):
                self.hmm_searches_dict[e['split']][source].append((hmm_hit['gene_name'], hmm_hit['gene_unique_identifier']),)

        self.progress.end()


    def get_nt_position_info(self, contig_name, pos_in_contig):
        """This function returns a tuple with three items for each nucleotide position.

            (in_partial_gene_call, in_complete_gene_call, base_pos_in_codon)

        See `init_nt_position_info_dict` for more info."""

        if not self.nt_positions_info:
            raise ConfigError, "get_nt_position_info: I am asked to return stuff, but self.nt_position_info is None!\
                                This may happen if you don't have the '.h5' file for your contigs database in the same\
                                directory with your contigs database. But if you do have it there, then anvi'o really\
                                needs an adult :("

        if not self.nt_positions_info.is_known_contig(contig_name):
            return (0, 0, 0)

        position_info = self.nt_positions_info.get(contig_name)[pos_in_contig]

        if not position_info:
            return (0, 0, 0)
        if position_info == 8:
            return (1, 0, 0)
        if position_info == 4:
            return (0, 1, 1)
        if position_info == 2:
            return (0, 1, 2)
        if position_info == 1:
            return (0, 1, 3)


    def init_functions(self, requested_sources=[], dont_panic=False):
        if not self.contigs_db_path:
            return

        self.progress.new('Initializing functions class')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path)

        gene_function_sources_in_db = set(contigs_db.meta['gene_function_sources'] or [])

        if requested_sources:
            missing_sources = [s for s in requested_sources if s not in gene_function_sources_in_db]
            if len(missing_sources):
                if dont_panic:
                    requested_sources = [s for s in requested_sources if s in gene_function_sources_in_db]
                else:
                    self.progress.end()
                    raise ConfigError, "Some of the functional sources you requested are missing from the contigs database '%s'. Here\
                                        they are (or here it is, whatever): %s." % \
                                                    (self.contigs_db_path, ', '.join(["'%s'" % s for s in missing_sources]))

            hits = contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name,
                                                                  '''source IN (%s)''' % (', '.join(["'%s'" % s for s in requested_sources])),
                                                                  error_if_no_data=False).values()
            self.gene_function_call_sources = requested_sources
        else:
            hits = contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name).values()
            self.gene_function_call_sources = gene_function_sources_in_db

        for hit in hits:
            gene_callers_id = hit['gene_callers_id']
            source = hit['source']
            accession = hit['accession']
            function = hit['function']
            e_value = hit['e_value']

            if gene_callers_id not in self.gene_function_calls_dict:
                self.gene_function_calls_dict[gene_callers_id] = dict([(s, None) for s in self.gene_function_call_sources])

            if self.gene_function_calls_dict[gene_callers_id][source]:
                if self.gene_function_calls_dict[gene_callers_id][source][1] < e_value:
                    # 'what we have:', self.gene_function_calls_dict[gene_callers_id][source]
                    # 'rejected    :', ('%s :: %s' % (function if function else 'unknown', accession), e_value)
                    continue

            entry = (accession, '%s' % (function if function else 'unknown'), e_value)
            self.gene_function_calls_dict[gene_callers_id][source] = entry

        contigs_db.disconnect()

        self.progress.end()

        self.gene_function_calls_initiated = True


    def search_splits_for_gene_functions(self, search_terms, verbose=False, full_report=False):
        if not isinstance(search_terms, list):
            raise ConfigError, "Search terms must be of type 'list'"

        search_terms = [s.strip() for s in search_terms]

        if len([s.strip().lower() for s in search_terms]) != len(set([s.strip().lower() for s in search_terms])):
            raise ConfigError, "Please do not use the same search term twice :/ Becasue, reasons. You know."

        for search_term in search_terms:
            if not len(search_term) >= 3:
                raise ConfigError, "A search term cannot be less than three characters"

        self.run.info('Search terms', '%d found' % (len(search_terms)))
        matching_gene_caller_ids = dict([(search_term, {}) for search_term in search_terms])
        matching_function_calls = dict([(search_term, {}) for search_term in search_terms])
        split_names = dict([(search_term, {}) for search_term in search_terms])
        full_report = []

        if not self.gene_function_calls_initiated:
            self.init_functions()

        contigs_db = ContigsDatabase(self.contigs_db_path)

        for search_term in search_terms:
            self.progress.new('Search function')
            self.progress.update('Searching for term "%s"' % search_term)
            response = contigs_db.db._exec('''select gene_callers_id, source, function from gene_functions where function LIKE "%%''' + search_term + '''%%";''').fetchall()

            full_report.extend([(r[0], r[1], r[2], search_term, self.gene_callers_id_to_split_name_dict[r[0]]) for r in response])

            matching_gene_caller_ids[search_term] = set([m[0] for m in response])
            matching_function_calls[search_term] = list(set([m[2] for m in response]))
            split_names[search_term] = [self.gene_callers_id_to_split_name_dict[gene_callers_id] for gene_callers_id in matching_gene_caller_ids[search_term]]

            self.progress.end()

            if len(matching_gene_caller_ids):
                self.run.info('Matches', '%d unique gene calls contained the search term "%s"' % (len(matching_gene_caller_ids[search_term]), search_term))
                if verbose:
                    self.run.warning('\n'.join(matching_function_calls[search_term][0:25]), header="Matching names for '%s' (up to 25)" % search_term, raw=True, lc='cyan')
            else:
                self.run.info('Matches', 'No matches found the search term "%s"' % (search_term), mc='red')

        contigs_db.disconnect()
        self.progress.end()

        return split_names, full_report


    def get_corresponding_codon_order_in_gene(self, gene_caller_id, contig_name, pos_in_contig):
        """Takes a gene caller id, a contig name, and a nucleotide position in that contig,
           and returns the order of codon the nucleotide matches to."""

        if not isinstance(pos_in_contig, int):
            raise ConfigError, "get_corresponding_codon_order_in_gene :: pos_in_contig must be of type 'int'"

        if not isinstance(gene_caller_id, int):
            raise ConfigError, "get_corresponding_codon_order_in_gene :: gene_caller_id must be of type 'int'"

        gene_call = self.genes_in_contigs_dict[gene_caller_id]

        if contig_name != gene_call['contig']:
            raise ConfigError, 'get_corresponding_codon_order_in_gene :: well, the gene call %d and the contig %s\
                                do not seem to have anything to do with each other :/ This is not a user-level error\
                                something must have gone very wrong somewhere in the code ...' % (gene_caller_id, contig_name)

        if not pos_in_contig >= gene_call['start'] or not pos_in_contig < gene_call['stop']:
            raise ConfigError, "get_corresponding_codon_order_in_gene :: position %d does not occur in gene call %d :(" \
                                                        % (pos_in_contig, gene_caller_id)

        start, stop = gene_call['start'], gene_call['stop']

        gene_length = stop - start
        num_codons_in_gene = int(gene_length / 3)

        if gene_call['direction'] == 'r':
            corresponding_codon_order_in_gene = num_codons_in_gene - int((pos_in_contig - start) / 3) - 1
        else:
            corresponding_codon_order_in_gene = int((pos_in_contig - start) / 3)

        return corresponding_codon_order_in_gene


    def get_gene_start_stops_in_contig(self, contig_name):
        """Return a list of (gene_callers_id, start, stop) tuples for each gene occurring
           in contig_name"""
        return self.contig_name_to_genes[contig_name]


    def get_AA_counts_dict(self, split_names=set([]), contig_names=set([]), gene_caller_ids=set([])):
        """Returns a dictionary of AA counts.

           The dict can be returned for a given collection of split names, contigs names,
           or gene calls. If none of these variables are specified, the dict will contain
           counts for all gene calls in the contigs database"""

        AA_counts_dict = {}

        # nothing to do here if the genes were not called:
        if not self.a_meta['genes_are_called']:
            return AA_counts_dict

        if len([True for v in [split_names, contig_names, gene_caller_ids] if v]) > 1:
            raise ConfigError, "get_AA_counts_dict :: If you want to get AA counts for a specific\
                                set of split names, contig names, or gene call ids, that is totally\
                                fine. But you can't request more than one at a time."

        # we need to understand what genes we're interested in first. it could be genes in
        # a collection, or it could be everything in the contigs database, etc
        gene_calls_of_interest = set([])

        if split_names:
            for split_name in split_names:
                for gene_call_id in self.split_name_to_gene_caller_ids_dict[split_name]:
                    gene_calls_of_interest.add(gene_call_id)
        elif contig_names:
            for contig_name in contig_names:
                for gene_call_id in [t[0] for t in self.contig_name_to_genes[contig_name]]:
                    gene_calls_of_interest.add(gene_call_id)
        elif gene_caller_ids:
            gene_calls_of_interest = set(gene_caller_ids)
        else:
            gene_calls_of_interest = set(self.genes_in_contigs_dict.keys())

        if not len(self.contig_sequences):
            self.init_contig_sequences()

        AAs = []
        for gene_call_id in gene_calls_of_interest:
            gene_call = self.genes_in_contigs_dict[gene_call_id]

            if gene_call['partial']:
                continue

            AAs.extend(utils.get_list_of_AAs_for_gene_call(gene_call, self.contig_sequences))

        AA_counts_dict['AA_counts'] = Counter(AAs)
        AA_counts_dict['total_AAs'] = sum(Counter(AAs).values())
        AA_counts_dict['total_gene_calls'] = len(gene_calls_of_interest)

        # add missing AAs into the dict .. if there are any
        for AA in codon_to_AA.values():
            if AA not in AA_counts_dict['AA_counts']:
                AA_counts_dict['AA_counts'][AA] = 0

        return AA_counts_dict


    def get_corresponding_gene_caller_ids_for_base_position(self, contig_name, pos_in_contig):
        """For a given nucleotide position and contig name, returns all matching gene caller ids"""
        gene_start_stops_in_contig = self.get_gene_start_stops_in_contig(contig_name)

        if not gene_start_stops_in_contig:
            return []

        corresponding_gene_calls = [gene_callers_id for (gene_callers_id, start, stop) in gene_start_stops_in_contig if pos_in_contig >= start and pos_in_contig < stop]

        return corresponding_gene_calls


    def get_sequences_for_gene_callers_ids(self, gene_caller_ids_list, reverse_complement_if_necessary=True):
        if not isinstance(gene_caller_ids_list, list):
            raise ConfigError, "Gene caller's ids must be of type 'list'"

        try:
            gene_caller_ids_list = [int(gene_callers_id) for gene_callers_id in gene_caller_ids_list]
        except:
            raise ConfigError, "List of IDs for gene calls contains non-integer values :/"

        if not len(self.contig_sequences):
            self.init_contig_sequences()

        sequences_dict = {}

        self.progress.new('Getting sequences')
        self.progress.update('...')
        for gene_callers_id in gene_caller_ids_list:
            gene_call = self.genes_in_contigs_dict[gene_callers_id]

            contig_name = gene_call['contig']
            start, stop = gene_call['start'], gene_call['stop']
            direction = gene_call['direction']
            sequence = self.contig_sequences[contig_name]['sequence'][start:stop]

            if direction == 'r' and reverse_complement_if_necessary:
                sequence = utils.rev_comp(sequence)
                rev_compd = "True"
            else:
                rev_compd = "False"

            sequences_dict[gene_callers_id] = {'sequence': sequence,
                                               'contig': contig_name,
                                               'start': start,
                                               'stop': stop,
                                               'direction': direction,
                                               'rev_compd': rev_compd,
                                               'length': stop - start}

        self.progress.end()

        return (gene_caller_ids_list, sequences_dict)


    def gen_FASTA_file_of_sequences_for_gene_caller_ids(self, gene_caller_ids_list=[], output_file_path=None, wrap=120, simple_headers=False, rna_alphabet=False):
        if not output_file_path:
            raise ConfigError, "gen_FASTA_file_of_sequences_for_gene_caller_ids function requires an explicit output file path.\
                                Anvi'o does not know how you managed to come here, but please go back and come again."

        filesnpaths.is_output_file_writable(output_file_path)

        if not isinstance(wrap, int):
            raise ConfigError, '"wrap" has to be an integer instance'
        if wrap == 0:
            wrap = None
        if wrap and wrap <= 20:
            raise ConfigError, 'Value for wrap must be larger than 20. Yes. Rules.'

        if not gene_caller_ids_list:
            gene_caller_ids_list = self.genes_in_contigs_dict.keys()
            self.run.warning("You did not provide any gene caller ids. As a result, anvi'o will give you back sequences for every\
                              %d gene call stored in the contigs database. %s" % (len(gene_caller_ids_list), ' Brace yourself.' if len(gene_caller_ids_list) > 10000 else ''))

        gene_caller_ids_list, sequences_dict = self.get_sequences_for_gene_callers_ids(gene_caller_ids_list)

        output = open(output_file_path, 'w')

        self.progress.new('Storing sequences')
        self.progress.update('...')
        for gene_callers_id in gene_caller_ids_list:
            entry = sequences_dict[gene_callers_id]

            if simple_headers:
                header = '%d' % (gene_callers_id)
            else:
                header = '%d|' % (gene_callers_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['contig', 'start', 'stop', 'direction', 'rev_compd', 'length']])

            if rna_alphabet:
                sequence = entry['sequence'].replace('T', 'U')
            else:
                sequence = entry['sequence']

            if wrap:
                sequence = textwrap.fill(sequence, wrap, break_on_hyphens=False)

            output.write('>%s\n' % header)
            output.write('%s\n' % sequence)

        output.close()

        self.progress.end()
        self.run.info('Output', output_file_path)


    def gen_TAB_delimited_file_for_split_taxonomies(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        if not self.a_meta['taxonomy_source']:
            raise ConfigError, "There is no taxonomy source in the contigs database :/"

        if not len(self.splits_taxonomy_dict):
            self.init_splits_taxonomy()

        if not len(self.splits_taxonomy_dict):
            raise ConfigError, "The splits taxonomy is empty. There is nothing to report. Could it be\
                                possible the taxonomy caller you used did not assign any taxonomy to\
                                anything?"

        self.run.info("Taxonomy", "Annotations for %d of %d total splits are recovered" % (len(self.splits_taxonomy_dict), len(self.splits_basic_info)))

        output = open(output_file_path, 'w')
        for split_name in sorted(self.splits_taxonomy_dict.keys()):
            output.write('{0}\t{1}\n'.format(split_name, self.splits_taxonomy_dict[split_name]))
        output.close()

        self.run.info("Output", output_file_path)


class PanSuperclass(object):
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        self.genome_names = []
        self.protein_clusters = {}
        self.protein_cluster_names = set([])
        self.protein_clusters_gene_alignments = {}
        self.protein_clusters_gene_alignments_available = False
        self.protein_clusters_function_sources = []
        self.protein_clusters_functions_dict = {}
        self.clusterings = {}
        self.views = {}
        self.collection_profile = {}

        self.num_protein_clusters = None
        self.num_genes_in_protein_clusters = None

        self.genomes_storage_is_available = False
        self.genomes_storage_has_functions = False
        self.functions_initialized = False

        try:
            self.pan_db_path = args.pan_db
        except:
            self.run.warning('PanSuperclass class called with args without pan_db_path member! Returning prematurely.')
            return

        filesnpaths.is_file_exists(self.pan_db_path)

        self.progress.new('Initializing the pan database superclass')

        self.progress.update('Creating an instance of the pan database')
        pan_db = PanDatabase(self.pan_db_path)

        self.progress.update('Setting profile self data dict')
        self.p_meta = pan_db.meta

        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if 'creation_date' in self.p_meta else 'unknown'
        self.p_meta['genome_names'] = sorted([s.strip() for s in self.p_meta['external_genome_names'].split(',') + self.p_meta['internal_genome_names'].split(',') if s])
        self.p_meta['num_genomes'] = len(self.p_meta['genome_names'])
        self.genome_names = self.p_meta['genome_names']
        self.protein_clusters_gene_alignments_available = self.p_meta['gene_alignments_computed']

        if self.p_meta['PCs_clustered']:
            self.p_meta['available_clusterings'] = sorted([s.strip() for s in self.p_meta['available_clusterings'].split(',')])
            self.clusterings = pan_db.db.get_table_as_dict(t.clusterings_table_name)
        else:
            self.p_meta['available_clusterings'] = None
            self.p_meta['default_clustering'] = None
            self.clusterings = None

        # recover all protein cluster names so others can access to this information
        # without having to initialize anything
        self.protein_cluster_names = set(pan_db.db.get_single_column_from_table(t.pan_protein_clusters_table_name, 'protein_cluster_id'))

        pan_db.disconnect()

        # create an instance of states table
        self.states_table = TablesForStates(self.pan_db_path)

        self.progress.end()

        if 'genomes_storage' in args.__dict__ and args.genomes_storage:
            self.genomes_storage = auxiliarydataops.GenomesDataStorage(args.genomes_storage,
                                                                       self.p_meta['genomes_storage_hash'],
                                                                       genome_names_to_focus=self.p_meta['genome_names'])
            self.genomes_storage_is_available = True
            self.genomes_storage_has_functions = self.genomes_storage.functions_are_available

        self.run.info('Pan DB', 'Initialized: %s (v. %s)' % (self.pan_db_path, anvio.__pan__version__))


    def init_protein_clusters_functions(self):
        self.progress.new('Initializing functions for protein clusters')
        self.progress.update('...')
        if not self.protein_clusters:
            raise ConfigError, "init_protein_clusters_functions is speaking! You called this function before you initialized\
                                protein clusters :/ One of us does not know what they're doing :("

        if not self.genomes_storage_has_functions:
            self.progress.end()
            self.run.warning("Genomes storage does not have any info about gene functions. Certain parts of the pangenomic\
                              workflow will not be accessible.")
            return

        # FIXME WE HAVE TO STORE AVAILABLE FUNCTIONS IN GENOMES STORAGE ATTRs!!!! THIS IS RIDICULOUS
        self.protein_clusters_function_sources = set([])
        for protein_cluster_id in self.protein_clusters:
            self.protein_clusters_functions_dict[protein_cluster_id] = {}
            for genome_name in self.genome_names:
                self.protein_clusters_functions_dict[protein_cluster_id][genome_name] = {}
                for gene_callers_id in self.protein_clusters[protein_cluster_id][genome_name]:
                    functions = self.genomes_storage.get_gene_functions(genome_name, gene_callers_id)
                    self.protein_clusters_functions_dict[protein_cluster_id][genome_name][gene_callers_id] = functions

                    if functions:
                        self.protein_clusters_function_sources.update(functions.keys())

        self.functions_initialized = True

        self.progress.end()


    def init_additional_layer_data(self):
        """Recover additional data stored in the pan database `additional data` table."""

        self.progress.new('Initializing additional layer data')
        self.progress.update('...')
        pan_db = PanDatabase(self.pan_db_path)
        self.additional_layers_dict = pan_db.db.get_table_as_dict('additional_data')
        self.additional_layers_headers = pan_db.db.get_meta_value('additional_data_headers').split(',')
        pan_db.disconnect()

        if len([h for h in self.additional_layers_headers if h not in self.additional_layers_dict.values()[0].keys()]):
            self.progress.end()
            raise ConfigError, "Something that should never happen happened :( At least one additional data header that\
                                appears in the self table of your pan database is not in the dictionary recovered for this\
                                data from another table. Anvi'o needs an adult :("

        # In fact we are done here since we have our `additional_layers_dict` all filled up with sweet data.
        # But if functions are initialized, we can also get a summary of protein clusters based on whether most
        # genes in them were annotated with known functions or not for a given annotation source. Of course,
        # for this to happen, we need to check whther functions were initialied prior to the call to
        # `init_additional_layer_data`.
        if not self.functions_initialized:
            # no? k. bye.
            self.progress.end()
            return

        # too many shitty nested loops here, but it is quite efficient since we work only with a dict
        # in memory
        for annotation_source in self.protein_clusters_function_sources:
            if annotation_source == 'COG_CATEGORY':
                # we don't need this one
                continue

            self.progress.update('Computing known/unknown dict for %s' % annotation_source)
            for protein_cluster_id in self.additional_layers_dict:
                hits = Counter({})
                for genome_id in self.protein_clusters_functions_dict[protein_cluster_id]:
                    for gene_callers_id in self.protein_clusters_functions_dict[protein_cluster_id][genome_id]:
                        if annotation_source in self.protein_clusters_functions_dict[protein_cluster_id][genome_id][gene_callers_id]:
                            hits[self.protein_clusters_functions_dict[protein_cluster_id][genome_id][gene_callers_id][annotation_source][0]] += 1
                        else:
                            hits['UNKNOWN'] += 1

                if not hits or hits.most_common()[0][0] == 'UNKNOWN':
                    self.additional_layers_dict[protein_cluster_id][annotation_source] = 'UNKNOWN'
                else:
                    self.additional_layers_dict[protein_cluster_id][annotation_source] = 'KNOWN'

            self.additional_layers_headers.append(annotation_source)

        self.progress.end()


    def init_protein_clusters(self):
        """Initializes the protein_clusters dictionary.

           At the end, the structure of this dictionary looks like this:

               {
                'PC_1': {'Genome_1': [gene_1, gene_2, (...)],
                         'Genome_2': [],
                         'Genome_3': [gene_1, gene_2],
                         (...)}
                'PC_2': {(...)},
                (...)
               }

          This function also initializes alignment summaries for each gene
          in each protein cluster. That information is stored in the dict
          `self.protein_clusters_gene_alignments`.
        """

        self.progress.new('Initializing protein clusters')
        self.progress.update('...')

        pan_db = PanDatabase(self.pan_db_path)

        protein_clusters_long_list = pan_db.db.get_table_as_dict(t.pan_protein_clusters_table_name)

        for entry in protein_clusters_long_list.values():
            genome_name = entry['genome_name']
            gene_callers_id = entry['gene_caller_id']
            protein_cluster_id = entry['protein_cluster_id']

            if protein_cluster_id not in self.protein_clusters:
                self.protein_clusters[protein_cluster_id] = {}
                for g in self.genome_names:
                    self.protein_clusters[protein_cluster_id][g] = []

            if self.protein_clusters_gene_alignments_available:
                if genome_name not in self.protein_clusters_gene_alignments:
                    self.protein_clusters_gene_alignments[genome_name] = {}

                self.protein_clusters_gene_alignments[genome_name][gene_callers_id] = entry['alignment_summary']

            self.protein_clusters[protein_cluster_id][genome_name].append(gene_callers_id)

        pan_db.disconnect()
        self.progress.end()


    def load_pan_views(self, splits_of_interest=None):
        pan_db = PanDatabase(self.pan_db_path)

        views_table = pan_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            self.views[view] = {'table_name': table_name,
                                'header': pan_db.db.get_table_structure(table_name)[1:],
                                'dict': pan_db.db.get_table_as_dict(table_name, keys_of_interest=splits_of_interest)}

        pan_db.disconnect()


    def get_summary_for_PCs_list(self, protein_cluster_ids):
        summary = {'genomes_contributing': set([]), 'num_gene_calls': 0, 'num_PCs': 0, 'functions': {}}

        if self.functions_initialized:
            for source in self.protein_clusters_function_sources:
                summary['functions'].update({source: Counter({})})

        for protein_cluster_id in protein_cluster_ids:
            single_summary = self.get_summary_for_PC_id(protein_cluster_id)
            summary['num_PCs'] += 1
            summary['genomes_contributing'] = summary['genomes_contributing'].union(single_summary['genomes_contributing'])
            summary['num_gene_calls'] += single_summary['num_gene_calls']

            if self.functions_initialized:
                for source in self.protein_clusters_function_sources:
                    for function in single_summary['functions'][source]:
                        summary['functions'][source][function] += single_summary['functions'][source][function]

        summary['genomes_contributing'] = sorted(list(summary['genomes_contributing']))

        return summary


    def get_summary_for_PC_id(self, protein_cluster_id):

        summary = {'genomes_contributing': set([]), 'num_gene_calls': 0, 'functions': {}}

        for genome_name in self.protein_clusters[protein_cluster_id]:
            num_gene_calls = len(self.protein_clusters[protein_cluster_id][genome_name])
            if num_gene_calls:
                summary['genomes_contributing'].add(genome_name)
                summary['num_gene_calls'] += num_gene_calls

        if self.functions_initialized:
            for source in self.protein_clusters_function_sources:
                summary['functions'].update({source: Counter({})})

            functions_dict = self.protein_clusters_functions_dict[protein_cluster_id]
            for genome_name in functions_dict:
                for gene_callers_id in functions_dict[genome_name]:
                    for source in functions_dict[genome_name][gene_callers_id]:
                        for function in functions_dict[genome_name][gene_callers_id][source]:
                            summary['functions'][source][function] += 1

        return summary


    def init_collection_profile(self, collection_name):
        pan_db = PanDatabase(self.pan_db_path)

        if not self.protein_clusters:
            raise ConfigError, "init_collection_profile wants to initialize the PC collection profile for '%s', but the\
                                the protein clusters dict is kinda empty. Someone forgot to initialize something maybe?" % \
                                        collection_name

        # get trimmed collection and bins_info dictionaries
        collection, bins_info, self.protein_clusters_in_pan_db_but_not_binned \
                    = self.collections.get_trimmed_dicts(collection_name, set(self.protein_clusters.keys()))

        # currenlty we are not really doing anything with this one, but we will be filling this up with
        # all sorts of amazing later.
        for bin_id in collection:
            self.collection_profile[bin_id] = {}

        self.progress.end()
        pan_db.disconnect()

        return collection, bins_info


class ProfileSuperclass(object):
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        self.gene_coverages_dict = None
        self.split_coverage_values = None

        self.split_names = set([])
        self.clusterings = {}
        self.views = {}
        self.collection_profile = {}

        try:
            self.profile_db_path = args.profile_db
        except:
            self.run.warning('ProfileSuperclass class called with args without profile_db member')
            return

        if not self.profile_db_path:
            return

        filesnpaths.is_file_exists(self.profile_db_path)

        self.progress.new('Initializing the profile database superclass')

        self.progress.update('Loading split names')
        self.split_names = get_split_names_in_profile_db(self.profile_db_path)

        self.progress.update('Creating an instance of the profile database')
        profile_db = ProfileDatabase(self.profile_db_path)

        self.progress.update('Setting profile self data dict')
        self.p_meta = profile_db.meta

        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if 'creation_date' in self.p_meta else 'unknown'
        self.p_meta['samples'] = sorted([s.strip() for s in self.p_meta['samples'].split(',')])
        self.p_meta['num_samples'] = len(self.p_meta['samples'])

        if self.p_meta['contigs_clustered'] and 'available_clusterings' in self.p_meta:
            self.p_meta['available_clusterings'] = sorted([s.strip() for s in self.p_meta['available_clusterings'].split(',')])
            self.clusterings = profile_db.db.get_table_as_dict(t.clusterings_table_name)
        elif self.p_meta['contigs_clustered'] and 'available_clusterings' not in self.p_meta:
            self.progress.end()
            self.run.warning("Your profile database thinks the hierarchical clustering was done, yet it contains no entries\
                              for any hierarchical clustering results. This is not good. Something must have gone wrong\
                              somewhere :/ To be on the safe side, anvi'o will assume this profile database has no\
                              clustering (which is literally the case, by the way, it is just the database itself is\
                              confused about that fact --it happens to the best of us).")
            self.progress.new('Initializing the profile database superclass')

            self.p_meta['contigs_clustered'] = False
            self.p_meta['available_clusterings'] = None
            self.p_meta['default_clustering'] = None
            self.clusterings = None
        else:
            self.p_meta['available_clusterings'] = None
            self.p_meta['default_clustering'] = None
            self.clusterings = None

        profile_db.disconnect()

        self.progress.update('Accessing the auxiliary data file')
        auxiliary_data_path = os.path.join(os.path.dirname(self.profile_db_path), 'AUXILIARY-DATA.h5')
        if not os.path.exists(auxiliary_data_path):
            self.auxiliary_profile_data_available = False
        else:
            self.auxiliary_profile_data_available = True
            self.split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(auxiliary_data_path, self.p_meta['contigs_db_hash'])

        self.progress.end()

        if self.auxiliary_profile_data_available:
            self.run.info('Auxiliary Data', 'Found: %s (v. %s)' % (auxiliary_data_path, anvio.__hdf5__version__))
        self.run.info('Profile DB', 'Initialized: %s (v. %s)' % (self.profile_db_path, anvio.__profile__version__))


    def init_gene_coverages_dict(self):
        if self.p_meta['blank']:
            # this is a blank profile, there is nothing to init here
            return

        if not self.a_meta['genes_are_called']:
            # genes were not identified/annotated
            return

        profile_db = ProfileDatabase(self.profile_db_path, quiet=True)

        self.progress.new('Reading gene coverages table')
        self.progress.update('...')
        gene_coverages_table = profile_db.db.get_table_as_dict(t.gene_coverages_table_name)
        profile_db.disconnect()

        if not len(gene_coverages_table):
            self.progress.end()
            self.run.warning('Something came up, please read this carefuly: your contigs database does\
                              contain information for open reading frames in your contigs. however, the\
                              gene coverages table in the profile database is empty. This happens when you\
                              annotate your contigs database with gene/function calls *after* you have\
                              profiled your samples. If you are OK with this situation, you can simply\
                              ignore this message. If you need to access to this information, you must\
                              re-profile your samples (and merge them) using your most update contigs\
                              database :/ Sorry!')
            return

        self.gene_coverages_dict = {}
        for gene_coverage_entry in gene_coverages_table.values():
            gene_callers_id = gene_coverage_entry['gene_callers_id']

            if gene_callers_id not in self.gene_coverages_dict:
                self.gene_coverages_dict[gene_callers_id] = {}

            self.gene_coverages_dict[gene_callers_id][gene_coverage_entry['sample_id']] = gene_coverage_entry['mean_coverage']

        self.progress.end()


    def get_variability_information_for_split(self, split_name, return_outliers=False, return_raw_results=False):
        if not split_name in self.split_names:
            raise ConfigError, "get_variability_information_for_split: The split name '%s' does not seem to be\
                                represented in this profile database. Are you sure you are looking for it\
                                in the right database?" % split_name

        profile_db = ProfileDatabase(self.profile_db_path)
        split_variability_information = profile_db.db.get_some_rows_from_table_as_dict(t.variable_nts_table_name, '''split_name = "%s"''' % split_name, error_if_no_data=False).values()
        profile_db.disconnect()

        if return_raw_results:
            return split_variability_information

        # they want pretty stuff...
        d = {}

        for sample_name in self.p_meta['samples']:
            d[sample_name] = {'variability': {0: {}, 1: {}, 2: {}, 3: {}}, 'competing_nucleotides': {}}

        for e in split_variability_information:
            if not return_outliers and e['cov_outlier_in_contig']:
                continue

            d[e['sample_id']]['variability'][e['base_pos_in_codon']][e['pos']] = e['departure_from_reference']
            d[e['sample_id']]['competing_nucleotides'][e['pos']] = e['competing_nts']

        return d


    def init_collection_profile(self, collection_name):
        profile_db = ProfileDatabase(self.profile_db_path, quiet=True)

        # get trimmed collection and bins_info dictionaries
        collection, bins_info, self.split_names_in_profile_db_but_not_binned \
                    = self.collections.get_trimmed_dicts(collection_name, self.split_names)

        for bin_id in collection:
            self.collection_profile[bin_id] = {}

        table_names = [] if self.p_meta['blank'] else [table_name for table_name in t.atomic_data_table_structure[1:-1]]

        samples_template = dict([(s, []) for s in self.p_meta['samples']])

        # anonymous function to convert single profile table dicts compatible with merged ones (#155):
        SINGLE_P = lambda d: dict([(s, dict([(self.p_meta['samples'][0], v) for v in d[s].values()])) for s in d])

        self.progress.new('Initializing the collection profile for "%s" ...' % collection_name)
        for table_name in table_names:
            self.progress.update('Populating collection profile for each "view" ... %s' % table_name)
            if self.p_meta['merged']:
                table_data = profile_db.db.get_table_as_dict('%s_splits' % table_name, omit_parent_column=True)
            else:
                table_data = SINGLE_P(profile_db.db.get_table_as_dict('atomic_data_splits', columns_of_interest=table_name, omit_parent_column=True))

            for bin_id in collection:
                # populate averages per bin
                averages = copy.deepcopy(samples_template)
                for split_name in collection[bin_id]:
                    if split_name not in table_data:
                        continue

                    for sample_name in samples_template:
                        averages[sample_name].append(table_data[split_name][sample_name])

                # finalize averages per bin:
                for sample_name in samples_template:
                    averages[sample_name] = numpy.mean(averages[sample_name])

                self.collection_profile[bin_id][table_name] = averages

        # generating precent recruitment of each bin plus __splits_not_binned__ in each sample:
        if self.p_meta['merged']:
            coverage_table_data = profile_db.db.get_table_as_dict('mean_coverage_splits', omit_parent_column=True)
        else:
            coverage_table_data = SINGLE_P(profile_db.db.get_table_as_dict('atomic_data_splits', columns_of_interest="mean_coverage", omit_parent_column=True))

        self.bin_percent_recruitment_per_sample = {}
        if self.p_meta['blank']:
            pass
        else:
            for sample in self.p_meta['samples']:
                percents = {}
                all_coverages_in_sample = sum([d[sample] for d in coverage_table_data.values()])

                for bin_id in collection:
                    bin_coverages_in_sample = sum([coverage_table_data[split_name][sample] for split_name in collection[bin_id]])
                    percents[bin_id] = bin_coverages_in_sample * 100 / all_coverages_in_sample

                splits_not_binned_coverages_in_sample = sum([coverage_table_data[split_name][sample] for split_name in self.split_names_in_profile_db_but_not_binned])
                percents['__splits_not_binned__'] = splits_not_binned_coverages_in_sample * 100 / all_coverages_in_sample
                self.bin_percent_recruitment_per_sample[sample] = percents

        self.progress.end()
        profile_db.disconnect()

        return collection, bins_info


    def load_views(self, splits_of_interest=None):
        profile_db = ProfileDatabase(self.profile_db_path)

        views_table = profile_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            self.views[view] = {'table_name': table_name,
                                'header': profile_db.db.get_table_structure(table_name)[1:],
                                'dict': profile_db.db.get_table_as_dict(table_name, keys_of_interest=splits_of_interest)}

        profile_db.disconnect()


class DatabasesMetaclass(ProfileSuperclass, ContigsSuperclass, object):
    """Essential data to load for a given run"""
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        filesnpaths.is_file_exists(args.contigs_db)
        filesnpaths.is_file_exists(args.profile_db)

        is_profile_db_and_contigs_db_compatible(args.profile_db, args.contigs_db)

        ContigsSuperclass.__init__(self, self.args, self.run, self.progress)
        ProfileSuperclass.__init__(self, self.args, self.run, self.progress)

        self.init_split_sequences()
        self.init_gene_coverages_dict()


####################################################################################################
#
#     DATABASES
#
####################################################################################################

class ProfileDatabase:
    """To create an empty profile database and/or access one."""
    def __init__(self, db_path, run=run, progress=progress, quiet=True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if os.path.exists(self.db_path):
            is_profile_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__profile__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            for key in ['min_contig_length', 'SNVs_profiled', 'AA_frequencies_profiled', 'min_coverage_for_variability', 'merged', 'blank', 'contigs_clustered', 'report_variability_full', 'num_contigs', 'num_splits', 'total_length', 'total_reads_mapped']:
                try:
                    self.meta[key] = int(self.meta[key])
                except:
                    pass

            self.samples = set([s.strip() for s in self.meta['samples'].split(',')])

            self.run.info('Profile database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
            self.run.info('Samples', self.meta['samples'], quiet=self.quiet)
        else:
            self.db = None


    def create(self, meta_values={}):
        is_db_ok_to_create(self.db_path, 'profile')

        self.db = db.DB(self.db_path, anvio.__profile__version__, new_database=True)

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

        # creating empty default tables
        self.db.create_table(t.clusterings_table_name, t.clusterings_table_structure, t.clusterings_table_types)
        self.db.create_table(t.gene_coverages_table_name, t.gene_coverages_table_structure, t.gene_coverages_table_types)
        self.db.create_table(t.variable_nts_table_name, t.variable_nts_table_structure, t.variable_nts_table_types)
        self.db.create_table(t.variable_aas_table_name, t.variable_aas_table_structure, t.variable_aas_table_types)
        self.db.create_table(t.views_table_name, t.views_table_structure, t.views_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)

        self.disconnect()

        self.run.info('Profile database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


class PanDatabase:
    """To create an empty pan database, and/or access to one."""
    def __init__(self, db_path, run=run, progress=progress, quiet=True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if os.path.exists(self.db_path):
            is_pan_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__pan__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            for key in ['num_genomes', 'pc_min_occurrence', 'use_ncbi_blast', 'diamond_sensitive', 'exclude_partial_gene_calls', 'num_protein_clusters', 'num_genes_in_protein_clusters', 'gene_alignments_computed']:
                try:
                    self.meta[key] = int(self.meta[key])
                except:
                    pass

            for key in ['min_percent_identity', 'maxbit', 'mcl_inflation']:
                try:
                    self.meta[key] = float(self.meta[key])
                except:
                    pass

            self.internal_genomes = [s.strip() for s in self.meta['internal_genome_names'].split(',')]
            self.external_genomes = [s.strip() for s in self.meta['external_genome_names'].split(',')]
            self.genomes = self.internal_genomes + self.external_genomes

            self.run.info('Pan database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
            self.run.info('Genomes', '%d found' % len(self.genomes), quiet=self.quiet)
        else:
            self.db = None


    def create(self, meta_values={}):
        is_db_ok_to_create(self.db_path, 'pan')

        self.db = db.DB(self.db_path, anvio.__pan__version__, new_database=True)

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

        # know thyself
        self.db.set_meta_value('db_type', 'pan')

        # creating empty default tables for pan specific operations:
        self.db.create_table(t.pan_protein_clusters_table_name, t.pan_protein_clusters_table_structure, t.pan_protein_clusters_table_types)

        # creating empty default tables for standard anvi'o profiles
        self.db.create_table(t.clusterings_table_name, t.clusterings_table_structure, t.clusterings_table_types)
        self.db.create_table(t.views_table_name, t.views_table_structure, t.views_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)

        self.disconnect()

        self.run.info('Pan database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


class ContigsDatabase:
    """To create an empty contigs database and/or access one."""
    def __init__(self, db_path, run=run, progress=progress, quiet=True, skip_init=False):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.meta = {}

        if not skip_init:
            self.init()


    def init(self):
        if os.path.exists(self.db_path):
            is_contigs_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__contigs__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            for key in ['split_length', 'kmer_size', 'total_length', 'num_splits', 'num_contigs', 'genes_are_called']:
                self.meta[key] = int(self.meta[key])

            self.meta['gene_function_sources'] = [s.strip() for s in self.meta['gene_function_sources'].split(',')] if self.meta['gene_function_sources'] else None

            if 'creation_date' not in self.meta:
                raise ConfigError, "The contigs database ('%s') seems to be corrupted :/ This happens if the process that\
                                    that generates the database ends prematurely. Most probably, you will need to generate\
                                    the contigs database from scratch. Sorry!" % (self.db_path)

            self.run.info('Contigs database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
            self.run.info('Number of contigs', self.meta['num_contigs'], quiet=self.quiet)
            self.run.info('Number of splits', self.meta['num_splits'], quiet=self.quiet)
            self.run.info('Total number of nucleotides', self.meta['total_length'], quiet=self.quiet)
            self.run.info('Split length', self.meta['split_length'], quiet=self.quiet)
        else:
            self.db = None


    def create(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        contigs_fasta = A('contigs_fasta')
        split_length = A('split_length')
        kmer_size = A('kmer_size')
        skip_gene_calling = A('skip_gene_calling')
        external_gene_calls = A('external_gene_calls')
        skip_mindful_splitting = A('skip_mindful_splitting')
        debug = A('debug')

        is_db_ok_to_create(self.db_path, 'contigs')

        if external_gene_calls:
            filesnpaths.is_file_exists(external_gene_calls)

        if external_gene_calls and skip_gene_calling:
            raise ConfigError, "You provided a file for external gene calls, and used requested gene calling to be\
                                skipped. Please make up your mind."

        filesnpaths.is_file_fasta_formatted(contigs_fasta)

        # go throught he FASTA file to make sure there are no surprises with deflines and sequence lengths.
        self.progress.new('Checking deflines and contig lengths')
        self.progress.update('tick tock ...')
        fasta = u.SequenceSource(contigs_fasta)
        while fasta.next():
            if not utils.check_contig_names(fasta.id, dont_raise=True):
                self.progress.end()
                raise ConfigError, "At least one of the deflines in your FASTA File does not comply with the 'simple deflines'\
                                    requirement of anvi'o. You can either use the script `anvi-script-reformat-fasta` to take\
                                    care of this issue, or read this section in the tutorial to understand the reason behind\
                                    this requirement (anvi'o is very upset for making you do this): %s" % \
                                        ('http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file')

            if len(fasta.seq) < kmer_size:
                self.progress.end()
                raise ConfigError, "At least one of the contigs in your input FASTA '%s' is shorter than the k-mer size. The k\
                                    is %d, and your contig is like %d :/ Anvi'o will not judge you for whatever you are doing\
                                    with such short contigs, but the length of each contig must be at least as long as your `k` for\
                                    k-mer analyis. You can use the script `anvi-script-reformat-fasta` to get rid of very short\
                                    contigs if you like." % (contigs_fasta, kmer_size, len(fasta.seq))
        fasta.close()
        self.progress.end()

        all_ids_in_FASTA = utils.get_all_ids_from_fasta(contigs_fasta)
        if len(all_ids_in_FASTA) != len(set(all_ids_in_FASTA)):
            raise ConfigError, "Every contig in the input FASTA file must have a unique ID. You know..."

        if not split_length:
            raise ConfigError, "Creating a new contigs database requires split length information to be\
                                provided. But the ContigsDatabase class was called to create one without this\
                                bit of information. Not cool."

        if not os.path.exists(contigs_fasta):
            raise ConfigError, "Creating a new contigs database requires a FASTA file with contigs to be provided."

        try:
            split_length = int(split_length)
        except:
            raise ConfigError, "Split size must be an integer."

        if split_length <= 0:
            split_length = sys.maxsize

        try:
            kmer_size = int(kmer_size)
        except:
            raise ConfigError, "K-mer size must be an integer."
        if kmer_size < 2 or kmer_size > 8:
            raise ConfigError, "We like our k-mer sizes between 2 and 8, sorry! (but then you can always change the\
                                source code if you are not happy to be told what you can't do, let us know how it goes!)."

        if skip_gene_calling:
            skip_mindful_splitting = True

        self.db = db.DB(self.db_path, anvio.__contigs__version__, new_database=True)

        # creating empty default tables
        self.db.create_table(t.hmm_hits_table_name, t.hmm_hits_table_structure, t.hmm_hits_table_types)
        self.db.create_table(t.hmm_hits_info_table_name, t.hmm_hits_info_table_structure, t.hmm_hits_info_table_types)
        self.db.create_table(t.hmm_hits_splits_table_name, t.hmm_hits_splits_table_structure, t.hmm_hits_splits_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.genes_in_contigs_table_name, t.genes_in_contigs_table_structure, t.genes_in_contigs_table_types)
        self.db.create_table(t.genes_in_splits_table_name, t.genes_in_splits_table_structure, t.genes_in_splits_table_types)
        self.db.create_table(t.splits_taxonomy_table_name, t.splits_taxonomy_table_structure, t.splits_taxonomy_table_types)
        self.db.create_table(t.taxon_names_table_name, t.taxon_names_table_structure, t.taxon_names_table_types)
        self.db.create_table(t.genes_taxonomy_table_name, t.genes_taxonomy_table_structure, t.genes_taxonomy_table_types)
        self.db.create_table(t.contig_sequences_table_name, t.contig_sequences_table_structure, t.contig_sequences_table_types)
        self.db.create_table(t.gene_function_calls_table_name, t.gene_function_calls_table_structure, t.gene_function_calls_table_types)
        self.db.create_table(t.gene_protein_sequences_table_name, t.gene_protein_sequences_table_structure, t.gene_protein_sequences_table_types)
        self.db.create_table(t.genes_in_splits_summary_table_name, t.genes_in_splits_summary_table_structure, t.genes_in_splits_summary_table_types)

        # know thyself
        self.db.set_meta_value('db_type', 'contigs')
        # this will be the unique information that will be passed downstream whenever this db is used:
        contigs_db_hash = '%08x' % random.randrange(16**8)
        self.db.set_meta_value('contigs_db_hash', contigs_db_hash)
        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)

        # first things first: do the gene calling on contigs. this part is important. we are doing the
        # gene calling first. so we understand wher genes start and end. this information will guide the
        # arrangement of the breakpoint of splits
        genes_in_contigs_dict = {}
        contig_name_to_gene_start_stops = {}
        if not skip_gene_calling:
            # temporarily disconnect to perform gene calls
            self.db.disconnect()

            gene_calls_tables = TablesForGeneCalls(self.db_path, contigs_fasta, debug=debug)

            # if the user provided a file for external gene calls, use it. otherwise do the gene calling yourself.
            if external_gene_calls:
                gene_calls_tables.use_external_gene_calls_to_populate_genes_in_contigs_table(external_gene_calls)
            else:
                gene_calls_tables.call_genes_and_populate_genes_in_contigs_table()

            # reconnect and learn about what's done
            self.db = db.DB(self.db_path, anvio.__contigs__version__)

            genes_in_contigs_dict = self.db.get_table_as_dict(t.genes_in_contigs_table_name)

            for gene_unique_id in genes_in_contigs_dict:
                e = genes_in_contigs_dict[gene_unique_id]
                if e['contig'] not in contig_name_to_gene_start_stops:
                    contig_name_to_gene_start_stops[e['contig']] = set([])

                contig_name_to_gene_start_stops[e['contig']].add((gene_unique_id, e['start'], e['stop']), )


        # here we will process each item in the contigs fasta file.
        fasta = u.SequenceSource(contigs_fasta)
        db_entries_contig_sequences = []

        contigs_kmer_table = KMerTablesForContigsAndSplits('kmer_contigs', k=kmer_size)
        splits_kmer_table = KMerTablesForContigsAndSplits('kmer_splits', k=kmer_size)

        # create an instance from AuxiliaryDataForNtPositions to store information
        # about each nt position while looping over contigs
        auxiliary_contigs_data_path = ''.join(self.db_path[:-3]) + '.h5'
        nt_positions_auxiliary = auxiliarydataops.AuxiliaryDataForNtPositions(auxiliary_contigs_data_path, contigs_db_hash, create_new=True)

        contigs_info_table = InfoTableForContigs(split_length)
        splits_info_table = InfoTableForSplits()

        recovered_split_lengths = []

        # THE INFAMOUS GEN CONTGS DB LOOP (because it is so costly, we call it South Loop)
        self.progress.new('South Loop')
        fasta.reset()
        while fasta.next():
            contig_name = fasta.id
            contig_sequence = fasta.seq

            self.progress.update('Contig "%d" ' % fasta.pos)

            genes_in_contig = contig_name_to_gene_start_stops[contig_name] if contig_name in contig_name_to_gene_start_stops else set([])

            self.progress.append('has %d genes, ' % len(genes_in_contig))
            if skip_mindful_splitting:
                contig_length, split_start_stops, contig_gc_content = contigs_info_table.append(contig_name, contig_sequence, set([]))
            else:
                contig_length, split_start_stops, contig_gc_content = contigs_info_table.append(contig_name, contig_sequence, genes_in_contig)

            # let's keep an eye on the returned split lengths
            if len(split_start_stops) > 1:
                recovered_split_lengths.extend([s[1] - s[0] for s in split_start_stops])

            self.progress.append('and %d nts. Now computing: auxiliary ... ' % contig_length)
            if genes_in_contig:
                nt_position_info_list = self.compress_nt_position_info(contig_length, genes_in_contig, genes_in_contigs_dict)
                nt_positions_auxiliary.append(contig_name, nt_position_info_list)

            contig_kmer_freq = contigs_kmer_table.get_kmer_freq(contig_sequence)

            self.progress.append('k-mers ...')
            for order in range(0, len(split_start_stops)):
                start, end = split_start_stops[order]
                split_name = contigops.gen_split_name(contig_name, order)

                # this is very confusing, because both contigs_kmer_table and splits_kmer_able in fact
                # holds kmer values for splits only. in one table, each split has a kmer value of their
                # contigs (to not lose the genomic context while clustering based on kmers), in the other
                # one each split holds its own kmer value.
                contigs_kmer_table.append(split_name, contig_sequence[start:end], kmer_freq=contig_kmer_freq)
                splits_kmer_table.append(split_name, contig_sequence[start:end])

                splits_info_table.append(split_name, contig_sequence[start:end], order, start, end, contig_gc_content, contig_name)

            db_entries_contig_sequences.append((contig_name, contig_sequence), )
        nt_positions_auxiliary.close()
        self.progress.end()

        self.db.set_meta_value('kmer_size', kmer_size)
        contigs_kmer_table.store(self.db)
        splits_kmer_table.store(self.db)
        contigs_info_table.store(self.db)
        splits_info_table.store(self.db)

        self.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.contig_sequences_table_name, db_entries_contig_sequences)

        # set some useful meta values:
        self.db.set_meta_value('num_contigs', contigs_info_table.total_contigs)
        self.db.set_meta_value('total_length', contigs_info_table.total_nts)
        self.db.set_meta_value('num_splits', splits_info_table.total_splits)
        self.db.set_meta_value('taxonomy_source', None)
        self.db.set_meta_value('gene_function_sources', None)
        self.db.set_meta_value('genes_are_called', (not skip_gene_calling))
        self.db.set_meta_value('splits_consider_gene_calls', (not skip_mindful_splitting))
        self.db.set_meta_value('creation_date', time.time())
        self.disconnect()

        if not skip_gene_calling:
            gene_calls_tables.populate_genes_in_splits_tables()

        self.run.info('Contigs database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of contigs', contigs_info_table.total_contigs, quiet=self.quiet)
        self.run.info('Number of splits', splits_info_table.total_splits, quiet=self.quiet)
        self.run.info('Total number of nucleotides', contigs_info_table.total_nts, quiet=self.quiet)
        self.run.info('Gene calling step skipped', skip_gene_calling, quiet=self.quiet)
        self.run.info("Splits broke genes (non-mindful mode)", skip_mindful_splitting, quiet=self.quiet)
        self.run.info('Desired split length (what the user wanted)', split_length, quiet=self.quiet)
        self.run.info("Average split length (wnat anvi'o gave back)", (int(round(numpy.mean(recovered_split_lengths)))) \
                                                                        if recovered_split_lengths \
                                                                            else "(Anvi'o did not create any splits)", quiet=self.quiet)


    def compress_nt_position_info(self, contig_length, genes_in_contig, genes_in_contigs_dict):
        """This function compresses information regarding each nucleotide position in a given contig
           into a small int. Every nucleotide position is represented by four bits depending on whether
           they occur in a complete opoen reading frame, and which base they correspond to in a codon.

                0000
                ||||
                ||| \
                |||   Third codon?
                || \
                ||   Second codon?
                | \
                |   First codon?
                 \
                   Whether the position in a partial gene call

           8: int('1000', 2); nt position is in a partial gene call
           4: int('0100', 2); nt position is in a complete gene call, and is at the 1st position in the codon
           2: int('0010', 2); nt position is in a complete gene call, and is at the 2nd position in the codon
           1: int('0001', 2); nt position is in a complete gene call, and is at the 3rd position in the codon
        """

        # first we create a list of zeros for each position of the contig
        nt_position_info_list = [0] * contig_length

        for gene_unique_id, start, stop in genes_in_contig:
            gene_call = genes_in_contigs_dict[gene_unique_id]

            # if the gene call is a partial one, meaning the gene was cut at the beginning or
            # at the end of the contig, we are not going to try to make sense of synonymous /
            # non-synonmous bases in that. the clients who wish to use these variables must also
            # be careful about the difference
            if gene_call['partial']:
                for nt_position in range(start, stop):
                    nt_position_info_list[nt_position] = 8
                continue

            if gene_call['direction'] == 'f':
                for nt_position in range(start, stop, 3):
                    nt_position_info_list[nt_position] = 4
                    nt_position_info_list[nt_position + 1] = 2
                    nt_position_info_list[nt_position + 2] = 1
            elif gene_call['direction'] == 'r':
                for nt_position in range(stop - 1, start - 1, -3):
                    nt_position_info_list[nt_position] = 4
                    nt_position_info_list[nt_position - 1] = 2
                    nt_position_info_list[nt_position - 2] = 1

        return nt_position_info_list


    def disconnect(self):
        self.db.disconnect()


class SamplesInformationDatabase:
    """To create an empty samples information database and/or access one.

       The purpose of this database is to deal with sample-specific information. Such as
       how should samples be organized in the interactive interface, or what environmental
       data available about them?
    """
    def __init__(self, db_path, run=run, progress=progress, quiet=True):
        self.db = None
        self.db_path = db_path

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.meta = {}
        self.init()


    def init(self):
        if not self.db_path:
            raise ConfigError, "When SamplesInformationDatabase is called, the db_path parameter cannot be\
                                'None' type :/"

        if os.path.exists(self.db_path):
            is_samples_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__samples__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])
            self.samples = set([s.strip() for s in self.meta['samples'].split(',')])
            self.sample_names_for_order = set([s.strip() for s in self.meta['sample_names_for_order'].split(',')]) \
                                                if self.meta['sample_names_for_order'] else self.samples
            self.samples_information_default_layer_order = self.meta['samples_information_default_layer_order'].split(',')

            self.run.info('Samples information database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
        else:
            self.db = None


    def get_samples_information_and_order_dicts(self):
        if not self.db:
            raise ConfigError, "The samples database has not been initialized. You are doing something wrong :/"

        samples = samplesops.SamplesInformation(run=self.run, progress=self.progress, quiet=self.quiet)

        samples_information_dict = samples.recover_samples_information_dict(self.db.get_table_as_dict(t.samples_information_table_name, error_if_no_data=False),
                                                                            self.db.get_table_as_dict(t.samples_attribute_aliases_table_name, error_if_no_data=False))
        samples_order_dict = self.db.get_table_as_dict(t.samples_order_table_name)

        return samples_information_dict, samples_order_dict


    def get_samples_information_default_layer_order(self):
        if not self.db:
            raise ConfigError, "The samples database has not been initialized. You are doing something wrong :/"

        return self.samples_information_default_layer_order


    def create(self, samples_information_path=None, samples_order_path=None):
        if not samples_information_path and not samples_order_path:
            raise ConfigError, "You must declare at least one of the input files to create a samples information\
                                database. Neither samples information, nor samples order file has been passed to\
                                the class :("

        is_db_ok_to_create(self.db_path, 'samples')

        samples = samplesops.SamplesInformation(run=self.run, progress=self.progress, quiet=self.quiet)
        samples.populate_from_input_files(samples_information_path, samples_order_path)

        self.db = db.DB(self.db_path, anvio.__samples__version__, new_database=True)

        # know thyself
        self.db.set_meta_value('db_type', 'samples_information')

        # set some useful meta values:
        self.db.set_meta_value('creation_date', time.time())

        # first create the easy one: the samples_order table.
        available_orders = samples.samples_order_dict.keys()
        db_entries = [(attribute, samples.samples_order_dict[attribute]['basic'], samples.samples_order_dict[attribute]['newick']) for attribute in samples.samples_order_dict]
        self.db.create_table(t.samples_order_table_name, t.samples_order_table_structure, t.samples_order_table_types)
        self.db._exec_many('''INSERT INTO %s VALUES (?,?,?)''' % t.samples_order_table_name, db_entries)
        self.db.set_meta_value('available_orders', ','.join(available_orders))

        # then create the table that holds aliases for sample attributes:
        self.db.create_table(t.samples_attribute_aliases_table_name, t.samples_attribute_aliases_table_structure, t.samples_attribute_aliases_table_types)
        db_entries = sorted([(alias, samples.aliases_to_attributes_dict[alias]) for alias in samples.aliases_to_attributes_dict])
        self.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.samples_attribute_aliases_table_name, db_entries)

        # then, create the harder one: the samples_information table.
        aliases = sorted(samples.aliases_to_attributes_dict.keys())
        samples_information_table_structure = ['samples'] + sorted(aliases)
        samples_information_table_types = ['str'] + ['str'] * len(aliases)
        self.db.create_table(t.samples_information_table_name, samples_information_table_structure, samples_information_table_types)
        db_entries = [tuple([sample] + [samples.samples_information_dict[sample][h] for h in samples_information_table_structure[1:]]) for sample in samples.samples_information_dict]
        self.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.samples_information_table_name, ','.join(['?'] * len(samples_information_table_structure))), db_entries)

        # store samples described into the self table
        self.db.set_meta_value('samples', ','.join(samples.sample_names) if samples.sample_names else None)
        self.db.set_meta_value('sample_names_for_order', ','.join(samples.sample_names_in_samples_order_file) if samples.sample_names_in_samples_order_file else None)
        self.db.set_meta_value('samples_information_default_layer_order', ','.join(samples.samples_information_default_layer_order) if hasattr(samples, 'samples_information_default_layer_order') else None)

        self.disconnect()

        self.run.info('Samples information database', 'A new samples information database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of samples', len(samples.sample_names), quiet=self.quiet)
        self.run.info('Number of organizations', len(available_orders), quiet=self.quiet)

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


    def append(self, seq_id, sequence, gene_start_stops=None):
        sequence_length = len(sequence)
        gc_content = utils.get_GC_content_for_sequence(sequence)

        # how many splits will there be?
        split_start_stops = utils.get_split_start_stops(sequence_length, self.split_length, gene_start_stops)

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
    def __init__(self, table_name, k=4):
        self.table_name = table_name
        self.kmers_class = kmers.KMers(k)
        self.kmers = sorted(list(self.kmers_class.kmers[k]))

        self.kmer_dict = {}
        self.db_entries = []

        self.kmers_table_structure = ['contig'] + self.kmers
        self.kmers_table_types = ['text'] + ['numeric'] * len(self.kmers)


    def get_kmer_freq(self, sequence):
        return self.kmers_class.get_kmer_frequency(sequence, dist_metric_safe=True)


    def append(self, seq_id, sequence, kmer_freq=None):
        if not kmer_freq:
            kmer_freq = self.kmers_class.get_kmer_frequency(sequence, dist_metric_safe=True)

        db_entry = tuple([seq_id] + [kmer_freq[kmer] for kmer in self.kmers])
        self.db_entries.append(db_entry)


    def store(self, db):
        db.create_table(self.table_name, self.kmers_table_structure, self.kmers_table_types)
        db._exec_many('''INSERT INTO %s VALUES (%s)''' % (self.table_name, (','.join(['?'] * len(self.kmers_table_structure)))), self.db_entries)


class TablesForViews(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run, progress)


    def create_new_view(self, data_dict, table_name, table_structure, table_types, view_name=None):
        """Creates a new view table, and adds an entry for it into the 'views' table.

        Entries in 'views' table appear in various places in the interface. However, we also generate
        view tables to store the type of data we do not wish to display on interfaces, but be able
        access from various other modules. A good example to this is the clustering recipes. When we
        profile a sample, we treat every stplit as their own entity with respect to their mean coverage.
        Although it is great for visualization purposes, it is not useful for clustering purposes since in
        most cases we wish splits to stay together in clustering output. Hence, we create a mean_coverage_splits
        table, where each split holds their own coverage, and we create a mean_coverage_contigs table where each
        split has the coverage of their parent. Clearly the second table is not useful to display. When a table
        is not added as an entry to the 'views' table, then it only exists in the database for other purposes
        than displaying it.

        If a new view does not have a 'view_id', it is not added the 'views' table to provide that flexibility.
        """

        anvio_db = DBClassFactory().get_db_object(self.db_path)

        views_in_db = anvio_db.db.get_table_as_dict(t.views_table_name)

        if view_name and view_name in views_in_db:
            raise ConfigError, "TablesForViews speaking: Yo yo yo. You already have a view in the db called '%s'.\
                                You can't create another one before you get rid of the existing one, because rules."\
                                                                        % view_name

        # first create the data table:
        anvio_db.db.drop_table(table_name)
        anvio_db.db.create_table(table_name, table_structure, table_types)
        db_entries = [tuple([item] + [data_dict[item][h] for h in table_structure[1:]]) for item in data_dict]
        anvio_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (table_name, ','.join(['?'] * len(table_structure))), db_entries)

        if view_name:
            anvio_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % t.views_table_name, (view_name, table_name))

        anvio_db.disconnect()


    def remove(self, view_name, table_names_to_blank=[]):
        anvio_db = DBClassFactory().get_db_object(self.db_path)
        anvio_db.db._exec('''DELETE FROM %s WHERE view_id = "%s"''' % (t.views_table_name, view_name))
        for table_name in table_names_to_blank:
            anvio_db.db._exec('''DELETE FROM %s''' % table_name)
        anvio_db.disconnect()


class TableForVariability(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run=self.run, progress=self.progress)

        self.num_entries = 0
        self.db_entries = []
        self.set_next_available_id(t.variable_nts_table_name)


    def append(self, profile):
        db_entry = tuple([self.next_id(t.variable_nts_table_name)] + [profile[h] for h in t.variable_nts_table_structure[1:]])
        self.db_entries.append(db_entry)
        self.num_entries += 1
        if self.num_entries % 100 == 0:
            self.progress.update('Information for %d SNV sites have been added ...' % self.num_entries)


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        profile_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''' % t.variable_nts_table_name, self.db_entries)
        profile_db.disconnect()


class AA_counts(ContigsSuperclass):
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.output_file_path = A('output_file')
        self.collection_name = A('collection_name')
        self.bin_ids_file_path = A('bin_ids_file')
        self.contigs_of_interest_file_path = A('contigs_of_interest')
        self.genes_of_interest_file_path = A('gene_caller_ids')

        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)

        self.counts_dict = {}

        # init contigs bro
        ContigsSuperclass.__init__(self, self.args, self.run, self.progress)

        error_msg = "You mixed up optional stuff :/ Please read the help."
        if self.profile_db_path:
            if self.contigs_of_interest_file_path or self.genes_of_interest_file_path:
                raise ConfigError, error_msg
            self.__AA_counts_for_bins()
        elif self.contigs_of_interest_file_path:
            if self.profile_db_path or self.genes_of_interest_file_path:
                raise ConfigError, error_msg
            self.__AA_counts_for_contigs()
        elif self.genes_of_interest_file_path:
            if self.profile_db_path or self.contigs_of_interest_file_path:
                raise ConfigError, error_msg
            self.__AA_counts_for_genes()
        else:
            self.__AA_counts_for_the_contigs_db()


    def __AA_counts_for_bins(self):
        if not self.collection_name:
            raise ConfigError, "You must declare a collection name along with the profile database."

        profile_db = ProfileDatabase(self.profile_db_path)
        collections_info_table = profile_db.db.get_table_as_dict(t.collections_info_table_name)
        collections_splits_table = profile_db.db.get_table_as_dict(t.collections_splits_table_name)
        profile_db.disconnect()

        if not len(collections_info_table):
            raise ConfigError, "There are no collections stored in the profile database :/"

        if not self.collection_name in collections_info_table:
            valid_collections = ', '.join(collections_info_table.keys())
            raise ConfigError, "'%s' is not a valid collection name. But %s: '%s'." \
                                    % (self.collection_name,
                                       'these are' if len(valid_collections) > 1 else 'this is',
                                       valid_collections)

        bin_names_in_collection = collections_info_table[self.collection_name]['bin_names'].split(',')

        if self.bin_ids_file_path:
            filesnpaths.is_file_exists(self.bin_ids_file_path)
            bin_names_of_interest = [line.strip() for line in open(self.bin_ids_file_path).readlines()]

            missing_bins = [b for b in bin_names_of_interest if b not in bin_names_in_collection]
            if len(missing_bins):
                raise ConfigError, "Some bin names you declared do not appear to be in the collection %s." \
                                            % self.collection_name
        else:
            bin_names_of_interest = bin_names_in_collection

        collection_dict = utils.get_filtered_dict(collections_splits_table, 'collection_name', set([self.collection_name]))
        collection_dict = utils.get_filtered_dict(collection_dict, 'bin_name', set(bin_names_of_interest))

        split_name_per_bin_dict = {}
        for bin_name in bin_names_of_interest:
            split_name_per_bin_dict[bin_name] = set([])

        for e in collection_dict.values():
            split_name_per_bin_dict[e['bin_name']].add(e['split'])

        for bin_name in bin_names_of_interest:
            self.counts_dict[bin_name] = self.get_AA_counts_dict(split_names=set(split_name_per_bin_dict[bin_name]))['AA_counts']


    def __AA_counts_for_contigs(self):
        filesnpaths.is_file_exists(self.contigs_of_interest_file_path)

        contigs_of_interest = [line.strip() for line in open(self.contigs_of_interest_file_path).readlines()]

        missing_contigs = [True for c in contigs_of_interest if c not in self.contigs_basic_info]
        if missing_contigs:
            raise ConfigError, "Some contig names you declared do not seem to be present in the contigs\
                                database :("

        for contig_name in contigs_of_interest:
            self.counts_dict[contig_name] = self.get_AA_counts_dict(contig_names=set([contig_name]))['AA_counts']


    def __AA_counts_for_genes(self):
        filesnpaths.is_file_exists(self.genes_of_interest_file_path)

        try:
            genes_of_interest = [int(line.strip()) for line in open(self.genes_of_interest_file_path).readlines()]
        except:
            raise ConfigError, "Gene call ids in your genes of interest file does not resemble anvi'o gene\
                                call ids (I tried to int them, and it didn't work!)"

        for gene_call in genes_of_interest:
            self.counts_dict[gene_call] = self.get_AA_counts_dict(gene_caller_ids=set([gene_call]))['AA_counts']


    def __AA_counts_for_the_contigs_db(self):
        self.counts_dict[self.args.contigs_db] = self.get_AA_counts_dict()['AA_counts']


    def report(self):
        if self.args.output_file:
            header = ['source'] + sorted(self.counts_dict.values()[0].keys())
            utils.store_dict_as_TAB_delimited_file(self.counts_dict, self.args.output_file, header)
            self.run.info('Output', self.args.output_file)

        return self.counts_dict


class TableForAAFrequencies(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run=self.run, progress=self.progress)

        self.num_entries = 0
        self.db_entries = []
        self.set_next_available_id(t.variable_aas_table_name)


    def append(self, profile):
        db_entry = tuple([self.next_id(t.variable_aas_table_name)] + [profile[h] for h in t.variable_aas_table_structure[1:]])
        self.db_entries.append(db_entry)
        self.num_entries += 1
        if self.num_entries % 100 == 0:
            self.progress.update('Information for %d codons have been added ...' % self.num_entries)


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (t.variable_aas_table_name, ','.join(['?'] * len(t.variable_aas_table_structure))), self.db_entries)
        profile_db.disconnect()


class TableForGeneCoverages(Table):
    '''The purpose of this class is to keep coverage values for each gene in contigs for found in a sample.
       Simply, you create an instance from it, keep sending contig instances from contig.py::Contig class along with
       a list of inferred start/stop locations for each reading frame. Once you are done, you call create_gene_coverages_table.'''
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run, progress)

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

        for gene_callers_id, start, stop in start_stop_pos_list:
            gene_coverage = numpy.mean(self.contig_coverages[contig.name][start:stop])
            self.add_gene_entry(gene_callers_id, sample_id, gene_coverage)


    def add_gene_entry(self, gene_callers_id, sample_id, coverage):
        self.genes.append({'gene_callers_id': gene_callers_id, 'sample_id': sample_id, 'mean_coverage': coverage})


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        db_entries = [tuple([self.next_id(t.gene_coverages_table_name)] + [gene[h] for h in t.gene_coverages_table_structure[1:]]) for gene in self.genes]
        profile_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.gene_coverages_table_name, db_entries)
        profile_db.disconnect()


class TablesForGeneCalls(Table):
    def __init__(self, db_path, contigs_fasta=None, run=run, progress=progress, debug=False):
        self.db_path = db_path
        self.contigs_fasta = contigs_fasta
        self.debug = debug

        is_contigs_db(self.db_path)

        if self.contigs_fasta:
            filesnpaths.is_file_exists(self.contigs_fasta)
            filesnpaths.is_file_fasta_formatted(self.contigs_fasta)


    def check_gene_calls_dict(self, gene_calls_dict):
        if not isinstance(gene_calls_dict, type({})):
            raise ConfigError, "Gene calls dict must be a dict instance :/"

        try:
            [int(g) for g in gene_calls_dict.keys()]
        except ValueError:
            raise ConfigError, "Keys of a gene calls dict must be integers!"

        if False in map(lambda x: x['direction'] in ['f', 'r'], gene_calls_dict.values()):
            raise ConfigError, "The values in 'direction' column can't be anything but 'f' (for forward)\
                                or 'r' (for reverse). You have other stuff, and it is not cool."

        if False in map(lambda x: x['stop'] > x['start'], gene_calls_dict.values()):
            raise ConfigError, "For each gene call, the stop position must be bigger than the start position.\
                                Your gene calls dict does not conform to that. If you have reverse gene calls\
                                you must use the 'direction' column to declare that."

        if False in map(lambda x: (x['stop'] - float(x['start'])) % 3.0 == 0, gene_calls_dict.values()):
            raise ConfigError, "Something is wrong with your gene calls. For every gene call, the (stop - start)\
                                should be multiply of 3. It is not the case for all, which is a deal breaker."


    def use_external_gene_calls_to_populate_genes_in_contigs_table(self, input_file_path):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress, simple=True)

        # take care of gene calls dict
        gene_calls_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_path,
                                                                     expected_fields=t.genes_in_contigs_table_structure,
                                                                     only_expected_fields=True,
                                                                     column_mapping=[int, str, int, int, str, int, str, str])

        # recover protein sequences. during this operation we are going to have to read all contig sequences
        # into the damn memory. anvi'o is doing a pretty bad job with memory management :(
        protein_sequences = {}

        contig_sequences = {}
        fasta = u.SequenceSource(self.contigs_fasta)
        while fasta.next():
            contig_sequences[fasta.id] = fasta.seq
        fasta.close()

        number_of_impartial_gene_calls = 0
        for gene_callers_id in gene_calls_dict:
            gene_call = gene_calls_dict[gene_callers_id]
            contig_name = gene_call['contig']

            if contig_name not in contig_sequences:
                raise ConfigError, "You are in big trouble :( The contig name '%s' in your external gene callers file\
                                    does not appear to be in the contigs FASTA file. How did this happen?" % contig_name

            if gene_call['partial']:
                protein_sequences[gene_callers_id] = ''
                number_of_impartial_gene_calls += 1
                continue

            sequence = contig_sequences[contig_name][gene_call['start']:gene_call['stop']]
            if gene_call['direction'] == 'r':
                sequence = utils.rev_comp(sequence)

            protein_sequences[gene_callers_id] = utils.get_DNA_sequence_translated(sequence, gene_callers_id)

        # populate genes_in_contigs, and gene_protein_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, protein_sequences)

        if number_of_impartial_gene_calls:
            self.run.warning('%d of your %d gene calls were impartial, hence the translated protein sequences for those\
                              were not stored in the database.' % (number_of_impartial_gene_calls, len(gene_calls_dict)))


    def call_genes_and_populate_genes_in_contigs_table(self, gene_caller='prodigal'):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress, simple=True)

        # get gene calls and protein sequences
        gene_calls_dict, protein_sequences = self.run_gene_caller(gene_caller)

        # make sure the returning gene calls dict is proper
        self.check_gene_calls_dict(gene_calls_dict)

        # populate genes_in_contigs, and gene_protein_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, protein_sequences)


    def run_gene_caller(self, gene_caller='prodigal'):
        """Runs gene caller, and returns gene_calls_dict, and protein sequences."""
        remove_fasta_after_processing = False

        if not self.contigs_fasta:
            self.contigs_fasta = self.export_sequences_table_in_db_into_FASTA_file()
            remove_fasta_after_processing = True

        if self.debug:
            self.run.info_single('--debug flag is [ON], which means temporary directories generated by\
                                 this run will not be removed', nl_after=2)

        gene_caller = genecalling.GeneCaller(self.contigs_fasta, gene_caller=gene_caller, debug=self.debug)

        gene_calls_dict, protein_sequences = gene_caller.process()

        if not self.debug and remove_fasta_after_processing:
            os.remove(self.contigs_fasta)

        return gene_calls_dict, protein_sequences


    def populate_genes_in_contigs_table(self, gene_calls_dict, protein_sequences):
        contigs_db = db.DB(self.db_path, anvio.__contigs__version__)

        # we know there is nothing there, but lets make double sure
        contigs_db._exec('''DELETE FROM %s''' % (t.genes_in_contigs_table_name))
        contigs_db._exec('''DELETE FROM %s''' % (t.gene_protein_sequences_table_name))

        self.progress.new('Processing')
        self.progress.update('Entering %d gene calls into the db ...' % (len(gene_calls_dict)))

        db_entries = [tuple([entry_id] + [gene_calls_dict[entry_id][h] for h in t.genes_in_contigs_table_structure[1:]]) for entry_id in gene_calls_dict]
        contigs_db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % t.genes_in_contigs_table_name, db_entries)

        db_entries = [tuple([entry_id] + [protein_sequences[entry_id]]) for entry_id in gene_calls_dict]
        contigs_db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.gene_protein_sequences_table_name, db_entries)

        self.progress.end()

        contigs_db.disconnect()


    def populate_genes_in_splits_tables(self):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)
        self.init_gene_calls_dict()

        genes_in_splits = GenesInSplits()
        # build a dictionary for fast access to all proteins identified within a contig
        gene_calls_in_contigs_dict = {}
        for gene_callers_id in self.gene_calls_dict:
            contig = self.gene_calls_dict[gene_callers_id]['contig']
            if contig in gene_calls_in_contigs_dict:
                gene_calls_in_contigs_dict[contig].add(gene_callers_id)
            else:
                gene_calls_in_contigs_dict[contig] = set([gene_callers_id])

        contigs_without_any_gene_calls = list(set(self.contigs_info.keys()) - set(gene_calls_in_contigs_dict.keys()))
        run.info('Contigs with at least one gene call', '%d of %d (%.1f%%)' % (len(gene_calls_in_contigs_dict),
                                                                               len(self.contigs_info),
                                                                               len(gene_calls_in_contigs_dict) * 100.0 / len(self.contigs_info)))

        for contig in contigs_without_any_gene_calls:
            gene_calls_in_contigs_dict[contig] = set([])

        splits_dict = {}
        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                gene_start_stops = []
                # here we go through all genes in the contig and identify the all the ones that happen to be in
                # this particular split to generate summarized info for each split. BUT one important that is done
                # in the following loop is genes_in_splits.add call, which populates GenesInSplits class.
                for gene_callers_id in gene_calls_in_contigs_dict[contig]:
                    if self.gene_calls_dict[gene_callers_id]['stop'] > start and self.gene_calls_dict[gene_callers_id]['start'] < stop:
                        gene_start_stops.append((self.gene_calls_dict[gene_callers_id]['start'], self.gene_calls_dict[gene_callers_id]['stop']), )
                        genes_in_splits.add(split_name, start, stop, gene_callers_id, self.gene_calls_dict[gene_callers_id]['start'], self.gene_calls_dict[gene_callers_id]['stop'])

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

                splits_dict[split_name] = {'num_genes': len(gene_start_stops),
                                           'avg_gene_length': numpy.mean([(l[1] - l[0]) for l in gene_start_stops]) if len(gene_start_stops) else 0.0,
                                           'ratio_coding': total_coding_nts * 1.0 / (stop - start),
                                           }

        # open connection
        contigs_db = ContigsDatabase(self.db_path)
        # push raw entries for splits table
        db_entries = [tuple([split] + [splits_dict[split][h] for h in t.genes_in_splits_summary_table_structure[1:]]) for split in splits_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.genes_in_splits_summary_table_name, db_entries)

        # push entries for genes in splits table
        db_entries = [tuple([entry_id] + [genes_in_splits.splits_to_prots[entry_id][h] for h in t.genes_in_splits_table_structure[1:]]) for entry_id in genes_in_splits.splits_to_prots]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.genes_in_splits_table_name, db_entries)
        # disconnect
        contigs_db.disconnect()


class TablesForHMMHits(Table):
    def __init__(self, db_path, num_threads_to_use=1, run=run, progress=progress):
        self.num_threads_to_use = num_threads_to_use
        self.db_path = db_path

        self.debug = False

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        if not self.genes_are_called:
            raise ConfigError, "It seems the contigs database '%s' was created with '--skip-gene-calling' flag.\
                                Nothing to do here :/" % (self.db_path)

        self.init_gene_calls_dict()

        if not len(self.gene_calls_dict):
            raise ConfigError, "Tables that should contain gene calls are empty. Which probably means the gene\
                                caller reported no genes for your contigs."

        self.set_next_available_id(t.hmm_hits_table_name)
        self.set_next_available_id(t.hmm_hits_splits_table_name)


    def populate_search_tables(self, sources={}):
        # if we end up generating a temporary file for protein sequences:
        if not len(sources):
            import anvio.data.hmm
            sources = anvio.data.hmm.sources

        if not sources:
            return

        target_files_dict = {}

        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # here we will go through targets and populate target_files_dict based on what we find among them.
        targets = set([s['target'] for s in sources.values()])
        for target in targets:

            alphabet, context = utils.anvio_hmm_target_term_to_alphabet_and_context(target)

            self.run.info('Target found', '%s:%s' % (alphabet, context))

            class Args: pass
            args = Args()
            args.contigs_db = self.db_path
            contigs_db = ContigsSuperclass(args)

            if context == 'GENE':
                if alphabet == 'AA':
                    target_files_dict['AA:GENE'] = os.path.join(tmp_directory_path, 'aa_gene_sequences.fa')
                    self.export_sequences_table_in_db_into_FASTA_file(t.gene_protein_sequences_table_name, output_file_path=target_files_dict['AA:GENE'])
                else:
                    target_files_dict['%s:GENE' % alphabet] = os.path.join(tmp_directory_path, '%s_gene_sequences.fa' % alphabet)
                    contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['%s:GENE' % alphabet],
                                                                               simple_headers=True,
                                                                               rna_alphabet=True if alphabet=='RNA' else False)
            elif context == 'CONTIG':
                if alphabet == 'AA':
                    pass # because you can't be here.
                else:
                    target_files_dict['%s:CONTIG' % alphabet] = os.path.join(tmp_directory_path, '%s_contig_sequences.fa' % alphabet)
                    utils.export_contigs_from_contigs_db(self.db_path,
                                                         target_files_dict['%s:CONTIG' % alphabet],
                                                         rna_alphabet=True if alphabet=='RNA' else False)

        commander = HMMer(target_files_dict, num_threads_to_use=self.num_threads_to_use)

        for source in sources:
            alphabet, context = utils.anvio_hmm_target_term_to_alphabet_and_context(sources[source]['target'])

            kind_of_search = sources[source]['kind']
            domain = sources[source]['domain']
            all_genes_searched_against = sources[source]['genes']
            hmm_model = sources[source]['model']
            reference = sources[source]['ref']

            hmm_scan_hits_txt = commander.run_hmmscan(source,
                                                      '%s:%s' % (alphabet, context),
                                                      kind_of_search,
                                                      domain,
                                                      all_genes_searched_against,
                                                      hmm_model,
                                                      reference)

            if not hmm_scan_hits_txt:
                search_results_dict = {}
            else:
                parser = parser_modules['search']['hmmscan'](hmm_scan_hits_txt)
                search_results_dict = parser.get_search_results()

            self.append(source, reference, kind_of_search, domain, all_genes_searched_against, search_results_dict)

        if not self.debug:
            commander.clean_tmp_dirs()
            for v in target_files_dict.values():
                os.remove(v)


    def append(self, source, reference, kind_of_search, domain, all_genes, search_results_dict):
        # we want to define unique identifiers for each gene first. this information will be used to track genes that will
        # break into multiple pieces due to arbitrary split boundaries. while doing that, we will add the 'source' info
        # into the dictionary, so it perfectly matches to the table structure

        if not len(search_results_dict):
            return

        for entry_id in search_results_dict:
            hit = search_results_dict[entry_id]

            gene_call = self.gene_calls_dict[hit['gene_callers_id']]

            hit['gene_unique_identifier'] = hashlib.sha224('_'.join([gene_call['contig'], hit['gene_name'], str(gene_call['start']), str(gene_call['stop'])])).hexdigest()
            hit['source'] = source

        self.delete_entries_for_key('source', source, [t.hmm_hits_info_table_name, t.hmm_hits_table_name, t.hmm_hits_splits_table_name])

        contigs_db = ContigsDatabase(self.db_path)

        # push information about this search result into serach_info table.
        db_entries = [source, reference, kind_of_search, domain, ', '.join(all_genes)]
        contigs_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.hmm_hits_info_table_name, db_entries)

        # then populate serach_data table for each contig.
        db_entries = []
        for hit in search_results_dict.values():
            entry_id = self.next_id(t.hmm_hits_table_name)
            db_entries.append(tuple([entry_id] + [hit[h] for h in t.hmm_hits_table_structure[1:]]))
            # tiny hack here: for each hit, we are generating a unique id (`entry_id`), and feeding that information
            #                 back into the dictionary to pass it to processing of splits, so each split-level
            #                 entry knows who is their parent.
            hit['hmm_hit_entry_id'] = entry_id

        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.hmm_hits_table_name, db_entries)

        db_entries = self.process_splits(search_results_dict)
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.hmm_hits_splits_table_name, db_entries)

        contigs_db.disconnect()


    def process_splits(self, search_results_dict):
        hits_per_contig = {}
        for hit in search_results_dict.values():
            contig_name = self.gene_calls_dict[hit['gene_callers_id']]['contig']

            if contig_name in hits_per_contig:
                hits_per_contig[contig_name].append(hit)
            else:
                hits_per_contig[contig_name] = [hit]

        db_entries_for_splits = []

        for contig in self.contigs_info:
            if contig not in hits_per_contig:
                # no hits for this contig. pity!
                continue

            for split_name in self.contig_name_to_splits[contig]:
                split_start = self.splits_info[split_name]['start']
                split_stop = self.splits_info[split_name]['end']

                # FIXME: this really needs some explanation.
                for hit in hits_per_contig[contig]:
                    hit_start = self.gene_calls_dict[hit['gene_callers_id']]['start']
                    hit_stop = self.gene_calls_dict[hit['gene_callers_id']]['stop']

                    if hit_stop > split_start and hit_start < split_stop:
                        gene_length = hit_stop - hit_start
                        # if only a part of the gene is in the split:
                        start_in_split = (split_start if hit_start < split_start else hit_start) - split_start
                        stop_in_split = (split_stop if hit_stop > split_stop else hit_stop) - split_start
                        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length

                        db_entry = tuple([self.next_id(t.hmm_hits_splits_table_name), hit['hmm_hit_entry_id'], split_name, percentage_in_split, hit['source']])
                        db_entries_for_splits.append(db_entry)

        return db_entries_for_splits


class TablesForCollections(Table):
    """Populates the collections_* tables, where collections of bins of contigs and splits are kept"""
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run, progress)

        # set these dudes so we have access to unique IDs:
        self.set_next_available_id(t.collections_bins_info_table_name)
        self.set_next_available_id(t.collections_contigs_table_name)
        self.set_next_available_id(t.collections_splits_table_name)


    def delete(self, collection_name):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        # remove any pre-existing information for 'collection_name'
        self.delete_entries_for_key('collection_name', collection_name, [t.collections_info_table_name, t.collections_contigs_table_name, t.collections_splits_table_name, t.collections_bins_info_table_name])


    def append(self, collection_name, collection_dict, bins_info_dict=None):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        if bins_info_dict:
            if set(collection_dict.keys()) - set(bins_info_dict.keys()):
                raise ConfigError, 'Bins in the collection dict do not match to the ones in the bins info dict.\
                                    They do not have to be identical, but for each bin id, there must be a unique\
                                    entry in the bins informaiton dict. There is something wrong with your input :/'

        # remove any pre-existing information for 'collection_name'
        self.delete(collection_name)

        num_splits_in_collection_dict = sum([len(splits) for splits in collection_dict.values()])
        splits_in_collection_dict = set(list(chain.from_iterable(collection_dict.values())))
        if len(splits_in_collection_dict) != num_splits_in_collection_dict:
            raise ConfigError, "TablesForCollections::append: %d of the split or contig IDs appear more than once in\
                                your collections input. It is unclear to anvi'o how did you manage to do this, but we\
                                cannot go anywhere with this :/" % (num_splits_in_collection_dict - len(splits_in_collection_dict))

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))

        # how many clusters are defined in 'collection_dict'?
        bin_names = collection_dict.keys()

        # push information about this search result into serach_info table.
        db_entries = tuple([collection_name, num_splits_in_collection_dict, len(bin_names), ','.join(bin_names)])
        database._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_info_table_name, db_entries)

        if not bins_info_dict:
            colors = utils.get_random_colors_dict(bin_names)
            for bin_name in bin_names:
                bins_info_dict[bin_name] = {'html_color': colors[bin_name], 'source': 'UNKNOWN'}

        # populate bins info table.
        db_entries = [(self.next_id(t.collections_bins_info_table_name), collection_name, b, bins_info_dict[b]['source'], bins_info_dict[b]['html_color']) for b in bin_names]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.collections_bins_info_table_name, db_entries)

        # populate splits table
        db_entries = []
        for bin_name in collection_dict:
            for split_name in collection_dict[bin_name]:
                db_entries.append(tuple([self.next_id(t.collections_splits_table_name), collection_name, split_name, bin_name]))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_splits_table_name, db_entries)
        num_splits = len(db_entries)


        # FIXME: This function can be called to populate the contigs database (via anvi-populate-collections), or
        # the profile database. when it is contigs database, the superclass Table has the self.splits_info variable
        # set when it is initialized. however, the Table instance is missing self.splis when it is initialized with
        # the profile database. hence some special controls for contigs db (note that collections_contigs_table is
        # only populated in the contigs database):
        if self.db_type == 'contigs':
            splits_only_in_collection_dict = [c for c in splits_in_collection_dict if c not in self.splits_info]
            splits_only_in_db = [c for c in self.splits_info if c not in splits_in_collection_dict]

            if len(splits_only_in_collection_dict):
                self.run.warning('%d of %d splits found in "%s" results are not in the database. This may be OK,\
                                          but you must be the judge of it. If this is somewhat surprising, please use caution\
                                          and make sure all is fine before going forward with you analysis.'\
                                                % (len(splits_only_in_collection_dict), len(splits_in_collection_dict), collection_name))

            if len(splits_only_in_db):
                self.run.warning('%d of %d splits found in the database were missing from the "%s" results. If this\
                                          does not make any sense, please make sure you know why before going any further.'\
                                                % (len(splits_only_in_db), len(self.splits_info), collection_name))

            # then populate contigs table.
            db_entries = self.process_contigs(collection_name, collection_dict)
            database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_contigs_table_name, db_entries)

        database.disconnect()

        self.run.info('Collections', '%s annotations for %d splits have been successfully added to the database at "%s".'\
                                        % (collection_name, num_splits, self.db_path), mc='green')


    def process_contigs(self, collection_name, collection_dict):
        db_entries_for_contigs = []

        split_to_bin_name = {}
        for bin_name in collection_dict:
            for split_name in collection_dict[bin_name]:
                split_to_bin_name[split_name] = bin_name

        contigs_processed = set([])
        for split_name in split_to_bin_name:
            if split_name not in self.splits_info:
                # which means this split only appears in the input file, but not in the database.
                continue

            contig_name = self.splits_info[split_name]['parent']

            if contig_name in contigs_processed:
                continue
            else:
                contigs_processed.add(contig_name)

            db_entry = tuple([self.next_id(t.collections_contigs_table_name), collection_name, contig_name, split_to_bin_name[split_name]])
            db_entries_for_contigs.append(db_entry)

        return db_entries_for_contigs


class TablesForStates(Table):
    def __init__(self, db_path):
        self.db_path = db_path
        self.states = {}

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run, progress)

        self.init()


    def init(self):
        anvio_db = DBClassFactory().get_db_object(self.db_path)
        self.states = anvio_db.db.get_table_as_dict(t.states_table_name)
        anvio_db.disconnect()


    def get_state(self, state_id):
        if state_id not in self.states:
            return None

        return self.states[state_id]


    def store_state(self, state_id, content, last_modified=None):
        self.remove_state(state_id)

        last_modified = datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S") if not last_modified else last_modified

        anvio_db = DBClassFactory().get_db_object(self.db_path)
        anvio_db.db._exec('''INSERT INTO %s VALUES (?,?,?)''' % t.states_table_name, (state_id, content, last_modified))
        self.states = anvio_db.db.get_table_as_dict(t.states_table_name)

        anvio_db.disconnect()


    def remove_state(self, state_id):
        self.delete_entries_for_key('name', state_id, [t.states_table_name])


class TablesForTaxonomy(Table):
    """Populate all tables with taxonomy information.

       Essentially takes in a dictionary of genes and taxon calls, populates three tables: taxon_names,
       splits_taxonomy, and gene_taxonomy when 'create' is called."""
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path
        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress)

        # this class keeps track of genes that occur in splits, and responsible
        # for generating the necessary table in the contigs database
        self.genes_in_splits = GenesInSplits()


    def create(self, genes_taxonomy_dict, taxon_names_dict, source='unkown source'):
        self.source = source

        if not self.genes_are_called:
            raise ConfigError, "Something is wrong. The contigs database says that genes were now called, and here\
                                you are trying to populate taxonomy tables for genes. No, thanks."

        self.init_gene_calls_dict()

        self.genes_taxonomy_dict = genes_taxonomy_dict
        self.taxon_names_dict = taxon_names_dict

        self.sanity_check() 

        # oepn connection
        contigs_db = ContigsDatabase(self.db_path)

        self.splits_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)

        # test whether there are already genes tables populated
        taxonomy_source = contigs_db.meta['taxonomy_source']
        if taxonomy_source:
            self.run.warning('Previous taxonomy information from "%s" is being replaced with the incoming data\
                              through "%s".' % (taxonomy_source, self.source))
            contigs_db.db._exec('''DELETE FROM %s''' % (t.splits_taxonomy_table_name))
            contigs_db.db._exec('''DELETE FROM %s''' % (t.taxon_names_table_name))
            contigs_db.db._exec('''DELETE FROM %s''' % (t.genes_taxonomy_table_name))

        # populate taxon mames table
        self.populate_taxon_names_table()

        # populate genes taxonomy table
        self.populate_genes_taxonomy_table()

        # compute and push split taxonomy information.
        self.populate_splits_taxonomy_table()

        # set the source
        contigs_db.db.remove_meta_key_value_pair('taxonomy_source')
        contigs_db.db.set_meta_value('taxonomy_source', self.source)

        # disconnect like a pro.
        contigs_db.disconnect()


    def populate_genes_taxonomy_table(self):
        # open connection
        contigs_db = ContigsDatabase(self.db_path)

        # push taxonomy data
        db_entries = [(gene_callers_id, self.genes_taxonomy_dict[gene_callers_id]) for gene_callers_id in self.genes_taxonomy_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.genes_taxonomy_table_name, db_entries)

        # disconnect
        contigs_db.disconnect()

        self.run.info('Genes taxonomy table', 'Taxonomy stored for %d gene calls' % len(db_entries))


    def populate_taxon_names_table(self):
        # open connection
        contigs_db = ContigsDatabase(self.db_path)

        db_entries = [tuple([t_name_id] + [self.taxon_names_dict[t_name_id][t_level] for t_level in t.taxon_names_table_structure[1:]]) for t_name_id in self.taxon_names_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.taxon_names_table_name, db_entries)

        contigs_db.disconnect()
        self.run.info('Taxon names table', 'Updated with %d unique taxon names' % len(db_entries))


    def sanity_check(self):
        """Basic checks to make sure things are in order at least to a minimum extent."""

        self.progress.new('Sanity checking ...')
        self.progress.update('Comparing gene caller ids in the incoming data and in the contigs database ..')
        gene_caller_ids_in_database = set(self.gene_calls_dict.keys())
        gene_caller_ids_in_taxonomy_dict = set(self.genes_taxonomy_dict.keys())
        gene_caller_ids_missing_in_db = gene_caller_ids_in_taxonomy_dict.difference(gene_caller_ids_in_database)
        self.progress.end()

        run.info("Num gene caller ids in the db", pp(len(gene_caller_ids_in_database)))
        run.info("Num gene caller ids in the incoming data", pp(len(gene_caller_ids_in_taxonomy_dict)))

        if gene_caller_ids_missing_in_db:
            raise ConfigError, "Taxonomy information for genes you are trying to import into the database contains\
                                %s gene caller ids that do not appear to be in the database. This is a step you must\
                                be very careful to make sure you are not importing annotations for genes that have\
                                nothing to do with your contigs database. To make sure of that, you should always work\
                                with `anvi-get-dna-sequences-for-gene-calls` or `anvi-get-aa-sequences-for-gene-calls` programs\
                                to get the data to annotate. For instance one of the gene caller ids you have in your\
                                input data that does not appear in the database is this one: '%s'. Anvi'o hopes it makes\
                                sense to you, because it definitely does not make any sense to anvi'o :("\
                                                        % (len(gene_caller_ids_missing_in_db), str(gene_caller_ids_missing_in_db.pop()))

        # check whether input matrix dict
        keys_found =  self.taxon_names_dict.values()[0].keys()
        missing_keys = [key for key in t.taxon_names_table_structure[1:] if key not in keys_found]
        if len(missing_keys):
            raise ConfigError, "Anvi'o is trying to get ready to create tables for taxonomy, but there is something\
                                wrong :( The taxonomy names dict (one of the required input dictionaries to the class\
                                seems to be missing a one or more keys that are necessary to finish the job. Here is \
                                a list of missing keys: %s. The complete list of input keys should contain these: %s."\
                                        % (', '.join(missing_keys), ', '.join(t.taxon_names_table_structure[1:]))

        if not len(self.taxon_names_dict):
            raise ConfigError, "Anvi'o is trying to get ready to create tables for taxonomy, but taxonomy names dict\
                                (one of the required input dictionaries to the class responsible for this task) seems\
                                to be empty."


    def populate_splits_taxonomy_table(self):
        """Populate the taxonomy information per split"""

        # build a dictionary for fast access to all proteins identified within a contig
        gene_caller_ids_in_contigs = {}
        for gene_callers_id in self.genes_taxonomy_dict:
            contig = self.gene_calls_dict[gene_callers_id]['contig']
            if contig in gene_caller_ids_in_contigs:
                gene_caller_ids_in_contigs[contig].add(gene_callers_id)
            else:
                gene_caller_ids_in_contigs[contig] = set([gene_callers_id])

        contigs_without_annotation = list(set(self.contigs_info.keys()) - set(gene_caller_ids_in_contigs.keys()))

        for contig in contigs_without_annotation:
            gene_caller_ids_in_contigs[contig] = set([])

        splits_dict = {}

        num_splits_processed = 0
        num_splits_with_taxonomy = 0

        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                num_splits_processed += 1

                splits_dict[split_name] = None
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                taxon_name_ids = []
                for gene_callers_id in gene_caller_ids_in_contigs[contig]:
                    if self.gene_calls_dict[gene_callers_id]['stop'] > start and self.gene_calls_dict[gene_callers_id]['start'] < stop:
                        taxon_name_ids.append(self.genes_taxonomy_dict[gene_callers_id])

                if not taxon_name_ids:
                    continue

                if len(set(taxon_name_ids)) == 1:
                    splits_dict[split_name] = taxon_name_ids[0]
                else:
                    d = Counter()
                    for taxon_name_id in taxon_name_ids:
                        d[taxon_name_id] += 1

                    most_frequent_taxon_name_id, occurrence = d.most_common()[0]
                    splits_dict[split_name] = most_frequent_taxon_name_id

                num_splits_with_taxonomy += 1

        # open connection
        contigs_db = ContigsDatabase(self.db_path)

        # push taxonomy data
        db_entries = [(split, splits_dict[split]) for split in splits_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.splits_taxonomy_table_name, db_entries)

        # disconnect
        contigs_db.disconnect()

        self.run.info('Splits taxonomy', 'Input data from "%s" annotated %d of %d splits (%.1f%%) with taxonomy.'\
                                            % (self.source, num_splits_with_taxonomy, num_splits_processed,
                                               num_splits_with_taxonomy * 100.0 / num_splits_processed))


class TableForProteinClusters(Table):
    """A class to populte  protein clusters tabel in a given pan db.

      Here is an example:

        >>> table = TableForProteinClusters(db_path)
        >>> for ...:
        >>>     table.add({'gene_caller_id': gene_caller_id, 'protein_cluster_id': protein_cluster_id, 'genome_name': genome_name})
        >>> table.store()
    """

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        is_pan_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pan__version__, run, progress)

        self.set_next_available_id(t.pan_protein_clusters_table_name)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in t.pan_protein_clusters_table_structure[1:]])


    def store(self):
        self.delete_contents_of_table(t.pan_protein_clusters_table_name, warning=False)

        db_entries = [tuple([self.next_id(t.pan_protein_clusters_table_name)] + entry) for entry in self.entries]
        pan_db = PanDatabase(self.db_path)
        pan_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.pan_protein_clusters_table_name, db_entries)
        pan_db.disconnect()


class TableForGeneFunctions(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        self.set_next_available_id(t.gene_function_calls_table_name)


    def create(self, functions_dict, drop_previous_annotations_first = False):
        self.sanity_check()

        # incoming stuff:
        gene_function_sources = set([v['source'] for v in functions_dict.values()])
        unique_num_genes = len(set([v['gene_callers_id'] for v in functions_dict.values()]))

        # oepn connection
        contigs_db = ContigsDatabase(self.db_path)

        # are there any previous annotations in the db:
        gene_function_sources_in_db = set(contigs_db.meta['gene_function_sources'] or [])

        # difference between sources in the db, and incoming sources:
        gene_function_sources_both_in_db_and_incoming_dict = gene_function_sources.intersection(gene_function_sources_in_db)

        # here we will do some magic. there are mulitple scenarios to consider here based on whether there
        # are functions already in the database, whether some of them matches to the incoming functions, etc.
        # let's go case-by-case:
        if not gene_function_sources_in_db:
            # set the sources and continue
            contigs_db.db.remove_meta_key_value_pair('gene_function_sources')
            contigs_db.db.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources)))

        elif gene_function_sources_in_db and drop_previous_annotations_first:
            # there are gene calls, but the user wants everything to be dropeped.
            self.run.warning("As per your request, anvi'o is DROPPING all previous function calls from %d sources\
                              before adding the incoming data, which contains %d entries originating from %d sources: %s" \
                                    % (len(gene_function_sources_in_db), len(functions_dict),
                                       len(gene_function_sources), ', '.join(gene_function_sources)))

            # clean the table and reset the next available ids
            contigs_db.db._exec('''DELETE FROM %s''' % (t.gene_function_calls_table_name))
            self.reset_next_available_id_for_table(t.gene_function_calls_table_name)

            # set the sources
            contigs_db.db.remove_meta_key_value_pair('gene_function_sources')
            contigs_db.db.set_meta_value('gene_function_sources', ','.join(gene_function_sources))

        elif gene_function_sources_in_db and gene_function_sources_both_in_db_and_incoming_dict:
            # some of the functions in the incoming dict match to what is already in the db. remove
            self.run.warning("Some of the annotaiton sources you want to add into the database are already in the db. So\
                              anvi'o will REPLACE those with the incoming data from these sources: %s" % \
                                            ', '.join(gene_function_sources_both_in_db_and_incoming_dict))

            # remove those entries for matching sources:
            for source in gene_function_sources_both_in_db_and_incoming_dict:
                contigs_db.db._exec('''DELETE FROM %s WHERE source = "%s"''' % (t.gene_function_calls_table_name, source))

            # set the sources
            contigs_db.db.remove_meta_key_value_pair('gene_function_sources')
            contigs_db.db.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources_in_db.union(gene_function_sources))))

        else:
            # fuctions in the db, but none of them match with the incoming annotation sources. totally new stuff.
            # good then. update sources
            contigs_db.db.remove_meta_key_value_pair('gene_function_sources')
            contigs_db.db.set_meta_value('gene_function_sources', ','.join(list(gene_function_sources_in_db.union(gene_function_sources))))

        # push the data
        db_entries = [tuple([self.next_id(t.gene_function_calls_table_name)] + [functions_dict[v][h] for h in t.gene_function_calls_table_structure[1:]]) for v in functions_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?)''' % t.gene_function_calls_table_name, db_entries)

        # disconnect like a pro.
        contigs_db.disconnect()

        self.run.info('Gene functions', '%d function calls from %d sources for %d unique gene calls has\
                                        been added to the contigs database.' % \
                                            (len(functions_dict), len(gene_function_sources), unique_num_genes))


    def sanity_check(self):
        pass


class GenesInSplits:
    def __init__(self):
        self.entry_id = 0
        self.splits_to_prots = {}

    def add(self, split_name, split_start, split_end, gene_callers_id, prot_start, prot_end):

        gene_length = prot_end - prot_start

        if gene_length <= 0:
            raise ConfigError, "dbops.py/GeneInSplits: OK. There is something wrong. We have this gene, '%s',\
                                which starts at position %d and ends at position %d. Well, it doesn't look right,\
                                does it?" % (gene_callers_id, prot_start, prot_end)

        # if only a part of the gene is in the split:
        start_in_split = (split_start if prot_start < split_start else prot_start) - split_start
        stop_in_split = (split_end if prot_end > split_end else prot_end) - split_start
        percentage_in_split = (stop_in_split - start_in_split) * 100.0 / gene_length

        self.splits_to_prots[self.entry_id] = {'split': split_name,
                                               'gene_callers_id': gene_callers_id,
                                               'start_in_split': start_in_split,
                                               'stop_in_split': stop_in_split,
                                               'percentage_in_split': percentage_in_split}
        self.entry_id += 1


####################################################################################################
#
#     HELPER FUNCTIONS
#
####################################################################################################

def is_contigs_db(db_path):
    filesnpaths.is_file_exists(db_path)
    if get_db_type(db_path) != 'contigs':
        raise ConfigError, "'%s' is not an anvi'o contigs database." % db_path


def is_pan_or_profile_db(db_path):
    if get_db_type(db_path) not in ['pan', 'profile']:
        raise ConfigError, "'%s' is neither a pan nor a profile database :/ Someone is in trouble."


def is_profile_db(db_path):
    if get_db_type(db_path) != 'profile':
        raise ConfigError, "'%s' is not an anvi'o profile database." % db_path


def is_pan_db(db_path):
    if get_db_type(db_path) != 'pan':
        raise ConfigError, "'%s' is not an anvi'o pan database." % db_path


def is_samples_db(db_path):
    if get_db_type(db_path) != 'samples_information':
        raise ConfigError, "'%s' is not an anvi'o samples database." % db_path


def is_db_ok_to_create(db_path, db_type):
    if os.path.exists(db_path):
        raise ConfigError, "Anvi'o will not overwrite an existing %s database. Please choose a different name\
                            or remove the existing database ('%s') first." % (db_type, db_path)

    if not db_path.lower().endswith('.db'):
        raise ConfigError, "Please make sure the file name for your new %s db has a '.db' extension. Anvi'o developers\
                            apologize for imposing their views on how anvi'o databases should be named, and are\
                            humbled by your cooperation." % db_type


def get_required_version_for_db(db_path):
    db_type = get_db_type(db_path)

    if db_type not in t.versions_for_db_types:
        raise ConfigError, "Anvi'o was trying to get the version of the -alleged- anvi'o database '%s', but it failed\
                            because it turns out it doesn't know anything about this '%s' type." % (db_path, db_type)

    return t.versions_for_db_types[db_type]


def get_db_type(db_path):
    filesnpaths.is_file_exists(db_path)

    try:
        database = db.DB(db_path, None, ignore_version=True)
    except:
        raise ConfigError, 'Are you sure "%s" is a database file? Because, you know, probably\
                            it is not at all..' % db_path

    tables = database.get_table_names()
    if 'self' not in tables:
        database.disconnect()
        raise ConfigError, "'%s' does not seem to be a anvi'o database..." % db_path

    db_type = database.get_meta_value('db_type')
    database.disconnect()

    return db_type


def is_profile_db_and_contigs_db_compatible(profile_db_path, contigs_db_path):
    is_contigs_db(contigs_db_path)
    is_profile_db(profile_db_path)

    contigs_db = ContigsDatabase(contigs_db_path)
    profile_db = ProfileDatabase(profile_db_path)

    a_hash = contigs_db.meta['contigs_db_hash']
    p_hash = profile_db.meta['contigs_db_hash']
    merged = profile_db.meta['merged']

    contigs_db.disconnect()
    profile_db.disconnect()

    if a_hash != p_hash:
        raise ConfigError, 'The contigs database and the profile database does not\
                            seem to be compatible. More specifically, this contigs\
                            database is not the one that was used when %s generated\
                            this profile database.'\
                                % 'anvi-merge' if merged else 'anvi-profile'

    return True


def is_profile_db_and_samples_db_compatible(profile_db_path, samples_db_path, manual_mode_exception=False):
    """Check whether every sample name in the profile database is represented in the samples information database"""
    is_profile_db(profile_db_path)
    is_samples_db(samples_db_path)

    profile_db = ProfileDatabase(profile_db_path)
    samples_db = SamplesInformationDatabase(samples_db_path)

    if 'merged' in profile_db.meta and not int(profile_db.meta['merged']):
        raise ConfigError, "Samples databases are only useful if you are working on a merged profile."

    if manual_mode_exception:
        # manual mode exception is a funny need. when the user wants to use --manual flag with anvi-interactive,
        # and provides a blank name for profile.db, sample names for that profile db are automaticalaly populated
        # from the data matrix file. so far so good. if the data matrix file is a mixed one, i.e., contains both
        # sample names in the conventional sense with numerical view data, and additional data, then sample names
        # coming from the first row of the file does not match with the entries in samples database. the best
        # solution for this is to enforce the use of view_data and additional_data inputs separately. unfortunately
        # that complicates the anvi-server interface unnecessarily.. so here we will simply pass the following
        # important tests, and hope that the user did a good job making sure sample names in samples database do match
        # the samples bit that appears in the data file :( I know, I know...
        samples_in_samples_db_but_not_in_profile_db = samples_db.samples - profile_db.samples
        if len(samples_in_samples_db_but_not_in_profile_db):
            raise ConfigError, "Anvi'o is upset with you :/ Please make sure your samples information files (or your\
                                samples database) contain sample names from your data file. These sample names are in\
                                your samples information, but not in your data file: '%s'. If this error does not make\
                                any sense to you, please contact an anvi'o developer." % ', '.join(samples_in_samples_db_but_not_in_profile_db)
        return


    missing_samples = profile_db.samples - samples_db.samples
    num_represented_samples = len(profile_db.samples) - len(missing_samples)

    if len(missing_samples):
        how_much_of_the_samples_are_represented_txt = 'none' if len(missing_samples) == len(profile_db.samples) else\
                                                      'only %d of %d' % (num_represented_samples, len(profile_db.samples))

        raise ConfigError, "The samples information database you provided ('%s') does not seem to agree well with the profile\
                            database ('%s'). More specifically, %s of the samples in the profile database are repesented in\
                            the samples information database. Names for these missing samples go like this: %s ...,\
                            while the sample names in the samples information database go like this: %s ... This could be due to\
                            a simple typo, or you may be using the wrong or outdated samples information database. You may need to\
                            regenerate the samples information database to fix this problem :/"\
                                                % (samples_db_path, profile_db_path, how_much_of_the_samples_are_represented_txt,
                                                   ', '.join(list(missing_samples)[0:3]), ', '.join(list(samples_db.samples)[0:3]))

    if samples_db.sample_names_for_order:
        missing_samples = samples_db.sample_names_for_order - profile_db.samples

        if len(missing_samples):
            raise ConfigError, "The samples order information in the samples database do not match with the sample names in\
                                the profile database (or the input data). To be precise, %d sample(s) occur(s) only in the\
                                samples database, and not found in the profile database (or in the input data). Here is some of\
                                them: %s ..." % (len(missing_samples), ', '.join(list(missing_samples)[0:3]))


def get_split_names_in_profile_db(profile_db_path):
    is_profile_db(profile_db_path)

    profile_db = ProfileDatabase(profile_db_path)

    if int(profile_db.meta['merged']):
        split_names = set(profile_db.db.get_single_column_from_table('mean_coverage_Q2Q3_splits', 'contig'))
    else:
        split_names = set(profile_db.db.get_single_column_from_table('atomic_data_splits', 'contig'))

    profile_db.disconnect()

    return split_names


def add_hierarchical_clustering_to_db(anvio_db_path, clustering_name, clustering_newick, distance, linkage, make_default=False, run=run):
    """Adds a new clustering into an anvi'o db"""

    # let's learn who we are dealing with:
    db_type = get_db_type(anvio_db_path)

    utils.is_this_name_OK_for_database('clustering_name parameter', clustering_name, stringent=False)

    # replace clustering id with a text that contains distance and linkage information
    clustering_id = ':'.join([clustering_name, distance, linkage])

    anvio_db = DBClassFactory().get_db_object(anvio_db_path)

    if t.clusterings_table_name not in anvio_db.db.get_table_names():
        raise ConfigError, "You can't a new clustering result into this %s database (%s). You know why? Becasue it doesn't\
                            have a table for 'clusterings' :(" % (db_type, anvio_db_path)

    try:
        available_clusterings = anvio_db.db.get_meta_value('available_clusterings').split(',')
    except:
        available_clusterings = []

    if clustering_id in available_clusterings:
        run.warning('Clustering for "%s" (with %s distance and %s linkage) is already in the database. Its content will\
                     be replaced with the new one.' % (clustering_name, distance, linkage))

        anvio_db.db._exec('''DELETE FROM %s where clustering = "%s"''' % (t.clusterings_table_name, clustering_id))
    else:
        available_clusterings.append(clustering_id)

    anvio_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % t.clusterings_table_name, tuple([clustering_id, clustering_newick]))

    try:
        anvio_db.db.remove_meta_key_value_pair('available_clusterings')
    except:
        pass
    anvio_db.db.set_meta_value('available_clusterings', ','.join(available_clusterings))

    try:
        anvio_db.db.remove_meta_key_value_pair('PCs_clustered' if db_type == 'pan' else 'contigs_clustered')
    except:
        pass
    anvio_db.db.set_meta_value('PCs_clustered' if db_type == 'pan' else 'contigs_clustered', True)

    try:
        anvio_db.db.get_meta_value('default_clustering')
        default_clustering_is_set = True
    except:
        default_clustering_is_set = False

    if make_default or not default_clustering_is_set:
        try:
            anvio_db.db.remove_meta_key_value_pair('default_clustering')
        except:
            pass
        anvio_db.db.set_meta_value('default_clustering', clustering_id)

    anvio_db.disconnect()

    run.info('New hierarchical clusetring', '"%s" has been added to the database...' % clustering_id)


def get_default_clustering_id(default_clustering_requested, clusterings_dict, progress=progress, run=run):
    """Get the proper default clustering given the desired default with respect to available clusterings.
    
       This is tricky. We have some deault clusterings defined in the constants. For instance, for the
       merged profiles we want the default to be 'tnf-cov', for single profiles we want it to be 'tnf',
       etc. The problem is that these defaults do not indicate any distance metric or linkages,
       even though anvi'o allows users to define those variables freely in cluster configurations.

       A clustering dict can contain multiple clustrings. The purpose of this function is to take the
       desired default into consideration, but then find a working one if it is not available, or there
       are multiple ones in the dict.
    """

    if not clusterings_dict:
        raise ConfigError, "You requested to get the default clustering given the clustering dictionary,\
                            but the clustering dict is empty :/ "

    matching_clustering_ids = [clustering for clustering in clusterings_dict if clustering.lower().startswith(default_clustering_requested.lower())]

    if not len(matching_clustering_ids):
        default_clustering = clusterings_dict.keys()[0]
        run.warning('`get_default_clustering_id` function is concerned, because nothing in the clusterings\
                     dict matched to the desired default clustring class "%s". So it literally set "%s"\
                     (a class of "%s") randomly as the default. Good luck :/' % (default_clustering_requested,
                                                                                 default_clustering,
                                                                                 default_clustering.split(':')[0]))
        return default_clustering
    elif len(matching_clustering_ids) == 1:
        return matching_clustering_ids[0]
    else:
        default_clustering = matching_clustering_ids[0]
        run.warning('`get_default_clustering_id` function is concerned, because there were multiple entries\
                     in the clusterings dict matched to the desired default clustring class "%s". So it set\
                     the first of all %d matching clusterings, which happened to be the "%s", as the\
                     default. We hope that will not screw up your mojo :/' % (default_clustering_requested,
                                                                              len(matching_clustering_ids),
                                                                              default_clustering))
        return default_clustering


def export_aa_sequences_from_contigs_db(contigs_db_path, output_file_path):
    filesnpaths.is_file_exists(contigs_db_path)
    filesnpaths.is_output_file_writable(output_file_path)

    class T(Table):
        def __init__(self, db_path, version, run=run, progress=progress):
            Table.__init__(self, db_path, version, run, progress)

    h = T(contigs_db_path, anvio.__contigs__version__)
    h.export_sequences_table_in_db_into_FASTA_file(t.gene_protein_sequences_table_name, output_file_path = output_file_path)

    return output_file_path


def get_all_item_names_from_the_database(db_path):
    """Return all split names or PC names in a given database"""

    all_items = set([])

    database = db.DB(db_path, get_required_version_for_db(db_path))
    db_type = database.get_meta_value('db_type')

    class Args: pass
    args = Args()

    if db_type == 'profile':
        args.profile_db = db_path
        all_items = set(ProfileSuperclass(args).split_names)
    elif db_type == 'pan':
        args.pan_db = db_path
        all_items = set(PanSuperclass(args).protein_cluster_names)
    elif db_type == 'contigs':
        args.contigs_db = db_path
        all_items = set(ContigsSuperclass(args).splits_basic_info.keys())
    else:
        raise ConfigError, "You wanted to get all items in the database %s, but no one here knows aobut its type. Seriously,\
                            what is '%s'?" % (db_path, db_type)

    if not len(all_items):
        raise ConfigError, "dbops::get_all_item_names_from_the_database speaking. Something that should never happen happened :/\
                            There seems to be nothing in this %s database. Anvi'o is as confused as you are. Please get in touch\
                            with a developer. They will love this story."

    return all_items
