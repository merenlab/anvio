# -*- coding: utf-8
"""
    Classes to create, access, and/or populate contigs and profile databases.
"""

import os
import time
import copy
import numpy
import random
import hashlib
import datetime
import operator
import textwrap
from itertools import chain
from collections import Counter

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
from anvio.hmmops import HMMSearch
from anvio.parsers import parser_modules
from anvio.tableops import Table


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


run = terminal.Run()
progress = terminal.Progress()


class ContigsSuperclass(object):
    def __init__(self, args, r = run, p = progress):
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
        self.split_to_genes_in_splits_ids = {} # for fast access to all self.genes_in_splits entries for a given split
        self.gene_callers_id_to_split_name_dict = {} # for fast access to a split name that contains a given gene callers id

        # the following two dicts require initialization:
        self.nt_positions_info_dicts_initiated = False
        self.nt_positions_in_complete_genes = {} # each nt position that falls into a complete gene call
        self.nt_positions_in_partial_genes = {} # nt positions in contigs that are in impartial gene calls
        self.nt_positions_first_base_in_codons = {}
        self.nt_positions_second_base_in_codons = {}
        self.nt_positions_third_base_in_codons = {}

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
        contigs_db = ContigsDatabase(self.contigs_db_path)

        self.progress.update('Setting contigs self data dict')
        self.a_meta = contigs_db.meta

        self.a_meta['creation_date'] = utils.get_time_to_date(self.a_meta['creation_date']) if self.a_meta.has_key('creation_date') else 'unknown'
        for key in ['split_length', 'kmer_size', 'total_length', 'num_splits', 'num_contigs', 'genes_are_called']:
            self.a_meta[key] = int(self.a_meta[key])

        self.progress.update('Reading contigs basic info')
        self.contigs_basic_info = contigs_db.db.get_table_as_dict(t.contigs_info_table_name, string_the_key = True)

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

        self.progress.update('Reading splits taxonomy')
        self.splits_taxonomy_dict = contigs_db.db.get_table_as_dict(t.splits_taxonomy_table_name)

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
            if split_name in self.split_to_genes_in_splits_ids:
                self.split_to_genes_in_splits_ids[split_name].add(entry_id)
            else:
                self.split_to_genes_in_splits_ids[split_name] = set([entry_id])

        for split_name in self.splits_basic_info:
            if not self.split_to_genes_in_splits_ids.has_key(split_name):
                self.split_to_genes_in_splits_ids[split_name] = set([])

        self.progress.update('Generating "gene caller id" to "split name" mapping dict')
        for entry in self.genes_in_splits.values():
            self.gene_callers_id_to_split_name_dict[entry['gene_callers_id']] = entry['split']

        self.progress.end()

        contigs_db.disconnect()
        self.run.info('Contigs DB', 'Initialized: %s (v. %s)' % (self.contigs_db_path, anvio.__contigs__version__))


    def init_contig_sequences(self, min_contig_length = 0):
        self.progress.new('Loading contig sequences')

        self.progress.update('Identifying contigs shorter than M')
        contigs_shorter_than_M = set([c for c in self.contigs_basic_info if self.contigs_basic_info[c]['length'] < min_contig_length])

        self.progress.update('Reading contig sequences')
        contigs_db = ContigsDatabase(self.contigs_db_path)
        self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name, string_the_key = True)
        contigs_db.disconnect()

        self.progress.update('Filtering out shorter contigs')
        for contig_name in contigs_shorter_than_M:
            self.contig_sequences.pop(contig_name)

        self.progress.end()

        return contigs_shorter_than_M


    def init_split_sequences(self, min_contig_length = 0):
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


    def init_non_singlecopy_gene_hmm_sources(self, split_names_of_interest = None, return_each_gene_as_a_layer = False):
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
        if not self.nt_positions_info_dicts_initiated:
            self.init_nt_position_info_dicts()

        in_partial_gene_call = 1 if pos_in_contig in self.nt_positions_in_partial_genes[contig_name] else 0
        in_complete_gene_call = 1 if not in_partial_gene_call and pos_in_contig in self.nt_positions_in_complete_genes[contig_name] else 0

        if in_complete_gene_call:
            if pos_in_contig in self.nt_positions_first_base_in_codons[contig_name]:
                pos_in_codon = 1
            elif pos_in_contig in self.nt_positions_second_base_in_codons[contig_name]:
                pos_in_codon = 2
            else:
                pos_in_codon = 3
        else:
            pos_in_codon = 0

        return (in_partial_gene_call, in_complete_gene_call, pos_in_codon)


    def init_nt_position_info_dicts(self):
        """This function initializes following three variables:

                - self.nt_positions_in_partial_genes
                - self.nt_positions_in_complete_genes
                - self.nt_positions_first_base_in_codons
                - self.nt_positions_second_base_in_codons
                - self.nt_positions_third_base_in_codons

        """

        self.progress.new('Initializing nt positional info (whether a nt position is within a gene and which base is it in a codon)')
        self.progress.update('...')

        for contig_name in self.contigs_basic_info.keys():
            self.nt_positions_in_partial_genes[contig_name] = set([])
            self.nt_positions_in_complete_genes[contig_name] = set([])

            self.nt_positions_first_base_in_codons[contig_name] = set([])
            self.nt_positions_second_base_in_codons[contig_name] = set([])
            self.nt_positions_third_base_in_codons[contig_name] = set([])

            for gene_unique_id, start, stop in self.contig_name_to_genes[contig_name]:
                gene_call = self.genes_in_contigs_dict[gene_unique_id]

                RANGE = range(start, stop)

                # if the gene call is a partial one, meaning the gene was cut at the beginning or
                # at the end of the contig, we are not going to try to make sense of synonymous /
                # non-synonmous bases in that. the clients who wish to use these variables must also
                # be careful about the difference
                if gene_call['partial']:
                    self.nt_positions_in_partial_genes[contig_name].update(RANGE)
                    continue

                self.nt_positions_in_complete_genes[contig_name].update(RANGE)

                if gene_call['direction'] == 'f':
                    self.nt_positions_first_base_in_codons[contig_name].update(range(start + 0, stop, 3))
                    self.nt_positions_second_base_in_codons[contig_name].update(range(start + 1, stop, 3))
                    self.nt_positions_third_base_in_codons[contig_name].update(range(start + 2, stop, 3))
                elif gene_call['direction'] == 'r':
                    self.nt_positions_first_base_in_codons[contig_name].update(range(stop - 1, start - 1, -3))
                    self.nt_positions_second_base_in_codons[contig_name].update(range(stop - 2, start - 1, -3))
                    self.nt_positions_third_base_in_codons[contig_name].update(range(stop - 3, start - 1, -3))

        self.progress.end()

        self.nt_positions_info_dicts_initiated = True


    def init_functions(self):
        if not self.contigs_db_path:
            return

        self.progress.new('Initializing functions class')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        self.gene_function_call_sources = contigs_db.meta['gene_function_sources'].split(',') if contigs_db.meta['gene_function_sources'] else []
        for hit in contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name).values():
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

            entry = ('%s :: %s' % (function if function else 'unknown', accession), e_value)
            self.gene_function_calls_dict[gene_callers_id][source] = entry

        contigs_db.disconnect()

        self.progress.end()

        self.gene_function_calls_initiated = True


    def search_splits_for_gene_functions(self, search_terms, verbose = False, full_report = False):
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
                    self.run.warning('\n'.join(matching_function_calls[search_term][0:25]), header="Matching names for '%s' (up to 25)" % search_term, raw = True, lc = 'cyan')
            else:
                self.run.info('Matches', 'No matches found the search term "%s"' % (search_term), mc = 'red')

        contigs_db.disconnect()
        self.progress.end()

        return split_names, full_report


    def get_sequences_for_gene_callers_ids(self, gene_caller_ids_list, reverse_complement_if_necessary = True):
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


    def gen_FASTA_file_of_sequences_for_gene_caller_ids(self, gene_caller_ids_list, output_file_path, wrap = 120):
        filesnpaths.is_output_file_writable(output_file_path)

        if type(wrap) != int:
            raise ConfigError, '"wrap" has to be an integer instance'
        if wrap == 0:
            wrap = None
        if wrap and wrap <= 20:
            raise ConfigError, 'Value for wrap must be larger than 20. Yes. Rules.'

        gene_caller_ids_list, sequences_dict = self.get_sequences_for_gene_callers_ids(gene_caller_ids_list)

        output = open(output_file_path, 'w')

        self.progress.new('Storing sequences')
        self.progress.update('...')
        for gene_callers_id in gene_caller_ids_list:
            entry = sequences_dict[gene_callers_id]
            header = '%d|' % (gene_callers_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['contig', 'start', 'stop', 'direction', 'rev_compd', 'length']])
            sequence = entry['sequence']

            if wrap:
                sequence = textwrap.fill(sequence, wrap, break_on_hyphens = False)

            output.write('>%s\n' % header)
            output.write('%s\n' % sequence)

        output.close()

        self.progress.end()
        self.run.info('Output', output_file_path)


class ProfileSuperclass(object):
    def __init__(self, args, r = run, p = progress):
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

        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if self.p_meta.has_key('creation_date') else 'unknown'
        self.p_meta['samples'] = sorted([s.strip() for s in self.p_meta['samples'].split(',')])
        self.p_meta['num_samples'] = len(self.p_meta['samples'])

        for key in ['merged', 'contigs_clustered', 'min_contig_length', 'total_length', 'num_splits', 'num_contigs']:
            self.p_meta[key] = int(self.p_meta[key])

        if self.p_meta['contigs_clustered']:
            self.p_meta['available_clusterings'] = sorted([s.strip() for s in self.p_meta['available_clusterings'].split(',')])
            self.clusterings = profile_db.db.get_table_as_dict(t.clusterings_table_name)
        else:
            self.p_meta['available_clusterings'] = None
            self.p_meta['default_clustering'] = None
            self.clusterings = None

        profile_db.disconnect()

        self.progress.update('Accessing the auxiliary data file')
        auxiliary_data_path = os.path.join(os.path.dirname(self.profile_db_path), 'AUXILIARY-DATA.h5')
        if not os.path.exists(auxiliary_data_path):
            self.auxiliary_data_available = False
        else:
            self.auxiliary_data_available = True
            self.split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(auxiliary_data_path, self.p_meta['contigs_db_hash'])

        self.progress.end()

        if self.auxiliary_data_available:
            self.run.info('Auxiliary Data', 'Found: %s (v. %s)' % (auxiliary_data_path, anvio.__hdf5__version__))
        self.run.info('Profile DB', 'Initialized: %s (v. %s)' % (self.profile_db_path, anvio.__profile__version__))


    def init_gene_coverages_dict(self):
        if not self.a_meta['genes_are_called']:
            # genes were not identified/annotated
            return

        profile_db = ProfileDatabase(self.profile_db_path, quiet = True)

        self.progress.update('Reading gene coverages table ...')
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

            if not self.gene_coverages_dict.has_key(gene_callers_id):
                self.gene_coverages_dict[gene_callers_id] = {}

            self.gene_coverages_dict[gene_callers_id][gene_coverage_entry['sample_id']] = gene_coverage_entry['mean_coverage']

        self.progress.end()


    def get_variability_information_for_split(self, split_name, return_outliers = False, return_raw_results = False):
        if not split_name in self.split_names:
            raise ConfigError, "get_variability_information_for_split: The split name '%s' does not seem to be\
                                represented in this profile database. Are you sure you are looking for it\
                                in the right database?" % split_name

        profile_db = ProfileDatabase(self.profile_db_path)
        split_variability_information = profile_db.db.get_some_rows_from_table_as_dict(t.variable_nts_table_name, '''split_name = "%s"''' % split_name, error_if_no_data = False).values()
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

            d[e['sample_id']]['variability'][e['pos_in_codon']][e['pos']] = e['departure_from_consensus']
            d[e['sample_id']]['competing_nucleotides'][e['pos']] = e['competing_nts']

        return d


    def init_collection_profile(self, collection):
        profile_db = ProfileDatabase(self.profile_db_path, quiet = True)

        table_names = [table_name for table_name in t.atomic_data_table_structure[1:-1]]

        samples_template = dict([(s, []) for s in self.p_meta['samples']])

        # anonymous function to convert single profile table dicts compatible with merged ones (#155):
        SINGLE_P = lambda d: dict([(s, dict([(self.p_meta['samples'][0], v) for v in d[s].values()])) for s in d])

        # get coverage_table_data temporarily to learn split names. it sucks to do it this way, but before init_split_sequences()
        # is called, there is no other way really. and calling init_split_sequences() before this is also not feasible. so, here
        # we are.
        if self.p_meta['merged']:
            coverage_table_data = profile_db.db.get_table_as_dict('mean_coverage_splits', omit_parent_column = True)
        else:
            coverage_table_data = SINGLE_P(profile_db.db.get_table_as_dict('atomic_data_splits', columns_of_interest = "mean_coverage", omit_parent_column = True))

        self.split_names_not_binned = set(coverage_table_data.keys())

        for bin_id in collection:
            self.split_names_not_binned -= set(collection[bin_id])
            self.collection_profile[bin_id] = {}

        for table_name in table_names:
            if self.p_meta['merged']:
                table_data = profile_db.db.get_table_as_dict('%s_splits' % table_name, omit_parent_column = True)
            else:
                table_data = SINGLE_P(profile_db.db.get_table_as_dict('atomic_data_splits', columns_of_interest = table_name, omit_parent_column = True))

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
        self.bin_percent_recruitment_per_sample = {}
        for sample in self.p_meta['samples']:
            percents = {}
            all_coverages_in_sample = sum([d[sample] for d in coverage_table_data.values()])

            for bin_id in collection:
                bin_coverages_in_sample = sum([coverage_table_data[split_name][sample] for split_name in collection[bin_id]])
                percents[bin_id] = bin_coverages_in_sample * 100 / all_coverages_in_sample

            splits_not_binned_coverages_in_sample = sum([coverage_table_data[split_name][sample] for split_name in self.split_names_not_binned])
            percents['__splits_not_binned__'] = splits_not_binned_coverages_in_sample * 100 / all_coverages_in_sample
            self.bin_percent_recruitment_per_sample[sample] = percents

        profile_db.disconnect()


    def load_views(self, splits_of_interest = None):
        profile_db = ProfileDatabase(self.profile_db_path)

        views_table = profile_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            self.views[view] = {'table_name': table_name,
                                'header': profile_db.db.get_table_structure(table_name)[1:],
                                'dict': profile_db.db.get_table_as_dict(table_name, keys_of_interest = splits_of_interest)}

        profile_db.disconnect()


class DatabasesMetaclass(ProfileSuperclass, ContigsSuperclass, object):
    """Essential data to load for a given run"""
    def __init__(self, args, r = run, p = progress):
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
            self.db = db.DB(self.db_path, anvio.__profile__version__)
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

        self.db = db.DB(self.db_path, anvio.__profile__version__, new_database = True)

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

        # creating empty default tables
        self.db.create_table(t.clusterings_table_name, t.clusterings_table_structure, t.clusterings_table_types)
        self.db.create_table(t.gene_coverages_table_name, t.gene_coverages_table_structure, t.gene_coverages_table_types)
        self.db.create_table(t.variable_nts_table_name, t.variable_nts_table_structure, t.variable_nts_table_types)
        self.db.create_table(t.views_table_name, t.views_table_structure, t.views_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_colors_table_name, t.collections_colors_table_structure, t.collections_colors_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)

        self.disconnect()

        self.run.info('Contigs database', 'A new database, %s, has been created.' % (self.db_path), quiet = self.quiet)


    def disconnect(self):
        self.db.disconnect()


class ContigsDatabase:
    """To create an empty contigs database and/or access one."""
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
            is_contigs_db(self.db_path)
            self.db = db.DB(self.db_path, anvio.__contigs__version__)
            meta_table = self.db.get_table_as_dict('self')
            self.meta = dict([(k, meta_table[k]['value']) for k in meta_table])

            if 'creation_date' not in self.meta:
                raise ConfigError, "The contigs database ('%s') seems to be corrupted :/ This happens if the process that\
                                    that generates the database ends prematurely. Most probably, you will need to generate\
                                    the contigs database from scratch. Sorry!" % (self.db_path)

            self.run.info('Contigs database', 'An existing database, %s, has been initiated.' % self.db_path, quiet = self.quiet)
            self.run.info('Number of contigs', self.meta['num_contigs'], quiet = self.quiet)
            self.run.info('Number of splits', self.meta['num_splits'], quiet = self.quiet)
            self.run.info('Total number of nucleotides', self.meta['total_length'], quiet = self.quiet)
            self.run.info('Split length', self.meta['split_length'], quiet = self.quiet)
        else:
            self.db = None


    def create(self, args):
        A = lambda x: args.__dict__[x] if args.__dict__.has_key(x) else None
        contigs_fasta = A('contigs_fasta')
        split_length = A('split_length')
        kmer_size = A('kmer_size')
        skip_gene_calling = A('skip_gene_calling')
        skip_mindful_splitting = A('skip_mindful_splitting')
        debug = A('debug')
 
        filesnpaths.is_file_fasta_formatted(contigs_fasta)

        # just take a quick look at the first defline to make sure this FASTA file complies with anvi'o's
        # "simple defline" rule.
        fasta = u.SequenceSource(contigs_fasta)
        fasta.next()
        defline = fasta.id
        fasta.close()

        if not utils.check_contig_names(defline, dont_raise = True):
            raise ConfigError, "The FASTA file you provided does not comply with the 'simple deflines' requirement of\
                                anvi'o. Please read this section in the tutorial to understand the reason behind this\
                                requirement (anvi'o is very upset for making you do this): %s" \
                                                      % ('http://merenlab.org/2015/05/02/anvio-tutorial/#preparation')

        if os.path.exists(self.db_path):
            raise ConfigError, "Anvi'o will not overwrite an existing contigs database. Please choose a different name\
                                or remove the existing database ('%s') first." % (self.db_path)

        if not split_length:
            raise ConfigError, "Creating a new contigs database requires split length information to be\
                                provided. But the ContigsDatabase class was called to create one without this\
                                bit of information. Not cool."

        if not os.path.exists(contigs_fasta):
            raise ConfigError, "Creating a new contigs database requires a FASTA file with contigs to be provided."


        if not self.db_path.lower().endswith('.db'):
            raise ConfigError, "Please make sure your output file name has a '.db' extension. anvio developers apologize\
                                for imposing their views on how local databases should be named, and are humbled by your\
                                cooperation."

        try:
            split_length = int(split_length)
        except:
            raise ConfigError, "Split size must be an integer."

        if split_length <= 0:
            raise ConfigError, "For obvious reasons, split length must be a positive number. Nice try, though..."

        try:
            kmer_size = int(kmer_size)
        except:
            raise ConfigError, "K-mer size must be an integer."
        if kmer_size < 2 or kmer_size > 8:
            raise ConfigError, "We like our k-mer sizes between 2 and 8, sorry! (but then you can always change the\
                                source code if you are not happy to be told what you can't do, let us know how it goes!)."

        if skip_gene_calling:
            skip_mindful_splitting = True

        self.db = db.DB(self.db_path, anvio.__contigs__version__, new_database = True)

        # creating empty default tables
        self.db.create_table(t.hmm_hits_table_name, t.hmm_hits_table_structure, t.hmm_hits_table_types)
        self.db.create_table(t.hmm_hits_info_table_name, t.hmm_hits_info_table_structure, t.hmm_hits_info_table_types)
        self.db.create_table(t.hmm_hits_splits_table_name, t.hmm_hits_splits_table_structure, t.hmm_hits_splits_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_colors_table_name, t.collections_colors_table_structure, t.collections_colors_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.genes_in_contigs_table_name, t.genes_in_contigs_table_structure, t.genes_in_contigs_table_types)
        self.db.create_table(t.genes_in_splits_table_name, t.genes_in_splits_table_structure, t.genes_in_splits_table_types)
        self.db.create_table(t.splits_taxonomy_table_name, t.splits_taxonomy_table_structure, t.splits_taxonomy_table_types)
        self.db.create_table(t.contig_sequences_table_name, t.contig_sequences_table_structure, t.contig_sequences_table_types)
        self.db.create_table(t.gene_function_calls_table_name, t.gene_function_calls_table_structure, t.gene_function_calls_table_types)
        self.db.create_table(t.gene_protein_sequences_table_name, t.gene_protein_sequences_table_structure, t.gene_protein_sequences_table_types)
        self.db.create_table(t.genes_in_splits_summary_table_name, t.genes_in_splits_summary_table_structure, t.genes_in_splits_summary_table_types)

        # know thyself
        self.db.set_meta_value('db_type', 'contigs')
        # this will be the unique information that will be passed downstream whenever this db is used:
        self.db.set_meta_value('contigs_db_hash', '%08x' % random.randrange(16**8))
        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)

        # first things first: do the gene calling on contigs. this part is important. we are doing the
        # gene calling first. so we understand wher genes start and end. this information will guide the
        # arrangement of the breakpoint of splits
        contig_name_to_gene_start_stops = {}
        if not skip_gene_calling:
            # temporarily disconnect to perform gene calls
            self.db.disconnect()

            gene_calls_tables = TablesForGeneCalls(self.db_path, contigs_fasta, debug = debug)
            gene_calls_tables.call_genes_and_populate_genes_in_contigs_table()

            # reconnect and learn about what's done
            self.db = db.DB(self.db_path, anvio.__contigs__version__)

            genes_in_contigs_dict = self.db.get_table_as_dict(t.genes_in_contigs_table_name)

            for gene_unique_id in genes_in_contigs_dict:
                e = genes_in_contigs_dict[gene_unique_id]
                if not contig_name_to_gene_start_stops.has_key(e['contig']):
                    contig_name_to_gene_start_stops[e['contig']] = set([])

                contig_name_to_gene_start_stops[e['contig']].add((e['start'], e['stop']), )

        # lets process and store the FASTA file.
        fasta = u.SequenceSource(contigs_fasta)
        db_entries_contig_sequences = []

        contigs_kmer_table = KMerTablesForContigsAndSplits('kmer_contigs', k=kmer_size)
        splits_kmer_table = KMerTablesForContigsAndSplits('kmer_splits', k=kmer_size)

        contigs_info_table = InfoTableForContigs(split_length)
        splits_info_table = InfoTableForSplits()

        recovered_split_lengths = []

        self.progress.new('Identifying splits, and computing k-mer frequencies')
        while fasta.next():
            self.progress.update('Contig %d ...' % fasta.pos)

            contig_name = fasta.id
            contig_sequence = fasta.seq

            gene_start_stops = contig_name_to_gene_start_stops[contig_name] if contig_name in contig_name_to_gene_start_stops else set([])

            if skip_mindful_splitting:
                contig_length, split_start_stops, contig_gc_content = contigs_info_table.append(contig_name, contig_sequence, set([]))
            else:
                contig_length, split_start_stops, contig_gc_content = contigs_info_table.append(contig_name, contig_sequence, gene_start_stops)

            if len(split_start_stops) > 1:
                recovered_split_lengths.extend([s[1] - s[0] for s in split_start_stops])

            contig_kmer_freq = contigs_kmer_table.get_kmer_freq(contig_sequence)

            for order in range(0, len(split_start_stops)):
                start, end = split_start_stops[order]
                split_name = contigops.gen_split_name(contig_name, order)

                # this is very confusing, because both contigs_kmer_table and splits_kmer_able in fact
                # holds kmer values for splits only. in one table, each split has a kmer value of their
                # contigs (to not lose the genomic context while clustering based on kmers), in the other
                # one each split holds its own kmer value.
                contigs_kmer_table.append(split_name, contig_sequence[start:end], kmer_freq = contig_kmer_freq)
                splits_kmer_table.append(split_name, contig_sequence[start:end])

                splits_info_table.append(split_name, contig_sequence[start:end], order, start, end, contig_gc_content, contig_name)

            db_entries_contig_sequences.append((contig_name, contig_sequence), )
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

        self.run.info('Contigs database', 'A new database, %s, has been created.' % (self.db_path), quiet = self.quiet)
        self.run.info('Number of contigs', contigs_info_table.total_contigs, quiet = self.quiet)
        self.run.info('Number of splits', splits_info_table.total_splits, quiet = self.quiet)
        self.run.info('Total number of nucleotides', contigs_info_table.total_nts, quiet = self.quiet)
        self.run.info('Gene calling step skipped', skip_gene_calling, quiet = self.quiet)
        self.run.info("Splits broke genes (non-mindful mode)", skip_mindful_splitting, quiet = self.quiet)
        self.run.info('Desired split length (what the user wanted)', split_length, quiet = self.quiet)
        self.run.info("Average split length (wnat anvi'o gave back)", (int(round(numpy.mean(recovered_split_lengths)))) \
                                                                        if recovered_split_lengths \
                                                                            else "(Anvi'o did not create any splits)", quiet = self.quiet)


    def disconnect(self):
        self.db.disconnect()


class SamplesInformationDatabase:
    """To create an empty samples information database and/or access one.
    
       The purpose of this database is to deal with sample-specific information. Such as
       how should samples be organized in the interactive interface, or what environmental
       data available about them?
    """
    def __init__(self, db_path, run=run, progress=progress, quiet = True):
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
            self.sample_names_for_order = set([s.strip() for s in self.meta['sample_names_for_order'].split(',')])
            self.samples_information_default_layer_order = self.meta['samples_information_default_layer_order'].split(',')

            self.run.info('Samples information database', 'An existing database, %s, has been initiated.' % self.db_path, quiet = self.quiet)
        else:
            self.db = None


    def get_samples_information_and_order_dicts(self):
        if not self.db:
            raise ConfigError, "The samples database has not been initialized. You are doing something wrong :/"

        samples = samplesops.SamplesInformation()

        samples_information_dict = samples.recover_samples_information_dict(self.db.get_table_as_dict(t.samples_information_table_name, error_if_no_data = False),
                                                                            self.db.get_table_as_dict(t.samples_attribute_aliases_table_name, error_if_no_data = False))
        samples_order_dict = self.db.get_table_as_dict(t.samples_order_table_name)

        return samples_information_dict, samples_order_dict


    def get_samples_information_default_layer_order(self):
        if not self.db:
            raise ConfigError, "The samples database has not been initialized. You are doing something wrong :/"

        return self.samples_information_default_layer_order


    def create(self, samples_information_path = None, samples_order_path = None):
        if not samples_information_path and not samples_order_path:
            raise ConfigError, "You must declare at least one of the input files to create a samples information\
                                database. Neither samples information, nor samples order file has been passed to\
                                the class :("

        if os.path.exists(self.db_path):
            raise ConfigError, "Anvi'o will not overwrite an existing samples information database. Please choose a\
                                different name or remove the existing database ('%s') first." % (self.db_path)

        samples = samplesops.SamplesInformation()
        samples.populate_from_input_files(samples_information_path, samples_order_path)

        if not self.db_path.lower().endswith('.db'):
            raise ConfigError, "Please make sure your output file name has a '.db' extension. anvio developers apologize\
                                for imposing their views on how local databases should be named, and are humbled by your\
                                cooperation."

        self.db = db.DB(self.db_path, anvio.__samples__version__, new_database = True)

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

        self.run.info('Samples information database', 'A new samples information database, %s, has been created.' % (self.db_path), quiet = self.quiet)
        self.run.info('Number of samples', len(samples.samples_information_dict), quiet = self.quiet)
        self.run.info('Number of organizations', len(available_orders), quiet = self.quiet)

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


    def append(self, seq_id, sequence, gene_start_stops = None):
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
        self.set_next_available_id(t.variable_nts_table_name)


    def append(self, profile):
        db_entry = tuple([self.next_id(t.variable_nts_table_name)] + [profile[h] for h in t.variable_nts_table_structure[1:]])
        self.db_entries.append(db_entry)
        self.num_entries += 1
        if self.num_entries % 100 == 0:
            progress.update('Information for %d SNP sites have been added ...' % self.num_entries)


    def store(self):
        profile_db = ProfileDatabase(self.db_path)
        profile_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''' % t.variable_nts_table_name, self.db_entries)
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
    def __init__(self, db_path, contigs_fasta = None, run=run, progress=progress, debug = False):
        self.db_path = db_path
        self.contigs_fasta = contigs_fasta
        self.debug = debug

        self.gene_calls_dict_id_to_db_unique_id = {}

        is_contigs_db(self.db_path)

        if self.contigs_fasta:
            filesnpaths.is_file_exists(self.contigs_fasta)
            filesnpaths.is_file_fasta_formatted(self.contigs_fasta)


    def call_genes_and_populate_genes_in_contigs_table(self, gene_caller = 'prodigal'):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress, simple = True)

        self.set_next_available_id(t.genes_in_contigs_table_name)

        # get gene calls and protein sequences
        gene_calls_dict, protein_sequences = self.run_gene_caller(gene_caller)

        # get a unique_id conversion dict
        self.set_gene_calls_dict_id_to_db_unique_id(gene_calls_dict)

        # populate genes_in_contigs, and gene_protein_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, protein_sequences)


    def set_gene_calls_dict_id_to_db_unique_id(self, gene_calls_dict):
        self.gene_calls_dict_id_to_db_unique_id = dict([(entry_id, self.next_id(t.genes_in_contigs_table_name)) for entry_id in gene_calls_dict])


    def run_gene_caller(self, gene_caller = 'prodigal'):
        """Runs gene caller, and returns gene_calls_dict, and protein sequences."""
        remove_fasta_after_processing = False

        if not self.contigs_fasta:
            self.contigs_fasta = self.export_sequences_table_in_db_into_FASTA_file()
            remove_fasta_after_processing = True

        if self.debug:
            self.run.info_single('--debug flag is [ON], which means temporary directories generated by\
                                 this run will not be removed', nl_after = 2)

        gene_caller = genecalling.GeneCaller(self.contigs_fasta, gene_caller = gene_caller, debug = self.debug)

        gene_calls_dict, protein_sequences = gene_caller.process()

        if not self.debug and remove_fasta_after_processing:
            os.remove(self.contigs_fasta)

        return gene_calls_dict, protein_sequences


    def populate_genes_in_contigs_table(self, gene_calls_dict, protein_sequences):
        contigs_db = db.DB(self.db_path, anvio.__contigs__version__)

        self.progress.new('Entering gene calls into the database')
        counter = 0
        for entry_id in gene_calls_dict:
            if counter % 10 == 0:
                self.progress.update('%d ...' % (counter))
            contigs_db._exec('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % t.genes_in_contigs_table_name, tuple([self.gene_calls_dict_id_to_db_unique_id[entry_id]] + [gene_calls_dict[entry_id][h] for h in t.genes_in_contigs_table_structure[1:]]))
            contigs_db._exec('''INSERT INTO %s VALUES (?,?)''' % t.gene_protein_sequences_table_name, tuple([self.gene_calls_dict_id_to_db_unique_id[entry_id]] + [protein_sequences[entry_id]]))
            counter += 1
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
            if gene_calls_in_contigs_dict.has_key(contig):
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
                        genes_in_splits.add(split_name, start, stop, self.gene_calls_dict_id_to_db_unique_id[gene_callers_id], self.gene_calls_dict[gene_callers_id]['start'], self.gene_calls_dict[gene_callers_id]['stop'])


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
    def __init__(self, db_path, run=run, progress=progress):
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


    def populate_search_tables(self, sources = {}, protein_sequences_fasta = None):
        # if we end up generating a temporary file for protein sequences:
        remove_fasta_file_upon_finish = False

        if not len(sources):
            import anvio.data.hmm
            sources = anvio.data.hmm.sources

        if not sources:
            return

        if not protein_sequences_fasta:
            # creating a temporary fasta file for protein sequences:
            protein_sequences_fasta = self.export_sequences_table_in_db_into_FASTA_file(t.gene_protein_sequences_table_name)
            remove_fasta_file_upon_finish = True

        commander = HMMSearch(protein_sequences_fasta)

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
                parser = parser_modules['search']['hmmscan'](hmm_scan_hits_txt)
                search_results_dict = parser.get_search_results()

            self.append(source, reference, kind_of_search, all_genes_searched_against, search_results_dict)

        if not self.debug:
            commander.clean_tmp_dirs()

        if remove_fasta_file_upon_finish and not self.debug:
            os.remove(protein_sequences_fasta)


    def append(self, source, reference, kind_of_search, all_genes, search_results_dict):
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
        db_entries = [source, reference, kind_of_search, ', '.join(all_genes)]
        contigs_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.hmm_hits_info_table_name, db_entries)

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

            if hits_per_contig.has_key(contig_name):
                hits_per_contig[contig_name].append(hit)
            else:
                hits_per_contig[contig_name] = [hit]

        db_entries_for_splits = []

        for contig in self.contigs_info:
            if not hits_per_contig.has_key(contig):
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
    """Populates the collections_* tables, where clusters of contigs and splits are kept"""
    def __init__(self, db_path, version, run=run, progress=progress):
        self.db_path = db_path
        self.version = version

        Table.__init__(self, self.db_path, version, run, progress)

        # set these dudes so we have access to unique IDs:
        self.set_next_available_id(t.collections_colors_table_name)
        self.set_next_available_id(t.collections_contigs_table_name)
        self.set_next_available_id(t.collections_splits_table_name)


    def append(self, source, clusters_dict, cluster_colors = None):
        if not len(source):
            raise ConfigError, 'Source identifier cannot be empty.'

        if cluster_colors:
            if set(clusters_dict.keys()) - set(cluster_colors.keys()):
                raise ConfigError, 'Entries in the cluster dict do not match to entries in the colors dict.\
                                    They do not have to be identical, but for each cluster id, there must be a color\
                                    in the colors dict).'

        # remove any pre-existing information for 'source'
        self.delete_entries_for_key('source', source, [t.collections_info_table_name, t.collections_contigs_table_name, t.collections_splits_table_name, t.collections_colors_table_name])

        num_splits_in_clusters_dict = sum([len(splits) for splits in clusters_dict.values()])
        splits_in_clusters_dict = set(list(chain.from_iterable(clusters_dict.values())))
        if len(splits_in_clusters_dict) != num_splits_in_clusters_dict:
            raise ConfigError, "TablesForCollections::append: %d of the split or contig IDs appear more than once in\
                                your collections input. It is unclear to anvio how did you manage to do this, but we\
                                cannot go anywhere with this :/" % (num_splits_in_clusters_dict - len(splits_in_clusters_dict))

        database = db.DB(self.db_path, self.version)

        # how many clusters are defined in 'clusters_dict'?
        cluster_ids = clusters_dict.keys()

        # push information about this search result into serach_info table.
        db_entries = tuple([source, num_splits_in_clusters_dict, len(cluster_ids), ','.join(cluster_ids)])
        database._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_info_table_name, db_entries)

        if not cluster_colors:
            cluster_colors = utils.get_random_colors_dict(cluster_ids)

        # populate colors table.
        db_entries = [(self.next_id(t.collections_colors_table_name), source, cid, cluster_colors[cid]) for cid in cluster_ids]
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_colors_table_name, db_entries)

        # populate splits table
        db_entries = []
        for cluster_id in clusters_dict:
            for split_name in clusters_dict[cluster_id]:
                db_entries.append(tuple([self.next_id(t.collections_splits_table_name), source, split_name, cluster_id]))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_splits_table_name, db_entries)
        num_splits = len(db_entries)


        # FIXME: This function can be called to populate the contigs database (via anvi-populate-collections), or
        # the profile database. when it is contigs database, the superclass Table has the self.splits_info variable
        # set when it is initialized. however, the Table instance is missing self.splis when it is initialized with
        # the profile database. hence some special controls for contigs db (note that collections_contigs_table is
        # only populated in the contigs database):
        if self.db_type == 'contigs':
            splits_only_in_clusters_dict = [c for c in splits_in_clusters_dict if c not in self.splits_info]
            splits_only_in_db = [c for c in self.splits_info if c not in splits_in_clusters_dict]

            if len(splits_only_in_clusters_dict):
                self.run.warning('%d of %d splits found in "%s" results are not in the database. This may be OK,\
                                          but you must be the judge of it. If this is somewhat surprising, please use caution\
                                          and make sure all is fine before going forward with you analysis.'\
                                                % (len(splits_only_in_clusters_dict), len(splits_in_clusters_dict), source))

            if len(splits_only_in_db):
                self.run.warning('%d of %d splits found in the database were missing from the "%s" results. If this\
                                          does not make any sense, please make sure you know why before going any further.'\
                                                % (len(splits_only_in_db), len(self.splits_info), source))

            # then populate contigs table.
            db_entries = self.process_contigs(source, clusters_dict)
            database._exec_many('''INSERT INTO %s VALUES (?,?,?,?)''' % t.collections_contigs_table_name, db_entries)

        database.disconnect()

        self.run.info('Collections', '%s annotations for %d splits have been successfully added to the database at "%s".'\
                                        % (source, num_splits, self.db_path), mc='green')


    def process_contigs(self, source, clusters_dict):
        db_entries_for_contigs = []

        split_to_cluster_id = {}
        for cluster_id in clusters_dict:
            for split_name in clusters_dict[cluster_id]:
                split_to_cluster_id[split_name] = cluster_id

        contigs_processed = set([])
        for split_name in split_to_cluster_id:
            if split_name not in self.splits_info:
                # which means this split only appears in the input file, but not in the database.
                continue

            contig_name = self.splits_info[split_name]['parent']

            if contig_name in contigs_processed:
                continue
            else:
                contigs_processed.add(contig_name)

            db_entry = tuple([self.next_id(t.collections_contigs_table_name), source, contig_name, split_to_cluster_id[split_name]])
            db_entries_for_contigs.append(db_entry)

        return db_entries_for_contigs


class TablesForStates(Table):
    def __init__(self, db_path, version):
        self.db_path = db_path
        self.version = version
        self.states = {}

        Table.__init__(self, self.db_path, self.version, run, progress)

        self.init()


    def init(self):
        is_profile_db(self.db_path)

        profile_db = db.DB(self.db_path, self.version)
        self.states = profile_db.get_table_as_dict(t.states_table_name)
        profile_db.disconnect()


    def get_state(self, state_id):
        if state_id not in self.states:
            return None

        return self.states[state_id]


    def store_state(self, state_id, content, last_modified = None):
        self.remove_state(state_id)

        last_modified = datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S") if not last_modified else last_modified

        profile_db = db.DB(self.db_path, self.version)
        profile_db._exec('''INSERT INTO %s VALUES (?,?,?)''' % t.states_table_name, (state_id, content, last_modified))
        self.states = profile_db.get_table_as_dict(t.states_table_name)

        profile_db.disconnect()


    def remove_state(self, state_id):
        self.delete_entries_for_key('name', state_id, [t.states_table_name])



class TableForSplitsTaxonomy(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        # this class keeps track of genes that occur in splits, and responsible
        # for generating the necessary table in the contigs database
        self.genes_in_splits = GenesInSplits()


    def create(self, genes_dict, parser):
        self.genes_dict = genes_dict

        if parser:
            self.parser = parser
        else:
            self.parser = "user_input"

        self.sanity_check()

        # oepn connection
        contigs_db = ContigsDatabase(self.db_path)

        self.splits_info = contigs_db.db.get_table_as_dict(t.splits_info_table_name)

        # test whether there are already genes tables populated
        taxonomy_source = contigs_db.meta['taxonomy_source']
        if taxonomy_source:
            self.run.warning('Previous taxonomy information from "%s" will be replaced with the incoming data\
                              through "%s"' % (taxonomy_source, self.parser))
            contigs_db.db._exec('''DELETE FROM %s''' % (t.splits_taxonomy_table_name))

        # compute and push split taxonomy information.
        self.populate_split_taxonomy_information()

        # set the parser
        contigs_db.db.remove_meta_key_value_pair('taxonomy_source')
        contigs_db.db.set_meta_value('taxonomy_source', self.parser)

        # disconnect like a pro.
        contigs_db.disconnect()


    def sanity_check(self):
        # check whether input matrix dict 
        keys_found = ['prot'] + self.genes_dict.values()[0].keys()
        missing_keys = [key for key in t.splits_taxonomy_table_structure[1:] if key not in keys_found]
        if len(missing_keys):
            raise ConfigError, "Your dictionary is missing one or more keys that are necessary to populate the\
                                taxonomy table. Here is a list of missing keys: %s. The complete list of input\
                                keys should contain these: %s." % (', '.join(missing_keys),
                                                                   ', '.join(t.splits_taxonomy_table_structure))


    def populate_split_taxonomy_information(self):
        # build a dictionary for fast access to all proteins identified within a contig
        prots_in_contig = {}
        for prot in self.genes_dict:
            contig = self.genes_dict[prot]['contig']
            if prots_in_contig.has_key(contig):
                prots_in_contig[contig].add(prot)
            else:
                prots_in_contig[contig] = set([prot])

        contigs_without_annotation = list(set(self.contigs_info.keys()) - set(prots_in_contig.keys()))

        for contig in contigs_without_annotation:
            prots_in_contig[contig] = set([])

        splits_dict = {}

        num_splits_processed = 0
        num_splits_with_taxonomy = 0

        for contig in self.contigs_info:
            for split_name in self.contig_name_to_splits[contig]:
                num_splits_processed += 1
                start = self.splits_info[split_name]['start']
                stop = self.splits_info[split_name]['end']

                taxa = []
                for prot in prots_in_contig[contig]:
                    if self.genes_dict[prot]['stop'] > start and self.genes_dict[prot]['start'] < stop:
                        taxa.append(self.genes_dict[prot]['t_species'])

                taxonomy_strings = [tt for tt in taxa if tt]

                splits_dict[split_name] = dict([(f, None) for f in t.splits_taxonomy_table_structure])
                distinct_taxa = set(taxonomy_strings)

                if not len(distinct_taxa):
                    continue

                if len(distinct_taxa) == 1:
                    splits_dict[split_name]['t_species'] = distinct_taxa.pop()
                else:
                    d = Counter()
                    for taxon in taxonomy_strings:
                        d[taxon] += 1
                    consensus, occurrence = sorted(d.items(), key=operator.itemgetter(1))[-1]
                    splits_dict[split_name]['t_species'] = consensus

                num_splits_with_taxonomy += 1

        # open connection
        contigs_db = ContigsDatabase(self.db_path)

        # push taxonomy data
        db_entries = [tuple([split] + [splits_dict[split][h] for h in t.splits_taxonomy_table_structure[1:]]) for split in splits_dict]
        contigs_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?)''' % t.splits_taxonomy_table_name, db_entries)

        # disconnect
        contigs_db.disconnect()

        self.run.info('Splits taxonomy', 'Input data from parser "%s" annotated %d of %d splits (%.1f%%) with taxonomy.'\
                                            % (self.parser, num_splits_with_taxonomy, num_splits_processed,
                                               num_splits_with_taxonomy * 100.0 / num_splits_processed))

class TableForGeneFunctions(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)

        self.set_next_available_id(t.gene_function_calls_table_name)

    def create(self, functions_dict):
        self.sanity_check()

        # oepn connection
        contigs_db = ContigsDatabase(self.db_path)

        gene_function_sources = set([v['source'] for v in functions_dict.values()])
        unique_num_genes = len(set([v['gene_callers_id'] for v in functions_dict.values()]))

        # test whether there are already genes tables populated
        gene_function_sources_in_db = contigs_db.meta['gene_function_sources']
        if gene_function_sources_in_db:
            self.run.warning('Previous function calls for genes will be replaced with incoming data\
                             which contains %d hits originating %d sources: %s' % \
                             (len(functions_dict), len(gene_function_sources), ', '.join(gene_function_sources)))
            contigs_db.db._exec('''DELETE FROM %s''' % (t.gene_function_calls_table_name))

        # set the sources
        contigs_db.db.remove_meta_key_value_pair('gene_function_sources')
        contigs_db.db.set_meta_value('gene_function_sources', ','.join(gene_function_sources))

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


def is_profile_db(db_path):
    filesnpaths.is_file_exists(db_path)
    if get_db_type(db_path) != 'profile':
        raise ConfigError, "'%s' is not an anvi'o profile database." % db_path


def is_samples_db(db_path):
    filesnpaths.is_file_exists(db_path)
    if get_db_type(db_path) != 'samples_information':
        raise ConfigError, "'%s' is not an anvi'o samples database." % db_path


def get_db_type(db_path):
    try:
        database = db.DB(db_path, None, ignore_version = True)
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


def is_profile_db_and_samples_db_compatible(profile_db_path, samples_db_path):
    """Check whether every sample name in the profile database is represented in the samples information database"""
    is_profile_db(profile_db_path)
    is_samples_db(samples_db_path)

    profile_db = ProfileDatabase(profile_db_path)
    samples_db = SamplesInformationDatabase(samples_db_path)

    if profile_db.meta.has_key('merged') and not int(profile_db.meta['merged']):
        raise ConfigError, "Samples databases are only useful if you are working on a merged profile."

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
        split_names = set(profile_db.db.get_single_column_from_table('mean_coverage_Q1Q3_splits', 'contig'))
    else:
        split_names = set(profile_db.db.get_single_column_from_table('atomic_data_splits', 'contig'))

    profile_db.disconnect()

    return split_names


def add_hierarchical_clustering_to_db(profile_db_path, clustering_id, clustering_newick, make_default = False, run = run):
    is_profile_db(profile_db_path)
    utils.is_this_name_OK_for_database('clustering_id', clustering_id)

    profile_db = ProfileDatabase(profile_db_path)

    try:
        available_clusterings = profile_db.db.get_meta_value('available_clusterings').split(',')
    except:
        available_clusterings = []

    if clustering_id in available_clusterings:
        run.warning('Clustering for the ID "%s" is already in the database. Its content will be replaced with\
                     the new result.' % clustering_id)

        profile_db.db._exec('''DELETE FROM %s where clustering = "%s"''' % (t.clusterings_table_name, clustering_id))
    else:
        available_clusterings.append(clustering_id)

    profile_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % t.clusterings_table_name, tuple([clustering_id, clustering_newick]))

    try:
        profile_db.db.remove_meta_key_value_pair('available_clusterings')
    except:
        pass
    profile_db.db.set_meta_value('available_clusterings', ','.join(available_clusterings))

    try:
        profile_db.db.remove_meta_key_value_pair('contigs_clustered')
    except:
        pass
    profile_db.db.set_meta_value('contigs_clustered', True)

    try:
        profile_db.db.get_meta_value('default_clustering')
        default_clustering_is_set = True
    except:
        default_clustering_is_set = False

    if make_default or not default_clustering_is_set:
        try:
            profile_db.db.remove_meta_key_value_pair('default_clustering')
        except:
            pass
        profile_db.db.set_meta_value('default_clustering', clustering_id)

    profile_db.disconnect()

    run.info('New hierarchical clusetring', '"%s" has been added to the database...' % clustering_id)

