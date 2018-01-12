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

from io import StringIO
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
import anvio.constants as constants
import anvio.contigops as contigops
import anvio.samplesops as samplesops
import anvio.filesnpaths as filesnpaths
import anvio.genecalling as genecalling
import anvio.auxiliarydataops as auxiliarydataops
import anvio.genomestorage as genomestorage

from anvio.tableops import Table
from anvio.drivers import Aligners
from anvio.errors import ConfigError
from anvio.drivers.hmmer import HMMer
from anvio.parsers import parser_modules
from anvio.sequence import get_list_of_outliers


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
aligners = Aligners()


class DBClassFactory:
    """Factory pattern to get the appropriate class for a given anvi'o db type"""
    def __init__(self):
        self.DB_CLASSES = {'profile': ProfileDatabase,
                           'contigs': ContigsDatabase,
                           'pan': PanDatabase}

    def get_db_class(self, db_path):
        db_type = get_db_type(db_path)

        if db_type not in self.DB_CLASSES:
            raise ConfigError("DBClassFactory speaking. I do not know a class for database type\
                                %s :/ I can deal with these though: '%s'" % (', '.join(self.DB_CLASSES)))

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
        self.gene_lengths = {}
        self.contig_name_to_genes = {}
        self.genes_in_splits = {} # keys of this dict are NOT gene caller ids. they are ids for each entry.
        self.genes_in_splits_summary_dict = {}
        self.genes_in_splits_summary_headers = []
        self.split_name_to_genes_in_splits_entry_ids = {} # for fast access to all self.genes_in_splits entries for a given split
        self.gene_callers_id_to_split_name_dict = {} # for fast access to a split name that contains a given gene callers id

        self.nt_positions_info = None

        self.gene_function_call_sources = []
        self.gene_function_calls_dict = {}
        self.gene_function_calls_initiated = False

        self.hmm_sources_info = {}
        self.hmm_searches_dict = {}   # <--- upon initiation, this dict only keeps hmm hits for non-singlecopy
        self.hmm_searches_header = [] #      gene searches... single-copy gene info is accessed through completeness.py

        self.singlecopy_gene_hmm_sources = set([])
        self.non_singlecopy_gene_hmm_sources = set([])

        # now all items are initialized, we will check whether we are being initialized from within
        # an object that is in `pan` or `manual` mode, neither of which will have a contigs database
        # associated with the call. so having done our part, we will quietly return from here hoping
        # that we are not driving a developer crazy somewhere by doing so.
        A = lambda x: self.__dict__[x] if x in self.__dict__ else None
        if A('mode') == 'pan' or A('mode') == 'manual':
            return

        self.contigs_db_path = args.contigs_db

        if not self.contigs_db_path:
            raise ConfigError("Someone (hopefully, you) is trying to initialize the Contigs Super Class without a contigs database path.\
                               There are many ways this can happen, but .. do you think you were trying to run anvi-interactive in\
                               manual mode but without a --manual flag right before this? Just a gut feeling... No? Well, then we \
                               are really in a big trouble. Please run what you did before seeing this again with a `--debug` flag,\
                               and send us an e-mail :(")

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

        self.progress.update('Populating gene lengths dict')
        self.gene_lengths = dict([(g, (self.genes_in_contigs_dict[g]['stop'] - self.genes_in_contigs_dict[g]['start'])) for g in self.genes_in_contigs_dict])

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

        self.singlecopy_gene_hmm_sources = set([s for s in list(self.hmm_sources_info.keys()) if self.hmm_sources_info[s]['search_type'] == 'singlecopy'])
        self.non_singlecopy_gene_hmm_sources = set([s for s in list(self.hmm_sources_info.keys()) if self.hmm_sources_info[s]['search_type'] != 'singlecopy'])

        self.progress.update('Generating "split name" to "gene entry ids" mapping dict')
        for entry_id in self.genes_in_splits:
            split_name = self.genes_in_splits[entry_id]['split']
            if split_name in self.split_name_to_genes_in_splits_entry_ids:
                self.split_name_to_genes_in_splits_entry_ids[split_name].add(entry_id)
            else:
                self.split_name_to_genes_in_splits_entry_ids[split_name] = set([entry_id])

        for split_name in self.splits_basic_info:
            if split_name not in self.split_name_to_genes_in_splits_entry_ids:
                self.split_name_to_genes_in_splits_entry_ids[split_name] = set([])

        self.progress.update('Generating "gene caller id" to "split name" mapping dict')
        for entry in list(self.genes_in_splits.values()):
            self.gene_callers_id_to_split_name_dict[entry['gene_callers_id']] = entry['split']

        self.progress.update('Accessing the auxiliary data file')

        self.nt_positions_info = {}
        for contig_name, nt_position_row in contigs_db.db.get_table_as_dict(t.nt_position_info_table_name).items():
            self.nt_positions_info[contig_name] = utils.convert_binary_blob_to_numpy_array(nt_position_row['position_info'], 'uint8')

        self.progress.end()

        contigs_db.disconnect()

        self.run.info('Contigs DB', 'Initialized: %s (v. %s)' % (self.contigs_db_path, anvio.__contigs__version__))


    def init_splits_taxonomy(self, t_level = 't_genus'):
        if not self.contigs_db_path:
            return

        if t_level not in t.taxon_names_table_structure[1:]:
            raise ConfigError("Pretty close. But the taxonomic level '%s' is not known to anvi'o. How about\
                                one of these: %s." % (t_level, ','.join(t.taxon_names_table_structure[1:])))

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

            for e in list(non_singlecopy_gene_hmm_results_dict.values()):
                hmm_hit = hmm_hits_table[e['hmm_hit_entry_id']]
                search_term = 'hmmx_%s_%s' % (self.hmm_sources_info[e['source']]['search_type'], hmm_hit['gene_name'])
                self.hmm_searches_dict[e['split']][search_term] = 1
        else:
            for source in self.non_singlecopy_gene_hmm_sources:
                search_type = 'hmms_%s' % self.hmm_sources_info[source]['search_type']
                sources_tmpl[source] = []
                self.hmm_searches_header.append((search_type, source),)

            for e in list(non_singlecopy_gene_hmm_results_dict.values()):
                hmm_hit = hmm_hits_table[e['hmm_hit_entry_id']]
                hmm_source = e['source']

                if not e['split'] in self.hmm_searches_dict:
                    self.hmm_searches_dict[e['split']] = copy.deepcopy(sources_tmpl)

                search_type = 'hmms_%s' % self.hmm_sources_info[e['source']]['search_type']

                # populate hmm_searches_dict with hmm_hit and unique identifier (see #180):
                self.hmm_searches_dict[e['split']][hmm_source].append((hmm_hit['gene_name'], hmm_hit['gene_unique_identifier']),)

        self.progress.end()


    def get_nt_position_info(self, contig_name, pos_in_contig):
        """This function returns a tuple with three items for each nucleotide position.

            (in_partial_gene_call, in_complete_gene_call, base_pos_in_codon)

        See `init_nt_position_info_dict` for more info."""

        if not self.nt_positions_info:
            raise ConfigError("get_nt_position_info: I am asked to return stuff, but self.nt_position_info is None!")

        if not contig_name in self.nt_positions_info:
            return (0, 0, 0)

        position_info = self.nt_positions_info[contig_name][pos_in_contig]

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
                    raise ConfigError("Some of the functional sources you requested are missing from the contigs database '%s'. Here\
                                        they are (or here it is, whatever): %s." % \
                                                    (self.contigs_db_path, ', '.join(["'%s'" % s for s in missing_sources])))

            hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name,
                                                                  '''source IN (%s)''' % (', '.join(["'%s'" % s for s in requested_sources])),
                                                                  error_if_no_data=False).values())
            self.gene_function_call_sources = requested_sources
        else:
            hits = list(contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name).values())
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
                if self.gene_function_calls_dict[gene_callers_id][source][2] < e_value:
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
            raise ConfigError("Search terms must be of type 'list'")

        search_terms = [s.strip() for s in search_terms]

        if len([s.strip().lower() for s in search_terms]) != len(set([s.strip().lower() for s in search_terms])):
            raise ConfigError("Please do not use the same search term twice :/ Becasue, reasons. You know.")

        for search_term in search_terms:
            if not len(search_term) >= 3:
                raise ConfigError("A search term cannot be less than three characters")

        self.run.info('Search terms', '%d found' % (len(search_terms)))
        matching_gene_caller_ids = dict([(search_term, {}) for search_term in search_terms])
        matching_accession_calls = dict([(search_term, {}) for search_term in search_terms])
        matching_function_calls = dict([(search_term, {}) for search_term in search_terms])
        split_names = dict([(search_term, {}) for search_term in search_terms])
        full_report = []

        if not self.gene_function_calls_initiated:
            self.init_functions()

        contigs_db = ContigsDatabase(self.contigs_db_path)

        for search_term in search_terms:
            self.progress.new('Search function')
            self.progress.update('Searching for term "%s"' % search_term)
            response = contigs_db.db._exec('''select gene_callers_id, source, accession, function from gene_functions where function LIKE "%%''' + search_term + '''%%" OR accession LIKE "%%''' + search_term + '''%%";''').fetchall()

            full_report.extend([(r[0], r[1], r[2], r[3], search_term, self.gene_callers_id_to_split_name_dict[r[0]]) for r in response])

            matching_gene_caller_ids[search_term] = set([m[0] for m in response])
            matching_accession_calls[search_term] = list(set([m[2] for m in response]))
            matching_function_calls[search_term] = list(set([m[3] for m in response]))
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
            raise ConfigError("get_corresponding_codon_order_in_gene :: pos_in_contig must be of type 'int'")

        if not isinstance(gene_caller_id, int):
            raise ConfigError("get_corresponding_codon_order_in_gene :: gene_caller_id must be of type 'int'")

        gene_call = self.genes_in_contigs_dict[gene_caller_id]

        if contig_name != gene_call['contig']:
            raise ConfigError('get_corresponding_codon_order_in_gene :: well, the gene call %d and the contig %s\
                                do not seem to have anything to do with each other :/ This is not a user-level error\
                                something must have gone very wrong somewhere in the code ...' % (gene_caller_id, contig_name))

        if not pos_in_contig >= gene_call['start'] or not pos_in_contig < gene_call['stop']:
            raise ConfigError("get_corresponding_codon_order_in_gene :: position %d does not occur in gene call %d :(" \
                                                        % (pos_in_contig, gene_caller_id))

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


    def get_AA_counts_dict(self, split_names=set([]), contig_names=set([]), gene_caller_ids=set([]), return_codons_instead=False):
        """Returns a dictionary of AA counts.

           The dict can be returned for a given collection of split names, contigs names,
           or gene calls. If none of these variables are specified, the dict will contain
           counts for all gene calls in the contigs database"""

        counts_dict = {}

        # nothing to do here if the genes were not called:
        if not self.a_meta['genes_are_called']:
            return counts_dict

        if len([True for v in [split_names, contig_names, gene_caller_ids] if v]) > 1:
            raise ConfigError("get_AA_counts_dict :: If you want to get AA counts for a specific\
                                set of split names, contig names, or gene call ids, that is totally\
                                fine. But you can't request more than one at a time.")

        # we need to understand what genes we're interested in first. it could be genes in
        # a collection, or it could be everything in the contigs database, etc
        gene_calls_of_interest = set([])

        if split_names:
            for split_name in split_names:
                for entry_id in self.split_name_to_genes_in_splits_entry_ids[split_name]:
                    gene_calls_of_interest.add(self.genes_in_splits[entry_id]['gene_callers_id'])
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

        # sequences is a list of gene sequences, where each gene sequence is itself a list
        # of either AAs, e.g. ["Ala","Ala","Trp",...] or codons, e.g. ["AAA","AAT","CCG",...]
        sequences = []
        for gene_call_id in gene_calls_of_interest:
            gene_call = self.genes_in_contigs_dict[gene_call_id]

            if gene_call['partial']:
                continue

            if return_codons_instead:
                sequences.extend(utils.get_list_of_codons_for_gene_call(gene_call, self.contig_sequences))
            else:
                sequences.extend(utils.get_list_of_AAs_for_gene_call(gene_call, self.contig_sequences))

        counts_dict['counts'] = Counter(sequences)
        counts_dict['total'] = sum(Counter(sequences).values())
        counts_dict['total_gene_calls'] = len(gene_calls_of_interest)

        # add missing AAs or codons into the dict ... if there are any
        items = constants.codons if return_codons_instead else constants.AA_to_codons.keys()
        for item in items:
            if item not in counts_dict['counts']:
                counts_dict['counts'][item] = 0

        return counts_dict


    def get_corresponding_gene_caller_ids_for_base_position(self, contig_name, pos_in_contig):
        """For a given nucleotide position and contig name, returns all matching gene caller ids"""
        gene_start_stops_in_contig = self.get_gene_start_stops_in_contig(contig_name)

        if not gene_start_stops_in_contig:
            return []

        corresponding_gene_calls = [gene_callers_id for (gene_callers_id, start, stop) in gene_start_stops_in_contig if pos_in_contig >= start and pos_in_contig < stop]

        return corresponding_gene_calls


    def get_sequences_for_gene_callers_ids(self, gene_caller_ids_list, reverse_complement_if_necessary=True):
        if not isinstance(gene_caller_ids_list, list):
            raise ConfigError("Gene caller's ids must be of type 'list'")

        if not len(gene_caller_ids_list):
            gene_caller_ids_list = list(self.genes_in_contigs_dict.keys())
            self.run.warning("You did not provide any gene caller ids. As a result, anvi'o will give you back sequences for every\
                              %d gene call stored in the contigs database. %s" % (len(gene_caller_ids_list), ' Brace yourself.' if len(gene_caller_ids_list) > 10000 else ''))

        try:
            gene_caller_ids_list = [int(gene_callers_id) for gene_callers_id in gene_caller_ids_list]
        except:
            raise ConfigError("List of IDs for gene calls contains non-integer values :/")

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
            raise ConfigError("We need an explicit output file path. Anvi'o does not know how you managed to come \
                               here, but please go back and come again.")

        filesnpaths.is_output_file_writable(output_file_path)

        if not isinstance(wrap, int):
            raise ConfigError('"wrap" has to be an integer instance')
        if wrap == 0:
            wrap = None
        if wrap and wrap <= 20:
            raise ConfigError('Value for wrap must be larger than 20. Yes. Rules.')

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


    def gen_GFF3_file_of_sequences_for_gene_caller_ids(self, gene_caller_ids_list=[], output_file_path=None, wrap=120, simple_headers=False, rna_alphabet=False):
        gene_caller_ids_list, sequences_dict = self.get_sequences_for_gene_callers_ids(gene_caller_ids_list)

        name_template = '' if simple_headers else ';Name={contig} {start} {stop} {direction} {rev_compd} {length}'

        self.progress.new('Storing sequences')
        self.progress.update('...')
        with open(output_file_path, 'wt') as output:
            output.write('##gff-version 3\n')
            for gene_callers_id in gene_caller_ids_list:
                entry = sequences_dict[gene_callers_id]
                output.write('{id}\t{source}\t{contig}\t1\t{length}\t.\t.\t.\tID={id}'.format(
                    id=gene_callers_id, source='IGS', contig=entry['contig'], length=entry['length']))
                output.write(name_template.format(entry))
                output.write('\n')

        self.progress.end()
        self.run.info('Output', output_file_path)


    def gen_TAB_delimited_file_for_split_taxonomies(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        if not self.a_meta['taxonomy_source']:
            raise ConfigError("There is no taxonomy source in the contigs database :/")

        if not len(self.splits_taxonomy_dict):
            self.init_splits_taxonomy()

        if not len(self.splits_taxonomy_dict):
            raise ConfigError("The splits taxonomy is empty. There is nothing to report. Could it be\
                                possible the taxonomy caller you used did not assign any taxonomy to\
                                anything?")

        self.run.info("Taxonomy", "Annotations for %d of %d total splits are recovered" % (len(self.splits_taxonomy_dict), len(self.splits_basic_info)))

        output = open(output_file_path, 'w')
        for split_name in self.splits_basic_info:
            if split_name in self.splits_taxonomy_dict:
                output.write('{0}\t{1}\n'.format(split_name, self.splits_taxonomy_dict[split_name]))
            else:
                output.write('{0}\t\n'.format(split_name))
        output.close()

        self.run.info("Output", output_file_path)


class PanSuperclass(object):
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        self.genome_names = []
        self.gene_clusters = {}
        self.gene_clusters_initialized = False
        self.gene_cluster_names = set([])
        self.gene_clusters_gene_alignments = {}
        self.gene_clusters_gene_alignments_available = False
        self.gene_clusters_function_sources = []
        self.gene_clusters_functions_dict = {}
        self.item_orders = {}
        self.views = {}
        self.collection_profile = {}

        # the following two are initialized via `init_items_additional_data()` and use information
        # stored in item additional data tables in the pan database
        self.items_additional_data_dict = None
        self.items_additional_data_keys = None

        self.num_gene_clusters = None
        self.num_genes_in_gene_clusters = None

        self.genomes_storage_is_available = False
        self.genomes_storage_has_functions = False
        self.functions_initialized = False

        try:
            self.pan_db_path = self.args.pan_db
        except:
            self.run.warning('PanSuperclass class called with args without pan_db_path member! Returning prematurely.')
            return

        filesnpaths.is_file_exists(self.pan_db_path)

        self.progress.new('Initializing the pan database superclass')

        self.progress.update('Creating an instance of the pan database')
        pan_db = PanDatabase(self.pan_db_path, run=self.run, progress=self.progress)

        self.progress.update('Setting profile self data dict')
        self.p_meta = pan_db.meta

        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if 'creation_date' in self.p_meta else 'unknown'
        self.p_meta['genome_names'] = sorted([s.strip() for s in self.p_meta['external_genome_names'].split(',') + self.p_meta['internal_genome_names'].split(',') if s])
        self.p_meta['num_genomes'] = len(self.p_meta['genome_names'])
        self.genome_names = self.p_meta['genome_names']
        self.gene_clusters_gene_alignments_available = self.p_meta['gene_alignments_computed']

        # FIXME: Is this the future where the pan db version is > 6? Great. Then the if statement here no longer
        # needs to check whether 'gene_clusters_ordered' is a valid key in self.p_meta:
        if 'gene_clusters_ordered' in self.p_meta and self.p_meta['gene_clusters_ordered']:
            self.p_meta['available_item_orders'] = sorted([s.strip() for s in self.p_meta['available_item_orders'].split(',')])
            self.item_orders = pan_db.db.get_table_as_dict(t.item_orders_table_name)

            # we need to convert data for 'basic' item orders to array in order to avoid compatibility issues with
            # other additional item orders in pan and full mode (otherwise interactive class gets complicated
            # unnecessarily).
            for item_order in self.item_orders:
                if self.item_orders[item_order]['type'] == 'basic':
                    try:
                        self.item_orders[item_order]['data'] = self.item_orders[item_order]['data'].split(',')
                    except:
                        raise ConfigError("Something is wrong with the basic order `%s` in this pan database :(" % (item_order))
        else:
            self.p_meta['available_item_orders'] = None
            self.p_meta['default_item_order'] = None
            self.item_orders = None

        # recover all gene cluster names so others can access to this information
        # without having to initialize anything
        self.gene_cluster_names = set(pan_db.db.get_single_column_from_table(t.pan_gene_clusters_table_name, 'gene_cluster_id'))

        pan_db.disconnect()

        # create an instance of states table
        self.states_table = TablesForStates(self.pan_db_path)

        self.progress.end()

        if 'genomes_storage' in args.__dict__ and args.genomes_storage:
            self.genomes_storage = genomestorage.GenomeStorage(args.genomes_storage,
                                                                       self.p_meta['genomes_storage_hash'],
                                                                       genome_names_to_focus=self.p_meta['genome_names'],
                                                                       run=self.run,
                                                                       progress=self.progress)
            self.genomes_storage_is_available = True
            self.genomes_storage_has_functions = self.genomes_storage.functions_are_available

        self.run.info('Pan DB', 'Initialized: %s (v. %s)' % (self.pan_db_path, anvio.__pan__version__))


    def get_sequences_for_gene_clusters(self, gene_cluster_names=set([]), skip_alignments=False, report_DNA_sequences=False):
        """Returns a dictionary of sequences (aligned or not) in a given gene cluster:

        {
            'GENE_CLUSTER_NAME_01': {
                            'GENOME_NAME_A': [('gene_callers_id_x', 'sequence_x'),
                                              ('gene_callers_id_y', 'sequence_y')],
                            'GENOME_NAME_B': [('gene_callers_id_z', 'sequence_z')],
                            (...)
                          },
            'GENE_CLUSTER_NAME_02': {
                            (...)
                          },
            (...)
        }

        By default, it will return amino acid sequences. You can ask for DNA sequences if setting
        the flag `report_DNA_sequences` True.

        """

        sequences = {}

        if not isinstance(gene_cluster_names, type(set([]))) or not gene_cluster_names:
            raise ConfigError("gene_cluster_names for get_sequences_for_gene_clusters must be a non-empty `list`.")

        if not self.genomes_storage_is_available:
            raise ConfigError("The pan anvi'o super class for is upset. You are attempting to get AA seqeunces for %s,\
                               but there is not genomes storage is available to get it." \
                                    % 'a gene cluster' if len(gene_cluster_names) > 1 else '%d gene_clusters' % len(gene_cluster_names))

        if not self.gene_clusters_initialized:
            self.init_gene_clusters()

        missing_gene_cluster_names = [p for p in gene_cluster_names if p not in self.gene_clusters]
        if len(missing_gene_cluster_names[0:5]):
            raise ConfigError("get_sequences_for_gene_clusters: %d of %d gene clusters are missing in the pan database. Not good :/\
                               Here are some of the missing ones; %s" \
                                        % (len(missing_gene_cluster_names), len(gene_cluster_names), ', '.join(missing_gene_cluster_names[0:5])))

        self.progress.new('Accessing gene cluster seqeunces')

        for gene_cluster_name in gene_cluster_names:
            self.progress.update("processing '%s' ..." % gene_cluster_name )
            sequences[gene_cluster_name] = {}
            for genome_name in self.gene_clusters[gene_cluster_name]:
                sequences[gene_cluster_name][genome_name] = {}
                for gene_callers_id in self.gene_clusters[gene_cluster_name][genome_name]:
                    sequence = self.genomes_storage.get_gene_sequence(genome_name, gene_callers_id, report_DNA_sequences=report_DNA_sequences)

                    if not skip_alignments and self.gene_clusters_gene_alignments_available:
                        alignment_summary = self.gene_clusters_gene_alignments[genome_name][gene_callers_id]
                        sequence = utils.restore_alignment(sequence, alignment_summary, from_aa_alignment_summary_to_dna=report_DNA_sequences)

                    sequences[gene_cluster_name][genome_name][gene_callers_id] = sequence

        self.progress.end()

        return sequences


    def write_sequences_in_gene_clusters_to_file(self, gene_cluster_names=set([]), skip_alignments=False, output_file_path=None, report_DNA_sequences=False):
        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        output_file = open(output_file_path, 'w')
        sequences = self.get_sequences_for_gene_clusters(gene_cluster_names=gene_cluster_names, skip_alignments=skip_alignments, report_DNA_sequences=report_DNA_sequences)

        self.progress.new('Writing gene cluster seqeunces to file')
        sequence_counter = 0
        for gene_cluster_name in gene_cluster_names:
            for genome_name in sequences[gene_cluster_name]:
                for gene_callers_id in sequences[gene_cluster_name][genome_name]:
                        output_file.write('>%08d|gene_cluster:%s|genome_name:%s|gene_callers_id:%d\n' % (sequence_counter,
                                                                                                         gene_cluster_name,
                                                                                                         genome_name,
                                                                                                         gene_callers_id))
                        output_file.write('%s\n' % sequences[gene_cluster_name][genome_name][gene_callers_id])
                        sequence_counter += 1
                        self.progress.update("processing '%s' ..." % gene_cluster_name)

        self.progress.end()
        output_file.close()

        self.run.info('Sequence type', 'DNA' if report_DNA_sequences else 'Amino acid', mc='green')
        self.run.info('Output FASTA file', output_file_path, mc='green')


    def write_sequences_in_gene_clusters_for_phylogenomics(self, gene_cluster_names=set([]), skip_alignments=False, output_file_path=None, skip_multiple_gene_calls=False, report_DNA_sequences=False, align_with=None):
        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        output_file = open(output_file_path, 'w')
        sequences = self.get_sequences_for_gene_clusters(gene_cluster_names=gene_cluster_names, skip_alignments=skip_alignments, report_DNA_sequences=report_DNA_sequences)


        if not self.gene_clusters_gene_alignments_available:
            aligner = aligners.select(align_with)

            run.warning("It seems sequences in gene clusters were not aligned during the pangenomic analysis, so we\
                         are going to have do it now .. which may take some time .. and it is totally your fault :/")
            progress.new("Aligning sequences")

        get_first_value = lambda x: next(iter(x.values()))
        get_first_key = lambda x: next(iter(x.keys()))

        silent_run = terminal.Run()
        silent_run.verbose = False

        output_buffer = dict({})
        for genome_name in self.genome_names:
            output_buffer[genome_name] = StringIO()

        skipped_gene_clusters = []
        for gene_cluster_name in gene_cluster_names:
            multiple_gene_calls = False
            multiple_gene_call_genome = None
            sequence_length = None

            for genome_name in self.genome_names:
                if len(sequences[gene_cluster_name][genome_name]) > 1:
                    multiple_gene_calls = True
                    multiple_gene_call_genome = genome_name
                elif self.gene_clusters_gene_alignments_available and len(sequences[gene_cluster_name][genome_name]) == 1:
                    sequence_length = len(get_first_value(sequences[gene_cluster_name][genome_name]))

            if multiple_gene_calls:
                if skip_multiple_gene_calls:
                    skipped_gene_clusters.append(gene_cluster_name)
                    continue
                else:
                    raise ConfigError("There are multiple gene calls in '%s' and sample '%s', if you want to continue use flag \
                                        --skip-multiple-gene-calls" % (gene_cluster_name, multiple_gene_call_genome))

            if not self.gene_clusters_gene_alignments_available:
                sequences_to_align = []
                for genome_name in self.genome_names:
                    if len(sequences[gene_cluster_name][genome_name]) == 1:
                        sequences_to_align.append((genome_name, get_first_value(sequences[gene_cluster_name][genome_name])))

                progress.update("Processing '" + gene_cluster_name + "'")

                aligned_sequences = aligner(run=silent_run).run_stdin(sequences_list=sequences_to_align)

                for genome_name in aligned_sequences:
                    gene_caller_id = get_first_key(sequences[gene_cluster_name][genome_name])
                    sequences[gene_cluster_name][genome_name][gene_caller_id] = aligned_sequences[genome_name]

                    if not sequence_length:
                        sequence_length = len(aligned_sequences[genome_name])

            for genome_name in self.genome_names:
                if len(sequences[gene_cluster_name][genome_name]) == 1:
                    output_buffer[genome_name].write(get_first_value(sequences[gene_cluster_name][genome_name]))
                else:
                    output_buffer[genome_name].write("-" * sequence_length)

        if not self.gene_clusters_gene_alignments_available:
            progress.end()

        if len(skipped_gene_clusters):
            self.run.warning("%s of %s gene_clusters contained multiple gene calls, and skipped during concatenation.\n '%s'" \
                                                        % (pp(len(skipped_gene_clusters)), pp(len(gene_cluster_names)), ", ".join(skipped_gene_clusters)))

        if len(skipped_gene_clusters) == len(gene_cluster_names):
            raise ConfigError("Well. You have no gene_clusters left.. Bye :/")

        for genome_name in self.genome_names:
            output_file.write('>%s\n' % genome_name)
            output_file.write(output_buffer[genome_name].getvalue())
            output_file.write('\n')
            output_buffer[genome_name].close()

        output_file.close()

        self.run.info('Sequence type', 'DNA' if report_DNA_sequences else 'Amino acid', mc='green')
        self.run.info('Output file for phylogenomics', output_file_path, mc='green')


    def init_gene_clusters_functions(self):
        self.progress.new('Initializing functions for gene clusters')
        self.progress.update('...')
        if not self.gene_clusters:
            raise ConfigError("init_gene_clusters_functions is speaking! You called this function before you initialized\
                                gene clusters :/ One of us does not know what they're doing :(")

        if not self.genomes_storage_has_functions:
            self.progress.end()
            self.run.warning("Genomes storage does not have any info about gene functions. Certain parts of the pangenomic\
                              workflow will not be accessible.")
            return

        # FIXME WE HAVE TO STORE AVAILABLE FUNCTIONS IN GENOMES STORAGE ATTRs!!!! THIS IS RIDICULOUS
        self.gene_clusters_function_sources = set([])
        for gene_cluster_id in self.gene_clusters:
            self.gene_clusters_functions_dict[gene_cluster_id] = {}
            for genome_name in self.genome_names:
                self.gene_clusters_functions_dict[gene_cluster_id][genome_name] = {}
                for gene_callers_id in self.gene_clusters[gene_cluster_id][genome_name]:
                    functions = self.genomes_storage.get_gene_functions(genome_name, gene_callers_id)
                    self.gene_clusters_functions_dict[gene_cluster_id][genome_name][gene_callers_id] = functions

                    if functions:
                        self.gene_clusters_function_sources.update(list(functions.keys()))

        self.functions_initialized = True

        self.progress.end()


    def init_items_additional_data(self):
        """Recover additional data stored in the pan database."""

        items_additional_data = TableForItemAdditionalData(self.args)
        self.items_additional_data_keys, self.items_additional_data_dict = items_additional_data.get()

        # In fact we are done here since we have our `items_additional_data_dict` all filled up with sweet data.
        # But if functions are initialized, we can also get a summary of gene clusters based on whether most
        # genes in them were annotated with known functions or not for a given annotation source. Of course,
        # for this to happen, we need to check whther functions were initialied prior to the call to
        # `init_items_additional_data`.
        if not self.functions_initialized:
            # no? k. bye.
            self.progress.end()
            return

        self.progress.new('Recovering functions')
        # too many shitty nested loops here, but it is quite efficient since we work only with a dict
        # in memory
        for annotation_source in self.gene_clusters_function_sources:
            if annotation_source == 'COG_CATEGORY':
                # we don't need this one
                continue

            self.progress.update('Computing known/unknown dict for %s' % annotation_source)
            for gene_cluster_id in self.items_additional_data_dict:
                hits = Counter({})
                for genome_id in self.gene_clusters_functions_dict[gene_cluster_id]:
                    for gene_callers_id in self.gene_clusters_functions_dict[gene_cluster_id][genome_id]:
                        if annotation_source in self.gene_clusters_functions_dict[gene_cluster_id][genome_id][gene_callers_id]:
                            hits[self.gene_clusters_functions_dict[gene_cluster_id][genome_id][gene_callers_id][annotation_source][0]] += 1
                        else:
                            hits['UNKNOWN'] += 1

                if not hits or hits.most_common()[0][0] == 'UNKNOWN':
                    self.items_additional_data_dict[gene_cluster_id][annotation_source] = 'UNKNOWN'
                else:
                    self.items_additional_data_dict[gene_cluster_id][annotation_source] = 'KNOWN'

            self.items_additional_data_keys.append(annotation_source)

        self.progress.end()


    def init_gene_clusters(self):
        """Initializes the gene_clusters dictionary.

           At the end, the structure of this dictionary looks like this:

               {
                'gene_cluster_1': {'Genome_1': [gene_1, gene_2, (...)],
                                   'Genome_2': [],
                                   'Genome_3': [gene_1, gene_2],
                                   (...)}
                'gene_cluster_2': {(...)},
                (...)
               }

          This function also initializes alignment summaries for each gene
          in each gene cluster. That information is stored in the dict
          `self.gene_clusters_gene_alignments`.
        """

        self.progress.new('Initializing gene clusters')
        self.progress.update('...')

        pan_db = PanDatabase(self.pan_db_path)

        gene_clusters_long_list = pan_db.db.get_table_as_dict(t.pan_gene_clusters_table_name)

        for entry in list(gene_clusters_long_list.values()):
            genome_name = entry['genome_name']
            gene_callers_id = entry['gene_caller_id']
            gene_cluster_id = entry['gene_cluster_id']

            if gene_cluster_id not in self.gene_clusters:
                self.gene_clusters[gene_cluster_id] = {}
                for g in self.genome_names:
                    self.gene_clusters[gene_cluster_id][g] = []

            if self.gene_clusters_gene_alignments_available:
                if genome_name not in self.gene_clusters_gene_alignments:
                    self.gene_clusters_gene_alignments[genome_name] = {}

                self.gene_clusters_gene_alignments[genome_name][gene_callers_id] = entry['alignment_summary']

            self.gene_clusters[gene_cluster_id][genome_name].append(gene_callers_id)

        pan_db.disconnect()
        self.progress.end()

        self.gene_clusters_initialized = True


    def load_pan_views(self, splits_of_interest=None):
        pan_db = PanDatabase(self.pan_db_path)

        views_table = pan_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            self.views[view] = {'table_name': table_name,
                                'header': pan_db.db.get_table_structure(table_name)[1:],
                                'dict': pan_db.db.get_table_as_dict(table_name, keys_of_interest=splits_of_interest)}

        pan_db.disconnect()


    def get_summary_for_gene_clusters_list(self, gene_cluster_ids):
        summary = {'genomes_contributing': set([]), 'num_gene_calls': 0, 'num_gene_clusters': 0, 'functions': {}}

        if self.functions_initialized:
            for source in self.gene_clusters_function_sources:
                summary['functions'].update({source: Counter({})})

        for gene_cluster_id in gene_cluster_ids:
            single_summary = self.get_summary_for_gene_cluster_id(gene_cluster_id)
            summary['num_gene_clusters'] += 1
            summary['genomes_contributing'] = summary['genomes_contributing'].union(single_summary['genomes_contributing'])
            summary['num_gene_calls'] += single_summary['num_gene_calls']

            if self.functions_initialized:
                for source in self.gene_clusters_function_sources:
                    for function in single_summary['functions'][source]:
                        summary['functions'][source][function] += single_summary['functions'][source][function]

        summary['genomes_contributing'] = sorted(list(summary['genomes_contributing']))

        return summary


    def get_summary_for_gene_cluster_id(self, gene_cluster_id):
        summary = {'genomes_contributing': set([]), 'num_gene_calls': 0, 'functions': {}}

        for genome_name in self.gene_clusters[gene_cluster_id]:
            num_gene_calls = len(self.gene_clusters[gene_cluster_id][genome_name])
            if num_gene_calls:
                summary['genomes_contributing'].add(genome_name)
                summary['num_gene_calls'] += num_gene_calls

        if self.functions_initialized:
            for source in self.gene_clusters_function_sources:
                summary['functions'].update({source: Counter({})})

            functions_dict = self.gene_clusters_functions_dict[gene_cluster_id]
            for genome_name in functions_dict:
                for gene_callers_id in functions_dict[genome_name]:
                    for source in functions_dict[genome_name][gene_callers_id]:
                        for function in functions_dict[genome_name][gene_callers_id][source]:
                            summary['functions'][source][function] += 1

        return summary


    def init_collection_profile(self, collection_name):
        pan_db = PanDatabase(self.pan_db_path)

        if not self.gene_clusters:
            raise ConfigError("init_collection_profile wants to initialize the collection profile for '%s', but the\
                                the gene clusters dict is kinda empty :/ Someone forgot to initialize something maybe?" \
                                        % collection_name)

        # get trimmed collection and bins_info dictionaries
        collection, bins_info, self.gene_clusters_in_pan_db_but_not_binned \
                    = self.collections.get_trimmed_dicts(collection_name, set(self.gene_clusters.keys()))

        # currently we are not really doing anything with this one, but we will be filling this up with
        # all sorts of amazing later.
        for bin_id in collection:
            self.collection_profile[bin_id] = {}

        self.progress.end()
        pan_db.disconnect()

        return collection, bins_info


class ProfileSuperclass(object):
    """Fancy super class to deal with profile db stuff.

       if you want to make use of this class directly (i.e., not as a superclass), get an instance
       like this:

            >>> import anvio.dbops as d
            >>> import argparse
            >>> args = argparse.Namespace(profile_db="/path/to/profile.db")
            >>> p = ProfileSuperclass(args)

       Alternatively, you can include a contigs database path (contigs_db) in args so you have access
       to some functions that would require that.

       Alternatively, you can define a set of split names of interest:

            >>> args.split_names_of_interest = set([split_names])
            >>> p = ProfileSuperclass(args)

       in which case some functions will initialize data only for those splits. This is one way to minimize
       the resources necessary to initialize gene_coverages if only a subset of bins in a collection is
       requested.
       """

    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        # this one is a large dictionary with coverage values for every nucletoide positon in every sample for
        # every split and initialized by the member function `init_split_coverage_values_per_nt_dict` --unless the
        # member funciton `init_gene_level_coverage_stats_dicts` is not called first, in which case it is
        # automatically initialized from within that function.
        self.split_coverage_values_per_nt_dict = None

        # these are initialized by the member function `init_gene_level_coverage_stats_dicts`. but you knew
        # that already becasue you are a smart ass.
        self.gene_level_coverage_stats_dict = {}

        # this one becomes the object that gives access to the auxiliary data ops for split coverages
        # used heavily in interactive interface to show stuff (see bottle routes and all).
        self.split_coverage_values = None

        # the following two are initialized via `init_items_additional_data()` and use information
        # stored in item additional data tables
        self.items_additional_data_dict = None
        self.items_additional_data_keys = None

        self.auxiliary_profile_data_available = None
        self.auxiliary_data_path = None

        self.split_names = set([])
        self.item_orders = {}
        self.views = {}
        self.collection_profile = {}

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        self.split_names_of_interest = A('split_names_of_interest')
        init_gene_coverages = A('init_gene_coverages')
        populate_nt_level_coverage = A('populate_nt_level_coverage')
        outliers_threshold = A('outliers_threshold')

        if self.split_names_of_interest and not isinstance(self.split_names_of_interest, type(set([]))):
            raise ConfigError("ProfileSuper says the argument `splits_of_interest` must be of type set().\
                               Someone screwed up somewhere :/")

        if not self.profile_db_path:
            self.run.warning("ProfileSuper is called with args without member profile_db. Anvi'o will assume\
                              you are a programmer, and will not raise an error. But the init function is returning\
                              prematurely. Just so you know.")
            return

        is_profile_db(self.profile_db_path)

        # we have a contigs db? let's see if it's for real.
        if self.contigs_db_path:
            is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.progress.new('Initializing the profile database superclass')

        self.progress.update('Loading split names')
        self.split_names = get_split_names_in_profile_db(self.profile_db_path)

        if self.split_names == self.split_names_of_interest:
            # the user is being silly. nick that split_names_of_interest
            self.split_names_of_interest = None

        split_names_missing = (self.split_names_of_interest - self.split_names) if self.split_names_of_interest else None
        if self.split_names_of_interest and len(split_names_missing):
            self.progress.end()
            raise ConfigError("You called ProfileSuper with a `split_names_of_interest` argument, yet it contained\
                               %d split names that does not occur in the profile database. Here is an example: '%s'." % \
                                                                (len(split_names_missing), split_names_missing.pop()))

        self.progress.update('Creating an instance of the profile database')
        profile_db = ProfileDatabase(self.profile_db_path)

        self.progress.update('Setting profile self data dict')
        self.p_meta = profile_db.meta

        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if 'creation_date' in self.p_meta else 'unknown'
        self.p_meta['samples'] = sorted([s.strip() for s in self.p_meta['samples'].split(',')])
        self.p_meta['num_samples'] = len(self.p_meta['samples'])

        if self.p_meta['blank'] and not self.p_meta['contigs_db_hash']:
            self.progress.end()
            raise ConfigError("ProfileSuperclass is upset, because it seems you are tyring to initialize a blank anvi'o profile\
                               database that is not associated with a contigs database. This will not work for multiple reasons.\
                               The current technical limitation is that blank profile databases that are in this situation do not\
                               keep track of split names they are working with. Yes. We too know that this is a serious design\
                               flaw, but THANKS for reminding anyway... The best way to address this is to make sure all anvi'o\
                               profile and pan databases maintain a table with all item names they are supposed to be working with.")

        if self.p_meta['contigs_ordered'] and 'available_item_orders' in self.p_meta:
            self.p_meta['available_item_orders'] = sorted([s.strip() for s in self.p_meta['available_item_orders'].split(',')])
            self.item_orders = profile_db.db.get_table_as_dict(t.item_orders_table_name)

            for item_order in self.item_orders:
                if self.item_orders[item_order]['type'] == 'basic':
                    try:
                        self.item_orders[item_order]['data'] = self.item_orders[item_order]['data'].split(',')
                    except:
                        self.progress.end()
                        raise ConfigError("Something is wrong with the basic order `%s` in this profile database :(" % (item_order))

        elif self.p_meta['contigs_ordered'] and 'available_item_orders' not in self.p_meta:
            self.progress.end()
            self.run.warning("Your profile database thinks the hierarchical item_order was done, yet it contains no entries\
                              for any hierarchical item_order results. This is not good. Something must have gone wrong\
                              somewhere :/ To be on the safe side, anvi'o will assume this profile database has no\
                              item_order (which is literally the case, by the way, it is just the database itself is\
                              confused about that fact --it happens to the best of us).")
            self.progress.new('Initializing the profile database superclass')

            self.p_meta['contigs_ordered'] = False
            self.p_meta['available_item_orders'] = None
            self.p_meta['default_item_order'] = None
            self.item_orders = None
        else:
            self.p_meta['available_item_orders'] = None
            self.p_meta['default_item_order'] = None
            self.item_orders = None

        profile_db.disconnect()

        self.progress.update('Accessing the auxiliary data file')
        self.auxiliary_data_path = get_auxiliary_data_path_for_profile_db(self.profile_db_path)
        if not os.path.exists(self.auxiliary_data_path):
            self.auxiliary_profile_data_available = False
        else:
            self.auxiliary_profile_data_available = True
            self.split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.auxiliary_data_path,
                                                                                           self.p_meta['contigs_db_hash'])

        self.progress.end()

        if init_gene_coverages:
            self.init_gene_level_coverage_stats_dicts(outliers_threshold=outliers_threshold, populate_nt_level_coverage=populate_nt_level_coverage)

        if self.auxiliary_profile_data_available:
            self.run.info('Auxiliary Data', 'Found: %s (v. %s)' % (self.auxiliary_data_path, anvio.__auxiliary_data_version__))

        if self.split_names_of_interest:
            self.run.info('Profile Super', 'Initialized with %d of %d splits: %s (v. %s)' % (len(self.split_names),
                                                                                             len(self.split_names_of_interest),
                                                                                             self.profile_db_path,
                                                                                             anvio.__profile__version__))
        else:
            self.run.info('Profile Super', 'Initialized with all %d splits: %s (v. %s)' % (len(self.split_names),
                                                                                           self.profile_db_path,
                                                                                           anvio.__profile__version__))


    def init_split_coverage_values_per_nt_dict(self):
        """This function will fill process the auxiliary data and fill this dictionary:

            - self.split_coverage_values_per_nt_dict

           If this is taking forever and you want to kill Meren, everyone will understand you.
        """

        if self.p_meta['blank']:
            self.run.warning("Someone asked gene coverages to be initialized when working with a blank profile database.\
                              Anvi'o will pretend nothing happened, and will return nothing. If you don't know what this\
                              is warning you about, just carry on.")
            return

        if not self.auxiliary_profile_data_available:
            self.run.warning("Gene-level detection and coverage values are always recovered from the auxiliary data files\
                              associated with profile databases. You don't seem to have one around for this profile database,\
                              and you shall get NO GENE COVERAGES OR ANYTHING :(")
            return

        if not self.contigs_db_path:
            self.run.warning("Someone wants to populate gene coverages data, but they called the profile super class without\
                              a contigs database path. Anvi'o will pretend nothing happened, but will return nothing back.\
                              Good luck with your downstream endeavors.")
            return

        self.progress.new('Initializing split coverage values per nt')
        self.progress.update('...')

        if self.split_names_of_interest:
            self.split_coverage_values_per_nt_dict = self.split_coverage_values.get_coverage_for_multiple_splits(self.split_names_of_interest)
        else:
            self.split_coverage_values_per_nt_dict = self.split_coverage_values.get_all()

        self.progress.end()


    def init_gene_level_coverage_stats_dicts(self, min_cov_for_detection=0, outliers_threshold=1.5, populate_nt_level_coverage=False):
        """This function will process `self.split_coverage_values_per_nt_dict` to populate
           `self.gene_level_coverage_stats_dict`.

           Note: if a `split_names_of_interest` argument is declared at the class level,
           this function will operate on those splits found in that set.
           """

        if not self.auxiliary_profile_data_available:
            raise ConfigError("Someone is asking gene level coverage stats to be computed, but then there is no auxiliary profile\
                               data does not seem to be available for this project. Yeah. That's what happens if you don't\
                               download everything from the server :(")

        contigs_db = ContigsSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        if not contigs_db.a_meta['genes_are_called']:
            self.run.warning("Well, someone wants to populate the gene coverages data, when in fact genes were not called :/\
                              Instead of giving an error, anvi'o will return prematurely, without really doing anything.")
            return

        if not contigs_db.a_meta['splits_consider_gene_calls']:
            self.run.warning("PLEASE READ THIS VERY CAREFULLY (remember, anvi'o never talks to you in CAPS, so it must be important).\
                              It seems when you generated your contigs database, you have skipped 'mindful' splitting of contigs.\
                              This means, some of the genes may be soft-broken into two or more pieces. For most things, it doesn't\
                              really matter, but here this will cause an issue as your gene coverages will average one of those splits\
                              without any biologically relevant reason. We could have done much better here, but it would have affected\
                              the performance very negatively. If you are seeing this warning, and go like 'crap, this will ruin\
                              everything because I possibly can not recover from this situation', then send us an e-mail, and we will\
                              think about whether we can be less lazy about stuff, and do things better.")

        sample_names = self.p_meta['samples']

        if not self.split_coverage_values_per_nt_dict:
            self.init_split_coverage_values_per_nt_dict()

        if self.split_names_of_interest:
            split_names = self.split_names_of_interest

            self.run.warning('A subset of genes (%d of %d, to be precise) are requested to initiate gene-level coverage stats for.\
                              No need to worry, this is just a warning in case you are as obsessed as wanting to know everything\
                              there is to know.' % (len(self.split_names_of_interest), len(self.split_names)), overwrite_verbose=True)

        else:
            split_names = self.split_names

        self.progress.new('Computing gene-level coverage stats ...')
        self.progress.update('...')

        num_splits, counter = len(split_names), 1
        # go through all the split names
        for split_name in split_names:
            if num_splits > 100 and counter % 100 == 0:
                self.progress.update('%d of %d splits ...' % (counter, num_splits))

            # recover split coverage values from the auxiliary data file:
            split_coverage = self.split_coverage_values_per_nt_dict[split_name]

            # identify entry ids for genes in `split_name`
            genes_in_splits_entries = contigs_db.split_name_to_genes_in_splits_entry_ids[split_name]

            # we have to go back, Kate :(
            if not genes_in_splits_entries:
                continue

            # we will go through each gene entry in the split
            for genes_in_splits_entry in genes_in_splits_entries:
                e = contigs_db.genes_in_splits[genes_in_splits_entry]
                gene_callers_id, gene_start, gene_stop = e['gene_callers_id'], e['start_in_split'], e['stop_in_split']
                gene_length = gene_stop - gene_start

                if gene_length <= 0:
                    raise ConfigError("What? :( How! The gene with the caller id '%d' has a length of %d :/ We are done\
                                       here!" % (gene_callers_id, gene_length))

                self.gene_level_coverage_stats_dict[gene_callers_id] = dict([(sample_name, dict([('mean_coverage', 0), ('gene_detection', 0)])) for sample_name in sample_names])

                # the magic happens here:
                for sample_name in sample_names:
                    # and recover the gene coverage array per position for a given sample:
                    gene_coverage_per_position = split_coverage[sample_name][gene_start:gene_stop]

                    mean_coverage = numpy.mean(gene_coverage_per_position)
                    detection = numpy.count_nonzero(gene_coverage_per_position) / gene_length

                    # findout outlier psitions, and get non-outliers
                    outliers_bool = get_list_of_outliers(gene_coverage_per_position, outliers_threshold)
                    non_outlier_positions = numpy.invert(outliers_bool)
                    non_outliers = gene_coverage_per_position[non_outlier_positions]

                    if not(len(non_outliers)):
                        non_outlier_mean_coverage = 0.0
                        non_outlier_coverage_std = 0.0
                    else:
                        non_outlier_mean_coverage = numpy.mean(non_outliers)
                        non_outlier_coverage_std = numpy.std(non_outliers)

                    self.gene_level_coverage_stats_dict[gene_callers_id][sample_name] = {'mean_coverage': mean_coverage,
                                                                                          'detection': detection,
                                                                                          'non_outlier_mean_coverage': non_outlier_mean_coverage,
                                                                                          'non_outlier_coverage_std':  non_outlier_coverage_std}
                    if populate_nt_level_coverage == True:
                        self.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['gene_coverage_per_position'] = gene_coverage_per_position
                        self.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['non_outlier_positions'] = non_outlier_positions

            counter += 1

        self.progress.end()


    def get_variability_information_for_split(self, split_name, skip_outlier_SNVs=False, return_raw_results=False):
        if not split_name in self.split_names:
            raise ConfigError("get_variability_information_for_split: The split name '%s' does not seem to be\
                                represented in this profile database. Are you sure you are looking for it\
                                in the right database?" % split_name)

        self.progress.new('Recovering variabilit information for split')
        self.progress.update('...')

        profile_db = ProfileDatabase(self.profile_db_path)
        split_variability_information = list(profile_db.db.get_some_rows_from_table_as_dict(t.variable_nts_table_name, '''split_name = "%s"''' % split_name, error_if_no_data=False).values())
        profile_db.disconnect()

        if return_raw_results:
            return split_variability_information

        # they want pretty stuff...
        d = {}

        for sample_name in self.p_meta['samples']:
            d[sample_name] = {'variability': {0: {}, 1: {}, 2: {}, 3: {}}, 'competing_nucleotides': {}}

        for e in split_variability_information:
            frequencies = utils.get_variabile_item_frequencies(e, engine='NT')
            e['n2n1ratio'], e['consensus'], e['departure_from_consensus'] = utils.get_consensus_and_departure_data(frequencies)

            if skip_outlier_SNVs and e['cov_outlier_in_contig']:
                continue

            d[e['sample_id']]['variability'][e['base_pos_in_codon']][e['pos']] = e['departure_from_reference']
            d[e['sample_id']]['competing_nucleotides'][e['pos']] = e

        self.progress.end()

        return d


    def init_items_additional_data(self):
        items_additional_data = TableForItemAdditionalData(self.args)
        self.items_additional_data_keys, self.items_additional_data_dict = items_additional_data.get()


    def init_collection_profile(self, collection_name):
        profile_db = ProfileDatabase(self.profile_db_path, quiet=True)

        # we only have a self.collections instance if the profile super has been inherited by summary super class.
        # the initialization of a collection profile should only be done through that module anyway. so we are
        # being cruel here, and sending the programmer back.
        if not hasattr(self, 'collections'):
            raise ConfigError("You are lost :/ You can only call `init_collection_profile` through an instance of \
                               the `SummarizerSuperClass`. Go back and come another way.")

        # get trimmed collection and bins_info dictionaries
        collection, bins_info, self.split_names_in_profile_db_but_not_binned \
                    = self.collections.get_trimmed_dicts(collection_name, self.split_names)

        for bin_id in collection:
            self.collection_profile[bin_id] = {}

        table_names = [] if self.p_meta['blank'] else [table_name for table_name in t.atomic_data_table_structure[1:-1]]

        samples_template = dict([(s, []) for s in self.p_meta['samples']])

        # anonymous function to convert single profile table dicts compatible with merged ones (#155):
        SINGLE_P = lambda d: dict([(s, dict([(self.p_meta['samples'][0], v) for v in list(d[s].values())])) for s in d])

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
                    averages[sample_name] = numpy.mean([a or 0 for a in averages[sample_name]])

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
                all_coverages_in_sample = sum([d[sample] for d in list(coverage_table_data.values())])

                for bin_id in collection:
                    bin_coverages_in_sample = sum([coverage_table_data[split_name][sample] for split_name in collection[bin_id]])

                    if all_coverages_in_sample == 0:
                        percents[bin_id] = 0.0
                    else:
                        percents[bin_id] = bin_coverages_in_sample * 100 / all_coverages_in_sample

                splits_not_binned_coverages_in_sample = sum([coverage_table_data[split_name][sample] for split_name in self.split_names_in_profile_db_but_not_binned])

                if all_coverages_in_sample == 0:
                    percents['__splits_not_binned__'] = 0.0
                else:
                    percents['__splits_not_binned__'] = splits_not_binned_coverages_in_sample * 100 / all_coverages_in_sample

                self.bin_percent_recruitment_per_sample[sample] = percents

        self.progress.end()
        profile_db.disconnect()

        return collection, bins_info


    def load_views(self, splits_of_interest=None, omit_parent_column=False):
        profile_db = ProfileDatabase(self.profile_db_path)

        views_table = profile_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            self.views[view] = {'table_name': table_name,
                                'header': profile_db.db.get_table_structure(table_name)[1:],
                                'dict': profile_db.db.get_table_as_dict(table_name, keys_of_interest=splits_of_interest, omit_parent_column=omit_parent_column)}

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

            for key in ['min_contig_length', 'SNVs_profiled', 'AA_frequencies_profiled', 'min_coverage_for_variability', 'merged', 'blank', 'contigs_ordered', 'report_variability_full', 'num_contigs', 'num_splits', 'total_length']:
                try:
                    self.meta[key] = int(self.meta[key])
                except:
                    pass

            sample_ids_list = [s.strip() for s in self.meta['samples'].split(',')]
            if 'total_reads_mapped' in self.meta:
                total_reads_mapped_list = [int(n.strip()) for n in self.meta['total_reads_mapped'].split(',')]
                self.meta['total_reads_mapped'] = dict([(sample_ids_list[i], total_reads_mapped_list[i]) for i in range(0, len(sample_ids_list))])
            else:
                self.meta['total_reads_mapped'] = dict([(sample_ids_list[i], 0) for i in range(0, len(sample_ids_list))])

            self.samples = set([s.strip() for s in self.meta['samples'].split(',')])

            self.run.info('Profile database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
            self.run.info('Samples', self.meta['samples'], quiet=self.quiet)
        else:
            self.db = None


    def touch(self):
        """Creates an empty profile database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db."""

        is_db_ok_to_create(self.db_path, 'profile')

        self.db = db.DB(self.db_path, anvio.__profile__version__, new_database=True)

        # creating empty default tables
        self.db.create_table(t.item_additional_data_table_name, t.item_additional_data_table_structure, t.item_additional_data_table_types)
        self.db.create_table(t.item_orders_table_name, t.item_orders_table_structure, t.item_orders_table_types)
        self.db.create_table(t.variable_nts_table_name, t.variable_nts_table_structure, t.variable_nts_table_types)
        self.db.create_table(t.variable_aas_table_name, t.variable_aas_table_structure, t.variable_aas_table_types)
        self.db.create_table(t.views_table_name, t.views_table_structure, t.views_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)

        return self.db


    def create(self, meta_values={}):
        self.touch()

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())

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

            for key in ['num_genomes', 'gene_cluster_min_occurrence', 'use_ncbi_blast', 'diamond_sensitive', 'exclude_partial_gene_calls', \
                        'num_gene_clusters', 'num_genes_in_gene_clusters', 'gene_alignments_computed', 'gene_clusters_ordered']:
                try:
                    self.meta[key] = int(self.meta[key])
                except:
                    pass

            for key in ['min_percent_identity', 'minbit', 'mcl_inflation']:
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
        self.db.create_table(t.pan_gene_clusters_table_name, t.pan_gene_clusters_table_structure, t.pan_gene_clusters_table_types)

        # creating empty default tables for standard anvi'o pan dbs
        self.db.create_table(t.item_additional_data_table_name, t.item_additional_data_table_structure, t.item_additional_data_table_types)
        self.db.create_table(t.item_orders_table_name, t.item_orders_table_structure, t.item_orders_table_types)
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

            for key in ['split_length', 'kmer_size', 'total_length', 'num_splits', 'num_contigs', 'genes_are_called', 'splits_consider_gene_calls']:
                self.meta[key] = int(self.meta[key])

            self.meta['gene_function_sources'] = [s.strip() for s in self.meta['gene_function_sources'].split(',')] if self.meta['gene_function_sources'] else None

            if 'creation_date' not in self.meta:
                raise ConfigError("The contigs database ('%s') seems to be corrupted :/ This happens if the process that\
                                    that generates the database ends prematurely. Most probably, you will need to generate\
                                    the contigs database from scratch. Sorry!" % (self.db_path))

            self.run.info('Contigs database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
            self.run.info('Number of contigs', self.meta['num_contigs'], quiet=self.quiet)
            self.run.info('Number of splits', self.meta['num_splits'], quiet=self.quiet)
            self.run.info('Total number of nucleotides', self.meta['total_length'], quiet=self.quiet)
            self.run.info('Split length', self.meta['split_length'], quiet=self.quiet)
        else:
            self.db = None


    def get_date(self):
        return time.time()


    def get_hash(self):
        return '%08x' % random.randrange(16**8)


    def touch(self):
        """Creates an empty contigs database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db."""

        is_db_ok_to_create(self.db_path, 'contigs')

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
        self.db.create_table(t.gene_amino_acid_sequences_table_name, t.gene_amino_acid_sequences_table_structure, t.gene_amino_acid_sequences_table_types)
        self.db.create_table(t.genes_in_splits_summary_table_name, t.genes_in_splits_summary_table_structure, t.genes_in_splits_summary_table_types)
        self.db.create_table(t.splits_info_table_name, t.splits_info_table_structure, t.splits_info_table_types)
        self.db.create_table(t.contigs_info_table_name, t.contigs_info_table_structure, t.contigs_info_table_types)
        self.db.create_table(t.nt_position_info_table_name, t.nt_position_info_table_structure, t.nt_position_info_table_types)

        return self.db


    def create(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        contigs_fasta = A('contigs_fasta')
        project_name = A('project_name')
        description_file_path = A('description')
        split_length = A('split_length')
        kmer_size = A('kmer_size')
        skip_gene_calling = A('skip_gene_calling')
        external_gene_calls = A('external_gene_calls')
        skip_mindful_splitting = A('skip_mindful_splitting')
        ignore_internal_stop_codons = A('ignore_internal_stop_codons')
        debug = A('debug')

        if external_gene_calls:
            filesnpaths.is_file_exists(external_gene_calls)

        if external_gene_calls and skip_gene_calling:
            raise ConfigError("You provided a file for external gene calls, and used requested gene calling to be\
                                skipped. Please make up your mind.")

        filesnpaths.is_file_fasta_formatted(contigs_fasta)
        contigs_fasta = os.path.abspath(contigs_fasta)

        # let the user see what's up
        self.run.info('Input FASTA file', contigs_fasta)

        if not project_name:
            project_name = '.'.join(os.path.basename(os.path.abspath(contigs_fasta)).split('.')[:-1])

            if project_name:
                self.run.warning("You are generating a new anvi'o contigs database, but you are not specifying a\
                                  project name for it. FINE. Anvi'o, in desperation, will use the input file name\
                                  to set the project name for this contigs database (which is '%s'). If you are not\
                                  happy with that, feel free to kill and restart this process. If you are not happy\
                                  with this name, but you don't like killing things either, maybe next time you\
                                  should either name your FASTA files better, or use the `--project-name` parameter\
                                  to set your desired name." % project_name, "Anvi'o made things up for you")
            else:
                raise ConfigError("Sorry, you must provide a project name for your contigs database :/ Anvi'o tried\
                                   to make up one, but failed.")

        self.run.info('Name', project_name, mc='green')
        self.run.info('Description', os.path.abspath(description_file_path) if description_file_path else 'No description is given', mc='green')

        if description_file_path:
            filesnpaths.is_file_plain_text(description_file_path)
            description = open(os.path.abspath(description_file_path), 'rU').read()
        else:
            description = ''

        # go throught he FASTA file to make sure there are no surprises with deflines and sequence lengths.
        self.progress.new('Checking deflines and contig lengths')
        self.progress.update('tick tock ...')
        fasta = u.SequenceSource(contigs_fasta)
        while next(fasta):
            if not utils.check_contig_names(fasta.id, dont_raise=True):
                self.progress.end()
                raise ConfigError("At least one of the deflines in your FASTA File does not comply with the 'simple deflines'\
                                    requirement of anvi'o. You can either use the script `anvi-script-reformat-fasta` to take\
                                    care of this issue, or read this section in the tutorial to understand the reason behind\
                                    this requirement (anvi'o is very upset for making you do this): %s" % \
                                        ('http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file'))

            if len(fasta.seq) < kmer_size:
                self.progress.end()
                raise ConfigError("At least one of the contigs in your input FASTA '%s' is shorter than the k-mer size. The k\
                                    is %d, and your contig is like %d :/ Anvi'o will not judge you for whatever you are doing\
                                    with such short contigs, but the length of each contig must be at least as long as your `k` for\
                                    k-mer analyis. You can use the script `anvi-script-reformat-fasta` to get rid of very short\
                                    contigs if you like." % (contigs_fasta, kmer_size, len(fasta.seq)))
        fasta.close()
        self.progress.end()

        all_ids_in_FASTA = utils.get_all_ids_from_fasta(contigs_fasta)
        if len(all_ids_in_FASTA) != len(set(all_ids_in_FASTA)):
            raise ConfigError("Every contig in the input FASTA file must have a unique ID. You know...")

        if not split_length:
            raise ConfigError("Creating a new contigs database requires split length information to be\
                                provided. But the ContigsDatabase class was called to create one without this\
                                bit of information. Not cool.")

        if not os.path.exists(contigs_fasta):
            raise ConfigError("Creating a new contigs database requires a FASTA file with contigs to be provided.")

        try:
            split_length = int(split_length)
        except:
            raise ConfigError("Split size must be an integer.")

        if split_length <= 0:
            split_length = sys.maxsize

        try:
            kmer_size = int(kmer_size)
        except:
            raise ConfigError("K-mer size must be an integer.")
        if kmer_size < 2 or kmer_size > 8:
            raise ConfigError("We like our k-mer sizes between 2 and 8, sorry! (but then you can always change the\
                                source code if you are not happy to be told what you can't do, let us know how it goes!).")

        if skip_gene_calling:
            skip_mindful_splitting = True

        # create a blank contigs database on disk, and set the self.db
        self.touch()

        # know thyself
        self.db.set_meta_value('db_type', 'contigs')
        self.db.set_meta_value('project_name', project_name)
        self.db.set_meta_value('description', description)

        # this will be the unique information that will be passed downstream whenever this db is used:
        contigs_db_hash = self.get_hash()
        self.db.set_meta_value('contigs_db_hash', contigs_db_hash)

        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)

        self.run.info('Split Length', pp(split_length))
        self.run.info('K-mer size', kmer_size)
        self.run.info('Skip gene calling?', skip_gene_calling)
        self.run.info('External gene calls provided?', external_gene_calls)
        self.run.info('Ignoring internal stop codons?', ignore_internal_stop_codons)
        self.run.info('Splitting pays attention to gene calls?', (not skip_mindful_splitting))

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
                gene_calls_tables.use_external_gene_calls_to_populate_genes_in_contigs_table(input_file_path=external_gene_calls, ignore_internal_stop_codons=ignore_internal_stop_codons)
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
        nt_positions_table = TableForNtPositions()
        contigs_info_table = InfoTableForContigs(split_length)
        splits_info_table = InfoTableForSplits()

        recovered_split_lengths = []

        # THE INFAMOUS GEN CONTGS DB LOOP (because it is so costly, we call it South Loop)
        self.progress.new('The South Loop')
        fasta.reset()
        while next(fasta):
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
                nt_positions_table.append(contig_name, nt_position_info_list)

            contig_kmer_freq = contigs_kmer_table.get_kmer_freq(contig_sequence)

            self.progress.append('k-mers ...')
            for order in range(0, len(split_start_stops)):
                start, end = split_start_stops[order]
                split_name = contigops.gen_split_name(contig_name, order)

                # this is very confusing, because both contigs_kmer_table and splits_kmer_able in fact
                # holds kmer values for splits only. in one table, each split has a kmer value of their
                # contigs (to not lose the genomic context while item_order based on kmers), in the other
                # one each split holds its own kmer value.
                contigs_kmer_table.append(split_name, contig_sequence[start:end], kmer_freq=contig_kmer_freq)
                splits_kmer_table.append(split_name, contig_sequence[start:end])

                splits_info_table.append(split_name, contig_sequence[start:end], order, start, end, contig_gc_content, contig_name)

            db_entries_contig_sequences.append((contig_name, contig_sequence), )

        self.progress.end()

        self.db.set_meta_value('kmer_size', kmer_size)
        nt_positions_table.store(self.db)
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
        self.db.set_meta_value('creation_date', self.get_date())
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
            raise ConfigError("When SamplesInformationDatabase is called, the db_path parameter cannot be\
                                'None' type :/")

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
            raise ConfigError("The samples database has not been initialized. You are doing something wrong :/")

        samples = samplesops.SamplesInformation(run=self.run, progress=self.progress, quiet=self.quiet)

        samples_information_dict = samples.recover_samples_information_dict(self.db.get_table_as_dict(t.samples_information_table_name, error_if_no_data=False),
                                                                            self.db.get_table_as_dict(t.samples_attribute_aliases_table_name, error_if_no_data=False))
        samples_order_dict = self.db.get_table_as_dict(t.samples_order_table_name)

        return samples_information_dict, samples_order_dict


    def get_samples_information_default_layer_order(self):
        if not self.db:
            raise ConfigError("The samples database has not been initialized. You are doing something wrong :/")

        return self.samples_information_default_layer_order


    def create(self, samples_information_path=None, samples_order_path=None, single_order_path=None, single_order_name=None):
        is_db_ok_to_create(self.db_path, 'samples')

        samples = samplesops.SamplesInformation(run=self.run, progress=self.progress, quiet=self.quiet)
        samples.populate_from_input_files(samples_information_path, samples_order_path, single_order_path, single_order_name)

        self.db = db.DB(self.db_path, anvio.__samples__version__, new_database=True)

        self.write_samples_to_database(samples)

        self.run.info('Samples information database', 'A new samples information database, %s, has been created.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of samples', len(samples.sample_names), quiet=self.quiet)
        self.run.info('Number of organizations', len(list(samples.samples_order_dict.keys())), quiet=self.quiet)


    def export_samples_db_files(self, order_output_path='samples-order.txt', information_output_path='samples-information.txt', output_file_prefix=None):
        """Export whatever information is stored in a ginve anvi'o samples database"""

        samples_information_dict, samples_order_dict = self.get_samples_information_and_order_dicts()

        if output_file_prefix:
            order_output_path = output_file_prefix + '-' + order_output_path
            information_output_path = output_file_prefix + '-' + information_output_path

        filesnpaths.is_output_file_writable(order_output_path)
        filesnpaths.is_output_file_writable(information_output_path)

        utils.store_dict_as_TAB_delimited_file(samples_order_dict, order_output_path, headers=['attributes', 'basic', 'newick'])
        utils.store_dict_as_TAB_delimited_file(samples_information_dict, information_output_path, headers=['samples'] + sorted(list(list(samples_information_dict.values())[0].keys())))

        self.run.info('Samples information file', information_output_path, mc='green')
        self.run.info('Samples order file', order_output_path, mc='green')


    def update(self, samples_information_path=None, samples_order_path=None, single_order_path=None, single_order_name=None):
        # first recover what is already in the database
        samples_information_dict, samples_order_dict = self.get_samples_information_and_order_dicts()

        # inherit a samples object and update its member dicts:
        samples = samplesops.SamplesInformation(run=self.run, progress=self.progress, quiet=self.quiet)
        samples.samples_order_dict = samples_order_dict
        samples.samples_information_dict = samples_information_dict

        # add what we have now
        samples.populate_from_input_files(samples_information_path, samples_order_path, single_order_path, single_order_name)

        self.db = db.DB(self.db_path, anvio.__samples__version__, new_database=False)

        self.write_samples_to_database(samples, update=True)

        self.run.info('Samples information database', 'Samples information database, %s, has been updated.' % (self.db_path), quiet=self.quiet)
        self.run.info('Number of samples', len(samples.sample_names), quiet=self.quiet)
        self.run.info('Number of organizations', len(list(samples.samples_order_dict.keys())), quiet=self.quiet)


    def write_samples_to_database(self, samples, update=False):
        if update:
            self.db.drop_table(t.samples_order_table_name)
            self.db.drop_table(t.samples_attribute_aliases_table_name)
            self.db.drop_table(t.samples_information_table_name)

            self.db.remove_meta_key_value_pair('samples')
            self.db.remove_meta_key_value_pair('available_orders')
            self.db.remove_meta_key_value_pair('sample_names_for_order')
            self.db.remove_meta_key_value_pair('samples_information_default_layer_order')
        else:
            # know thyself
            self.db.set_meta_value('db_type', 'samples_information')

            # set some useful meta values:
            self.db.set_meta_value('creation_date', time.time())


        # first create the easy one: the samples_order table.
        available_orders = list(samples.samples_order_dict.keys())
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


class TableForItemAdditionalData(Table):
    """Implements the item additional data class.

       This is the class where we maintain the item additional data tables in anvi'o
       pan and profile databases. Related issue: https://github.com/merenlab/anvio/issues/662.
    """

    def __init__(self, args, r=run, p=progress, table_name=t.item_additional_data_table_name):
        self.run = r
        self.progress = p
        self.table_name = table_name

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.db_path = A('pan_or_profile_db') or A('profile_db') or A('pan_db')
        self.just_do_it = A('just_do_it')

        if not self.db_path:
            raise ConfigError("ItemAdditionalData class is inherited with args object that did not\
                               contain any database path :/ Even though any of the following would\
                               have worked: `pan_or_profile_db`, `profile_db`, `pan_db` :(")

        is_pan_or_profile_db(self.db_path)
        self.db_type = get_db_type(self.db_path)
        self.db_version = get_required_version_for_db(self.db_path)

        database = db.DB(self.db_path, self.db_version)
        self.item_additional_data_keys = database.get_single_column_from_table(self.table_name, 'key')
        database.disconnect()

        Table.__init__(self, self.db_path, self.db_version, self.run, self.progress)


    def export(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        keys, data = self.get()
        utils.store_dict_as_TAB_delimited_file(data, output_file_path, headers=['item'] + keys)

        self.run.info('Output file', output_file_path)


    def get(self):
        """Will return the additional data keys and the dict."""

        self.progress.new('Recovering item additional keys and data')
        self.progress.update('...')
        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        item_additional_data = database.get_table_as_dict(self.table_name)
        item_additional_data_keys = database.get_single_column_from_table(self.table_name, 'key', unique=True)
        item_names = database.get_single_column_from_table(self.table_name, 'item_name', unique=True)
        database.disconnect()

        if not len(item_names):
            self.progress.end()
            return [], {}

        d = {}
        for item_name in item_names:
            d[item_name] = {}

        for entry in item_additional_data.values():
            item_name = entry['item_name']
            key = entry['key']
            value = entry['value']

            if entry['type'] in ['int', 'float']:
                d[item_name][key] = eval(entry['type'])(value)
            else:
                d[item_name][key] = value

        for item_name in d:
            for key in item_additional_data_keys:
                if key not in d[item_name]:
                    d[item_name][key] = None

        self.progress.end()

        return item_additional_data_keys, d


    def remove(self, keys_list):
        '''Give this guy a list of key, and watch their demise.'''

        if not isinstance(keys_list, list):
            raise ConfigError("The remove function in ItemAdditionalData class wants you to watch\
                               yourself before you wreck yourself. In other words, can you please\
                               make sure the keys you send is of type `list` thankyouverymuch.")

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))

        item_additional_data_keys = sorted(database.get_single_column_from_table(self.table_name, 'key', unique=True))

        if not len(item_additional_data_keys):
            self.run.info_single('There is nothing to remove --the items additional data table is already empty :(')
            database.disconnect()

            return

        missing_keys = [k for k in keys_list if k not in item_additional_data_keys]
        if len(missing_keys) and not self.just_do_it:
            database.disconnect()
            raise ConfigError("The following keys you wanted to remove from the items additional data table are\
                               not really in the table: '%s'. Anvi'o is confused :/" % (', '.join(missing_keys)))

        if keys_list:
            for key in keys_list:
                if key not in item_additional_data_keys:
                    # what the hell, user?
                    return

                database._exec('''DELETE from %s WHERE key="%s"''' % (self.table_name, key))

            self.run.warning("Data for the following keys removed from the database: '%s'. #SAD." % (', '.join(keys_list)))
        else:
            database._exec('''DELETE from %s''' % (self.table_name))

            self.run.warning("All data from the items additional data table is removed (wow).")

        database.disconnect()


    def check_item_names(self, data_dict):
        """Compares item names found in the data dict to the ones in the db"""

        items_in_db = get_all_item_names_from_the_database(self.db_path)
        items_in_data = set(data_dict.keys())

        items_in_data_but_not_in_db = items_in_data.difference(items_in_db)
        if len(items_in_data_but_not_in_db):
            raise ConfigError("Well. %d of %d item names in your additional data are only in your data (which\
                               that they are not in the %s database you are working with (which means bad news)).\
                               Since there is no reason to add additional data for items that do not exist in your\
                               database, anvi'o will stop you right there. Please fix your data and come again. In\
                               case you want to see a random item that is only in your data, here is one: %s. Stuff\
                               in your db looks like this: %s." \
                                    % (len(items_in_data_but_not_in_db), len(items_in_data), self.db_type, \
                                       items_in_data_but_not_in_db.pop(), items_in_db.pop()))

        items_in_db_but_not_in_data = items_in_db.difference(items_in_data)
        if len(items_in_db_but_not_in_data):
            self.run.warning("Your input contains additional data for only %d of %d total number of items in your %s\
                              database. Just wanted to make sure you know what's up, but we cool." \
                                % (len(items_in_db) - len(items_in_db_but_not_in_data), len(items_in_db), self.db_type))


    def add(self, keys_list, data_dict, skip_check_names=False):
        """Main function to add data into the item additional data table.

           * `data_dict`: a dictionary that should follow this format:

                d = {
                        'item_name_01': {'key_01': value,
                                         'key_02': value,
                                         'key_03': value
                                         },
                        'item_name_02': {'key_01': value,
                                         'key_03': value,
                                         },
                        (...)
                    }

           * `keys_list`: is a list of keys one or more of which should appear for each item
                          in `data_dict`.
        """

        if not isinstance(keys_list, list):
            raise ConfigError("List of keys must be of type `list`. Go away.")

        if not isinstance(data_dict, dict):
            raise ConfigError("Nope. Your data must be of type `dict`.")

        self.run.warning(None, 'New additional data...', lc="yellow")
        key_types = {}
        for key in keys_list:
            if '!' in key:
                predicted_key_type = "stackedbar"
            else:
                type_class = utils.get_predicted_type_of_items_in_a_dict(data_dict, key)
                predicted_key_type = type_class.__name__ if type_class else None

            key_types[key] = predicted_key_type
            self.run.info('Key "%s"' % key, 'Predicted type: %s' % (key_types[key]), \
                                            nl_after = 1 if key == keys_list[-1] else 0)

        # we be responsible here.
        keys_already_in_db = [c for c in keys_list if c in self.item_additional_data_keys]
        if len(keys_already_in_db):
            if self.just_do_it:
                self.run.warning('The following keys in your data dict will replace the ones that are already\
                                  in your %s database: %s.' % (self.db_type, ', '.join(keys_already_in_db)))

                self.remove(keys_already_in_db)
            else:
                run.info('Keys already in the db', ', '.join(keys_already_in_db), nl_before=2, mc='red')

                raise ConfigError("Some of the keys in your new data appear to be in the database already. If you\
                                   want to replace those in the database with the ones in your new data use the\
                                   `--just-do-it` flag, and watch anvi'o make an exception just for you and complain\
                                   about nothin' for this once.")

        if skip_check_names:
            self.run.warning("You (or the programmer) asked anvi'o to NOT check the consistency of item names\
                              between your additional data and the %s database you are attempting to update. So be it.\
                              Anvi'o will not check anything, but if things don't look the way you expected them to look,\
                              you will not blame anvi'o for your poorly prepared data, but choose between yourself or\
                              Obama." % (self.db_type))
        else:
            self.check_item_names(data_dict)

        db_entries = []
        self.set_next_available_id(self.table_name)
        for item_name in data_dict:
            for key in data_dict[item_name]:
                db_entries.append(tuple([self.next_id(self.table_name),
                                         item_name,
                                         key,
                                         data_dict[item_name][key],
                                         key_types[key]]))

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        database._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % self.table_name, db_entries)
        database.disconnect()

        self.run.info('New data added to the db', '%s.' % (', '.join(keys_list)))


    def list_keys(self):
        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))
        item_additional_data_keys = sorted(database.get_single_column_from_table(self.table_name, 'key', unique=True))

        if not len(item_additional_data_keys):
            self.run.info_single('There are no item additional data in this database :/', nl_before=1, nl_after=1, mc='red')
        else:
            self.run.warning('', 'AVAILABLE DATA KEYS (%d FOUND)' % (len(item_additional_data_keys)), lc='yellow')
            for key in item_additional_data_keys:
                rows = database.get_some_rows_from_table_as_dict(self.table_name, 'key="%s"' % key)
                self.run.info_single('%s (%s, describes %d items)' % (key, list(rows.values())[0]['type'], len(rows)),
                                     nl_after = 1 if key == item_additional_data_keys[-1] else 0)

        database.disconnect()


    def populate_from_file(self, additional_data_file_path, skip_check_names=False):
        filesnpaths.is_file_tab_delimited(additional_data_file_path)

        keys = utils.get_columns_of_TAB_delim_file(additional_data_file_path)
        data = utils.get_TAB_delimited_file_as_dictionary(additional_data_file_path)

        if not len(keys):
            raise ConfigError("There is something wrong with the additional data file at %s.\
                               It does not seem to have any additional keys for data :/" \
                                            % (additional_data_file_path))

        self.add(keys, data, skip_check_names)


class TablesForViews(Table):
    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        Table.__init__(self, self.db_path, get_required_version_for_db(db_path), run, progress)


    def create_new_view(self, data_dict, table_name, table_structure, table_types, view_name=None, append_mode=False):
        """Creates a new view table, and adds an entry for it into the 'views' table.

        Entries in 'views' table appear in various places in the interface. However, we also generate
        view tables to store the type of data we do not wish to display on interfaces, but be able
        access from various other modules. A good example to this is the item_order recipes. When we
        profile a sample, we treat every stplit as their own entity with respect to their mean coverage.
        Although it is great for visualization purposes, it is not useful for item_order purposes since in
        most cases we wish splits to stay together in item_order output. Hence, we create a mean_coverage_splits
        table, where each split holds their own coverage, and we create a mean_coverage_contigs table where each
        split has the coverage of their parent. Clearly the second table is not useful to display. When a table
        is not added as an entry to the 'views' table, then it only exists in the database for other purposes
        than displaying it.

        If a new view does not have a 'view_id', it is not added the 'views' table to provide that flexibility.
        """

        anvio_db = DBClassFactory().get_db_object(self.db_path)

        views_in_db = anvio_db.db.get_table_as_dict(t.views_table_name)

        if not append_mode:
            if view_name and view_name in views_in_db:
                raise ConfigError("TablesForViews speaking: Yo yo yo. You already have a view in the db called '%s'.\
                                    You can't create another one before you get rid of the existing one, because rules."\
                                                                            % view_name)

            # first create the data table:
            anvio_db.db.drop_table(table_name)

        try:
            anvio_db.db.create_table(table_name, table_structure, table_types)
        except:
            if not append_mode:
                raise ConfigError("Table already exists")

        db_entries = [tuple([item] + [data_dict[item][h] for h in table_structure[1:]]) for item in data_dict]
        anvio_db.db._exec_many('''INSERT INTO %s VALUES (%s)''' % (table_name, ','.join(['?'] * len(table_structure))), db_entries)

        if view_name and view_name not in views_in_db:
            anvio_db.db._exec('''INSERT INTO %s VALUES (?,?)''' % t.views_table_name, (view_name, table_name))

        anvio_db.disconnect()


    def remove(self, view_name, table_names_to_blank=[]):
        anvio_db = DBClassFactory().get_db_object(self.db_path)
        anvio_db.db._exec('''DELETE FROM %s WHERE view_id = "%s"''' % (t.views_table_name, view_name))
        for table_name in table_names_to_blank:
            if table_name in anvio_db.db.get_table_names():
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


    def append(self, profile, quiet=False):
        db_entry = tuple([self.next_id(t.variable_nts_table_name)] + [profile[h] for h in t.variable_nts_table_structure[1:]])
        self.db_entries.append(db_entry)
        self.num_entries += 1
        if not quiet and self.num_entries % 100 == 0:
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
                raise ConfigError(error_msg)
            self.__AA_counts_for_bins()
        elif self.contigs_of_interest_file_path:
            if self.profile_db_path or self.genes_of_interest_file_path:
                raise ConfigError(error_msg)
            self.__AA_counts_for_contigs()
        elif self.genes_of_interest_file_path:
            if self.profile_db_path or self.contigs_of_interest_file_path:
                raise ConfigError(error_msg)
            self.__AA_counts_for_genes()
        else:
            self.__AA_counts_for_the_contigs_db()


    def __AA_counts_for_bins(self):
        if not self.collection_name:
            raise ConfigError("You must declare a collection name along with the profile database.")

        profile_db = ProfileDatabase(self.profile_db_path)
        collections_info_table = profile_db.db.get_table_as_dict(t.collections_info_table_name)
        collections_splits_table = profile_db.db.get_table_as_dict(t.collections_splits_table_name)
        profile_db.disconnect()

        if not len(collections_info_table):
            raise ConfigError("There are no collections stored in the profile database :/")

        if not self.collection_name in collections_info_table:
            valid_collections = ', '.join(list(collections_info_table.keys()))
            raise ConfigError("'%s' is not a valid collection name. But %s: '%s'." \
                                    % (self.collection_name,
                                       'these are' if len(valid_collections) > 1 else 'this is',
                                       valid_collections))

        bin_names_in_collection = collections_info_table[self.collection_name]['bin_names'].split(',')

        if self.bin_ids_file_path:
            filesnpaths.is_file_exists(self.bin_ids_file_path)
            bin_names_of_interest = [line.strip() for line in open(self.bin_ids_file_path).readlines()]

            missing_bins = [b for b in bin_names_of_interest if b not in bin_names_in_collection]
            if len(missing_bins):
                raise ConfigError("Some bin names you declared do not appear to be in the collection %s." \
                                            % self.collection_name)
        else:
            bin_names_of_interest = bin_names_in_collection

        collection_dict = utils.get_filtered_dict(collections_splits_table, 'collection_name', set([self.collection_name]))
        collection_dict = utils.get_filtered_dict(collection_dict, 'bin_name', set(bin_names_of_interest))

        split_name_per_bin_dict = {}
        for bin_name in bin_names_of_interest:
            split_name_per_bin_dict[bin_name] = set([])

        for e in list(collection_dict.values()):
            split_name_per_bin_dict[e['bin_name']].add(e['split'])

        for bin_name in bin_names_of_interest:
            self.counts_dict[bin_name] = self.get_AA_counts_dict(split_names=set(split_name_per_bin_dict[bin_name]))['counts']


    def __AA_counts_for_contigs(self):
        filesnpaths.is_file_exists(self.contigs_of_interest_file_path)

        contigs_of_interest = [line.strip() for line in open(self.contigs_of_interest_file_path).readlines()]

        missing_contigs = [True for c in contigs_of_interest if c not in self.contigs_basic_info]
        if missing_contigs:
            raise ConfigError("Some contig names you declared do not seem to be present in the contigs\
                                database :(")

        for contig_name in contigs_of_interest:
            self.counts_dict[contig_name] = self.get_AA_counts_dict(contig_names=set([contig_name]))['counts']


    def __AA_counts_for_genes(self):
        filesnpaths.is_file_exists(self.genes_of_interest_file_path)

        try:
            genes_of_interest = [int(line.strip()) for line in open(self.genes_of_interest_file_path).readlines()]
        except:
            raise ConfigError("Gene call ids in your genes of interest file does not resemble anvi'o gene\
                                call ids (I tried to int them, and it didn't work!)")

        for gene_call in genes_of_interest:
            self.counts_dict[gene_call] = self.get_AA_counts_dict(gene_caller_ids=set([gene_call]))['counts']


    def __AA_counts_for_the_contigs_db(self):
        self.counts_dict[self.args.contigs_db] = self.get_AA_counts_dict()['counts']


    def report(self):
        if self.args.output_file:
            header = ['source'] + sorted(list(self.counts_dict.values())[0].keys())
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


class TablesForGeneCalls(Table):
    def __init__(self, db_path, contigs_fasta=None, run=run, progress=progress, debug=False):
        self.run = run
        self.progress = progress
        self.db_path = db_path
        self.contigs_fasta = contigs_fasta
        self.debug = debug

        is_contigs_db(self.db_path)

        if self.contigs_fasta:
            filesnpaths.is_file_exists(self.contigs_fasta)
            filesnpaths.is_file_fasta_formatted(self.contigs_fasta)


    def check_gene_calls_dict(self, gene_calls_dict):
        if not isinstance(gene_calls_dict, type({})):
            raise ConfigError("Gene calls dict must be a dict instance :/")

        try:
            [int(g) for g in list(gene_calls_dict.keys())]
        except ValueError:
            raise ConfigError("Keys of a gene calls dict must be integers!")

        if False in [x['direction'] in ['f', 'r'] for x in list(gene_calls_dict.values())]:
            raise ConfigError("The values in 'direction' column can't be anything but 'f' (for forward)\
                                or 'r' (for reverse). You have other stuff, and it is not cool.")

        if False in [x['stop'] > x['start'] for x in list(gene_calls_dict.values())]:
            raise ConfigError("For each gene call, the stop position must be bigger than the start position.\
                                Your gene calls dict does not conform to that. If you have reverse gene calls\
                                you must use the 'direction' column to declare that.")

        if False in [(x['stop'] - float(x['start'])) % 3.0 == 0 for x in list(gene_calls_dict.values())]:
            raise ConfigError("Something is wrong with your gene calls. For every gene call, the (stop - start)\
                                should be multiply of 3. It is not the case for all, which is a deal breaker.")


    def use_external_gene_calls_to_populate_genes_in_contigs_table(self, input_file_path, gene_calls_dict=None, ignore_internal_stop_codons=False):
        """Add genes to the contigs database.

           Either provide an `input_file_path` for external gene calls, or provide an
           external gene calls dictionary. The format should follow this:

                {
                  "1": {
                      "contig": "contig_name",
                      "start": 20,
                      "stop": 1544,
                      "direction": "f",
                      "partial": 0,
                      "source": "source_name",
                      "version": "unknown"
                  },

                  "2": {
                    (...)
                  },

                (...)
                }

            If you provide a `gene_calls_dict`, they will be APPENDED to the database. So you
            need to make sure gene caller ids in your dict does not overlap with the ones in
            the database.

        """

        # by default we assume that this is a pristine run. but if the user sends a dictionary
        append_to_the_db = False

        gene_calls_found = False
        # let's do a rigorous check whether the user provided a gene_calls_dict.
        if (gene_calls_dict is not None and gene_calls_dict is not False):
            if not isinstance(gene_calls_dict, dict):
                raise ConfigError("'Use external gene calls' function received a non-empty gene_calls_dict object,\
                                    but it is of type '%s', and not '%s'" % (type(gene_calls_dict), type({})))

            # congrats, we have a dict.
            gene_calls_found = True

            if not len(gene_calls_dict):
                # but it is empty ... silly user.
                self.run.info_single("'Use external gene calls' function found an empty gene calls dict, returning\
                                      prematurely and assuming you know what's up. If you don't, stop here and try to\
                                      identify what decisions you've made might have led you to this weird point your\
                                      workflow (or 'life', totally up to you and your mood, but anvi'o thinks you've\
                                      done great so far.", nl_before=1, nl_after=1)
                return


        if (not input_file_path and not gene_calls_found) or (input_file_path and gene_calls_found):
            raise ConfigError("You must provide either an input file, or an gene calls dict to process external\
                               gene calls. You called `use_external_gene_calls_to_populate_genes_in_contigs_table`\
                               with wrong parameters.")

        Table.__init__(self, self.db_path, anvio.__contigs__version__, self.run, self.progress, simple=True)

        # take care of gene calls dict
        if not gene_calls_found:
            gene_calls_dict = utils.get_TAB_delimited_file_as_dictionary(input_file_path,
                                                                         expected_fields=t.genes_in_contigs_table_structure,
                                                                         only_expected_fields=True,
                                                                         column_mapping=[int, str, int, int, str, int, str, str])

            if not len(gene_calls_dict):
                raise ConfigError("You provided an external gene calls file, but it returned zero gene calls. Assuming that\
                                   this is an error, anvi'o will stop here and complain. If this is not an error and you\
                                   in fact expected this, the proper way of doing this is to use `--skip-gene-calls` flag,\
                                   instead of providing an emtpy external gene calls file. You don't agree? You need this\
                                   for some weird step for you weird pipeline? Let us know, and we will consider changing\
                                   this.")

            self.run.info("External gene calls", "%d gene calls recovered and will be processed." % len(gene_calls_dict))
        else:
            # FIXME: we need to make sure the gene caller ids in the incoming directory is not going to
            #        overwrite an existing gene call. Something like this would have returned the
            #        current max, which could be cross-checked with what's in the dict:
            #
            #            contigs_db = ContigsDatabase(self.db_path)
            #            next_id = contigs_db.db.get_max_value_in_column('genes_in_contigs', 'gene_callers_id') + 1
            #            contigs_db.disconnect()
            append_to_the_db = True

        # recover amino acid sequences. during this operation we are going to have to read all contig sequences
        # into the damn memory. anvi'o is doing a pretty bad job with memory management :(
        amino_acid_sequences = {}

        contig_sequences = {}
        if self.contigs_fasta:
            fasta = u.SequenceSource(self.contigs_fasta)
            while next(fasta):
                contig_sequences[fasta.id] = {'sequence': fasta.seq}
            fasta.close()
        else:
            class Args: None
            args = Args()
            args.contigs_db = self.db_path
            contigs_db = ContigsSuperclass(args, r=terminal.Run(verbose=False))
            contigs_db.init_contig_sequences()
            contig_sequences = contigs_db.contig_sequences

        num_genes_with_internal_stops = 0
        number_of_impartial_gene_calls = 0
        for gene_callers_id in gene_calls_dict:
            gene_call = gene_calls_dict[gene_callers_id]
            contig_name = gene_call['contig']

            if contig_name not in contig_sequences:
                # remove the partial contigs database so things don't get screwed later
                os.remove(self.db_path)
                raise ConfigError("You are in big trouble :( The contig name '%s' in your external gene callers file\
                                    does not appear to be in the contigs FASTA file. How did this happen?" % contig_name)

            if gene_call['partial']:
                amino_acid_sequences[gene_callers_id] = ''
                number_of_impartial_gene_calls += 1
                continue

            sequence = contig_sequences[contig_name]['sequence'][gene_call['start']:gene_call['stop']]
            if gene_call['direction'] == 'r':
                sequence = utils.rev_comp(sequence)

            amino_acid_sequence = utils.get_DNA_sequence_translated(sequence, gene_callers_id)

            # check if there are any internal stops:
            if amino_acid_sequence.find('*') > -1:
                if ignore_internal_stop_codons:
                    amino_acid_sequence = amino_acid_sequence.replace('*', 'X')
                    num_genes_with_internal_stops += 1
                else:
                    os.remove(self.db_path)
                    raise ConfigError("Oops. Anvi'o run into an amino acid seqeunce (that corresponds to the gene callers id '%s')\
                                       which had an internal stop codon :/ This usually indicates that your external gene calls\
                                       have problems. If you still want to continue, you can ask anvi'o to ignore internal stop\
                                       codons on your own risk. It will probably look very ugly on your screen, but here is the\
                                       DNA sequence for that gene in case you don't trust anvi'o (which only would be fair since\
                                       anvi'o does not trust you either): %s" % (str(gene_callers_id), sequence))

            amino_acid_sequences[gene_callers_id] = amino_acid_sequence

        # populate genes_in_contigs, and gene_amino_acid_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, amino_acid_sequences, append_to_the_db=append_to_the_db)

        if num_genes_with_internal_stops:
            percent_genes_with_internal_stops = num_genes_with_internal_stops * 100.0 / len(gene_calls_dict)
            self.run.warning("Please read this carefully: Your external gene calls contained open reading frames with internal\
                              stop codons, and you asked anvi'o to ignore those. Anvi'o replaced internal stop codons with 'X'\
                              characters, and stored them in the contigs database that way. %d of your genes, which corresponded\
                              to %.2f%% of the total %d genes, had internal stop codons. We hope you are happy." % \
                                        (num_genes_with_internal_stops, percent_genes_with_internal_stops, len(gene_calls_dict)))

        if number_of_impartial_gene_calls:
            self.run.warning('%d of your %d gene calls were impartial, hence the translated amino acid sequences for those\
                              were not stored in the database.' % (number_of_impartial_gene_calls, len(gene_calls_dict)))


    def call_genes_and_populate_genes_in_contigs_table(self, gene_caller='prodigal'):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress, simple=True)

        # get gene calls and amino acid sequences
        gene_calls_dict, amino_acid_sequences = self.run_gene_caller(gene_caller)

        # make sure the returning gene calls dict is proper
        self.check_gene_calls_dict(gene_calls_dict)

        # populate genes_in_contigs, and gene_amino_acid_sequences table in contigs db.
        self.populate_genes_in_contigs_table(gene_calls_dict, amino_acid_sequences)


    def run_gene_caller(self, gene_caller='prodigal'):
        """Runs gene caller, and returns gene_calls_dict, and amino acid sequences."""
        remove_fasta_after_processing = False

        if not self.contigs_fasta:
            self.contigs_fasta = self.export_sequences_table_in_db_into_FASTA_file()
            remove_fasta_after_processing = True

        if self.debug:
            self.run.info_single('--debug flag is [ON], which means temporary directories generated by\
                                 this run will not be removed', nl_after=2)

        gene_caller = genecalling.GeneCaller(self.contigs_fasta, gene_caller=gene_caller, debug=self.debug)

        gene_calls_dict, amino_acid_sequences = gene_caller.process()

        if not self.debug and remove_fasta_after_processing:
            os.remove(self.contigs_fasta)

        return gene_calls_dict, amino_acid_sequences


    def populate_genes_in_contigs_table(self, gene_calls_dict, amino_acid_sequences, append_to_the_db=False):
        contigs_db = db.DB(self.db_path, anvio.__contigs__version__)

        if not append_to_the_db:
            contigs_db._exec('''DELETE FROM %s''' % (t.genes_in_contigs_table_name))
            contigs_db._exec('''DELETE FROM %s''' % (t.gene_amino_acid_sequences_table_name))
        else:
            # so we are in the append mode. We must remove all the previous entries from genes in contigs
            # that matches to the incoming sources. otherwhise we may end up with many duplicates in the db.
            sources = set([v['source'] for v in gene_calls_dict.values()])
            for source in sources:
                contigs_db._exec('''DELETE FROM %s WHERE source = "%s"''' % (t.genes_in_contigs_table_name, source))

        self.progress.new('Processing')
        self.progress.update('Entering %d gene calls into the db ...' % (len(gene_calls_dict)))

        db_entries = [tuple([entry_id] + [gene_calls_dict[entry_id][h] for h in t.genes_in_contigs_table_structure[1:]]) for entry_id in gene_calls_dict]
        contigs_db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?,?,?,?)''' % t.genes_in_contigs_table_name, db_entries)

        db_entries = [tuple([entry_id] + [amino_acid_sequences[entry_id]]) for entry_id in gene_calls_dict]
        contigs_db._exec_many('''INSERT INTO %s VALUES (?,?)''' % t.gene_amino_acid_sequences_table_name, db_entries)

        self.progress.end()

        contigs_db.disconnect()


    def populate_genes_in_splits_tables(self):
        Table.__init__(self, self.db_path, anvio.__contigs__version__, run, progress)
        self.init_gene_calls_dict()

        genes_in_splits = GenesInSplits()
        # build a dictionary for fast access to all genes identified within a contig
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
            raise ConfigError("It seems the contigs database '%s' was created with '--skip-gene-calling' flag.\
                                Nothing to do here :/" % (self.db_path))

        self.init_gene_calls_dict()

        if not len(self.gene_calls_dict):
            raise ConfigError("Tables that should contain gene calls are empty. Which probably means the gene\
                                caller reported no genes for your contigs.")

        self.set_next_available_id(t.hmm_hits_table_name)
        self.set_next_available_id(t.hmm_hits_splits_table_name)


    def populate_search_tables(self, sources={}):
        # if we end up generating a temporary file for amino acid sequences:
        if not len(sources):
            import anvio.data.hmm
            sources = anvio.data.hmm.sources

        if not sources:
            return

        target_files_dict = {}

        tmp_directory_path = filesnpaths.get_temp_directory_path()

        # here we will go through targets and populate target_files_dict based on what we find among them.
        targets = set([s['target'] for s in list(sources.values())])
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
                    self.export_sequences_table_in_db_into_FASTA_file(t.gene_amino_acid_sequences_table_name, output_file_path=target_files_dict['AA:GENE'])
                else:
                    target_files_dict['%s:GENE' % alphabet] = os.path.join(tmp_directory_path, '%s_gene_sequences.fa' % alphabet)
                    contigs_db.gen_FASTA_file_of_sequences_for_gene_caller_ids(output_file_path=target_files_dict['%s:GENE' % alphabet],
                                                                               simple_headers=True,
                                                                               rna_alphabet=True if alphabet=='RNA' else False)
            elif context == 'CONTIG':
                if alphabet == 'AA':
                    raise ConfigError("You are somewhere you shouldn't be. You came here because you thought it would be OK\
                                       to ask for AA sequences in the CONTIG context. The answer to that is 'no, thanks'. If\
                                       you think this is dumb, please let us know.")
                else:
                    target_files_dict['%s:CONTIG' % alphabet] = os.path.join(tmp_directory_path, '%s_contig_sequences.fa' % alphabet)
                    utils.export_sequences_from_contigs_db(self.db_path,
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
                                                      alphabet,
                                                      context,
                                                      kind_of_search,
                                                      domain,
                                                      all_genes_searched_against,
                                                      hmm_model,
                                                      reference)

            if not hmm_scan_hits_txt:
                search_results_dict = {}
            else:
                parser = parser_modules['search']['hmmscan'](hmm_scan_hits_txt, alphabet=alphabet, context=context)
                search_results_dict = parser.get_search_results()

            if not len(search_results_dict):
                run.info_single("The HMM source '%s' returned 0 hits. SAD (but it's stil OK)." % source, nl_before=1)


            if context == 'CONTIG':
                # we are in trouble here. because our search results dictionary contains no gene calls, but contig
                # names that contain our hits. on the other hand, the rest of the code outside of this if statement
                # expects a `search_results_dict` with gene callers id in it. so there are two things we need to do
                # to do. one is to come up with some new gene calls and add them to the contigs database. so things
                # will go smoothly downstream. two, we will need to update our `search_results_dict` so it looks
                # like a a dictionary the rest of the code expects with `gene_callers_id` fields. both of these
                # steps are going to be taken care of in the following function. magic.

                self.run.warning("Alright! You just called an HMM profile that runs on contigs. Because it is not\
                                 working with anvi'o gene calls directly, the resulting hits will need to be added\
                                 as 'new gene calls' into the contigs database. This is a new feature, and if it\
                                 starts screwing things up for you please let us know. Other than that you're pretty\
                                 much golden. Carry on.",
                                 header="Psst. Your fancy HMM profile '%s' speaking" % source,
                                 lc="green")

                num_hits_before = len(search_results_dict)
                search_results_dict = utils.get_pruned_HMM_hits_dict(search_results_dict)
                num_hits_after = len(search_results_dict)

                if num_hits_before != num_hits_after:
                    self.run.info('Pruned', '%d out of %d hits were removed due to redundancy' % (num_hits_before - num_hits_after, num_hits_before))

                search_results_dict = self.add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(kind_of_search, search_results_dict)

            self.append(source, reference, kind_of_search, domain, all_genes_searched_against, search_results_dict)

        if not self.debug:
            commander.clean_tmp_dirs()
            for v in list(target_files_dict.values()):
                os.remove(v)


    def add_new_gene_calls_to_contigs_db_and_update_serach_results_dict(self, source, search_results_dict):
        """Add new gene calls to the contigs database and update the HMM `search_results_dict`.

           When we are looking for HMM hits in the context of CONTIGS, our hits do not
           related to the gene calls we already have in a given contigs database. One
           slution is to add additional gene calls for a given set of HMM hits to keep
           them in the database."""

        # we will first learn the next available id in the gene callers table
        contigs_db = ContigsDatabase(self.db_path)
        next_id = contigs_db.db.get_max_value_in_column('genes_in_contigs', 'gene_callers_id') + 1
        contigs_db.disconnect()

        additional_gene_calls = {}
        for e in search_results_dict.values():
            start = e['start']
            stop = e['stop']

            if stop > start:
                direction = 'f'
            else:
                direction = 'r'
                stop, start = start, stop

            partial = 0 if ((stop - start) % 3 == 0) else 1

            # add a new gene call in to the dictionary
            additional_gene_calls[next_id] = {'contig': e['contig_name'],
                                              'start': start,
                                              'stop': stop,
                                              'direction': direction,
                                              'partial': partial,
                                              'source': source,
                                              'version': 'unknown'
                                            }

            # update the search results dictionary with gene callers id:
            e['gene_callers_id'] = next_id

            # update the next available gene callers id:
            next_id += 1

        # update the contigs db with the gene calls in `additional_gene_calls` dict.
        gene_calls_table = TablesForGeneCalls(self.db_path, run=terminal.Run(verbose=False))
        gene_calls_table.use_external_gene_calls_to_populate_genes_in_contigs_table(input_file_path=None,
                                                                                    gene_calls_dict=additional_gene_calls,
                                                                                    ignore_internal_stop_codons=True)
        gene_calls_table.populate_genes_in_splits_tables()

        # refresh the gene calls dict
        self.init_gene_calls_dict()

        self.run.info('Gene calls added to db', '%d (from source "%s")' % (len(additional_gene_calls), source))

        return search_results_dict


    def append(self, source, reference, kind_of_search, domain, all_genes, search_results_dict):
        # we want to define unique identifiers for each gene first. this information will be used to track genes that will
        # break into multiple pieces due to arbitrary split boundaries. while doing that, we will add the 'source' info
        # into the dictionary, so it perfectly matches to the table structure

        for entry_id in search_results_dict:
            hit = search_results_dict[entry_id]

            gene_call = self.gene_calls_dict[hit['gene_callers_id']]

            hit['gene_unique_identifier'] = hashlib.sha224('_'.join([gene_call['contig'], hit['gene_name'], str(gene_call['start']), str(gene_call['stop'])]).encode('utf-8')).hexdigest()
            hit['source'] = source

        self.delete_entries_for_key('source', source, [t.hmm_hits_info_table_name, t.hmm_hits_table_name, t.hmm_hits_splits_table_name])

        contigs_db = ContigsDatabase(self.db_path)

        # push information about this search result into serach_info table.
        db_entries = [source, reference, kind_of_search, domain, ', '.join(all_genes)]
        contigs_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.hmm_hits_info_table_name, db_entries)

        # if our search results were empty, we can return from here.
        if not len(search_results_dict):
            contigs_db.disconnect()
            return

        # then populate serach_data table for each contig.
        db_entries = []
        for hit in list(search_results_dict.values()):
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
        for hit in list(search_results_dict.values()):
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


    def append(self, collection_name, collection_dict, bins_info_dict={}):
        utils.is_this_name_OK_for_database('collection name', collection_name, stringent=False)

        if bins_info_dict:
            if set(collection_dict.keys()) - set(bins_info_dict.keys()):
                raise ConfigError('Bins in the collection dict do not match to the ones in the bins info dict.\
                                    They do not have to be identical, but for each bin id, there must be a unique\
                                    entry in the bins informaiton dict. There is something wrong with your input :/')

        # remove any pre-existing information for 'collection_name'
        self.delete(collection_name)

        num_splits_in_collection_dict = sum([len(splits) for splits in list(collection_dict.values())])
        splits_in_collection_dict = set(list(chain.from_iterable(list(collection_dict.values()))))
        if len(splits_in_collection_dict) != num_splits_in_collection_dict:
            raise ConfigError("TablesForCollections::append: %d of the split or contig IDs appear more than once in\
                                your collections input. It is unclear to anvi'o how did you manage to do this, but we\
                                cannot go anywhere with this :/" % (num_splits_in_collection_dict - len(splits_in_collection_dict)))

        database = db.DB(self.db_path, get_required_version_for_db(self.db_path))

        # how many clusters are defined in 'collection_dict'?
        bin_names = list(collection_dict.keys())

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

        self.run.info('Collections', 'The collection "%s" that describes %s splits has been successfully added to the database at "%s".'\
                                        % (collection_name, pp(num_splits), self.db_path), mc='green')


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


    def list_states(self):
        state_names = sorted(list(self.states.keys()))

        self.run.warning('', 'AVAILABLE STATES (%d FOUND)' % (len(self.states)), lc='yellow')
        for state_name in state_names:
            self.run.info_single('%s (last modified on %s)' % (state_name, self.states[state_name]['last_modified']),
                                 nl_after = 1 if state_name == state_names[-1] else 0)


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
            raise ConfigError("Something is wrong. The contigs database says that genes were now called, and here\
                                you are trying to populate taxonomy tables for genes. No, thanks.")

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
            raise ConfigError("Taxonomy information for genes you are trying to import into the database contains\
                                %s gene caller ids that do not appear to be in the database. This is a step you must\
                                be very careful to make sure you are not importing annotations for genes that have\
                                nothing to do with your contigs database. To make sure of that, you should always work\
                                with `anvi-get-dna-sequences-for-gene-calls` or `anvi-get-aa-sequences-for-gene-calls` programs\
                                to get the data to annotate. For instance one of the gene caller ids you have in your\
                                input data that does not appear in the database is this one: '%s'. Anvi'o hopes it makes\
                                sense to you, because it definitely does not make any sense to anvi'o :("\
                                                        % (len(gene_caller_ids_missing_in_db), str(gene_caller_ids_missing_in_db.pop())))

        # check whether input matrix dict
        keys_found =  list(self.taxon_names_dict.values())[0].keys()
        missing_keys = [key for key in t.taxon_names_table_structure[1:] if key not in keys_found]
        if len(missing_keys):
            raise ConfigError("Anvi'o is trying to get ready to create tables for taxonomy, but there is something\
                                wrong :( The taxonomy names dict (one of the required input dictionaries to the class\
                                seems to be missing a one or more keys that are necessary to finish the job. Here is \
                                a list of missing keys: %s. The complete list of input keys should contain these: %s."\
                                        % (', '.join(missing_keys), ', '.join(t.taxon_names_table_structure[1:])))

        if not len(self.taxon_names_dict):
            raise ConfigError("Anvi'o is trying to get ready to create tables for taxonomy, but taxonomy names dict\
                                (one of the required input dictionaries to the class responsible for this task) seems\
                                to be empty.")


    def populate_splits_taxonomy_table(self):
        """Populate the taxonomy information per split"""

        # build a dictionary for fast access to all genes identified within a contig
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


class TableForGeneClusters(Table):
    """A class to populte gene clusters table in a given pan db.

      Here is an example:

        >>> table = TableForGeneClusters(db_path)
        >>> for ...:
        >>>     table.add({'gene_caller_id': gene_caller_id, 'gene_cluster_id': gene_cluster_id, 'genome_name': genome_name})
        >>> table.store()
    """

    def __init__(self, db_path, run=run, progress=progress):
        self.db_path = db_path

        is_pan_db(db_path)

        self.run = run
        self.progress = progress

        Table.__init__(self, self.db_path, anvio.__pan__version__, run, progress)

        self.set_next_available_id(t.pan_gene_clusters_table_name)

        self.entries = []


    def add(self, entry_dict):
        self.entries.append([entry_dict[key] for key in t.pan_gene_clusters_table_structure[1:]])


    def store(self):
        self.delete_contents_of_table(t.pan_gene_clusters_table_name, warning=False)

        db_entries = [tuple([self.next_id(t.pan_gene_clusters_table_name)] + entry) for entry in self.entries]
        pan_db = PanDatabase(self.db_path)
        pan_db.db._exec_many('''INSERT INTO %s VALUES (?,?,?,?,?)''' % t.pan_gene_clusters_table_name, db_entries)
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
        gene_function_sources = set([v['source'] for v in list(functions_dict.values())])
        unique_num_genes = len(set([v['gene_callers_id'] for v in list(functions_dict.values())]))

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


class TableForNtPositions(object):
    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.numpy_data_type = 'uint8'
        self.db_entries = []


    def append(self, contig_name, position_info_list):
        position_info_blob = utils.convert_numpy_array_to_binary_blob(numpy.array(position_info_list, dtype=self.numpy_data_type))
        self.db_entries.append((contig_name, position_info_blob, ))


    def store(self, db):
        db.insert_many(t.nt_position_info_table_name, entries=self.db_entries)
        self.db_entries = []


class GenesInSplits:
    def __init__(self):
        self.entry_id = 0
        self.splits_to_prots = {}

    def add(self, split_name, split_start, split_end, gene_callers_id, prot_start, prot_end):

        gene_length = prot_end - prot_start

        if gene_length <= 0:
            raise ConfigError("dbops.py/GeneInSplits: OK. There is something wrong. We have this gene, '%s',\
                                which starts at position %d and ends at position %d. Well, it doesn't look right,\
                                does it?" % (gene_callers_id, prot_start, prot_end))

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
        raise ConfigError("'%s' is not an anvi'o contigs database." % db_path)


def is_pan_or_profile_db(db_path):
    if get_db_type(db_path) not in ['pan', 'profile']:
        raise ConfigError("'%s' is neither a pan nor a profile database :/ Someone is in trouble." % db_path)


def is_profile_db(db_path):
    if get_db_type(db_path) != 'profile':
        raise ConfigError("'%s' is not an anvi'o profile database." % db_path)


def is_pan_db(db_path):
    if get_db_type(db_path) != 'pan':
        raise ConfigError("'%s' is not an anvi'o pan database." % db_path)


def is_samples_db(db_path):
    if get_db_type(db_path) != 'samples_information':
        raise ConfigError("'%s' is not an anvi'o samples database." % db_path)


def is_db_ok_to_create(db_path, db_type):
    if os.path.exists(db_path):
        raise ConfigError("Anvi'o will not overwrite an existing %s database. Please choose a different name\
                            or remove the existing database ('%s') first." % (db_type, db_path))

    if not db_path.lower().endswith('.db'):
        raise ConfigError("Please make sure the file name for your new %s db has a '.db' extension. Anvi'o developers\
                            apologize for imposing their views on how anvi'o databases should be named, and are\
                            humbled by your cooperation." % db_type)


def get_required_version_for_db(db_path):
    db_type = get_db_type(db_path)

    if db_type not in t.versions_for_db_types:
        raise ConfigError("Anvi'o was trying to get the version of the -alleged- anvi'o database '%s', but it failed\
                            because it turns out it doesn't know anything about this '%s' type." % (db_path, db_type))

    return t.versions_for_db_types[db_type]


def is_blank_profile(db_path):
    if get_db_type(db_path) != 'profile':
        return False

    database = db.DB(db_path, None, ignore_version=True)
    blank = database.get_meta_value('blank')
    database.disconnect()

    return blank


def get_db_type(db_path):
    filesnpaths.is_file_exists(db_path)

    try:
        database = db.DB(db_path, None, ignore_version=True)
    except:
        raise ConfigError('Are you sure "%s" is a database file? Because, you know, probably\
                            it is not at all..' % db_path)

    tables = database.get_table_names()
    if 'self' not in tables:
        database.disconnect()
        raise ConfigError("'%s' does not seem to be a anvi'o database..." % db_path)

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
        raise ConfigError('The contigs database and the profile database does not\
                           seem to be compatible. More specifically, this contigs\
                           database is not the one that was used when %s generated\
                           this profile database (%s != %s).'\
                               % ('anvi-merge' if merged else 'anvi-profile', a_hash, p_hash))

    return True


def is_profile_db_and_samples_db_compatible(profile_db_path, samples_db_path, manual_mode_exception=False):
    """Check whether every sample name in the profile database is represented in the samples information database"""
    is_profile_db(profile_db_path)
    is_samples_db(samples_db_path)

    profile_db = ProfileDatabase(profile_db_path)
    samples_db = SamplesInformationDatabase(samples_db_path)

    if 'merged' in profile_db.meta and not int(profile_db.meta['merged']):
        raise ConfigError("Samples databases are only useful if you are working on a merged profile.")

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
            raise ConfigError("Anvi'o is upset with you :/ Please make sure your samples information files (or your\
                                samples database) contain sample names from your data file. These sample names are in\
                                your samples information, but not in your data file: '%s'. If this error does not make\
                                any sense to you, please contact an anvi'o developer." % ', '.join(samples_in_samples_db_but_not_in_profile_db))
        return


    missing_samples = profile_db.samples - samples_db.samples
    num_represented_samples = len(profile_db.samples) - len(missing_samples)

    if len(missing_samples):
        how_much_of_the_samples_are_represented_txt = 'none' if len(missing_samples) == len(profile_db.samples) else\
                                                      'only %d of %d' % (num_represented_samples, len(profile_db.samples))

        raise ConfigError("The samples information database you provided ('%s') does not seem to agree well with the profile\
                            database ('%s'). More specifically, %s of the samples in the profile database are repesented in\
                            the samples information database. Names for these missing samples go like this: %s ...,\
                            while the sample names in the samples information database go like this: %s ... This could be due to\
                            a simple typo, or you may be using the wrong or outdated samples information database. You may need to\
                            regenerate the samples information database to fix this problem :/"\
                                                % (samples_db_path, profile_db_path, how_much_of_the_samples_are_represented_txt,
                                                   ', '.join(list(missing_samples)[0:3]), ', '.join(list(samples_db.samples)[0:3])))

    if samples_db.sample_names_for_order:
        missing_samples = samples_db.sample_names_for_order - profile_db.samples

        if len(missing_samples):
            raise ConfigError("The samples order information in the samples database do not match with the sample names in\
                                the profile database (or the input data). To be precise, %d sample(s) occur(s) only in the\
                                samples database, and not found in the profile database (or in the input data). Here is some of\
                                them: %s ..." % (len(missing_samples), ', '.join(list(missing_samples)[0:3])))


def get_split_names_in_profile_db(profile_db_path):
    is_profile_db(profile_db_path)

    profile_db = ProfileDatabase(profile_db_path)

    if int(profile_db.meta['blank']):
        run.warning("dbops::get_split_names_in_profile_db is speaking. Someone asked for the split names in a blank profile database.\
                     Sadly, anvi'o does not keep track of split names in blank profile databases. This function will return an\
                     empty set as split names to not kill your mojo, but whatever you were trying to do will not work :(")
        split_names = set([])
    elif int(profile_db.meta['merged']):
        split_names = set(profile_db.db.get_single_column_from_table('mean_coverage_Q2Q3_splits', 'contig'))
    else:
        split_names = set(profile_db.db.get_single_column_from_table('atomic_data_splits', 'contig'))

    profile_db.disconnect()

    return split_names


def get_auxiliary_data_path_for_profile_db(profile_db_path):
    return  os.path.join(os.path.dirname(profile_db_path), 'AUXILIARY-DATA.db')


def get_description_in_db(anvio_db_path, run=run):
    """Reads the description in an anvi'o database"""

    anvio_db = db.DB(anvio_db_path, None, ignore_version=True)
    description = None
    try:
        description = anvio_db.get_meta_value('description')
    except ConfigError:
        description = '_No description is available_'

    anvio_db.disconnect()
    return description


def update_description_in_db_from_file(anvio_db_path, description_file_path, run=run):
    filesnpaths.is_file_plain_text(description_file_path)
    description = open(os.path.abspath(description_file_path), 'rU').read()

    update_description_in_db(anvio_db_path, description, run=run)


def update_description_in_db(anvio_db_path, description, run=run):
    """Updates the description in an anvi'o database"""

    if not isinstance(description, str):
        raise ConfigError("Description parameter must be of type `string`.")

    db_type = get_db_type(anvio_db_path)

    anvio_db = db.DB(anvio_db_path, None, ignore_version=True)
    anvio_db.remove_meta_key_value_pair('description')
    anvio_db.set_meta_value('description', description)
    anvio_db.disconnect()

    run.info_single("The anvi'o %s database has just been updated with a description that contains %d words\
                     and %d characters." % (db_type, len(description.split()), len(description)))


def do_hierarchical_clustering_of_items(anvio_db_path, clustering_configs, split_names=[], database_paths={}, input_directory=None, default_clustering_config=None, \
                                distance=constants.distance_metric_default, linkage=constants.linkage_method_default, run=run, progress=progress):
    """This is just an orphan function that computes hierarchical clustering w results
       and calls the `add_hierarchical_clustering_to_db` function with correct input.

       Ugly but useful --yet another one of those moments in which we sacrifice
       important principles for simple conveniences."""

    from anvio.clusteringconfuguration import ClusteringConfiguration
    from anvio.clustering import order_contigs_simple

    for config_name in clustering_configs:
        config_path = clustering_configs[config_name]

        config = ClusteringConfiguration(config_path, input_directory, db_paths=database_paths, row_ids_of_interest=split_names)

        try:
            clustering_name, newick = order_contigs_simple(config, distance=distance, linkage=linkage, progress=progress)
        except Exception as e:
            progress.end()
            run.warning('Clustering has failed for "%s": "%s"' % (config_name, e))
            continue

        _, distance, linkage = clustering_name.split(':')

        add_hierarchical_clustering_to_db(anvio_db_path, config_name, newick, distance=distance, linkage=linkage, make_default=config_name == default_clustering_config, run=run)



def add_hierarchical_clustering_to_db(anvio_db_path, clustering_name, clustering_newick, distance, linkage, make_default=False, run=run):
    """Backwards compatibility function.

       We can fix all instances of `add_hierarchical_clustering_to_db` everywhere in the code to work
       with `add_items_order_to_db` function directly, and end this tyranny."""

    add_items_order_to_db(anvio_db_path,
                          clustering_name,
                          order_data=clustering_newick,
                          order_data_type_newick=True,
                          distance=distance,
                          linkage=linkage,
                          make_default=False,
                          run=run)


def add_items_order_to_db(anvio_db_path, order_name, order_data, order_data_type_newick=True, distance=None, linkage=None, make_default=False, run=run):
    """Adds a new clustering into an anvi'o db"""

    if order_data_type_newick and (not distance or not linkage):
        raise ConfigError("You are trying to add a newick-formatted clustering dendrogram to the database without providing\
                           distance and linkage data that generated this dendrogram :/")

    if not order_data_type_newick and (distance or linkage):
        raise ConfigError("Distance and linkage variables are only relevant if you are trying to add a newick-formatted\
                           clustering dendrogram. But your function call suggests you are not.")

    # let's learn who we are dealing with:
    db_type = get_db_type(anvio_db_path)

    # replace clustering id with a text that contains distance and linkage information
    if order_data_type_newick:
        order_name = ':'.join([order_name, distance, linkage])
    else:
        order_name = ':'.join([order_name, 'NA', 'NA'])

    anvio_db = DBClassFactory().get_db_object(anvio_db_path)

    if t.item_orders_table_name not in anvio_db.db.get_table_names():
        raise ConfigError("You can't add a new items order into this %s database (%s). You know why? Becasue it doesn't\
                           have a table for 'item_order' :(" % (db_type, anvio_db_path))

    try:
        available_item_orders = anvio_db.db.get_meta_value('available_item_orders').split(',')
    except:
        available_item_orders = []

    if order_name in available_item_orders:
        run.warning('Clustering for "%s" is already in the database. Its content will\
                     be replaced with the new one.' % (order_name))

        anvio_db.db._exec('''DELETE FROM %s where name = "%s"''' % (t.item_orders_table_name, order_name))
    else:
        available_item_orders.append(order_name)

    anvio_db.db._exec('''INSERT INTO %s VALUES (?,?,?)''' % t.item_orders_table_name, tuple([order_name, 'newick' if order_data_type_newick else 'basic', order_data]))

    anvio_db.db.set_meta_value('available_item_orders', ','.join(available_item_orders))

    # We don't consider basic orders as orders becasue we are rebels.
    if order_data_type_newick:
        anvio_db.db.set_meta_value('gene_clusters_ordered' if db_type == 'pan' else 'contigs_ordered', True)

    try:
        anvio_db.db.get_meta_value('default_item_order')
        default_item_order_is_set = True
    except:
        default_item_order_is_set = False

    if make_default or not default_item_order_is_set:
        anvio_db.db.set_meta_value('default_item_order', order_name)

    anvio_db.disconnect()

    run.info('New items order', '"%s" (type %s) has been added to the database...' % (order_name, 'newick' if order_data_type_newick else 'basic'))


def get_default_item_order_name(default_item_order_requested, item_orders_dict, progress=progress, run=run):
    """Get the proper default item_order given the desired default with respect to available item_orders.

       This is tricky. We have some deault item_orders defined in the constants. For instance, for the
       merged profiles we want the default to be 'tnf-cov', for single profiles we want it to be 'tnf',
       etc. The problem is that these defaults do not indicate any distance metric or linkages,
       even though anvi'o allows users to define those variables freely in cluster configurations.

       A item_order dict can contain multiple clustrings. The purpose of this function is to take the
       desired default into consideration, but then find a working one if it is not available, or there
       are multiple ones in the dict.
    """

    if not item_orders_dict:
        raise ConfigError("You requested to get the default item_order given the item_order dictionary,\
                            but the item_order dict is empty :/ ")

    matching_item_order_names = [item_order for item_order in item_orders_dict if item_order.lower().split(':')[0] == default_item_order_requested.lower()]

    if not len(matching_item_order_names):
        default_item_order = list(item_orders_dict.keys())[0]
        run.warning('`get_default_item_order_name` function is concerned, because nothing in the item_orders\
                     dict matched to the desired order class "%s". So the order literally set to "%s"\
                     (a class of "%s") randomly as the default order. Good luck :/' % (default_item_order_requested,
                                                                                 default_item_order,
                                                                                 default_item_order.split(':')[0]))
        return default_item_order
    elif len(matching_item_order_names) == 1:
        return matching_item_order_names[0]
    else:
        default_item_order = matching_item_order_names[0]
        run.warning('`get_default_item_order_name` function is concerned, because there were multiple entries\
                     in the item_orders dict matched to the desired default order class "%s". So it set\
                     the first of all %d matching item_orders, which happened to be the "%s", as the\
                     default. We hope that will not screw up your mojo :/' % (default_item_order_requested,
                                                                              len(matching_item_order_names),
                                                                              default_item_order))
        return default_item_order


def export_aa_sequences_from_contigs_db(contigs_db_path, output_file_path, gene_caller_ids=set([]), quiet=False):
    filesnpaths.is_file_exists(contigs_db_path)
    filesnpaths.is_output_file_writable(output_file_path)

    class T(Table):
        def __init__(self, db_path, version, run=run, progress=progress, quiet=False):
            Table.__init__(self, db_path, version, run, progress, quiet=quiet)

    h = T(contigs_db_path, anvio.__contigs__version__, quiet=quiet)
    h.export_sequences_table_in_db_into_FASTA_file(t.gene_amino_acid_sequences_table_name,
                                                   output_file_path=output_file_path,
                                                   item_names=gene_caller_ids)

    return output_file_path


def get_all_item_names_from_the_database(db_path):
    """Return all split names or gene cluster names in a given database"""

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
        all_items = set(PanSuperclass(args).gene_cluster_names)
    elif db_type == 'contigs':
        args.contigs_db = db_path
        all_items = set(ContigsSuperclass(args).splits_basic_info.keys())
    else:
        raise ConfigError("You wanted to get all items in the database %s, but no one here knows aobut its type. Seriously,\
                            what is '%s'?" % (db_path, db_type))

    if not len(all_items):
        raise ConfigError("dbops::get_all_item_names_from_the_database speaking. Something that should never happen happened :/\
                            There seems to be nothing in this %s database. Anvi'o is as confused as you are. Please get in touch\
                            with a developer. They will love this story.")

    return all_items
