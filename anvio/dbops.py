# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes to create, access, and/or populate contigs, tRNASeq, and profile databases.
"""

import os
import re
import sys
import time
import copy
import json
import numpy
import random
import argparse
import textwrap
import itertools
import scipy.signal

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from io import StringIO
from collections import Counter

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.fastalib as u
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.contigops as contigops
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccolections
import anvio.genomestorage as genomestorage
import anvio.auxiliarydataops as auxiliarydataops
import anvio.homogeneityindex as homogeneityindex

from anvio.drivers import Aligners
from anvio.errors import ConfigError
from anvio.dbinfo import DBInfo as dbi

from anvio.tables.states import TablesForStates
from anvio.tables.genecalls import TablesForGeneCalls
from anvio.tables.ntpositions import TableForNtPositions
from anvio.tables.miscdata import TableForItemAdditionalData
from anvio.tables.miscdata import TableForLayerAdditionalData
from anvio.tables.kmers import KMerTablesForContigsAndSplits
from anvio.tables.genelevelcoverages import TableForGeneLevelCoverages
from anvio.tables.contigsplitinfo import TableForContigsInfo, TableForSplitsInfo

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize
aligners = Aligners()


class DBClassFactory:
    """Factory pattern to get the appropriate class for a given anvi'o db type"""
    def __init__(self):
        self.DB_CLASSES = {'profile': ProfileDatabase,
                           'contigs': ContigsDatabase,
                           'trnaseq': TRNASeqDatabase,
                           'pan': PanDatabase,
                           'genes': GenesDatabase}

    def get_db_class(self, db_path):
        db_type = utils.get_db_type(db_path)

        if db_type not in self.DB_CLASSES:
            raise ConfigError("DBClassFactory speaking. I do not know a class for database type "
                               "%s :/ I can deal with these though: '%s'" % (db_type, ', '.join(self.DB_CLASSES)))
        return self.DB_CLASSES[db_type]


    def get_db_object(self, db_path):
        anvio_db_class = self.get_db_class(db_path)
        return anvio_db_class(db_path)


class ContigsSuperclass(object):
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        # if ContigsSuperclass is being used as a base class, then the class that inherits
        # it may have set a `split_names_of_interest` variable. in which case we don't want to
        # overwrite it, and use it for any task that requires a focus on a particular set of
        # splits. But if there is non set by the class that inherits ContigsSuperclass, then
        # we don't want functions that use/assume `split_names_of_interest` variable to throw
        # errors because the variable is not set. so what we are doing here ensures that.
        if not hasattr(self, 'split_names_of_interest'):
            self.split_names_of_interest = set([])

        if hasattr(self.args, 'split_names_of_interest'):
            self.split_names_of_interest = set(self.args.split_names_of_interest)

        self.a_meta = {}

        self.splits_basic_info = {}
        self.splits_taxonomy_dict = {}
        self.split_sequences = {}
        self.contigs_basic_info = {}
        self.nt_positions_info = {}

        self.contig_sequences = {}
        self.gene_caller_ids_included_in_contig_sequences_initialized = set([])

        self.genes_in_contigs_dict = {}
        self.gene_lengths = {}
        self.contig_name_to_genes = {}
        self.genes_in_splits = {} # keys of this dict are NOT gene caller ids. they are ids for each entry.
        self.split_name_to_genes_in_splits_entry_ids = {} # for fast access to all self.genes_in_splits entries for a given split
        self.gene_callers_id_to_split_name_dict = {} # for fast access to a split name that contains a given gene callers id

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
        D = lambda x: self.__dict__[x] if x in self.__dict__ else None
        if D('mode') == 'pan' or D('mode') == 'functional' or D('mode') == 'manual':
            return

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.contigs_db_path = A('contigs_db')

        if not self.contigs_db_path:
            raise ConfigError("Someone (hopefully, you) is trying to initialize the Contigs Super Class without a contigs database path. "
                              "There are many ways this can happen, but .. do you think you were trying to run anvio-interactive in "
                              "manual mode but without a --manual flag right before this? Just a gut feeling... No? Maybe you created "
                              "an instance of the profile superclass without a `contigs_db` argument in the namespace? No? Well, then we "
                              "may be in a really big trouble. Please run what you did before seeing this again with a `--debug` flag, "
                              "and send us an e-mail :(")

        filesnpaths.is_file_exists(self.contigs_db_path)

        self.progress.new('Loading the contigs DB')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path, run=self.run, progress=self.progress)

        self.a_meta = contigs_db.meta

        self.a_meta['creation_date'] = utils.get_time_to_date(self.a_meta['creation_date']) if 'creation_date' in self.a_meta else 'unknown'

        ####################################################################################
        # MINDFULLY READING STUFF FROM THE DATABASE
        ####################################################################################
        # read SPLITS and GENES basic information.
        self.splits_basic_info = contigs_db.db.smart_get(t.splits_info_table_name, 'split', self.split_names_of_interest, progress=self.progress)
        self.genes_in_splits = contigs_db.db.smart_get(t.genes_in_splits_table_name, 'split', self.split_names_of_interest, progress=self.progress)

        # if there are no splits names of interest, contig names of interest will be an empty set.
        # that's OK, because `smart_get` will take care of it.
        contig_names_of_interest = set([self.splits_basic_info[s]['parent'] for s in self.split_names_of_interest])

        # read CONTIGS and GENES basic information.
        self.contigs_basic_info = contigs_db.db.smart_get(t.contigs_info_table_name, 'contig', contig_names_of_interest, string_the_key=True, progress=self.progress)
        self.genes_in_contigs_dict = contigs_db.db.smart_get(t.genes_in_contigs_table_name, 'contig', contig_names_of_interest, progress=self.progress)

        # because this table is as dumb as Eric, it needs some special attention
        if self.split_names_of_interest:
            self.progress.update('Reading **SOME** entries in the nucleotide positions info table :)')
            where_clause = """contig_name IN (%s)""" % ','.join(['"%s"' % c for c in contig_names_of_interest])
            for contig_name, nt_positions_info in contigs_db.db.get_some_rows_from_table(t.nt_position_info_table_name, where_clause=where_clause):
                self.nt_positions_info[contig_name] = utils.convert_binary_blob_to_numpy_array(nt_positions_info, 'uint8')
        else:
            self.progress.update('Reading **ALL** entries in the nucleotide positions info table :(')
            for contig_name, nt_positions_info in contigs_db.db.get_all_rows_from_table(t.nt_position_info_table_name):
                self.nt_positions_info[contig_name] = utils.convert_binary_blob_to_numpy_array(nt_positions_info, 'uint8')
        ####################################################################################
        # /MINDFULLY READING STUFF FROM THE DATABASE
        ####################################################################################

        self.progress.update('Populating gene lengths dict')
        self.gene_lengths = dict([(g, (self.genes_in_contigs_dict[g]['stop'] - self.genes_in_contigs_dict[g]['start'])) for g in self.genes_in_contigs_dict])

        self.progress.update('Populating contig name to gene IDs dict')
        for contig_name in self.contigs_basic_info:
            self.contig_name_to_genes[contig_name] = set([])
        for gene_unique_id in self.genes_in_contigs_dict:
            e = self.genes_in_contigs_dict[gene_unique_id]
            self.contig_name_to_genes[e['contig']].add((gene_unique_id, e['start'], e['stop']), )


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

        self.progress.end()

        contigs_db.disconnect()

        self.run.info('Contigs DB', 'Initialized: %s (v. %s)' % (self.contigs_db_path, anvio.__contigs__version__))


    def init_splits_taxonomy(self, t_level='t_genus'):
        if not self.contigs_db_path:
            return

        if t_level not in t.taxon_names_table_structure[1:]:
            raise ConfigError("Pretty close. But the taxonomic level '%s' is not known to anvi'o. How about "
                               "one of these: %s." % (t_level, ','.join(t.taxon_names_table_structure[1:])))

        self.progress.new('Initializing splits taxonomy')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        splits_taxonomy_table = contigs_db.db.smart_get(t.splits_taxonomy_table_name, 'split', self.split_names_of_interest, string_the_key=True, error_if_no_data=False, progress=self.progress)
        taxon_names_table = contigs_db.db.get_table_as_dict(t.taxon_names_table_name)

        for split_name in splits_taxonomy_table:
            taxon_id = splits_taxonomy_table[split_name]['taxon_id']
            if taxon_id:
                if t_level in taxon_names_table[taxon_id] and taxon_names_table[taxon_id][t_level]:
                    self.splits_taxonomy_dict[split_name] = taxon_names_table[taxon_id][t_level]

        contigs_db.disconnect()
        self.progress.end()

        if len(splits_taxonomy_table):
            self.run.info('Splits taxonomy', 'Initiated for taxonomic level for "%s"' % t_level)


    def init_contig_sequences(self, min_contig_length=0, gene_caller_ids_of_interest=set([]), split_names_of_interest=set([]), contig_names_of_interest=set([])):
        contigs_db = ContigsDatabase(self.contigs_db_path)

        if not len(split_names_of_interest) and not len(gene_caller_ids_of_interest) and not len(contig_names_of_interest):
            split_names_of_interest = self.split_names_of_interest

        too_many_args = False
        if len(gene_caller_ids_of_interest):
            if len(split_names_of_interest):
                too_many_args = True
                opt1, opt2 = 'gene caller ids of interest', 'split names of interest'
            elif len(contig_names_of_interest):
                too_many_args = True
                opt1, opt2 = 'gene caller ids of interest', 'contig names of interest'
        elif len(split_names_of_interest):
            if len(contig_names_of_interest):
                too_many_args = True
                opt1, opt2 = 'split names of interest', 'contig names of interest'
        if too_many_args:
            raise ConfigError("Ehem. Someone just called `init_contig_sequences` with %s AND %s. "
                              "Someone should make up their mind and go for only one of those." % (opt1, opt2))

        # are we going to read everything, or only those that are of interest?
        if contig_names_of_interest:
            subset_provided = True
        elif gene_caller_ids_of_interest:
            subset_provided = True
            contig_names_of_interest = set([])
            for gene_callers_id in self.genes_in_contigs_dict:
                if gene_callers_id in gene_caller_ids_of_interest:
                    contig_names_of_interest.add(self.genes_in_contigs_dict[gene_callers_id]['contig'])
        elif split_names_of_interest:
            subset_provided = True
            contig_names_of_interest = set([self.splits_basic_info[s]['parent'] for s in split_names_of_interest])
        else:
            subset_provided = False

        if subset_provided and not len(contig_names_of_interest):
            raise ConfigError("Anvi'o was trying to identify the contig names of interest in `init_contig_sequences` "
                              "and then after a few steps there were no contig names of interest at all :( Something "
                              "fishy happened, and code is Jon Snow.")

        self.progress.new('Loading contig sequences')
        self.contig_sequences = contigs_db.db.smart_get(t.contig_sequences_table_name, 'contig', contig_names_of_interest, string_the_key=True, progress=self.progress)

        # now we have our contig sequences, and it is time to make sure we
        # know which gene caller ids they represent if only a subset of contigs
        # were requested:
        if subset_provided:
            self.gene_caller_ids_included_in_contig_sequences_initialized = set(contigs_db.db.smart_get(t.genes_in_contigs_table_name, 'contig', contig_names_of_interest, string_the_key=False, progress=self.progress).keys())
        else:
            self.gene_caller_ids_included_in_contig_sequences_initialized = set(self.genes_in_contigs_dict.keys())

        self.progress.end()

        if subset_provided:
            self.run.warning(f"Someone asked the Contigs Superclass to initialize only a subset of contig sequences. "
                             f"Usually this is a good thing and means that some good code somewhere is looking after "
                             f"you. Just FYI, this class will only know about {P('contig sequence', len(contig_names_of_interest))} "
                             f"instead of all the things in the database.", header="THE MORE YOU KNOW 🌈", lc='yellow')

        contigs_db.disconnect()

        self.progress.new('Filtering contig sequences')
        self.progress.update('Identifying contigs shorter than M')
        contigs_shorter_than_M = set([c for c in self.contigs_basic_info if self.contigs_basic_info[c]['length'] < min_contig_length])

        self.progress.update('Filtering out shorter contigs')
        for contig_name in contigs_shorter_than_M:
            self.contig_sequences.pop(contig_name) if contig_name in self.contig_sequences else None

        self.progress.end()

        return contigs_shorter_than_M


    def get_items_additional_data_for_functions_per_split_summary(self, source, split_names_of_interest, data_dict={}, keys_list=[]):
        """Get items additional data layers to display the frequency of function names
           for each split in a given contigs database so it can be shown as a stacked bar
           chart in the anvi'o interactive interface.

        Parameters
        ==========
        source : str
            A functional annotation source that is in the contigs database.
        split_names_of_interest : list
            Split names to be considered.
        data_dict : dict
            An optional `items_additional_data_dict` type dictionary to update.
        keys_list : list
            An optional `items_additional_data_keys` type list to update.

        Returns
        =======
        data_dict : dict
            An `items_additional_data_dict` type dictionary.
        keys_list : list
            An `items_additional_data_keys` type list.
        """

        if not self.gene_function_calls_initiated:
            raise ConfigError("For this to work, someone needs to initialize gene functions first :/")

        if source not in self.gene_function_call_sources:
            raise ConfigError(f"Nice try. Your '{source}' is not a valid function annotation source in this "
                              f"contigs database. You will need to choose one of these: {', '.join(self.gene_function_call_sources)}")

        if len(data_dict) or len(keys_list):
            if not len(data_dict) and len(keys_list):
                raise ConfigError("If you are sending a data dictionary to expand with function summaries "
                                  "per split, then you also need to send a keys dictionary to be expanded.")

        # learn the number of categories for the function source
        function_source_categories = set([])
        for entry in self.gene_function_calls_dict.values():
            if entry[source]:
                function_source_category = entry[source][1].split('!!!')[0]
                function_source_categories.add(function_source_category)

        function_source_categories = sorted(list(function_source_categories))

        # we can't use this strategy if there are many categories for a given
        # function source
        if len(function_source_categories) > 10:
            raise ConfigError(f"The functional annotation source '{source}' has {len(function_source_categories)} "
                              f"which is way too many to summarize into 'per split' data. If you think this is "
                              f"a dumb reason to not do this, please let us know and we will try to find a better "
                              f"solution to this. In the current implementation, any function annotation source "
                              f"that has up to 10 categories is good.")

        # create a template dictionaries to hold the category frequencies and default entries
        # for splits with no information
        _frequency_of_categories = dict([(cat, 0) for cat in function_source_categories])

        for split_name in split_names_of_interest:
            frequency_of_categories = copy.deepcopy(_frequency_of_categories)

            for entry_id in self.split_name_to_genes_in_splits_entry_ids[split_name]:
                gene_callers_id = self.genes_in_splits[entry_id]['gene_callers_id']
                if gene_callers_id in self.gene_function_calls_dict and self.gene_function_calls_dict[gene_callers_id][source]:
                    gene_function = self.gene_function_calls_dict[gene_callers_id][source][1].split('!!!')[0]
                    frequency_of_categories[gene_function] += 1

            if split_name not in data_dict:
                data_dict[split_name] = {}
                for key in keys_list:
                    data_dict[split_name][key] = None

            for category in function_source_categories:
                data_dict[split_name][f"{source}!{category}"] = frequency_of_categories[category]

        for category in function_source_categories:
            keys_list.append(f"{source}!{category}")

        return data_dict, keys_list


    def init_split_sequences(self, min_contig_length=0, split_names_of_interest=set([])):
        if not len(split_names_of_interest):
            split_names_of_interest = self.split_names_of_interest

        contigs_shorter_than_M = self.init_contig_sequences(min_contig_length)

        if not len(self.splits_basic_info):
            self.run.info_single("Anvi'o was attempting to initialize split sequences, but the splits basic info dictionary "
                                 "was mysteriously empty. So you are warned.", mc="red")
            return

        self.progress.new('Computing split sequences from contigs')

        self.progress.update('Discarding split names coming from short contigs')
        split_names_to_discard = set([])
        for split_name in self.splits_basic_info:
            if self.splits_basic_info[split_name]['parent'] in contigs_shorter_than_M:
                split_names_to_discard.add(split_name)

        for split_name in split_names_to_discard:
            self.splits_basic_info.pop(split_name)

        if not len(self.splits_basic_info):
            self.progress.end()
            raise ConfigError("Something bad happened :/ The minimum length criterion of %d matched %d split names, and "
                              "removed all splits from the splits basic info dict. How could this happen? What have "
                              "you done?" % (min_contig_length, len(contigs_shorter_than_M)))

        # user asks for a specific set of splits to be initialized? maybe a better idea
        # is to set those names at a higher level in the contigs super, but for now this
        # will do it.
        if len(split_names_of_interest):
            missing_split_names = [s for s in split_names_of_interest if s not in self.splits_basic_info]
            if len(missing_split_names):
                self.progress.end()
                raise ConfigError("The `init_split_sequences` function was called with a set of split names of interest "
                                  "but %d of %d of those split names were missing from the splits basic info dict, which "
                                  "contained %d split names. Note that if you have been using a `min_contig_length` cutoff "
                                  "that may have resulted in the removal of your splits from the primary dict of splits. "
                                  "Regardless, here is one of the split names that you requested and were missing: '%s'. "
                                  "And here is one found in the splits basic info dict: '%s'." % \
                                                (len(missing_split_names), len(split_names_of_interest), len(self.splits_basic_info),
                                                 missing_split_names[0], list(self.splits_basic_info.keys())[0]))

            self.progress.end()
            if len(split_names_of_interest) != len(self.splits_basic_info):
                self.run.info_single(f"FYI: A subset of split sequences are being initialized (to be precise, only "
                                     f"{len(split_names_of_interest)} of {len(self.splits_basic_info)} splits the contigs database "
                                     f"knows about). Nothing to worry about. Probably.", mc="cyan", nl_after=1, nl_before=1)
            self.progress.new('Computing split sequences from contigs')
        else:
            split_names_of_interest = list(self.splits_basic_info.keys())

        self.progress.update('Generating split sequences dict')
        for split_name in split_names_of_interest:
            split = self.splits_basic_info[split_name]

            if split['parent'] in contigs_shorter_than_M:
                contigs_shorter_than_M.remove(split['parent'])
                continue

            if self.contigs_basic_info[split['parent']]['num_splits'] == 1:
                self.split_sequences[split_name] = {'sequence': self.contig_sequences[split['parent']]['sequence']}
            else:
                self.split_sequences[split_name] = {'sequence': self.contig_sequences[split['parent']]['sequence'][split['start']:split['end']]}

        self.progress.end()


    def init_non_singlecopy_gene_hmm_sources(self, split_names_of_interest=set([]), return_each_gene_as_a_layer=False):
        if not len(split_names_of_interest):
            split_names_of_interest = self.split_names_of_interest

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
        """Returns a tuple with 3 pieces of information for a given nucleotide position.

        This function accesses the self.nt_positions_info dictionary of arrays (each key is a contig
        name) to return the tuple: (in_noncoding_gene_call, in_coding_gene_call, base_pos_in_codon).

        Notes
        =====
        - If you plan on calling this function many times, consider instead `self.get_gene_info_for_each_position`
        """

        if (not self.a_meta['genes_are_called']) or \
           (not contig_name in self.nt_positions_info) or \
           (not len(self.nt_positions_info[contig_name])):
            return (0, 0, 0)

        if not self.nt_positions_info:
            raise ConfigError("get_nt_position_info :: I am asked to return stuff, but self.nt_positions_info is None!")

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
        """This method initializes a dictionary of function calls.

        It establishes the following attributes:
            self.gene_function_call_sources     a list of functional annotation sources
            self.gene_function_calls_dict       a dictionary of structure (accession, function, evalue) = self.gene_function_calls_dict[gene_callers_id][source]
                                                (the tuple can be None if there is no annotation from that source for the gene call)

        If requested_sources are provided, the dictionary only includes gene calls from those sources.
        If self.split_names_of_interest has a value, the dictionary only includes gene calls from those splits.

        Afterwards, it sets self.gene_function_calls_initiated to True.

        Note: the global argument RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE affects the behavior of this function. If False, we get 
        the best hit per gene (lowest e-value) for a given annotation source. If True, we get all hits.
        """
        if not self.contigs_db_path:
            return

        self.progress.new('Initializing functions class')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path)

        gene_function_sources_in_db = set(contigs_db.meta['gene_function_sources'] or [])

        where_clauses = []
        if requested_sources:
            self.check_functional_annotation_sources(requested_sources, dont_panic=dont_panic)
            self.gene_function_call_sources = requested_sources
            where_clauses.append('''source IN (%s)''' % (', '.join(["'%s'" % s for s in requested_sources])))
        else:
            self.gene_function_call_sources = gene_function_sources_in_db

        if self.split_names_of_interest:
            gene_caller_ids_of_interest = set(self.genes_in_contigs_dict.keys())
            where_clauses.append('''gene_callers_id IN (%s)''' % (', '.join([f"{g}" for g in gene_caller_ids_of_interest])))
        else:
            gene_caller_ids_of_interest = set([])

        if len(where_clauses):
            where_clause = ' AND '.join(where_clauses)
            hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
        else:
            hits = list(contigs_db.db.get_table_as_dict(t.gene_function_calls_table_name, error_if_no_data=False).values())

        for hit in hits:
            gene_callers_id = hit['gene_callers_id']
            source = hit['source']
            accession = hit['accession']
            function = hit['function']
            e_value = hit['e_value']

            if gene_callers_id not in self.gene_function_calls_dict:
                self.gene_function_calls_dict[gene_callers_id] = dict([(s, None) for s in self.gene_function_call_sources])

            entry = (accession, '%s' % (function if function else 'unknown'), e_value)

            if self.gene_function_calls_dict[gene_callers_id][source]:
                if anvio.RETURN_ALL_FUNCTIONS_FROM_SOURCE_FOR_EACH_GENE:
                    previous_entry_acc, previous_entry_func, previous_entry_evalue = self.gene_function_calls_dict[gene_callers_id][source]
                    combined_acc = f"{previous_entry_acc}!!!{accession}"
                    combined_func = f"{previous_entry_func}!!!{entry[1]}"
                    combined_evalue = f"{previous_entry_evalue}!!!{e_value}"
                    entry = (combined_acc, combined_func, combined_evalue)
                else:
                    if e_value and self.gene_function_calls_dict[gene_callers_id][source][2] < e_value:
                        # 'what we have:', self.gene_function_calls_dict[gene_callers_id][source]
                        # 'rejected    :', ('%s :: %s' % (function if function else 'unknown', accession), e_value)
                        continue

            self.gene_function_calls_dict[gene_callers_id][source] = entry

        contigs_db.disconnect()

        self.progress.end()

        self.gene_function_calls_initiated = True


    def list_function_sources(self):
        ContigsDatabase(self.contigs_db_path).list_function_sources()


    def check_functional_annotation_sources(self, sources=None, dont_panic=False):
        """Checks whether a given list of sources for functional annotation is valid.

           When `dont_panic` is True, quietly returns a list of sources that are matching to the
           ones in the database.
        """

        if not sources:
            return

        if sources and not isinstance(sources, list):
            raise ConfigError("Sources for functional annotations must be of type `list`.")

        contigs_db = ContigsDatabase(self.contigs_db_path, run=terminal.Run(verbose=False))
        gene_function_sources_in_db = contigs_db.meta['gene_function_sources']
        contigs_db.disconnect()

        missing_sources = [s for s in sources if s not in gene_function_sources_in_db]

        if len(missing_sources):
            if dont_panic:
                # quietly return matching sources
                return [s for s in sources if s in gene_function_sources_in_db]
            else:
                self.progress.reset()
                raise ConfigError("Some of the functional sources you requested are missing from the contigs database '%s'. Here "
                                  "they are (or here it is, whatever): %s." % \
                                                 (self.contigs_db_path, ', '.join(["'%s'" % s for s in missing_sources])))


    def search_for_gene_functions(self, search_terms, requested_sources=None, verbose=False, full_report=False, delimiter=',', case_sensitive=False, exact_match=False):
        if not isinstance(search_terms, list):
            raise ConfigError("Search terms must be of type 'list'")

        if requested_sources:
            if isinstance(requested_sources, str):
                requested_sources = list(set([source.strip() for source in requested_sources.split(delimiter)]))
            elif isinstance(requested_sources, list):
                pass
            else:
                raise ConfigError("Requested sources for annotations must be of type 'list' or 'str'")

        self.check_functional_annotation_sources(requested_sources)

        search_terms = [s.strip() for s in search_terms]

        if len([s.strip().lower() for s in search_terms]) != len(set([s.strip().lower() for s in search_terms])):
            raise ConfigError("Please do not use the same search term twice :/ Because, reasons. You know.")

        for search_term in search_terms:
            if not len(search_term) >= 3:
                raise ConfigError("A search term cannot be less than three characters")

        self.run.info('Search terms', f"{len(search_terms)} found: '{'|'.join(search_terms)}'")
        self.run.info("Case sensitive search?", "True" if case_sensitive else "False")
        self.run.info("Exact match?", "True" if exact_match else "False")
        matching_gene_caller_ids = dict([(search_term, {}) for search_term in search_terms])
        split_names = dict([(search_term, {}) for search_term in search_terms])
        full_report = []

        contigs_db = ContigsDatabase(self.contigs_db_path)

        for search_term in search_terms:
            self.progress.new('Search functions')
            self.progress.update('Searching for term "%s"' % search_term)

            if case_sensitive:
                contigs_db.db._exec('PRAGMA case_sensitive_like=ON;')
            else:
                contigs_db.db._exec('PRAGMA case_sensitive_like=OFF;')

            if exact_match:
                query = '''select gene_callers_id, source, accession, function from gene_functions where (function = "''' \
                                + search_term + '''" OR accession = "''' + search_term + '''")'''
            else:
                query = '''select gene_callers_id, source, accession, function from gene_functions where (function LIKE "%%''' \
                                + search_term + '''%%" OR accession LIKE "%%''' + search_term + '''%%")'''

            if requested_sources:
                query += ''' AND source IN (%s);''' % (', '.join(["'%s'" % s for s in requested_sources]))
            else:
                query += ';'

            response = contigs_db.db._exec(query).fetchall()

            # the response now contains all matching gene calls found in the contigs database. this may cause an issue
            # (just like the one reported here: https://github.com/merenlab/anvio/issues/1515) if the user is working
            # with only a subset of splits in the contigs database (for instance through `anvi-refine`). here we will
            # remove gene calls for which we don't have a split name associated:
            response = [r for r in response if r[0] in self.gene_callers_id_to_split_name_dict]

            # now we are good to go with extending the report
            full_report.extend([(r[0], r[1], r[2], r[3], search_term, self.gene_callers_id_to_split_name_dict[r[0]]) for r in response])

            matching_gene_caller_ids[search_term] = set([m[0] for m in response])
            split_names[search_term] = [self.gene_callers_id_to_split_name_dict[gene_callers_id] for gene_callers_id in matching_gene_caller_ids[search_term]]

            self.progress.end()

            if len(full_report):
                self.run.info('Matches', '%d unique genes contained the search term "%s"' % (len(set([e[0] for e in full_report if e[4] == search_term])), search_term), nl_before=1)
                if verbose:
                    observed_functions = set([])
                    self.run.warning('', header="Sneak peak into matching functions (up to 25)", raw=True, lc='cyan')
                    matches = [e for e in full_report if e[4] == search_term]
                    for entry in random.sample(matches, len(matches)):
                        if entry[3] in observed_functions:
                            continue
                        else:
                            self.run.info_single(f"{'; '.join(entry[3].split('!!!'))} :: {'; '.join(entry[2].split('!!!'))} within '{entry[1]}'", mc='cyan', cut_after=0)

                        observed_functions.add(entry[3]) 

                        if len(observed_functions) == 25:
                            break
            else:
                self.run.info('Matches', 'No matches found the search term "%s"' % (search_term), mc='red')

        contigs_db.disconnect()
        self.progress.end()

        return split_names, full_report


    def nt_position_to_gene_caller_id(self, contig_name, position_in_contig):
        """Returns the gene caller id for a given nucleotide position.

        Parameters
        ==========
        contig_name : str
        position_in_contig : int

        Returns
        =======
        gene_callers_id : int
            The gene call that covers the `position_in_contig` in a given contig.
            If the `position_in_contig` does not occur in any gene context in
            `contig_name`, this function will return `-1`
        """

        if not isinstance(position_in_contig, int):
            raise ConfigError("get_gene_caller_id_for_position_in_contig :: position_in_contig must be of type 'int'")

        if contig_name not in self.contigs_basic_info:
            raise ConfigError("Contig name '{contig_name} does not occur in this contigs-db :/'")

        if not self.a_meta['genes_are_called']:
            raise ConfigError(f"`get_gene_caller_id_for_position_in_contig` is speaking: You wanted to get back the gene call "
                              f"that occurs at the genome position {position_in_contig}. But it seems genes were not called "
                              f"for this contigs-db file, therefore that action is not one anvi'o can help you with :/")

        for gene_caller_id in self.genes_in_contigs_dict:
            gene_call = self.genes_in_contigs_dict[gene_caller_id]
            if gene_call['contig'] == contig_name and position_in_contig >= gene_call['start'] and position_in_contig < gene_call['stop']:
                return gene_caller_id

        return -1


    def get_corresponding_codon_order_in_gene(self, gene_caller_id, contig_name, pos_in_contig):
        """Returns the order of codon a given nucleotide belongs to.

        Parameters
        ==========
        gene_caller_id : int

        contig_name : str

        pos_in_contig : int

        Notes
        =====
        - If calling many times, consider `self.get_gene_info_for_each_position`
        """

        if not isinstance(pos_in_contig, int):
            raise ConfigError("get_corresponding_codon_order_in_gene :: pos_in_contig must be of type 'int'")

        if not isinstance(gene_caller_id, int):
            raise ConfigError("get_corresponding_codon_order_in_gene :: gene_caller_id must be of type 'int'")

        gene_call = self.genes_in_contigs_dict[gene_caller_id]

        if contig_name != gene_call['contig']:
            raise ConfigError("get_corresponding_codon_order_in_gene :: well, the gene call %d and the contig %s "
                              "do not seem to have anything to do with each other :/ This is not a user-level error "
                              "something must have gone very wrong somewhere in the code ..." % (gene_caller_id, contig_name))

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
        """Return a list of (gene_callers_id, start, stop) tuples for each gene occurring in contig_name"""
        return self.contig_name_to_genes[contig_name]


    def get_AA_counts_dict(self, split_names=set([]), contig_names=set([]), gene_caller_ids=set([]), return_codons_instead=False):
        """Returns a dictionary of AA counts.

        The dict can be returned for a given collection of split names, contigs names,
        or gene calls. If none of these variables are specified, the dict will contain
        counts for all gene calls in the contigs database
        """

        counts_dict = {}

        # nothing to do here if the genes were not called:
        if not self.a_meta['genes_are_called']:
            return counts_dict

        if len([True for v in [split_names, contig_names, gene_caller_ids] if v]) > 1:
            raise ConfigError("get_AA_counts_dict :: If you want to get AA counts for a specific "
                               "set of split names, contig names, or gene call ids, that is totally "
                               "fine. But you can't request more than one at a time.")

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

            if gene_call['call_type'] != constants.gene_call_types['CODING']:
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
        """For a given nucleotide position and contig name, returns all matching gene caller ids

        Notes
        =====
        - If you're calling this function many times, consider using
          self.get_gene_info_for_each_position
        """

        gene_start_stops_in_contig = self.get_gene_start_stops_in_contig(contig_name)

        if not gene_start_stops_in_contig:
            return []

        corresponding_gene_calls = [gene_callers_id
                                    for (gene_callers_id, start, stop) in gene_start_stops_in_contig
                                    if pos_in_contig >= start and pos_in_contig < stop]

        return corresponding_gene_calls


    def get_gene_info_for_each_position(self, contig_name, info='all'):
        """For a given contig, calculate per-position gene info

        Returns a dictionary of arrays, each with length equal to the contig length. Each key in the
        dictionary describes a different piece of gene info. By default, the dictionary has the
        following keys:

            'corresponding_gene_call' : To what gene_caller_id does the nt belong (-1 if there are 0
                                        or multiple gene calls)?
            'codon_order_in_gene'     : To which codon does the nt belong (0-indexed, -1 if
                                        corresponding_gene_call is -1)?
            'in_noncoding_gene_call'  : Does this position lie in a gene that is noncoding (0 or 1)?
            'in_coding_gene_call'     : Does this position lie in a gene that is coding (0 or 1)?
            'base_pos_in_codon'       : To what codon position (1, 2, or 3) does the nt belong (0 if
                                        corresponding_gene_call is -1, or gene does not have codons,
                                        e.g. ribosomal proteins)?
            'forward'                 : 1 if gene direction is forward, 0 if it is reverse
            'gene_start'              : Where in the contig does the gene start?
            'gene_stop'               : Where in the contig does the gene end?

        Parameters
        ==========
        contig_name : str

        info : list, 'all'
            A list of desired info names. By default, 'all' corresponds to
            ['corresponding_gene_call', 'codon_order_in_gene', 'in_noncoding_gene_call',
            'in_coding_gene_call', 'base_pos_in_codon', 'forward', 'gene_start', 'gene_stop'].
            Please note that this is just a convenience for the programmer: _all_ keys are
            calculated, and then only the requested subset is returned.

        Notes
        =====
        - If you are interested in just a few nt positions, use instead the "per-nucleotide"
          functions `get_nt_position_info`, `get_corresponding_gene_caller_ids_for_base_position`,
          and `get_corresponding_codon_order_in_gene`
        - This function gives per-nucleotide info about things that are potentially not
          per-nucleotide, e.g. "corresponding_gene_call". That means the output of this function is
          inherently redundant, and sometimes that's okay.
        """

        available_info = [
            'in_noncoding_gene_call',
            'in_coding_gene_call',
            'base_pos_in_codon',
            'corresponding_gene_call',
            'codon_order_in_gene',
            'forward',
            'gene_start',
            'gene_stop',
        ]

        if info == 'all':
            column_names = available_info
        else:
            for i in info:
                if i not in available_info:
                    raise ConfigError("get_gene_info_for_each_position :: %s is not an available choice for info" % i)
            column_names = info

        output = {}
        contig_length = len(self.contig_sequences[contig_name]['sequence'])
        data_shape = (contig_length, len(available_info))

        # Init the array. First 3 columns have defaults of 0, last 5 have defaults of -1
        data = -numpy.ones(data_shape).astype(int)
        data[:, :3] = 0

        # First, we populate the first 3 columns of data, 'in_noncoding_gene_call',
        # 'in_coding_gene_call', and 'base_pos_in_codon'. This is done straightforwardly by
        # accessing self.nt_positions_info

        if (not self.a_meta['genes_are_called']) or (not contig_name in self.nt_positions_info) or (not len(self.nt_positions_info[contig_name])):
            # In these cases everything gets 0
            pass
        else:
            data[self.nt_positions_info[contig_name] == 8, :3] = numpy.array([1,0,0])
            data[self.nt_positions_info[contig_name] == 4, :3] = numpy.array([0,1,1])
            data[self.nt_positions_info[contig_name] == 2, :3] = numpy.array([0,1,2])
            data[self.nt_positions_info[contig_name] == 1, :3] = numpy.array([0,1,3])

        # Next, we calculte the next 5 columns. As a first pass, we populate the splice of `data`
        # corresponding to each gene call and set the "gene_caller_id" and "codon_order_in_gene"
        # columns. This first ignores the fact that gene calls may overlap.

        gene_calls = self.get_gene_start_stops_in_contig(contig_name)

        for gene_caller_id, start, stop in gene_calls:
            positions = numpy.arange(start, stop)

            direction = self.genes_in_contigs_dict[gene_caller_id]['direction']

            if direction == 'r':
                codon_order_in_gene = (stop - start) / 3 - numpy.floor((positions - start) / 3) - 1
            else:
                codon_order_in_gene = numpy.floor((positions - start) / 3)

            data[start:stop, 3] = gene_caller_id
            data[start:stop, 4] = codon_order_in_gene
            data[start:stop, 5] = direction == 'f'
            data[start:stop, 6] = start
            data[start:stop, 7] = stop

        # Next, we compare each gene call to every other gene call. If they overlap, find the
        # overlapping region and set all columns to their defaults. This conservatively says, "if
        # there are two gene calls corresponding to a nucleotide position, anvi'o will simply say it
        # does not belong to any gene."

        gene_calls_to_compare = gene_calls.copy()

        for gene_call_1 in gene_calls:
            _, start1, stop1 = gene_call_1
            gene_calls_to_compare.remove(gene_call_1)

            for _, start2, stop2 in gene_calls_to_compare:
                if ((start1 < stop2  and stop1 > start2) or (stop1  > start2 and stop2 > start1)):
                    # There is overlap
                    overlap_start, overlap_end = max(start1, start2), min(stop1, stop2)
                    data[overlap_start:overlap_end, :3] = 0
                    data[overlap_start:overlap_end, 3:] = -1

        # Finally, we look for genes that have base_pos_in_codon == 0. These do not have
        # codon_order_in_genes, and so we must set them to -1
        data[data[:, 2] == 0, 4] = -1

        # Recast the requested info into `output` and return
        for i, c in enumerate(available_info):
            if c in column_names:
                output[c] = data[:, i]

        return output


    def get_gene_amino_acid_sequence(self, gene_caller_ids):
        """A much faster way to get back amino acid sequences for genes.

        Paremeters
        ==========
        gene_caller_ids : list
            A list of one or more gene caller ids.
        """

        if not isinstance(gene_caller_ids, list):
            raise ConfigError("Anvi'o is disappoint. Gene caller ids sent to this function must "
                              "be of type `list`.")

        contigs_db = ContigsDatabase(self.contigs_db_path)
        d = contigs_db.db.smart_get(t.gene_amino_acid_sequences_table_name, 'gene_callers_id', gene_caller_ids)
        contigs_db.disconnect()

        sequences = {}
        for gene_callers_id in gene_caller_ids:
            if gene_callers_id in d:
                sequences[gene_callers_id] = d[gene_callers_id]['sequence']
            else:
                sequences[gene_callers_id] = None

        return sequences


    def get_sequences_for_gene_callers_ids(self, gene_caller_ids_list=[], output_file_path=None, reverse_complement_if_necessary=True,
                                           include_aa_sequences=False, flank_length=0, output_file_path_external_gene_calls=None,
                                           simple_headers=False, list_defline_variables=False, defline_format='{gene_caller_id}',
                                           report_aa_sequences=False, wrap=120, rna_alphabet=False):

        ##################################################################################################
        #
        # DEFLIINE FORMATTING REPORTING RELATED PRE-CHECKS
        #
        ##################################################################################################
        # available options to determine deflines through user-provided f-strings. the dictionary is
        # populated below, and if you make any changes here, please don't forget to update it there too:
        defline_data_dict = {'gene_caller_id': None,
                             'contig_name': None,
                             'start': None,
                             'stop': None,
                             'direction': None,
                             'length': None,
                             'contigs_db_project_name': None}

        # if the user needs to see the list, show the list and quit
        if list_defline_variables:
            self.run.warning("Here are the variables you can use to provide a user-defined defline template: ")
            for key in defline_data_dict.keys():
                self.run.info_single("{%s}" % key)
            self.run.info_single("Remember, by default, anvi'o will only use '{gene_caller_id}' to format the deflines of "
                                 "FASTA files it produces.", level=0, nl_before=1, nl_after=1, mc='red')

            sys.exit()

        ##################################################################################################
        #
        # BUNCH OF SANITY CHECKS BEFORE WE GET INTO BUSINESS
        #
        ##################################################################################################

        if not isinstance(gene_caller_ids_list, list):
            raise ConfigError("Gene caller's ids must be of type 'list'")

        if not len(gene_caller_ids_list):
            gene_caller_ids_list = list(self.genes_in_contigs_dict.keys())
            self.run.warning("You did not provide any gene caller ids. As a result, anvi'o will give you back sequences for every "
                             "%d gene call stored in the contigs database. %s" % (len(gene_caller_ids_list), ' Brace yourself.' if len(gene_caller_ids_list) > 10000 else ''))

        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        filesnpaths.is_output_file_writable(output_file_path_external_gene_calls) if output_file_path_external_gene_calls else None

        if not isinstance(wrap, int):
            raise ConfigError('"wrap" has to be an integer instance')
        if wrap == 0:
            wrap = None
        if wrap and wrap <= 20:
            raise ConfigError('Value for wrap must be larger than 20. Yes. Rules.')

        try:
            gene_caller_ids_list = [int(gene_callers_id) for gene_callers_id in gene_caller_ids_list]
        except:
            raise ConfigError("List of IDs for gene calls contains non-integer values :/")

        if flank_length:
            try:
                int(flank_length)
            except ValueError:
                raise ConfigError("flank_length must be an integer.")

            if int(flank_length) < 0:
                raise ConfigError("flank_length must be a positive integer.")

            if flank_length and include_aa_sequences:
                raise ConfigError("You indicated --flank-length and --get-aa-sequences. anvi'o can only add flanking regions to nucleotide sequences.")

        if output_file_path_external_gene_calls and not simple_headers:
            raise ConfigError("If you are asking anvi'o to create an external gene calls file for your gene sequences, you can't also "
                              "also ask FASTA file headers for gene sequences to be not simple. External gene calls file and the FASTA "
                              "file must match, and anvi'o will have to take care of it without your supervision.")

        # if we came all the way down here without a defline format, let's set one up:
        if not defline_format:
            defline_format = "{gene_caller_id}"

        # we will also check if the `defline_format` is composed of variables that are defined in
        # the  `defline_data_dict` which is filled later
        utils.get_f_string_evaluated_by_dict(defline_format, defline_data_dict)

        ##################################################################################################
        #
        # BUSINESS TIME
        #
        ##################################################################################################

        # finally getting our sequences initialized. please NOTE that we do it only if there are no
        # contig sequences available OR if the gene caller ids of interest is not represented among
        # those that were previously initialized.
        if not len(self.contig_sequences) or not set(gene_caller_ids_list).issubset(self.gene_caller_ids_included_in_contig_sequences_initialized):
            self.init_contig_sequences(gene_caller_ids_of_interest=set(gene_caller_ids_list))

        if include_aa_sequences or report_aa_sequences or output_file_path_external_gene_calls:
            contigs_db = ContigsDatabase(self.contigs_db_path)
            aa_sequences_dict = contigs_db.db.get_table_as_dict(t.gene_amino_acid_sequences_table_name)
            contigs_db.disconnect()
        else:
            aa_sequences_dict = None

        # we will put everything we learn into this data structure.
        sequences_dict = {}

        # this will track whether we successfully recovered desired flanks for genes.
        gene_flanks_truncated = {'start': set([]), 'end': set([])}

        # ACTION
        self.progress.new('Working on sequences data structure')
        self.progress.update('...')
        for gene_callers_id in gene_caller_ids_list:
            gene_call = copy.deepcopy(self.genes_in_contigs_dict[gene_callers_id])

            contig_name = gene_call['contig']
            start, stop = gene_call['start'], gene_call['stop']

            if flank_length:
                # if the user asked for flanking sequences, then the start / stop positions
                # to excise the gene sequence with flanks from the contig sequence will
                # have to be updated. PLEASE NOTE THAT THESE START / STOP ARE NOT FOR
                # THE GENE CALL ANYMORE.
                contig_length = len(self.contig_sequences[contig_name]['sequence'])

                start = start - int(flank_length)
                if start < 0:
                    gene_flanks_truncated['start'].add(gene_callers_id)
                    start = 0

                stop = stop + int(flank_length)
                if stop > contig_length:
                    gene_flanks_truncated['end'].add(gene_callers_id)
                    stop = contig_length

            direction = gene_call['direction']
            sequence = self.contig_sequences[contig_name]['sequence'][start:stop]

            # NOTE: this is a bit ugly here. if there is an external gene calls output file path,
            # we simply do not reverse-complement reverse genes EVEN when `reverse_complement_if_necessary`
            # variable is True. and we do it quietly.
            if not output_file_path_external_gene_calls:
                if direction == 'r' and reverse_complement_if_necessary:
                    sequence = utils.rev_comp(sequence)
                    rev_compd = "True"
                else:
                    rev_compd = "False"
            else:
                rev_compd = "False"

            ###################################################################################
            #
            # UPDATING THE GENE CALL
            #
            ###################################################################################

            # if the user asked for flanking sequences, in the next conditional WE WILL UPDATE
            # GENE START/STOP POSITIONS this is a special case of reporting for specific applications.
            if flank_length:
                # update start / stop positions of the gene GIVEN the sequence we are reporting
                # since they are currently showing start/stop positions GIVEN the contig
                # they were on.
                original_gene_call_start = gene_call['start']
                original_gene_call_stop = gene_call['stop']
                sequence_with_flank_start = start
                gene_call['start'] = original_gene_call_start - sequence_with_flank_start
                gene_call['stop'] = original_gene_call_stop - original_gene_call_start + gene_call['start']

            # update the sequence
            gene_call['sequence'] = sequence
            gene_call['length'] = gene_call['stop'] - gene_call['start']
            gene_call['rev_compd'] = rev_compd

            if include_aa_sequences or report_aa_sequences or output_file_path_external_gene_calls:
                if gene_callers_id in aa_sequences_dict:
                    gene_call['aa_sequence'] = aa_sequences_dict[gene_callers_id]['sequence']
                else:
                    gene_call['aa_sequence'] = None

            # let's populate the dictionary that holds all the information that could be used to report
            # gene FASTA files. if you change anything in this dictionary, please don't forget to
            # update the list of variables where it is first defined in this function.
            defline_data_dict = {'gene_caller_id': gene_callers_id,
                                 'contig_name': gene_call['contig'],
                                 'start': gene_call['start'],
                                 'stop': gene_call['stop'],
                                 'direction': gene_call['direction'],
                                 'length': gene_call['length'],
                                 'contigs_db_project_name': self.a_meta['project_name_str']} 

            if output_file_path_external_gene_calls:
                # if the user is asking for an external gene calls file, the FASTA file for sequences
                # should not start with digits and we also need to set the contig name in sequences
                # dict properly so the FASTA file and the external gene calls file correspond to each other
                gene_call['header'] = '%s_%d' % (gene_call['contig'], gene_callers_id)
                gene_call['contig'] = gene_call['header']

                if not flank_length:
                    gene_call['start'] = 0
                    gene_call['stop'] = gene_call['length']
            else:
                gene_call['header'] = utils.get_f_string_evaluated_by_dict(defline_format, defline_data_dict)
                if not simple_headers:
                    gene_call['header'] += ' ' + ';'.join(['%s:%s' % (k, str(gene_call[k])) for k in ['contig', 'start', 'stop', 'direction', 'rev_compd', 'length']])

            # adding the updated gene call to our sequences dict.
            sequences_dict[gene_callers_id] = gene_call

        self.progress.end()

        if len(gene_flanks_truncated['start']) or len(gene_flanks_truncated['end']):
            missing_starts = f"{len(gene_flanks_truncated['start'])} genes were too close to the contig start to have the entire flank length at their beginning (gene caller ids: {', '.join([str(g) for g in gene_flanks_truncated['start']])}). "
            missing_ends = f"{len(gene_flanks_truncated['end'])} genes were too close to the contig end to have the entire flank lenght at their end (gene caller ids: {', '.join([str(g) for g in gene_flanks_truncated['end']])}). "

            msg = ""
            msg += missing_starts if len(gene_flanks_truncated['start']) else ""
            msg += missing_ends if len(gene_flanks_truncated['end']) else ""

            self.run.warning(f"While anvi'o was trying to add flanking regions of {flank_length} nucleotides to your gene sequences "
                             f"but things didn't go as smoothly as you may have hoped :/ {msg}")

        skipped_gene_calls = []

        # Storing the FASTA file for seqeunces
        if output_file_path:
            with open(output_file_path, 'w') as output:
                self.progress.new('Storing sequences')
                self.progress.update('...')
                for gene_callers_id in gene_caller_ids_list:
                    entry = sequences_dict[gene_callers_id]

                    sequence = None

                    if report_aa_sequences and rna_alphabet:
                        raise ConfigError("You can not request AA sequences repored in RNA alphabet.")
                    elif rna_alphabet:
                        sequence = entry['sequence'].replace('T', 'U')
                    elif report_aa_sequences:
                        sequence = entry['aa_sequence']
                    else:
                        sequence = entry['sequence']

                    if not sequence:
                        skipped_gene_calls.append(gene_callers_id)
                        continue

                    if wrap:
                        sequence = textwrap.fill(sequence, wrap, break_on_hyphens=False)

                    output.write('>%s\n' % entry['header'])
                    output.write('%s\n' % sequence)

                self.progress.end()

            self.run.info('Output FASTA', output_file_path)

        # Storing the external gene calls file for sequences stored as FASTA
        if output_file_path_external_gene_calls:
            utils.store_dict_as_TAB_delimited_file(sequences_dict,
                                                   output_file_path_external_gene_calls,
                                                   headers=['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version', 'aa_sequence'])
            self.run.info('Output external gene calls', output_file_path_external_gene_calls)

        if len(skipped_gene_calls):
            self.run.warning("Gene caller IDs %s have empty AA sequences and skipped." % (", ".join(map(str, skipped_gene_calls))))

        return (gene_caller_ids_list, sequences_dict)


    def gen_GFF3_file_of_sequences_for_gene_caller_ids(self, gene_caller_ids_list=[], output_file_path=None, wrap=120, simple_headers=False, rna_alphabet=False, gene_annotation_source=None):
        gene_caller_ids_list, sequences_dict = self.get_sequences_for_gene_callers_ids(gene_caller_ids_list)

        name_template = '' if simple_headers else ';Name={contig} {start} {stop} {direction} {rev_compd} {length}'

        # let's see the situation with functions
        if gene_annotation_source and not self.a_meta['gene_function_sources']:
            raise ConfigError("The code came all the way here with a `gene_annotation_source` but without the proper contigs db configuration o_O "
                              "Anvi'o needs an adult :(")

        if not gene_annotation_source and self.a_meta['gene_function_sources']:
            self.run.warning(f"You didn't ask for functions to be included in your GFF output file, but anvi'o found "
                             f"in your contigs database all these in case you would like to include one of them via "
                             f"the parameter `--annotation-source`: {', '.join(self.a_meta['gene_function_sources'])}",
                             header="THE MORE YOU KNOW 🌈", lc='yellow')

        gene_functions_found = False
        if gene_annotation_source and self.a_meta['gene_function_sources']:
            if gene_annotation_source in self.a_meta['gene_function_sources']:
                self.init_functions(requested_sources=[gene_annotation_source])
                gene_functions_found = True
                self.run.warning(f"Anvi'o found gene function annotations by '{gene_annotation_source}' in your contigs database "
                                 f"and will include that information in the output GFF file.", header="FUNCTIONS FOUND 🎊", lc="green")
            else:
                raise ConfigError(f"The function annotation source you have requested ('{gene_annotation_source}') does not "
                                  f"seem to be in this contigs database :/ Here is what we have: {', '.join(self.a_meta['gene_function_sources'])}")

        self.progress.new('Storing sequences')
        self.progress.update('...')
        with open(output_file_path, 'wt') as output:
            output.write('##gff-version 3\n')
            for gene_callers_id in gene_caller_ids_list:
                entry = sequences_dict[gene_callers_id]
                strand = entry['direction'].replace('f','+').replace('r','-')

                if entry['source'] == 'Transfer_RNAs':
                    seq_type = "tRNA"
                elif entry['source'].startswith('Ribosomal_RNA'):
                    seq_type = "rRNA"
                else:
                    seq_type = "CDS"

                entry_id = '___'.join([self.a_meta['project_name_str'], str(gene_callers_id)])
                attributes = f"ID={entry_id}"
                if gene_functions_found:
                    if gene_callers_id in self.gene_function_calls_dict:
                        accession_str, function_str, evalue = self.gene_function_calls_dict[gene_callers_id][gene_annotation_source]
                        accessions = accession_str.split('!!!')
                        functions = function_str.split('!!!')

                        # Clean functions for GFF3 compliance
                        cleaned_functions = []
                        for function in functions:
                            # now we have the function text, but it may contain characters that are offensive to the GFF3
                            # specifications, so we will simply replace them with space character. for more information
                            # please see https://github.com/merenlab/anvio/issues/2070
                            for character in [':', ';', ',', '&', '\t']:
                                function = function.replace(character, ' ')
                            cleaned_functions.append(function)

                        # Build attribute fields
                        name_field = ",".join(accessions) # Not strictly correct for gff3 as only meant to have 1 entry for Name
                        dbxref_field = f"{gene_annotation_source}:" + ",".join(accessions)
                        product_field = ",".join(cleaned_functions)

                        attributes += f";Name={name_field};db_xref={dbxref_field};product={product_field}"

                output.write(f"{entry['contig']}\t.\t{seq_type}\t{entry['start'] + 1}\t{entry['stop']}\t.\t{strand}\t.\t{attributes}")
                output.write(name_template.format(entry))
                output.write('\n')

        self.progress.end()
        self.run.info('Output', output_file_path)


    def gen_TAB_delimited_file_for_split_taxonomies(self, output_file_path):
        filesnpaths.is_output_file_writable(output_file_path)

        if not self.a_meta['gene_level_taxonomy_source']:
            raise ConfigError("There is no taxonomy source for genes in the contigs database :/")


        self.progress.new('Initializing splits taxonomy')
        self.progress.update('...')

        contigs_db = ContigsDatabase(self.contigs_db_path)
        splits_taxonomy_table = contigs_db.db.smart_get(t.splits_taxonomy_table_name, 'split', self.split_names_of_interest, string_the_key=True, error_if_no_data=False, progress=self.progress)
        taxon_names_table = contigs_db.db.get_table_as_dict(t.taxon_names_table_name)

        output = open(output_file_path, 'w')
        column_list = ['split_name'] + [k for k in taxon_names_table[list(taxon_names_table.keys())[0]].keys()]
        header = "\t".join(column_list)
        output.write(f"{header}\n")
        for split_name in self.splits_basic_info:
            if split_name in splits_taxonomy_table:
                taxon_id = splits_taxonomy_table[split_name]['taxon_id']
                if taxon_id:
                    self.splits_taxonomy_dict[split_name] = taxon_names_table[taxon_id]
                    taxonomy_string = "\t".join(str(x) for x in self.splits_taxonomy_dict[split_name].values())
                    output.write(f"{split_name}\t{taxonomy_string}\n")
                else:
                    output.write(f"{split_name}\t\n")
        output.close()

        contigs_db.disconnect()
        self.progress.end()
        self.run.info("Taxonomy", "Annotations for %d of %d total splits are recovered" % (len(splits_taxonomy_table), len(self.splits_basic_info)))
        self.run.info("Output", output_file_path)


class PanSuperclass(object):
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        self.run.width = 60

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.pan_db_path = A('pan_db')
        self.genomes_storage_path = A('genomes_storage')
        self.skip_init_functions = A('skip_init_functions')
        self.just_do_it = A('just_do_it')
        self.include_gc_identity_as_function = A('include_gc_identity_as_function')
        # these attributes involve determination of gene cluster function summaries
        self.discard_ties = A('discard_ties')
        self.consensus_threshold = A('consensus_threshold')

        self.genome_names = []
        self.gene_clusters = {}
        self.gene_clusters_initialized = False
        self.gene_cluster_names = set([])
        self.gene_cluster_names_in_db = set([])
        self.gene_clusters_gene_alignments = {}
        self.gene_clusters_gene_alignments_available = False
        # these two are initialized by self.init_gene_clusters_functions():
        self.gene_clusters_function_sources = []
        self.gene_clusters_functions_dict = {}
        # this one below is initialized by self.init_gene_cluster_functions_summary():
        self.gene_clusters_functions_summary_dict = {}
        self.gene_callers_id_to_gene_cluster = {}
        self.item_orders = {}
        self.views = {}
        self.collection_profile = {}

        # the following two are initialized via `init_items_additional_data()` and use information
        # stored in item additional data tables in the pan database
        self.items_additional_data_dict = None
        self.items_additional_data_keys = None

        # let's figure out whether this pan database has gene cluster homogeneity data available
        k = TableForItemAdditionalData(self.args).get_available_data_keys()
        self.functional_homogeneity_info_is_available = 'functional_homogeneity_index' in k
        self.geometric_homogeneity_info_is_available = 'geometric_homogeneity_index' in k
        self.combined_homogeneity_info_is_available = 'combined_homogeneity_index' in k

        self.num_gene_clusters = None
        self.num_genes_in_gene_clusters = None

        self.genomes_storage_is_available = False
        self.genomes_storage_has_functions = False
        self.functions_initialized = False

        if not self.pan_db_path:
            self.run.warning('PanSuperclass class called with args without a pan_db variable! Returning prematurely.')
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

        self.p_meta['available_item_orders'], self.item_orders = get_item_orders_from_db(self.pan_db_path)

        # recover all gene cluster names so others can access to this information
        # without having to initialize anything
        self.gene_cluster_names = set(pan_db.db.get_single_column_from_table(t.pan_gene_clusters_table_name, 'gene_cluster_id'))

        if not self.gene_cluster_names:
            raise ConfigError("You seem to have no gene clusters in this pan database :/ This is weird,\
                               sad, and curious at the same time. Probably you will have to go back to\
                               previous outputs of your worklow to make sure everything worked out properly.")

        pan_db.disconnect()

        # create an instance of states table
        self.states_table = TablesForStates(self.pan_db_path)

        self.progress.end()

        if self.genomes_storage_path:
            self.genomes_storage = genomestorage.GenomeStorage(self.genomes_storage_path,
                                                               self.p_meta['genomes_storage_hash'],
                                                               genome_names_to_focus=self.p_meta['genome_names'],
                                                               skip_init_functions=self.skip_init_functions,
                                                               run=self.run,
                                                               progress=self.progress)
            self.genomes_storage_is_available = True
            self.genomes_storage_has_functions = self.genomes_storage.functions_are_available
        else:
            self.run.warning("The pan database is being initialized without a genomes storage.")

        F = lambda x: '[YES]' if x else '[NO]'
        self.run.info('Pan DB', 'Initialized: %s (v. %s)' % (self.pan_db_path, anvio.__pan__version__))
        self.run.info('Gene cluster homogeneity estimates', 'Functional: %s; Geometric: %s; Combined: %s' % \
                         (F(self.functional_homogeneity_info_is_available), F(self.geometric_homogeneity_info_is_available),
                          F(self.combined_homogeneity_info_is_available)),
                      mc="cyan")


    def get_sequences_for_gene_clusters(self, gene_clusters_dict=None, gene_cluster_names=set([]), skip_alignments=False, report_DNA_sequences=False):
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

        You can call this function either with a `gene_clusters_dict`, or with a `gene_clusters_names` set. If you call it
        with the dict, it will operate on it. If you call it with names, it will use self.gene_clusters to find the names
        you specified. This looks like a shitty design, but was required to support exploratory / ad hoc user wishes through
        both command line and interactive anvi'o interfaces.

        By default, this function will report amino acid sequences. You can ask for DNA sequences if setting
        the flag `report_DNA_sequences` True.

        """

        if gene_clusters_dict and gene_cluster_names:
            raise ConfigError("OK. get_sequences_for_gene_clusters is speaking: You can call this function either with a "
                              "`gene_clusters_dict`, or with a `gene_clusters_names` set. If you call it with the dict, "
                              "it will operate on it. If you call it with names, it will use self.gene_clusters to find the names "
                              "you specified. This looks like a shitty design, but was required to support exploratory / ad hoc user wishes through "
                              "both command line and interactive anvi'o interfaces.")

        if not skip_alignments and self.gene_clusters_gene_alignments_available and report_DNA_sequences:
            if self.just_do_it:
                self.run.warning("Since you are using the flag `--just-do-it`, anvi'o will attempt to first access to the DNA sequences of gene calls "
                                 "found in gene clusters of interst, and align them from scratch instead of using the previously computed amino acid "
                                 "sequence alignment summary data (related info: https://github.com/merenlab/anvio/issues/772). Anvi'o apologizes in "
                                 "advance if this explodes downstream.", header="MINI WARNING AS YOU ARE GETTING SEQUENCES FOR GENE CLUSTERS", lc="yellow")
            else:
                self.run.warning("Anvi'o is getting sequences for the gene clusters you are interested in. While anvi'o aligned amino acid sequences "
                                 "for these gene clusters, you have requested DNA sequences. Anvi'o stores alignment information in a compressed form "
                                 "(which we call amino acid sequence alignment summary) which can also align DNA sequences instantaneously. BUT, due to "
                                 "intricacies associated with gene calling, the amino acid sequence of a gene stored in the contigs database MAY differ "
                                 "from its DNA sequence (true story). For those VERY RARE instances, the alignment summary for the amino acid sequence "
                                 "may no longer be used to make sense of the DNA sequence (see https://github.com/merenlab/anvio/issues/772 for more "
                                 "information). One way to address that is to do another alignment on the fly for the DNA sequences for genes found in "
                                 "gene clusters. But you need ot ask for that specifically, which you can do by including `--just-do-it` flag in your "
                                 "command line and re-run this program.", header="MINI WARNING AS YOU ARE GETTING SEQUENCES FOR GENE CLUSTERS", lc="yellow")
                skip_alignments = True

        sequences = {}

        if not gene_cluster_names and not gene_clusters_dict:
            raise ConfigError("get_sequences_for_gene_clusters is speaking: You must call this function either with a `gene_clusters_dict` "
                              "or with a `gene_cluster_names` set.")

        if not gene_cluster_names:
            gene_cluster_names = set(list(gene_clusters_dict.keys()))

        if not isinstance(gene_cluster_names, type(set([]))) or not gene_cluster_names:
            raise ConfigError("gene_cluster_names for get_sequences_for_gene_clusters must be a non-empty `set`.")

        if not self.genomes_storage_is_available:
            raise ConfigError("The pan anvi'o super class for is upset. You are attempting to get AA sequences for %s,\
                               but there is not genomes storage is available to get it." \
                                    % 'a gene cluster' if len(gene_cluster_names) > 1 else '%d gene_clusters' % len(gene_cluster_names))

        if gene_clusters_dict is None:
            if not self.gene_clusters_initialized:
                self.init_gene_clusters()
            gene_clusters_dict = self.gene_clusters

        missing_gene_cluster_names = [p for p in gene_cluster_names if p not in gene_clusters_dict]
        if len(missing_gene_cluster_names[0:5]):
            raise ConfigError("get_sequences_for_gene_clusters: %d of %d gene clusters are missing in your data. Not good :/ "
                              "Here are some of the missing ones; %s" \
                                        % (len(missing_gene_cluster_names), len(gene_cluster_names), ', '.join(missing_gene_cluster_names[0:5])))

        self.progress.new('Accessing gene cluster sequences', progress_total_items=len(gene_cluster_names))

        for gene_cluster_name in gene_cluster_names:
            self.progress.increment()
            self.progress.update("processing '%s' ..." % gene_cluster_name )
            sequences[gene_cluster_name] = {}
            for genome_name in gene_clusters_dict[gene_cluster_name]:
                sequences[gene_cluster_name][genome_name] = {}
                for gene_callers_id in gene_clusters_dict[gene_cluster_name][genome_name]:
                    sequence = self.genomes_storage.get_gene_sequence(genome_name, gene_callers_id, report_DNA_sequences=report_DNA_sequences)

                    if not skip_alignments and self.gene_clusters_gene_alignments_available:
                        alignment_summary = self.gene_clusters_gene_alignments[genome_name][gene_callers_id]
                        sequence = utils.restore_alignment(sequence, alignment_summary, from_aa_alignment_summary_to_dna=report_DNA_sequences)

                    sequences[gene_cluster_name][genome_name][gene_callers_id] = sequence

        self.progress.end()

        return sequences


    def compute_AAI_for_gene_cluster(self, gene_cluster):
        """Simple AAI calculator for sequences in a gene cluster"""

        def calculate_sequence_identity(s1, s2):
            matches = sum(r1 == r2 for r1, r2 in zip(s1, s2))
            identity = matches / len(s1)
            return identity
    
        # turn the sequences dict into a more useful format for this
        # step
        d = {}
        for gene_cluster_name in gene_cluster:
            for genome_name in gene_cluster[gene_cluster_name]:
                for gene_call in gene_cluster[gene_cluster_name][genome_name]:
                    d[f"{genome_name}_{gene_call}"] = gene_cluster[gene_cluster_name][genome_name][gene_call]
    
        # if there is only one sequence, then there is nothing to do here.
        if len(d) == 1:
            return (0.0, 0.0, 0.0)
    
        # get all pairs of sequences
        pairs_of_sequences = list(itertools.combinations(d.keys(), 2))
    
        # FIXME: we need a smart strategy here to deal with gaps, but meren is old and it is 11pm.

        # calculate pairwise identities
        identities = []
        for s1, s2 in pairs_of_sequences:
            identity = calculate_sequence_identity(d[s1], d[s2])
            identities.append(identity)
    
        # calculate AAI values
        AAI_min = min(identities)
        AAI_max = max(identities)
        AAI_avg = numpy.mean(identities)

        return (AAI_min, AAI_max, AAI_avg)


    def compute_homogeneity_indices_for_gene_clusters(self, gene_cluster_names=set([]), gene_clusters_failed_to_align=set([]), num_threads=1):
        if gene_cluster_names is None:
            self.run.warning("The function `compute_homogeneity_indices_for_gene_clusters` did not receive any gene "
                             "cluster names to work with. If you are a programmer, you should know that you are "
                             "doing it wrong. If you are a user, please get in touch with a programmer because this "
                             "is not normal. This function will now return prematurely without computing anything :(")
            return None

        if self.args.quick_homogeneity:
            self.run.warning("Performing quick homogeneity calculations (skipping horizontal geometric calculations) "
                             "per the '--quick-homogeneity' flag")

        sequences = self.get_sequences_for_gene_clusters(gene_cluster_names=gene_cluster_names, skip_alignments=False)

        homogeneity_calculator = homogeneityindex.HomogeneityCalculator(quick_homogeneity=self.args.quick_homogeneity)

        AAI_calculator = self.compute_AAI_for_gene_cluster

        self.progress.new('Computing gene cluster homogeneity indices', progress_total_items=len(gene_cluster_names))
        self.progress.update('Initializing %d threads...' % num_threads)

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        for gene_cluster_name in gene_cluster_names:
            input_queue.put(gene_cluster_name)

        workers = []
        for i in range(num_threads):
            worker = multiprocessing.Process(target=PanSuperclass.homogeneity_worker,
                                             args=(input_queue, output_queue, sequences, gene_clusters_failed_to_align, homogeneity_calculator, AAI_calculator, self.run))
            workers.append(worker)
            worker.start()

        results_dict = {}
        received_gene_clusters = 0
        while received_gene_clusters < len(sequences):
            try:
                homogeneity_dict = output_queue.get()
                if homogeneity_dict:
                    results_dict[homogeneity_dict['gene cluster']] = {'functional_homogeneity_index': homogeneity_dict['functional'],
                                                                      'geometric_homogeneity_index': homogeneity_dict['geometric'],
                                                                      'combined_homogeneity_index': homogeneity_dict['combined'],
                                                                      'AAI_min': homogeneity_dict['AAI_min'],
                                                                      'AAI_max': homogeneity_dict['AAI_max'],
                                                                      'AAI_avg': homogeneity_dict['AAI_avg']}

                received_gene_clusters += 1
                self.progress.increment(increment_to=received_gene_clusters)
                self.progress.update('Processed %d gene clusters using %d threads' % (received_gene_clusters, num_threads))
            except KeyboardInterrupt:
                print("Recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        self.progress.end()

        return results_dict


    @staticmethod
    def homogeneity_worker(input_queue, output_queue, gene_clusters_dict, gene_clusters_failed_to_align, homogeneity_calculator, AAI_calculator, run):
        r = terminal.Run()
        r.verbose = False

        while True:
            gene_cluster_name = input_queue.get(True)

            funct_index = {}
            geo_index = {}

            gene_cluster = {gene_cluster_name: gene_clusters_dict[gene_cluster_name]}
            indices_dict = {'gene cluster': gene_cluster_name}

            try:
                funct_index, geo_index, combined_index = homogeneity_calculator.get_homogeneity_dicts(gene_cluster)
                indices_dict['functional'] = funct_index[gene_cluster_name]
                indices_dict['geometric'] = geo_index[gene_cluster_name]
                indices_dict['combined'] = combined_index[gene_cluster_name]

                indices_dict['AAI_min'], indices_dict['AAI_max'], indices_dict['AAI_avg'] = AAI_calculator(gene_cluster)

            except:
                if gene_cluster_name not in str(gene_clusters_failed_to_align):
                    progress.reset()
                    run.warning(f"Homogeneity indices computation for gene cluster '{gene_cluster_name}' failed. This can happen due to one of "
                                f"three reasons: (1) this gene cluster is named incorrectly, does not exist in the database, or is formatted "
                                f"into the input dictionary incorrectly, (2) there is an alignment mistake in the gene cluster, and not all "
                                f" genes are alignedto be the same lenght; or (3) the homogeneity calculator was initialized incorrectly. As "
                                f"you can see, this is a rare circumstance, and anvi'o will set this gene cluster's homogeneity indices to "
                                f"`-1` so things can move on, but we highly recommend you to take a look at your data to make sure you are "
                                f"satisfied with your analysis.", overwrite_verbose=True)

                indices_dict['functional'] = -1
                indices_dict['geometric'] = -1
                indices_dict['combined'] = -1
                indices_dict['AAI_min'] = -1
                indices_dict['AAI_max'] = -1
                indices_dict['AAI_avg'] = -1

            output_queue.put(indices_dict)


    def write_sequences_in_gene_clusters_to_file(self, gene_clusters_dict=None, gene_cluster_names=set([]), \
                                                  skip_alignments=False, output_file_path=None, report_DNA_sequences=False):
        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        output_file = open(output_file_path, 'w')
        sequences_dict = self.get_sequences_for_gene_clusters(gene_clusters_dict=gene_clusters_dict,
                                                              gene_cluster_names=gene_cluster_names,
                                                              skip_alignments=skip_alignments,
                                                              report_DNA_sequences=report_DNA_sequences)

        self.progress.new('Writing gene cluster sequences to file')
        sequence_counter = 0
        for gene_cluster_name in sequences_dict:
            for genome_name in sequences_dict[gene_cluster_name]:
                for gene_callers_id in sequences_dict[gene_cluster_name][genome_name]:
                        output_file.write('>%08d|gene_cluster:%s|genome_name:%s|gene_callers_id:%d\n' % (sequence_counter,
                                                                                                         gene_cluster_name,
                                                                                                         genome_name,
                                                                                                         gene_callers_id))
                        output_file.write('%s\n' % sequences_dict[gene_cluster_name][genome_name][gene_callers_id])
                        sequence_counter += 1
                        self.progress.update("processing '%s' ..." % gene_cluster_name)

        self.progress.end()
        output_file.close()

        if len(gene_clusters_dict) == 1:
            gene_cluster_name = list(gene_clusters_dict.keys())[0]
            self.run.info('Gene cluster name', gene_cluster_name)
        self.run.info('Sequence type', 'DNA' if report_DNA_sequences else 'Amino acid')
        self.run.info('Num sequences reported', sequence_counter)
        self.run.info('Output FASTA file', output_file_path, mc='green', nl_after=1)


    def write_sequences_in_gene_clusters_for_phylogenomics(self, gene_clusters_dict=None, skip_alignments=False, \
                                                output_file_path=None, report_DNA_sequences=False, align_with=None, \
                                                separator=None, partition_file_path=None):
        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        if partition_file_path:
            filesnpaths.is_output_file_writable(partition_file_path)

        if not separator:
            separator = 'NNN' if report_DNA_sequences else 'XXX'

        output_file = open(output_file_path, 'w')
        sequences_dict = self.get_sequences_for_gene_clusters(gene_clusters_dict=gene_clusters_dict,
                                                              skip_alignments=skip_alignments,
                                                              report_DNA_sequences=report_DNA_sequences)

        if not self.gene_clusters_gene_alignments_available:
            aligner = aligners.select(align_with)

            run.warning("It seems sequences in gene clusters were not aligned during the pangenomic analysis, so we "
                        "are going to have do it now .. which may take some time .. and it is totally your fault :/")
            progress.new("Aligning sequences")
        elif align_with:
            run.warning("Your gene clusters are already aligned, yet you are asking for them to be aligned with "
                        "'%s' :( If you know what's going on (i.e. you are here because you run a command and "
                        "used the '--align-with' parameter or something), here anvi'o lets you know that it will "
                        "not use '%s' becase things are already aligned." % (align_with, align_with))

        get_first_value = lambda x: next(iter(x.values()))
        get_first_key = lambda x: next(iter(x.keys()))

        silent_run = terminal.Run()
        silent_run.verbose = False

        output_buffer = dict({})
        for genome_name in self.genome_names:
            output_buffer[genome_name] = StringIO()

        gene_cluster_representative_seq_lengths = []
        gene_cluster_names = list(sequences_dict.keys())
        for gene_cluster_name in gene_cluster_names:
            multiple_gene_calls = False
            multiple_gene_call_genome = None
            sequence_length = None

            for genome_name in self.genome_names:
                if len(sequences_dict[gene_cluster_name][genome_name]) > 1:
                    multiple_gene_calls = True
                    multiple_gene_call_genome = genome_name
                elif self.gene_clusters_gene_alignments_available and len(sequences_dict[gene_cluster_name][genome_name]) == 1:
                    sequence_length = len(get_first_value(sequences_dict[gene_cluster_name][genome_name]))

            if multiple_gene_calls:
                raise ConfigError("There are multiple gene calls in '%s' and sample '%s', which is not appropriate for phylogenomic "
                                  "analyses. Please use advanced filters (see help if you are not sure what this means) to remove "
                                  "gene clusters from your analysis if they contain multiple gene calls from any "
                                  "of the genomes in your pan (not to tell you what to do, but '--max-num-genes-from-each-genome 1' "
                                  "would make sure gene clusters that contain multiple genes from a given genome would be removed "
                                  "from your final list)." % (gene_cluster_name, multiple_gene_call_genome))

            if not self.gene_clusters_gene_alignments_available:
                sequences_to_align = []
                for genome_name in self.genome_names:
                    if len(sequences_dict[gene_cluster_name][genome_name]) == 1:
                        sequences_to_align.append((genome_name, get_first_value(sequences_dict[gene_cluster_name][genome_name])))

                progress.update("Processing '" + gene_cluster_name + "'")

                aligned_sequences = aligner(run=silent_run).run_stdin(sequences_list=sequences_to_align)

                for genome_name in aligned_sequences:
                    gene_caller_id = get_first_key(sequences_dict[gene_cluster_name][genome_name])
                    sequences_dict[gene_cluster_name][genome_name][gene_caller_id] = aligned_sequences[genome_name]

                    if not sequence_length:
                        sequence_length = len(aligned_sequences[genome_name])

            gene_cluster_representative_seq_lengths.append((gene_cluster_name, sequence_length),)

            for genome_name in self.genome_names:
                if len(sequences_dict[gene_cluster_name][genome_name]) == 1:
                    sequence = get_first_value(sequences_dict[gene_cluster_name][genome_name])
                else:
                    sequence = "-" * sequence_length

                output_buffer[genome_name].write(sequence)

                if not gene_cluster_name == gene_cluster_names[-1]:
                    output_buffer[genome_name].write(separator)

        # see https://github.com/merenlab/anvio/issues/1333
        if partition_file_path:
            utils.gen_NEXUS_format_partition_file_for_phylogenomics(partition_file_path,
                                                                    gene_cluster_representative_seq_lengths,
                                                                    separator,
                                                                    run=self.run,
                                                                    progress=self.progress)

        if not self.gene_clusters_gene_alignments_available:
            progress.end()

        for genome_name in self.genome_names:
            output_file.write('>%s gene_clusters:%s|separator:%s\n' % (genome_name, ','.join(gene_cluster_names), separator))
            output_file.write(output_buffer[genome_name].getvalue())
            output_file.write('\n')
            output_buffer[genome_name].close()

        output_file.close()

        self.run.info('Sequence type', 'DNA' if report_DNA_sequences else 'Amino acid', mc='green')
        self.run.info('Output file for phylogenomics', output_file_path, mc='green')


    def get_gene_cluster_function_summary(
        self,
        gene_cluster_id,
        functional_annotation_source,
        discard_ties: bool = False,
        consensus_threshold: float = None
    ):
        """
        Returns the most frequently occurring functional annotation across genes in a gene cluster
        for a given functional annotation source.

        A single function and accession number is determined purely based on frequencies of
        occurrence. If multiple most frequent functions are tied, then, by default, one is chosen
        arbitrarily.

        A threshold can be set for the proportion of the genes that must contain the most frequent
        annotation.

        Parameters
        ==========
        gene_cluster_id : str
            A gene cluster id that is defined in the `gene_clusters_functions_dict`

        functional_annotation_source : str
            A known functional annotation source

        discard_ties : bool
            If multiple annotations are most frequent among genes in a cluster, then do not assign
            an annotation to the cluster itself when this argument is True. By default, this
            argument is False, so one of the most frequent annotations would be arbitrarily chosen.

        consensus_threshold : float
            Without this argument (default None), the annotation most frequent among genes in a
            cluster is assigned to the cluster itself. With this argument (a value on [0, 1]), at
            least this proportion of genes in the cluster must have the most frequent annotation for
            the cluster to be annotated.

        Returns
        =======
        (most_common_accession, most_common_function) : tuple
            The representative function and its accession id for `gene_cluster_id` and
            `functional_annotation_source`
        """

        if not self.functions_initialized:
            raise ConfigError("Although it is an expensive step, this function currently requires functions to be "
                              "initialized first :/ If you are a programmer and need this functionality to be much "
                              "more effective for ad hoc use, please let us know (the right solution in that case is "
                              "to work directly with the genome storage database to recover functions for a single gene "
                              "cluster).")

        functions_counter = Counter({})
        num_annotated_genes = 0
        for genome_name in self.gene_clusters_functions_dict[gene_cluster_id]:
            for gene_caller_id in self.gene_clusters_functions_dict[gene_cluster_id][genome_name]:
                if functional_annotation_source in self.gene_clusters_functions_dict[gene_cluster_id][genome_name][gene_caller_id]:
                    annotation_blob = self.gene_clusters_functions_dict[gene_cluster_id][genome_name][gene_caller_id][functional_annotation_source]
                    functions_counter[annotation_blob] += 1
                    num_annotated_genes += 1

        if not len(functions_counter):
            return (None, None)

        functions_counter = functions_counter.most_common()

        if consensus_threshold and functions_counter[0][1] / num_annotated_genes < consensus_threshold:
            return (None, None)

        max_freq = 0
        candidates = []
        for annotation_blob, freq in functions_counter:
            if freq < max_freq:
                break
            max_freq = freq
            candidates.append(annotation_blob)

        if discard_ties and len(candidates) > 1:
            return (None, None)

        most_common_accession, most_common_function = functions_counter[0][0].split('|||')

        return (most_common_accession, most_common_function)


    def get_gene_clusters_functions_summary_dict(self, functional_annotation_source):
        """Returns a dictionary where each gene cluster is associated with a single function.

           See `get_gene_cluster_function_summary` for details since this function is merely
           calling it for each gene cluster.
        """

        if functional_annotation_source not in self.gene_clusters_function_sources:
            raise ConfigError("Your favorite functional annotation source '%s' does not seem to be among one of the sources "
                              "that are available to you. Here are the ones you should choose from: %s." % (functional_annotation_source, ', '.join(self.gene_clusters_function_sources)))

        if not self.functions_initialized:
            self.init_gene_clusters_functions()

        if not len(self.gene_clusters_functions_dict):
            raise ConfigError("The gene clusters functions dict seems to be empty. We assume this error makes "
                              "zero sense to you, and it probably will not help you to know that it also makes "
                              "zero sense to anvi'o too :/ Maybe you forgot to provide a genomes storage?")

        gene_clusters_functions_summary_dict = {}

        self.progress.new(f'Summarizing "{functional_annotation_source}" for gene clusters')
        self.progress.update('Creating a dictionary')
        for gene_cluster in self.gene_clusters_functions_dict:
            accession, function = self.get_gene_cluster_function_summary(gene_cluster, functional_annotation_source, discard_ties=self.discard_ties, consensus_threshold=self.consensus_threshold)
            gene_clusters_functions_summary_dict[gene_cluster] = {'gene_cluster_function': function, 'gene_cluster_function_accession': accession}

        self.progress.end()

        return gene_clusters_functions_summary_dict


    def init_gene_clusters_functions_summary_dict(self, source_list = None, gene_clusters_of_interest = None):
        """ A function to initialize the `gene_clusters_functions_summary_dict` by calling
            the atomic function `get_gene_clusters_functions_summary_dict` for all the
            functional annotation sources.

            You can restrict which sources to use and which gene clusters to initialize using the `source_list` and
            `gene_clusters_of_interest` parameters, respectively. Use this ability with caution as it is possible that
            downstream code may expect all sources and/or all gene clusters to be initialized.
        """

        if not self.functions_initialized:
            self.init_gene_clusters_functions()

        if not len(self.gene_clusters_functions_dict):
            self.run.warning("Someone asked anvi'o to initialize a gene cluster functions summary dict, but it seems there "
                             "are no gene cluster functions even after initializing functions for the pangenome. So we move "
                             "on without any summary dict for functions and/or drama about it to let the downstream analyses "
                             "decide how to punish the unlucky.")
            return
        
        gene_clusters_to_init = self.gene_clusters_functions_dict.keys()
        if gene_clusters_of_interest:
            gene_clusters_to_init = gene_clusters_of_interest

            # make sure that all of the requested gene clusters exist
            for cluster_id in gene_clusters_of_interest:
                if cluster_id not in self.gene_clusters:
                    raise ConfigError(f"Someone requested the function `init_gene_clusters_functions_summary_dict()` "
                                      f"to work on a gene cluster called {cluster_id} but it does not exist in the "
                                      f"pangenome. :/")

        sources_to_summarize = self.gene_clusters_function_sources
        if source_list:
            sources_to_summarize = source_list

            # make sure that all of the requested sources are in the pangenome
            for s in source_list:
                if s not in self.gene_clusters_function_sources:
                    raise ConfigError(f"Someone requested the function `init_gene_clusters_functions_summary_dict()` "
                                      f"to use the annotation source {s}, but that source is not found in the pangenome.")

        self.progress.new('Generating a gene cluster functions summary dict', progress_total_items=len(gene_clusters_to_init))
        counter = 0
        for gene_cluster_id in gene_clusters_to_init:
            if counter % 100 == 0:
                self.progress.increment(increment_to=counter)
                self.progress.update(f'{gene_cluster_id} ...')

            self.gene_clusters_functions_summary_dict[gene_cluster_id] = {}

            for functional_annotation_source in sources_to_summarize:
                accession, function = self.get_gene_cluster_function_summary(gene_cluster_id, functional_annotation_source, discard_ties=self.discard_ties, consensus_threshold=self.consensus_threshold)
                self.gene_clusters_functions_summary_dict[gene_cluster_id][functional_annotation_source] = {'function': function, 'accession': accession}

            counter += 1

        self.progress.end()


    def init_gene_clusters_functions(self):
        """This function initializes a dictionary of functions for a pangenome.

        It establishes the following attributes:
            self.gene_clusters_function_sources     a set() of functional annotation sources
            self.gene_clusters_functions_dict       a 3-level dictionary with the structure
                                                    "accession|||function" = self.gene_clusters_functions_dict[gene_cluster_id][genome_name][gene_callers_id][source]

        Afterwards it sets self.functions_initialized to True.
        """

        if not self.genomes_storage_is_available:
            self.run.warning("Someone tried to initialize gene cluster functions, but it seems there is no genomes "
                             "storage available to this run. That's OK. But no gene clusters functions for you "
                             "obviously.")
            return

        self.progress.new('Initializing functions for gene clusters')
        self.progress.update('...')
        if not self.gene_clusters:
            raise ConfigError("init_gene_clusters_functions is speaking! You called this function before you initialized "
                               "gene clusters :/ One of us does not know what they're doing :(")

        if not self.genomes_storage_has_functions:
            self.progress.end()
            self.run.warning("Genomes storage does not have any info about gene functions. Certain parts of the pangenomic "
                             "workflow will not be accessible.")
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

        if self.include_gc_identity_as_function:
            self.progress.reset()
            self.run.info_single("Gene cluster identities are being added as functions into the functions dictionary. "
                                 "Functional annotation resources will include `IDENTITY` as an option. See here why "
                                 "(apart from the fact that you asked for it by using the flag `--include-gc-identity-as-function`): "
                                 "https://github.com/merenlab/anvio/issues/1196", nl_after=1, mc='green')
            for gene_cluster_id in self.gene_clusters:
                for genome_name in self.genome_names:
                    for gene_callers_id in self.gene_clusters[gene_cluster_id][genome_name]:
                        self.gene_clusters_functions_dict[gene_cluster_id][genome_name][gene_callers_id]['IDENTITY'] = '%s|||%s' % (gene_cluster_id, gene_cluster_id)
            self.gene_clusters_function_sources.update(['IDENTITY'])

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


    def get_all_genome_names_in_gene_clusters_dict(self, gene_clusters_dict):
        """Returns all genome names found in a `gene_clusters_dict`"""

        all_genomes = set([])

        for entry in gene_clusters_dict.values():
            for genome_name in entry:
                all_genomes.add(genome_name)

        return all_genomes


    def get_gene_clusters_in_genomes_dict(self, gene_clusters_dict):
        """Goes through the `gene_clusters_dict and returns gene clusters in genomes dict
           as well as all genome names"""

        all_genomes = self.get_all_genome_names_in_gene_clusters_dict(gene_clusters_dict)

        gene_clusters_in_genomes_dict = dict([(genome_name, set([])) for genome_name in all_genomes])

        for genome_name in all_genomes:
            gene_clusters_in_genomes_dict[genome_name] = set([])

        for gene_cluster_name in gene_clusters_dict:
            for genome_name in gene_clusters_dict[gene_cluster_name]:
                if len(gene_clusters_dict[gene_cluster_name][genome_name]):
                    gene_clusters_in_genomes_dict[genome_name].add(gene_cluster_name)

        return gene_clusters_in_genomes_dict


    def get_num_genes_contributed_per_genome_dict(self, gene_clusters_dict):
        """Get a dictionary of gene cluster names and the number of genes contributed to each per genome"""


    def get_basic_gene_clusters_stats(self, gene_clusters_dict):
        """Returns two dictionaries: a dictionary of gene cluster names and their number of occurrences across all genomes,
           and a dictionary for number of genes contributed by each genome into each gene cluster"""

        all_gene_clusters = set(list(gene_clusters_dict.keys()))
        all_genomes = self.get_all_genome_names_in_gene_clusters_dict(gene_clusters_dict)

        num_occurrences_accross_genomes = dict([(gene_cluster_name, set([])) for gene_cluster_name in all_gene_clusters])
        num_genes_contributed_per_genome = dict([(gene_cluster_name, dict()) for gene_cluster_name in all_gene_clusters])

        for gene_cluster_name in num_genes_contributed_per_genome:
            for genome_name in all_genomes:
                num_genes_contributed_per_genome[gene_cluster_name][genome_name] = 0

        for gene_cluster_name in gene_clusters_dict:
            for genome_name in gene_clusters_dict[gene_cluster_name]:
                if len(gene_clusters_dict[gene_cluster_name][genome_name]):
                    num_occurrences_accross_genomes[gene_cluster_name].add(genome_name)
                    num_genes_contributed_per_genome[gene_cluster_name][genome_name] = len(gene_clusters_dict[gene_cluster_name][genome_name])

        for gene_cluster_name in num_occurrences_accross_genomes:
            num_occurrences_accross_genomes[gene_cluster_name] = len(num_occurrences_accross_genomes[gene_cluster_name])

        return num_occurrences_accross_genomes, num_genes_contributed_per_genome


    def get_num_gene_clusters_missing_per_genome_dict(self, gene_clusters_dict):
        """Get a dictionary of how many gene_clusters each genome is missing"""

        all_genomes = self.get_all_genome_names_in_gene_clusters_dict(gene_clusters_dict)
        gene_clusters_in_genomes_dict = self.get_gene_clusters_in_genomes_dict(gene_clusters_dict)

        gene_cluster_names = set(list(gene_clusters_dict.keys()))

        num_gene_clusters_missing_per_genome = dict([(genome_name, 0) for genome_name in all_genomes])

        for genome_name in all_genomes:
            for gene_cluster_name in gene_cluster_names:
                if gene_cluster_name not in gene_clusters_in_genomes_dict[genome_name]:
                    num_gene_clusters_missing_per_genome[genome_name] += 1

        return num_gene_clusters_missing_per_genome


    def filter_gene_clusters_from_gene_clusters_dict(self, gene_clusters_dict, min_num_genomes_gene_cluster_occurs=0,
             max_num_genomes_gene_cluster_occurs=sys.maxsize, min_num_genes_from_each_genome=0, max_num_genes_from_each_genome=sys.maxsize,
             min_functional_homogeneity_index=-1, max_functional_homogeneity_index=1, min_geometric_homogeneity_index=-1,
             max_geometric_homogeneity_index=1, min_combined_homogeneity_index=-1, max_combined_homogeneity_index=1):
        """This takes in your `gene_clusters_dict`, and removes gene_clusters based on their occurrences across genomes.

           The `min_num_genomes_gene_cluster_occurs` parameter defines what is the minimum number of genomes you want a gene to
           be present. It removes all the gene_clusters that do not fit into that criterion. In contrast, `max_num_genomes_gene_cluster_occurs`
           parameter will remove any gene cluster that occurs in more genomes than what the paramter asks for."""


        def check(param, param_pretty):
            if param is None:
                raise ConfigError("filter_gene_clusters_from_gene_clusters_dict is speaking: %s is is literally 'None'. It can't be." % param_pretty)

            if not isinstance(param, int):
                try:
                    param = int(param)
                except ValueError:
                    raise ConfigError("The parameter %s must occur should be of type int :/ It is of type %s" % (param_pretty,
                                                                                                                 type(param)))

            return param


        min_num_genomes_gene_cluster_occurs = check(min_num_genomes_gene_cluster_occurs, '--min-num-genomes-gene-cluster-occurs')
        max_num_genomes_gene_cluster_occurs = check(max_num_genomes_gene_cluster_occurs, '--max-num-genomes-gene-cluster-occurs')
        min_num_genes_from_each_genome = check(min_num_genes_from_each_genome, '--min-num-genes-from-each-genome')
        max_num_genes_from_each_genome = check(max_num_genes_from_each_genome, '--max-num-genes-from-each-genome')

        # behave depending on the availability of homogeneity indices in the database. if either or both of the homogeneity indices are
        # missing, we will force the parameters to remain as their defaults and keep the user posted.
        if not self.functional_homogeneity_info_is_available:
            if min_functional_homogeneity_index != -1 or max_functional_homogeneity_index != 1:
                self.run.warning("You are trying to filter your gene clusters by functional homogeneity, when your pan database does not "
                                 "include information about functional homogeneity. You can always compute this index for all of your "
                                 "gene clusters using 'anvi-compute-gene-cluster-homogeneity', but anvi'o will override your decision for now. "
                                 "You will not be able to filter your gene clusters by functional homogeneity at this time.")
                min_functional_homogeneity_index = -1
                max_functional_homogeneity_index = 1

        if not self.geometric_homogeneity_info_is_available:
            if min_geometric_homogeneity_index != -1 or max_geometric_homogeneity_index != 1:
                self.run.warning("You are trying to filter your gene clusters by geometric homogeneity, when your pan database does not "
                                 "include information about geometric homogeneity. You can always compute this index for all of your "
                                 "gene clusters using 'anvi-compute-gene-cluster-homogeneity', but anvi'o will override your decision for now. "
                                 "You will not be able to filter your gene clusters by geometric homogeneity at this time.")
                min_geometric_homogeneity_index = -1
                max_geometric_homogeneity_index = 1

        if not self.combined_homogeneity_info_is_available:
            if min_combined_homogeneity_index != -1 or max_combined_homogeneity_index != 1:
                self.run.warning("You are trying to filter your gene clusters by combined homogeneity, when your pan database does not "
                                 "include information about combined homogeneity. You can always compute this index for all of your "
                                 "gene clusters using 'anvi-compute-gene-cluster-homogeneity', but anvi'o will override your decision for now. "
                                 "You will not be able to filter your gene clusters by combined homogeneity at this time.")
                min_combined_homogeneity_index = -1
                max_combined_homogeneity_index = 1

        if min_num_genomes_gene_cluster_occurs < 0 or max_num_genomes_gene_cluster_occurs < 0:
            raise ConfigError("When you ask for a negative value for the minimum or maximum number of genomes a gene cluster is expected "
                              "to be found, you are pushing the boundaries of physics instead of biology. Let's focus on one field of science "
                              "at a time :(")

        if min_num_genes_from_each_genome < 0 or max_num_genes_from_each_genome < 0:
            raise ConfigError("Nice try. Min or max number of genes from each genome per gene cluster can't be a negative value.")

        if min_num_genomes_gene_cluster_occurs > max_num_genomes_gene_cluster_occurs:
            raise ConfigError("Min number of genomes a gene cluster should occur can't be larger than the max number of genomes a gene cluster "
                              "should occur. You're making anvi'o come up with the stupidest error messages.")

        if min_num_genes_from_each_genome > max_num_genes_from_each_genome:
            raise ConfigError("Min number of genes for each gene cluster can't be larger than the .. pfft. Anvi'o refuses to continue with this "
                              "error message. Check your parameters :(")

        if self.functional_homogeneity_info_is_available and self.geometric_homogeneity_info_is_available:
            if (min_functional_homogeneity_index < 0 and min_functional_homogeneity_index != -1) or (min_geometric_homogeneity_index < 0 and min_geometric_homogeneity_index != -1):
                raise ConfigError("Geometric and Functional homogeneity indices have a mininum value of 0, along with an error value of -1. You can either ask for "
                                  "values of 0 or greater, or put in '-1'. These are hard limits.")

            if max_functional_homogeneity_index > 1 or max_geometric_homogeneity_index > 1:
                raise ConfigError("Geometric and Functional homogeneity indices have a maximum possible value of 1. Your parameters exceed this hard upper limit. "
                                  "Please check your parameters.")

            if max_functional_homogeneity_index < min_functional_homogeneity_index or max_geometric_homogeneity_index < min_geometric_homogeneity_index:
                raise ConfigError("Please. Check your parameters. Make sure that minimum values are less than (or equal to) maximum values. We beg you")

        if self.combined_homogeneity_info_is_available:
            if min_combined_homogeneity_index < 0 and min_combined_homogeneity_index != -1:
                raise ConfigError("The combined homogeneity index have a mininum value of 0, along with an error value of -1. You can either ask for "
                                  "values of 0 or greater, or put in '-1'. These are hard limits.")

            if max_combined_homogeneity_index > 1:
                raise ConfigError("The combined homogeneity index has a maximum possible value of 1. Your parameters exceed this hard upper limit. "
                                  "Please check your parameters.")

            if max_combined_homogeneity_index < min_combined_homogeneity_index:
                raise ConfigError("Please. Check your parameters. Make sure that minimum values are less than (or equal to) maximum values. We beg you")

        all_genomes = self.get_all_genome_names_in_gene_clusters_dict(gene_clusters_dict)

        if max_num_genomes_gene_cluster_occurs == sys.maxsize:
            max_num_genomes_gene_cluster_occurs = len(all_genomes)

        if max_num_genes_from_each_genome == sys.maxsize:
            max_num_genes_from_each_genome = len(all_genomes)

        if min_num_genomes_gene_cluster_occurs > len(all_genomes):
            raise ConfigError("You have %d genomes, and you are asking anvi'o to remove any gene cluster that occurs in less than %d of them. "
                              "On the one hand, it is totally OK to make up a number like that. On the other, anvi'o would like to think that "
                              "that is not what you're doing." % (len(all_genomes), min_num_genomes_gene_cluster_occurs))

        gene_cluster_occurrences_accross_genomes, num_genes_contributed_per_genome = self.get_basic_gene_clusters_stats(gene_clusters_dict)
        if self.functional_homogeneity_info_is_available and self.geometric_homogeneity_info_is_available and not self.combined_homogeneity_info_is_available:
            homogeneity_keys, homogeneity_dict = TableForItemAdditionalData(self.args).get(['functional_homogeneity_index', 'geometric_homogeneity_index'])
        elif self.functional_homogeneity_info_is_available and self.geometric_homogeneity_info_is_available and self.combined_homogeneity_info_is_available:
            homogeneity_keys, homogeneity_dict = TableForItemAdditionalData(self.args).get(['functional_homogeneity_index', 'geometric_homogeneity_index', 'combined_homogeneity_index'])


        gene_clusters_to_remove = set([])
        all_gene_clusters = set(list(gene_cluster_occurrences_accross_genomes.keys()))
        for gene_cluster_name in all_gene_clusters:
            num_occurrence = gene_cluster_occurrences_accross_genomes[gene_cluster_name]
            num_genes_from_genomes = num_genes_contributed_per_genome[gene_cluster_name]

            if num_occurrence < min_num_genomes_gene_cluster_occurs or num_occurrence > max_num_genomes_gene_cluster_occurs:
                gene_clusters_to_remove.add(gene_cluster_name)
                continue

            if len([g for g in all_genomes if num_genes_from_genomes[g] < min_num_genes_from_each_genome or num_genes_from_genomes[g] > max_num_genes_from_each_genome]):
                gene_clusters_to_remove.add(gene_cluster_name)
                continue

            try:
                if homogeneity_dict[gene_cluster_name]['functional_homogeneity_index'] < min_functional_homogeneity_index or homogeneity_dict[gene_cluster_name]['functional_homogeneity_index'] > max_functional_homogeneity_index:
                    gene_clusters_to_remove.add(gene_cluster_name)
                    continue

                if homogeneity_dict[gene_cluster_name]['geometric_homogeneity_index'] < min_geometric_homogeneity_index or homogeneity_dict[gene_cluster_name]['geometric_homogeneity_index'] > max_geometric_homogeneity_index:
                    gene_clusters_to_remove.add(gene_cluster_name)
                    continue

            except:
                if min_functional_homogeneity_index == -1 and max_functional_homogeneity_index == 1 and min_geometric_homogeneity_index == -1 and max_geometric_homogeneity_index == 1:
                    continue #No need to raise an error if the parameters are default/all at their bounds

                raise ConfigError("Bad news: anvi'o was unable to retrieve functional and geometric homogeneity indices for gene cluster %s. This could be because "
                                  "functional and geometric homogeneity was not computed for this gene cluster when the pangenomic analysis was created. "
                                  "The good news is that you can fix that! Take a look at the anvi-compute-gene-cluster-homogeneity script" % gene_cluster_name)

            try:
                if homogeneity_dict[gene_cluster_name]['combined_homogeneity_index'] < min_combined_homogeneity_index or homogeneity_dict[gene_cluster_name]['combined_homogeneity_index'] > max_combined_homogeneity_index:
                        gene_clusters_to_remove.add(gene_cluster_name)
                        continue
            except:
                if min_combined_homogeneity_index == -1  and max_combined_homogeneity_index == 1:
                    continue
                raise ConfigError("Bad news: anvi'o was unable to retrieve the combined homogeneity index for gene cluster %s. This could be because combined homogeneity was not "
                                      "computed for this gene cluster when the pangenomic analysis was created. The good news is that you can fix that! "
                                      "Take a look at the anvi-compute-gene-cluster-homogeneity script" % gene_cluster_name)


        gene_clusters_to_keep = all_gene_clusters.difference(gene_clusters_to_remove)

        if not len(gene_clusters_to_keep):
            raise ConfigError("Bad news: the combination of your filters resulted in zero gene clusters :/ These are the filtesr anvi'o used: --min-num-genomes-gene-cluster-occurs %(min_oc)d,\
                               --max-num-genomes-gene-cluster-occurs %(max_oc)d, --min-num-genes-from-each-genome %(min_g)d, --max-num-genes-from-each-genome %(max_g)d, \
                               --min-functional-homogeneity-index %(min_fh)f, --max-functional-homogeneity-index %(max_fh)f, --min-geometric-homogeneity-index %(min_gh)f, \
                               --max-geometric-homogeneity-index %(max_gh)f, --min-combined-homogeneity-index %(min_ch)f, and --max-combined-homogeneity-index %(max_ch)f. \
                               None of your %(all_gcs)d gene clusters in your %(all_gs)d genomes that were included this analysis matched to this combination \
                               (please note that number of genomes may be smaller than the actual number of genomes in the original pan genome \
                               if other filters were applied to the gene clusters dictionary prior)." % \
                                            {'min_oc': min_num_genomes_gene_cluster_occurs, 'max_oc': max_num_genomes_gene_cluster_occurs,
                                             'min_g': min_num_genes_from_each_genome, 'max_g': max_num_genes_from_each_genome,
                                             'min_fh': min_functional_homogeneity_index, 'max_fh': max_functional_homogeneity_index,
                                             'min_gh': min_geometric_homogeneity_index, 'max_gh': max_geometric_homogeneity_index,
                                             'min_ch': min_combined_homogeneity_index, 'max_ch': max_combined_homogeneity_index,
                                             'all_gcs': len(all_gene_clusters), 'all_gs': len(all_genomes)})

        msg = "Based on --min-num-genomes-gene-cluster-occurs %d, --max-num-genomes-gene-cluster-occurs %d, \
               --min-num-genes-from-each-genome %d, --max-num-genes-from-each-genome %d, --min-functional-homogeneity-index %0.3f, \
               --max-functional-homogeneity-index %0.3f, --min-geometric-homogeneity-index %0.3f,  \
               --max-geometric-homogeneity-index %0.3f, --min-combined-homogeneity-index %0.3f,  \
               --max-combined-homogeneity-index %0.3f (some of these may be default values, no need to panic)." \
                            % (min_num_genomes_gene_cluster_occurs, max_num_genomes_gene_cluster_occurs,
                               min_num_genes_from_each_genome, max_num_genes_from_each_genome, min_functional_homogeneity_index,
                               max_functional_homogeneity_index, min_geometric_homogeneity_index, max_functional_homogeneity_index,
                               min_combined_homogeneity_index, max_combined_homogeneity_index)

        # Baris Metin: lambda functions are ugly.
        # Meren Urat : YOU'RE UGLY :(
        M = lambda l: ', '.join(l) if anvio.DEBUG \
                else ', '.join(list(l)[0:3] + ['(... %d more (`--debug` will show all))' % (len(l) - 3)]) if len(l) > 3 \
                   else ', '.join(l)

        self.run.warning(msg, "GENE CLUSTER FILTERS", lc="cyan")
        self.run.info('All gene clusters (%d)' % len(all_gene_clusters), M(all_gene_clusters))
        self.run.info('Gene clusters that passed the filter (%d)' % (len(gene_clusters_to_keep)), M(gene_clusters_to_keep), mc='green')
        self.run.info('Genes clusters that failed the filter (%d)' % (len(gene_clusters_to_remove)), M(gene_clusters_to_remove) if gene_clusters_to_remove else 'None.', nl_after=1, mc='red')

        if len(gene_clusters_to_remove):
            for gene_cluster_name in gene_clusters_to_remove:
                gene_clusters_dict.pop(gene_cluster_name)

            return (gene_clusters_dict, gene_clusters_to_remove)
        else:
            return (gene_clusters_dict, set([]))


    def filter_genomes_from_clusters_dict(self, gene_clusters_dict, max_num_gene_clusters_missing_from_genome=0):
        """This takes a `gene_clusters_dict`, and goes through every genome to identify genomes that lack more than \
           `max_num_gene_clusters_missing_from_genome` from a list of gene_clusters.

           Note that it returns a filtered dictionary, AND the genomes that are removed."""

        if not isinstance(max_num_gene_clusters_missing_from_genome, int):
            try:
                max_num_gene_clusters_missing_from_genome = int(max_num_gene_clusters_missing_from_genome)
            except ValueError:
                raise ConfigError("Well. The parameter max number of gene clusters missing from genome must be of type int.")

        if max_num_gene_clusters_missing_from_genome < 0:
            raise ConfigError("The parameter max number of gene clusters missing from genome can't be smaller than zero. "
                              "Well, it can be, as it is the case in this particlar instance, but maybe then you should "
                              "try a different platform to analyze your stuff.")

        all_genomes = self.get_all_genome_names_in_gene_clusters_dict(gene_clusters_dict)
        num_gene_clusters_missing_per_genome = self.get_num_gene_clusters_missing_per_genome_dict(gene_clusters_dict)

        genomes_to_remove = set([])
        for genome_name in num_gene_clusters_missing_per_genome:
            if num_gene_clusters_missing_per_genome[genome_name] > max_num_gene_clusters_missing_from_genome:
                genomes_to_remove.add(genome_name)

        genomes_to_keep = all_genomes.difference(genomes_to_remove)

        self.run.warning(None, "FILTER GENOMES (--max-num-gene-clusters-missing-from-genome)", lc="cyan")
        self.run.info('All genomes found (%d)' % len(all_genomes), ', '.join(all_genomes))
        self.run.info('Genomes that missed AT MOST %d of the %d gene_clusters (%d)' % (max_num_gene_clusters_missing_from_genome, len(gene_clusters_dict), len(genomes_to_keep)), ', '.join(genomes_to_keep), mc='green')
        self.run.info('Genomes that are no more in the analysis (%d)' % (len(genomes_to_remove)), ', '.join(genomes_to_remove) if genomes_to_remove else 'None. Lovely.', mc='red', nl_after=1)

        if len(genomes_to_remove) == len(all_genomes):
            raise ConfigError("Bad news: using --max-num-gene-clusters-missing-from-genome paramter with '%d' removed all of your "
                              "%d genomes from the analysis. This means every genome you have in your pangenome misses at least %d "
                              "of your %d gene clusters. Now you know :/" \
                                    % (max_num_gene_clusters_missing_from_genome,
                                       len(all_genomes),
                                       max_num_gene_clusters_missing_from_genome,
                                       len(gene_clusters_dict)))

        if len(genomes_to_remove):
            for gene_cluster_name in gene_clusters_dict:
                for genome_name in genomes_to_remove:
                    gene_clusters_dict[gene_cluster_name].pop(genome_name)
            return (gene_clusters_dict, genomes_to_remove)
        else:
            return (gene_clusters_dict, set([]))


    def filter_gene_clusters_dict(self, args, gene_clusters_dict=None):
        """Returns filtered gene clusters dicts, without editing the `self.gene_clusters` dict,
           UNLESS, it is provided a dictionary for gene clusters.

           WIHTOUT GENE CLUSTERS DICT PARAMTER
           ========================================================================================
           This function looks for variables max_num_gene_clusters_missing_from_genome and
           min_num_genomes_gene_cluster_occurs in args (see the example below), applies filters
           to a fresh copy of self.gene_clusters, and returns a new dictionary. The reason we
           avoid operating on the actual dictionary is to (1) allow testing of multiple parameters
           without having to re-initialize the actual gene clusters dictionary from the database.
           Although this seems to be a memory-intensive way of doing it, it will offer some
           flexibility to interface operations.

                >>>
                >>> pan = dbops.PanSuperclass(args)
                >>> pan.init_gene_clusters()
                >>>
                >>> args = argparse.Namespace(max_num_gene_clusters_missing_from_genome=5, min_num_genomes_gene_cluster_occurs=10)
                >>> filtered_gene_clusters = pan.filter_gene_clusters_dict(args)
                >>>

           WITH GENE CLUSTERS DICT PARAMTER
           ========================================================================================
           Does pretty much the same thing above, but it does not generate a fresh copy of the self.gene_clusters
           and only operates on the incoming dictionary.
        """

        if not gene_clusters_dict:
            if not self.gene_clusters_initialized:
                raise ConfigError("You need to initialize the gene clusters dictionary if you want to apply filters on it. "
                                  "See relevant memeber functions in your instance of PanSuperclass.")

            gene_clusters_dict = copy.deepcopy(self.gene_clusters)

        if not isinstance(gene_clusters_dict, dict):
            raise ConfigError("Houston, we have a problem. The gene clusters dict seems to be of type %s and not dict :/"\
                            % type(gene_clusters_dict))

        # let's see what the user wants.
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        max_num_gene_clusters_missing_from_genome = A('max_num_gene_clusters_missing_from_genome')
        min_num_genomes_gene_cluster_occurs = A('min_num_genomes_gene_cluster_occurs')
        max_num_genes_from_each_genome = A('max_num_genes_from_each_genome')
        min_num_genes_from_each_genome = A('min_num_genes_from_each_genome')
        max_num_genomes_gene_cluster_occurs = A('max_num_genomes_gene_cluster_occurs')
        min_functional_homogeneity_index = A('min_functional_homogeneity_index')
        max_functional_homogeneity_index = A('max_functional_homogeneity_index')
        min_geometric_homogeneity_index = A('min_geometric_homogeneity_index')
        max_geometric_homogeneity_index = A('max_geometric_homogeneity_index')
        min_combined_homogeneity_index = A('min_combined_homogeneity_index')
        max_combined_homogeneity_index = A('max_combined_homogeneity_index')
        add_into_items_additional_data_table = A('add_into_items_additional_data_table')
        gene_clusters_names_of_interest = A('gene_clusters_names_of_interest')
        just_do_it = A('just_do_it')

        # keep only the names we are interested in.
        if gene_clusters_names_of_interest:
            unwanted_keys = set(gene_clusters_dict.keys()) - set(gene_clusters_names_of_interest)
            for key in unwanted_keys:
                del gene_clusters_dict[key]


        # remove genomes from the dict if necessary.
        if max_num_gene_clusters_missing_from_genome:
            gene_clusters_dict, genomes_removed = self.filter_genomes_from_clusters_dict(gene_clusters_dict, max_num_gene_clusters_missing_from_genome)

        # remove gene clusters from the dict if necessary
        if min_num_genomes_gene_cluster_occurs or max_num_genomes_gene_cluster_occurs:
            gene_clusters_dict, gene_clusters_removed = \
                    self.filter_gene_clusters_from_gene_clusters_dict(gene_clusters_dict,
                                                                      min_num_genomes_gene_cluster_occurs,
                                                                      max_num_genomes_gene_cluster_occurs,
                                                                      min_num_genes_from_each_genome,
                                                                      max_num_genes_from_each_genome,
                                                                      min_functional_homogeneity_index,
                                                                      max_functional_homogeneity_index,
                                                                      min_geometric_homogeneity_index,
                                                                      max_geometric_homogeneity_index,
                                                                      min_combined_homogeneity_index,
                                                                      max_combined_homogeneity_index)

        # this is where we add the items in the resulting filtered dict into the items additonal data
        # table:
        if add_into_items_additional_data_table:
            data_key = add_into_items_additional_data_table

            if not self.gene_cluster_names_in_db:
                self.init_gene_clusters()

            items_additional_data_dict = {}

            for gene_cluster_name in self.gene_cluster_names_in_db:
                if gene_cluster_name in gene_clusters_dict:
                    items_additional_data_dict[gene_cluster_name] = {data_key: 'TRUE'}
                else:
                    items_additional_data_dict[gene_cluster_name] = {data_key: None}

            items_additional_data_table = TableForItemAdditionalData(argparse.Namespace(pan_db=self.pan_db_path, just_do_it=just_do_it))
            items_additional_data_table.add(items_additional_data_dict, [data_key])

        return gene_clusters_dict


    def init_gene_clusters(self, gene_cluster_ids_to_focus=set([])):
        """Initializes the gene_clusters dictionary (only for `gene_cluster_ids_to_focus` if necessary).

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
        self.gene_cluster_names_in_db = pan_db.db.get_single_column_from_table(t.pan_gene_clusters_table_name, 'gene_cluster_id', unique=True)

        for entry in list(gene_clusters_long_list.values()):
            genome_name = entry['genome_name']
            gene_callers_id = entry['gene_caller_id']
            gene_cluster_id = entry['gene_cluster_id']

            if genome_name not in self.gene_callers_id_to_gene_cluster:
                self.gene_callers_id_to_gene_cluster[genome_name] = {}

            self.gene_callers_id_to_gene_cluster[genome_name][gene_callers_id] = gene_cluster_id

            if gene_cluster_ids_to_focus and gene_cluster_id not in gene_cluster_ids_to_focus:
                continue

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

        if not len(gene_cluster_ids_to_focus):
            self.run.info_single("Gene clusters are initialized for all %d gene clusters in the database." % len(self.gene_clusters), nl_before=1, nl_after=1)
        else:
            self.run.info_single("A short announcement for the curious: anvi'o found %d gene clusters in the database, attempted to "
                                 "initialize a gene clusters dictionary for %d of them as requested by the user or the programmer, and "
                                 "managed to get back a gene clusters dictionary with %d items. We just hope all these make sense to you." \
                                % (len(self.gene_cluster_names_in_db), len(gene_cluster_ids_to_focus), len(self.gene_clusters)), nl_after=1, nl_before=1)

        # gene cluster names were set when we first initialized the class, but if we are here, it means the user may have
        # alrady updated the list of gene clusters. let's keep ourselves up-to-date:
        self.gene_cluster_names = set(list(self.gene_clusters.keys()))
        self.gene_clusters_initialized = True


    def load_pan_views(self, split_names_of_interest=None):
        pan_db = PanDatabase(self.pan_db_path)

        views_table = pan_db.db.get_table_as_dict(t.views_table_name)

        for view in views_table:
            table_name = views_table[view]['target_table']
            data, header = pan_db.db.get_view_data(table_name, split_names_of_interest=split_names_of_interest)
            self.views[view] = {'table_name': table_name,
                                'header': header,
                                'dict': data}

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
            raise ConfigError("init_collection_profile wants to initialize the collection profile for '%s', but the "
                               "the gene clusters dict is kinda empty :/ Someone forgot to initialize something maybe?" \
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


    def search_for_gene_functions(self, search_terms, requested_sources=None, verbose=False, full_report=False, case_sensitive=False, exact_match=False):
        if not isinstance(search_terms, list):
            raise ConfigError("Search terms must be of type 'list'")

        # FIXME: Even though the function header contaisn `case_sensitive` and `exact_match` variables,
        #        the function below pays no attention to these variables. FURTHERMORE, this is a pretty
        #        shitty way to handle function search operations -- there are two functions in this file
        #        with the same name, with lots of redundant code. We need to consolidate them into a single
        #        class that can deal with pan or contigs db files seamlessly :/

        search_terms = [s.strip() for s in search_terms]

        if len([s.strip().lower() for s in search_terms]) != len(set([s.strip().lower() for s in search_terms])):
            raise ConfigError("Please do not use the same search term twice :/ Becasue, reasons. You know.")

        for search_term in search_terms:
            if not len(search_term) >= 3:
                raise ConfigError("A search term cannot be less than three characters")

        self.run.info('Search terms', '%d found' % (len(search_terms)))
        gene_clusters = {}
        full_report = []

        genomes_db = db.DB(self.genomes_storage_path, anvio.__genomes_storage_version__)

        for search_term in search_terms:
            self.progress.new('Search functions')
            self.progress.update('Searching for term "%s"' % search_term)

            query = '''select gene_callers_id, source, accession, function, genome_name from ''' + t.genome_gene_function_calls_table_name + ''' where (function LIKE "%%''' \
                            + search_term + '''%%" OR accession LIKE "%%''' + search_term + '''%%")'''

            query += ''' AND genome_name IN (%s) ''' % (', '.join(["'%s'" % s for s in self.p_meta['genome_names']]))

            if requested_sources:
                query += ''' AND source IN (%s);''' % (', '.join(["'%s'" % s for s in requested_sources]))
            else:
                query += ';'

            results = genomes_db._exec(query).fetchall()
            gene_clusters[search_term] = []

            found_mismatch = False
            for result in results:
                gene_caller_id, source, accession, function, genome_name = result

                # we're finding gene caller ids in the genomes storage, but they may not end up in any
                # of the final gene clusters stored in the pan database due to various reasons. for
                # instance, if the user set a min occurrence parameter, a singleton will not be found
                # in the pan db yet it will return a functional hit.
                if not gene_caller_id in self.gene_callers_id_to_gene_cluster[genome_name]:
                    gene_cluster_id = 'n/a'
                    found_mismatch = True
                else:
                    gene_cluster_id = self.gene_callers_id_to_gene_cluster[genome_name][gene_caller_id]

                gene_dict = self.genomes_storage.gene_info[genome_name][gene_caller_id]

                full_report.extend([(gene_caller_id, genome_name, source, accession, function, search_term,
                    gene_cluster_id, gene_dict['dna_sequence'], gene_dict['aa_sequence'])])

                gene_clusters[search_term].append(gene_cluster_id)
            self.progress.end()

            if found_mismatch:
                self.run.warning("Some of the search results for the term '%s' found in your genomes storage do not seem to "
                                "belong any gene cluster in your pan database. This may be due to filtering parameters used (ex: --min-occurrence) "
                                "during the pangenome analysis. Gene cluster ids for these results will appear as 'n/a' in the report." % search_term)

        genomes_db.disconnect()
        self.progress.end()

        return gene_clusters, full_report


    def list_function_sources(self):
        genome_storage = genomestorage.GenomeStorage(self.genomes_storage_path, run=terminal.Run(verbose=False))
        gene_function_sources = genome_storage.db.get_meta_value('gene_function_sources').split(',')
        genome_storage.close()

        if not len(gene_function_sources):
            self.run.info_single('No functional annotations found in this genomes storage :/', nl_before=1, nl_after=1, mc='red')
        else:
            self.run.warning('', 'AVAILABLE FUNCTIONS (%d FOUND)' % (len(gene_function_sources)), lc='yellow')
            for source in gene_function_sources:
                self.run.info_single('%s' % (source), nl_after = 1 if source == gene_function_sources[-1] else 0)


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

    Alternatively, you can define a set of split names of interest to gain performance when it is
    needed. There are two ways to do that, which are mutually exclusive (so you have to grow up and
    pick one). One way is to explicitly mention which splits are of interest (for control freaks):

         >>> args.split_names_of_interest = set([split_names])
         >>> p = ProfileSuperclass(args)

    The second way to initialize ProfileSuper with a subset of splits a profile database contains
    is to use the collections framework (the elegant way of doing this). For which, you need to
    set collection name:

        >>> args.collection_name = 'collection_name'
        >>> args.bin_ids = 'bin_1,bin_2,bin_3' # if no bin_ids is provided, all bins will be used
        >>> p = ProfileSuperClass(args)

    The best practice is to set anvi'o programs to put together `args` objects with these variables.
    """

    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        # these are initialized by the member function `init_gene_level_coverage_stats_dicts`. but you knew
        # that already because you are a smart ass.
        self.gene_level_coverage_stats_dict = {}
        self.split_coverage_values_per_nt_dict = {}

        # this one becomes the object that gives access to the auxiliary data ops for split coverages
        # used heavily in interactive interface to show stuff (see bottle routes and all).
        self.split_coverage_values = None

        # the following two are initialized via `init_items_additional_data()` and use information
        # stored in item additional data tables UNLESS THEY ARE ALREADY INITIALIZED IN THE CONTEXT
        # FROM WITHIN PROFILE SUPERCLASS IS INHERITED (SUCH A THING IS HAPPENING AT THE INTERACTIVE
        # CLASS)
        if super() and 'items_additional_data_dict' in dir(self) and 'items_additional_data_keys' in dir(self):
            pass
        else:
            self.items_additional_data_dict = None
            self.items_additional_data_keys = None

        if super() and 'layers_additional_data_dict' in dir(self) and 'layers_additional_data_keys' in dir(self):
            pass
        else:
            self.layers_additional_data_keys, self.layers_additional_data_dict = TableForLayerAdditionalData(self.args, p=self.progress).get_all()

        self.auxiliary_profile_data_available = None
        self.auxiliary_data_path = None
        self.genes_db_available = None
        self.genes_db_path = None

        self.split_names = set([])
        self.item_orders = {}
        self.views = {}
        self.collection_profile = {}

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.profile_db_path = A('profile_db')
        self.contigs_db_path = A('contigs_db')
        init_gene_coverages = A('init_gene_coverages')
        skip_consider_gene_dbs = A('skip_consider_gene_dbs')
        init_split_coverage_values_per_nt = A('init_split_coverage_values_per_nt')
        outliers_threshold = A('outliers_threshold')
        zeros_are_outliers = A('zeros_are_outliers')
        self.inseq_stats = A('inseq_stats')

        # early on let's check some ground truth
        if not self.profile_db_path:
            self.run.warning("ProfileSuper is called with args without member profile_db. Anvi'o will assume "
                             "you are a programmer, and will not raise an error. But the init function is returning "
                             "prematurely. Just so you know.")
            return

        utils.is_profile_db(self.profile_db_path)

        # NOTE for programmers. The next few lines are quite critical for the flexibility of ProfileSuper.
        # Should we initialize the profile super for a specific list of splits? This is where we take care of that.
        # the user can initialize the profile super two ways: by providing split names of interest explicitly, or
        # by providing collection name and bin names in args.
        if not hasattr(self, 'split_names_of_interest'):
            self.split_names_of_interest = set([])
        elif hasattr(self, 'collection_name') and self.collection_name:
            # if self.split_names_of_interest is defined upstream somewhere,
            # but if we ALSO have a collection name here, we want to recover those
            # split names relevant to the collection name later. so in this case,
            # we will OVERWRITE this variable, which is kind of dangerous. IF this is
            # causing you headache and sets your split names of interest to null, the
            # best way to avoid this is to set the collection_name in your args to None,
            # initialize the ProfileSuper, and then set it back to its original value
            # after (which has been the case in the refine mode in anvio/interactive.py).
            if hasattr(self, 'split_names_of_interest'):
                self.run.warning("Note for programmers: even though the ProfileSuper was initialized with "
                                 "`self.split_names_of_interest`, the contents of that variable just got "
                                 "set to `set([])` becasue the context also had `self.collection_name` in "
                                 "it. If you are concerned about this, read the comments in the code where "
                                 "this warning is shown. If you are a user, there probably is nothing for "
                                 "you to be concerned about (unless things go haywire downstream, in which "
                                 "case our sincere condolences).")

            self.split_names_of_interest = set([])

        if A('split_names_of_interest'):
            self.split_names_of_interest = set(self.args.split_names_of_interest)

        self.collection_name = A('collection_name')

        # figure out bin names, if there is one to figure out
        if A('bin_id') and A('bin_names_list'):
            raise ConfigError("ProfileSuper says you can't use both `bin_id` and `bin_names_list` as argument. Pick "
                              "one, and stick with it. ProfileSuper is grumpy.")
        if A('bin_id'):
            self.bin_names = [A('bin_id')]
        elif A('bin_names_list'):
            if isinstance(A('bin_names_list'), list):
                self.bin_names = A('bin_names_list')
            elif isinstance(A('bin_names_list'), str):
                self.bin_names = A('bin_names_list').split(',')
            else:
                raise ConfigError("ProfileSuper says `bin_names_list` can either be a string of comma-separated bin "
                                  "names, or a proper Python `list` of bin names. But not %s." % (type(A('bin_names_list'))))
        else:
            self.bin_names = None

        if self.split_names_of_interest and not isinstance(self.split_names_of_interest, type(set([]))):
            raise ConfigError("ProfileSuper says the argument `split_names_of_interest` must be of type set(). "
                              "Someone screwed up somewhere :/")
        elif self.split_names_of_interest and self.collection_name:
            raise ConfigError("ProfileSuper is initialized with args that contain both `split_names_of_interest`,\
                               and `collection_name`. You can initialize the ProfileSuper with either of those. As\
                               a programmer if you have no control over incoming `args` and just passing things\
                               around, you might need to implement a workaround to set either of those params to None\
                               and then reset them back to their original in `args` once you are done with\
                               ProfileSuper.")

        if self.split_names_of_interest:
            self.run.warning("ProfileSuperClass is inherited with a set of split names of interest, which means it will be "
                             "initialized using only the %d split names specified" % (len(self.split_names_of_interest)))
        elif self.collection_name and not utils.is_blank_profile(self.profile_db_path):
            self.run.warning("ProfileSuperClass found a collection focus, which means it will be initialized using only "
                             "the splits in the profile database that are affiliated with the collection %s and "
                             "%s it describes." % (self.collection_name, \
                                                   'bins "%s" ' % ', '.join(self.bin_names) if self.bin_names else 'all bins'))
            self.split_names_of_interest = ccolections.GetSplitNamesInBins(self.args).get_split_names_only()


        # we have a contigs db? let's see if it's for real.
        if self.contigs_db_path:
            utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.progress.new('Initializing the profile database superclass')

        self.progress.update('Loading split names')
        if utils.is_blank_profile(self.profile_db_path) and self.contigs_db_path:
            self.split_names = utils.get_all_item_names_from_the_database(self.contigs_db_path)
        else:
            self.split_names = utils.get_all_item_names_from_the_database(self.profile_db_path)

        split_names_missing = (self.split_names_of_interest - self.split_names) if self.split_names_of_interest else None
        if self.split_names_of_interest and len(split_names_missing):
            self.progress.end()
            raise ConfigError("%d of the %d split names of interest does not occur in the profile database. Here is "
                              "an example: '%s'." % (len(split_names_missing), len(self.split_names_of_interest), split_names_missing.pop()))

        self.progress.update('Creating an instance of the profile database')
        profile_db = ProfileDatabase(self.profile_db_path)

        self.progress.update('Setting profile self data dict')
        self.p_meta = profile_db.meta

        self.p_meta['db_variant'] = str(utils.get_db_variant(self.profile_db_path))
        self.p_meta['creation_date'] = utils.get_time_to_date(self.p_meta['creation_date']) if 'creation_date' in self.p_meta else 'unknown'
        self.p_meta['samples'] = sorted([s.strip() for s in self.p_meta['samples'].split(',')])
        self.p_meta['num_samples'] = len(self.p_meta['samples'])

        if self.p_meta['blank'] and not self.p_meta['contigs_db_hash']:
            self.progress.end()
            raise ConfigError("ProfileSuperclass is upset, because it seems you are tyring to initialize a blank anvi'o profile "
                              "database that is not associated with a contigs database. This will not work for multiple reasons. "
                              "The current technical limitation is that blank profile databases that are in this situation do not "
                              "keep track of split names they are working with. Yes. We too know that this is a serious design "
                              "flaw, but THANKS for reminding anyway... The best way to address this is to make sure all anvi'o "
                              "profile and pan databases maintain a table with all item names they are supposed to be working with.")

        self.progress.end()

        # learn the number of mapped reads and set it in a nice variable VERY CAREFULLY (blank profiles don't have it,
        # and some ancient anvi'o databases may be lacking it).
        if self.p_meta['blank']:
            self.num_mapped_reads_per_sample = None
        else:
            if self.p_meta['merged']:
                if 'total_reads_mapped' not in self.p_meta:
                    self.num_mapped_reads_per_sample = None
                else:
                    total_reads_mapped = [int(num_reads) for num_reads in self.p_meta['total_reads_mapped'].split(',')]
                    self.num_mapped_reads_per_sample = {self.p_meta['samples'][i]: total_reads_mapped[i] for i in range(0, len(self.p_meta['samples']))}
            else:
                sample_name = self.p_meta['samples'][0]

                # create a copy of the args to avoid some misunderstanding downstream when THIS class is initiated with
                # a gene_mode = True statement. Here we are trying to learn the layer additional data from a profile
                # database. But if genes_mod = True, TableForLayerAdditionalData class will try to initialize for a genes
                # database rather than a profile database. this workaround addresses that issue:
                args = copy.deepcopy(self.args)
                args.gene_mode = False
                keys, data = TableForLayerAdditionalData(args, p=self.progress).get()

                if 'total_reads_mapped' not in data[sample_name]:
                    self.num_mapped_reads_per_sample = None
                else:
                    self.num_mapped_reads_per_sample = {sample_name: int(data[sample_name]['total_reads_mapped'])}

        profile_db.disconnect()

        self.p_meta['available_item_orders'], self.item_orders = get_item_orders_from_db(self.profile_db_path)

        if not self.item_orders:
            self.p_meta['default_item_order'] = None

        self.progress.new('Initializing the profile database superclass')
        self.progress.update('Accessing the auxiliary data file')
        self.auxiliary_data_path = get_auxiliary_data_path_for_profile_db(self.profile_db_path)
        if not os.path.exists(self.auxiliary_data_path):
            self.auxiliary_profile_data_available = False
        else:
            self.auxiliary_profile_data_available = True
            self.split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.auxiliary_data_path,
                                                                                         self.p_meta['contigs_db_hash'],
                                                                                         db_variant=self.p_meta['db_variant'])

        if self.collection_name and self.bin_names and len(self.bin_names) == 1 and not skip_consider_gene_dbs:
            self.progress.update('Accessing the genes database')
            self.genes_db_path = utils.get_genes_database_path_for_bin(self.profile_db_path,
                                                                       self.collection_name,
                                                                       self.bin_names[0])
            if not os.path.exists(self.genes_db_path):
                self.genes_db_available = False
            else:
                self.genes_db_available = True

        self.progress.end()

        if init_gene_coverages:
            self.init_gene_level_coverage_stats_dicts(outliers_threshold=outliers_threshold,
                                                      zeros_are_outliers=zeros_are_outliers,
                                                      init_split_coverage_values_per_nt=init_split_coverage_values_per_nt)

        if self.auxiliary_profile_data_available:
            self.run.info('Auxiliary Data', 'Found: %s (v. %s)' % (self.auxiliary_data_path, anvio.__auxiliary_data_version__))

        if self.split_names_of_interest:
            self.run.info('Profile Super', 'Initialized with %d of %d splits: %s (v. %s)' % (len(self.split_names_of_interest),
                                                                                             len(self.split_names),
                                                                                             self.profile_db_path,
                                                                                             anvio.__profile__version__))
        else:
            self.run.info('Profile Super', 'Initialized with all %d splits: %s (v. %s)' % (len(self.split_names),
                                                                                           self.profile_db_path,
                                                                                           anvio.__profile__version__))


    def create_blank_genes_database(self, split_names):
        if self.genes_db_available:
            raise ConfigError("You can't create a blank genes database when there is already one :/")

        splits_hash = utils.get_hash_for_list(split_names)

        meta_values = {'anvio': __version__,
                       'contigs_db_hash': self.p_meta['contigs_db_hash'],
                       'collection_name': self.collection_name,
                       'bin_name': self.bin_names[0],
                       'splits_hash': splits_hash}

        # make sure the GENES directory exists.
        filesnpaths.gen_output_directory(os.path.dirname(self.genes_db_path), progress=self.progress, run=self.run, delete_if_exists=False, dont_warn=False)

        # generate a blank genes database here:
        GenesDatabase(self.genes_db_path).create(meta_values=meta_values)

        # and immediately copy-paste the layers additional data and states table into
        # it for convenience:
        genes_database = db.DB(self.genes_db_path, None, ignore_version=True)
        genes_database.copy_paste(table_name=t.layer_additional_data_table_name, source_db_path=self.profile_db_path)
        genes_database.copy_paste(table_name=t.states_table_name, source_db_path=self.profile_db_path)
        genes_database.disconnect()

        # vamp the stage like a pro:
        self.genes_db_available = True
        self.genes_db_path = self.genes_db_path


    def store_gene_level_coverage_stats_into_genes_db(self, parameters, gene_caller_ids_to_exclude=None):
        table_for_gene_level_coverages = TableForGeneLevelCoverages(self.genes_db_path,
                                                                    parameters,
                                                                    split_names=self.split_names_of_interest,
                                                                    mode="INSEQ" if self.inseq_stats else "STANDARD",
                                                                    run=self.run)

        if gene_caller_ids_to_exclude:
            for gene_callers_id in gene_caller_ids_to_exclude:
                if gene_callers_id in self.gene_level_coverage_stats_dict:
                    self.gene_level_coverage_stats_dict.pop(gene_callers_id)

        table_for_gene_level_coverages.store(self.gene_level_coverage_stats_dict)


    def init_gene_level_coverage_stats_from_genes_db(self, parameters):
        if not (self.collection_name and len(self.bin_names) == 1):
            raise ConfigError("The function `get_gene_level_coverage_stats_dicts_for_a_bin` can only be called from an instance "
                              "of the profile super class that is initalized with a collection name and a single bin.")

        table_for_gene_level_coverages = TableForGeneLevelCoverages(self.genes_db_path,
                                                                    parameters,
                                                                    split_names=self.split_names_of_interest,
                                                                    mode="INSEQ" if self.inseq_stats else "STANDARD",
                                                                    run=self.run)

        self.gene_level_coverage_stats_dict = table_for_gene_level_coverages.read()


    def init_gene_level_coverage_stats_dicts(self, min_cov_for_detection=0, outliers_threshold=1.5, zeros_are_outliers=False, callback=None, callback_interval=100, init_split_coverage_values_per_nt=False, gene_caller_ids_of_interest=set([])):
        """This function will populate both `self.split_coverage_values_per_nt_dict` and
           `self.gene_level_coverage_stats_dict`.

           Note: if a `split_names_of_interest` argument is declared at the class level,
           this function will operate on those splits found in that set.

           If there is a collection and a single bin name, then this function will work with
           the genes database to read form or to write to.
           """

        # let's get these parameters set
        parameters = {
            'min_cov_for_detection': min_cov_for_detection,
            'outliers_threshold': outliers_threshold,
            'zeros_are_outliers': zeros_are_outliers,
            'mode': 'INSEQ' if self.inseq_stats else 'STANDARD'
        }

        if self.p_meta['blank']:
            self.run.warning("Someone asked gene coverages to be initialized when working with a blank profile database. "
                             "Anvi'o will pretend nothing happened, and will return nothing. If you don't know what this "
                             "is warning you about, just carry on.")
            return

        if not self.auxiliary_profile_data_available:
            raise ConfigError("Someone is asking gene-level coverage stats to be computed, but then there is no AUXILIARY-DATA.db "
                              "available for this project. Yeah. That's what happens if you don't download everything from the "
                              "server :(")

        contigs_db = ContigsSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        if not contigs_db.a_meta['genes_are_called']:
            self.run.warning("Well, someone wants to populate the gene coverages data, when in fact genes were not called :/ "
                             "Instead of giving an error, anvi'o will return prematurely, without really doing anything.")
            return

        if not contigs_db.a_meta['splits_consider_gene_calls']:
            self.run.warning("PLEASE READ THIS VERY CAREFULLY (remember, anvi'o never talks to you in CAPS, so it must be important). "
                             "It seems when you generated your contigs database, you have skipped 'mindful' splitting of contigs. "
                             "This means, some of the genes may be soft-broken into two or more pieces. For most things, it doesn't "
                             "really matter, but here this will cause an issue as your gene coverages will average one of those splits "
                             "without any biologically relevant reason. We could have done much better here, but it would have affected "
                             "the performance very negatively. If you are seeing this warning, and go like 'crap, this will ruin "
                             "everything because I possibly can not recover from this situation', then send us an e-mail, and we will "
                             "think about whether we can be less lazy about stuff, and do things better.")

        if self.split_names_of_interest:
            split_names = self.split_names_of_interest

            self.run.warning('A subset of splits (%d of %d, to be precise) are requested to initiate gene-level coverage stats for. '
                             'No need to worry, this is just a warning in case you are as obsessed as wanting to know everything '
                             'there is to know.' % (len(self.split_names_of_interest), len(self.split_names)))
        else:
            split_names = self.split_names

        if self.genes_db_path and self.genes_db_available:
            # THIS IS A SPECIAL CASE, where someone is initializing the gene-level coverage
            # stats for a single bin. In this case anvi'o will want to work with a genes
            # database to read from, or to populate one for later uses. But the proplem is,
            # we may be called from a part of the code that doesn't know what KIND of genes
            # database is being called. Therefore, the `parameters` dict we are about to send
            # to `init_gene_level_coverage_stats_from_genes_db` may contain the default
            # `mode` value set before. BUT WE DON'T WANT THAT.
            mode_set_in_db = db.DB(self.genes_db_path, client_version=None, ignore_version=True).get_meta_value('mode')

            self.run.warning("A gene stats database of type '%s' is found (anvi'o hopes that this is the type of stats you "
                             "were expecting to find)." % mode_set_in_db.upper())

            parameters['mode'] = mode_set_in_db

            # since we are here and learned the mode, we can also set the self.inseq_stats variable IF the table
            # is actually inseq stats table. if we don't do this, the interactive interface will never load the inseq
            # data because this variable is not set anywhere :/ the best practice would have been using a mode variable
            # rather than an operation specific boolean flag, but well .. apologies to future generations of developers:
            if mode_set_in_db == "INSEQ":
                self.inseq_stats = True

            self.init_gene_level_coverage_stats_from_genes_db(parameters)
        elif self.genes_db_path and not self.genes_db_available:
            self.run.warning("You don't seem to have a genes database associated with your profile database. "
                             "Genes database is an optional anvi'o database to store gene-level coverage and "
                             "stats dicts. Anvi'o will attempt to create one for you.", lc="cyan")
            self.create_blank_genes_database(split_names)

        if len(self.gene_level_coverage_stats_dict):
            # FIXME: the design here is a fucking mess and needs a fresh look. we need split coverage per nt values initiated to
            #        to calculate the gene-level coverage stats. but if we have a genes database, gene-level coverage stats are
            #        already computed, so we don't need it .. but if we are using MCGC, even if we have gene-level coverage stats
            #        available, we need split coverages per-nt recalculated... that's why we have a miserable variable now called
            #        `init_split_coverage_values_per_nt`, when populating this disctionary should somehow be handled intrinsically
            #        by our otherwise talented superclass. we really sort this out in a more beautiful way, in my opinion. otherwise
            #        we will get stuck somwhere that will take forever to dig ourselves out.
            if init_split_coverage_values_per_nt:
                self.init_split_coverage_values_per_nt_dict(split_names)

            # we have nothing to do here anymore. and can return.
            return

        # if we have not 'returned' yet it means we gotta go through this
        self.init_split_coverage_values_per_nt_dict(split_names)

        self.progress.new('Computing gene-level coverage stats in %s mode...' % ('INSEQ' if self.inseq_stats else 'STANDARD'))
        self.progress.update('...')

        num_splits, counter = len(split_names), 1
        failed_gene_caller_ids_set = set([])
        # go through all the split names
        for split_name in split_names:
            if num_splits > 10 and counter % 10 == 0:
                self.progress.update('%d of %d splits ...' % (counter, num_splits))

            if len(gene_caller_ids_of_interest):
                gene_level_coverage_stats, failed_gene_caller_ids = self.get_gene_level_coverage_stats(split_name, contigs_db, gene_caller_ids_of_interest=gene_caller_ids_of_interest, **parameters)
            else:
                gene_level_coverage_stats, failed_gene_caller_ids = self.get_gene_level_coverage_stats(split_name, contigs_db, **parameters)

            if len(failed_gene_caller_ids):
                failed_gene_caller_ids_set.update(failed_gene_caller_ids)

            self.gene_level_coverage_stats_dict.update(gene_level_coverage_stats)

            if callback and counter % callback_interval == 0:
                callback()

            counter += 1

        self.progress.end()

        if len(failed_gene_caller_ids_set):
            self.run.warning(f"Please read this carefully as something sad just happened. While anvi'o was trying to recover "
                             f"gene level coverage stats, it became clear that a few gene calls were not found in splits they "
                             f"were meant to be found. These genes will not be a part of any downstream reporting :/ It is "
                             f"extremely difficult to even entertain the idea why this might have happened, we suspect it "
                             f"is some sort of artifact left behind from the use of external gene calls. Regardless, here "
                             f"are the gene calls that cause you this headache, in case you would like to go after this and "
                             f"find out why is this happening: {', '.join([str(f) for f in failed_gene_caller_ids_set])}",
                             header="🚑 SOMETHING WEIRD HAPPENED 🚑")

        if callback:
            callback()
        else:
            if self.genes_db_path:
                # we computed all the stuff, and we can as well store them into the genes db.
                self.store_gene_level_coverage_stats_into_genes_db(parameters, gene_caller_ids_to_exclude=failed_gene_caller_ids_set)


    def init_split_coverage_values_per_nt_dict(self, split_names=None):
        if not self.auxiliary_profile_data_available:
            raise ConfigError("What you're trying to do requires the AUXILIARY-DATA.db file :/ Please make sure it is in the "
                              "same directory with the profile database you are working with.")

        self.progress.new('Computing split coverage values per nt ...')
        self.progress.update('...')

        if not split_names:
            split_names = self.split_names_of_interest

        num_splits, counter = len(split_names), 1
        for split_name in split_names:
            if num_splits > 10 and counter % 10 == 0:
                self.progress.update('%d of %d splits ...' % (counter, num_splits))

            self.split_coverage_values_per_nt_dict[split_name] = self.split_coverage_values.get(split_name)

        self.progress.end()


    def get_gene_level_coverage_stats_entry_for_default(self, gene_callers_id, split_coverage, sample_name, gene_start, gene_stop, gene_length, outliers_threshold=1.5):
        """Returns coverage stats for a single gene in default mode.

        The alternative to this mode is the INSEQ/Tn-SEQ mode that is handled in `get_gene_level_coverage_stats_entry_for_inseq`,
        where coverage statistics are computed differently.
        """
        # and recover the gene coverage array per position for a given sample:
        gene_coverage_values_per_nt = split_coverage[sample_name][gene_start:gene_stop]

        # if we fail the following, there is something absolutely very very wrong with
        # this gene call as it does not seem to be in the split in which it is supposed
        # to be.
        if not len(gene_coverage_values_per_nt):
            return None

        mean_coverage = numpy.mean(gene_coverage_values_per_nt)
        detection = numpy.count_nonzero(gene_coverage_values_per_nt) / gene_length

         # findout outlier positions, and get non-outliers
        outliers_bool = utils.get_list_of_outliers(gene_coverage_values_per_nt, outliers_threshold)
        non_outlier_positions = numpy.invert(outliers_bool)
        non_outliers = gene_coverage_values_per_nt[non_outlier_positions]


        if not (len(non_outliers)):
            non_outlier_mean_coverage = 0.0
            non_outlier_coverage_std = 0.0
        else:
            non_outlier_mean_coverage = numpy.mean(non_outliers)
            non_outlier_coverage_std = numpy.std(non_outliers)

        return {'gene_callers_id': gene_callers_id,
                'sample_name': sample_name,
                'mean_coverage': mean_coverage,
                'detection': detection,
                'non_outlier_mean_coverage': non_outlier_mean_coverage,
                'non_outlier_coverage_std': non_outlier_coverage_std,
                'gene_coverage_values_per_nt': gene_coverage_values_per_nt,
                'non_outlier_positions': non_outlier_positions}


    def get_gene_level_coverage_stats_entry_for_inseq(self, gene_callers_id, split_coverage, sample_name, gene_start, gene_stop, gene_length, outliers_threshold=0.9):
        """Returns coverage stats for a single gene in INSEQ/Tn-SEQ mode.

           The alternative to this mode is the default mode that is handled in `get_gene_level_coverage_stats_entry_for_default`,
           where coverage statistics are computed in most conventional ways.
        """

        # Lets ignore those pesky warnings...
        numpy.seterr(divide='ignore', over='ignore')

        if not len(self.num_mapped_reads_per_sample):
            raise ConfigError("Total read counts were not set for this database, without which INSEQ/Tn-SEQ coverage stats "
                              "can't be recovered :/ This number is automatically set by anvi'o during profiling given the "
                              "short read information in BAM files that match to contigs of interest. If you are working with "
                              "a 'blank' anvi'o profile database, there is no hope for you (since there is no BAM files involved "
                              "in that workflow), but if you are working with a legacy database there are other ways to set this "
                              "number (for instance, by using `anvi-db-info` program to set a `total_reads_mapped` variable). "
                              "If you want to do this but have no idea how this would work, please get in touch with the anvi'o "
                              "community and someone will help you :)")

        total_read_counts_in_sample = self.num_mapped_reads_per_sample[sample_name]
        gene_coverage_values_per_nt = split_coverage[sample_name][gene_start:gene_stop]

        # if we fail the following, there is something absolutely very very wrong with
        # this gene call as it does not seem to be in the split in which it is supposed
        # to be.
        if not len(gene_coverage_values_per_nt):
            return None

        # variables for INSEQ/Tn-SEQ views
        mean_coverage = 0
        total_counts_of_sites_in_gene = 0
        total_counts_of_sites_in_gene_normalized = 0
        mean_three_prime = 0
        below_threshold = 0

        # Split gene coverage values into splits that are nonzero
        insertion_splits_nonzero = [numpy.array(list(g)) for k, g in itertools.groupby(gene_coverage_values_per_nt, lambda x: x != 0) if k]
        if insertion_splits_nonzero:
            for split in insertion_splits_nonzero:
                ta_sites = scipy.signal.find_peaks(split)[0]

                gene_sites_stats = list()
                if ta_sites.any():
                    # If we found high peaks, then we check if they are TA sites.
                    for ta_index in ta_sites:
                        index_highest_peak = numpy.where(gene_coverage_values_per_nt == split[ta_index])[0][0]
                        value_highest_peak = gene_coverage_values_per_nt[index_highest_peak]

                        # Need to check left and right indexes to see it part of TA site
                        check_sites_dict = {gene_coverage_values_per_nt[index_highest_peak + 1]: index_highest_peak + 1,
                                            gene_coverage_values_per_nt[index_highest_peak - 1]: index_highest_peak - 1}

                        # Need to check left and right indexes to see it part of TA site
                        value_second_highest_possible_peak = min(check_sites_dict.keys(), key=lambda x: abs(x - value_highest_peak))
                        index_second_highest_possible_peak = check_sites_dict[value_second_highest_possible_peak]

                        # Decide which one is left and right.
                        if index_highest_peak > index_second_highest_possible_peak:
                            index_left_ta, index_right_ta = index_second_highest_possible_peak, index_highest_peak
                        elif index_highest_peak < index_second_highest_possible_peak:
                            index_left_ta, index_right_ta = index_highest_peak, index_second_highest_possible_peak

                        value_left_ta, value_right_ta = gene_coverage_values_per_nt[index_left_ta], gene_coverage_values_per_nt[index_right_ta]

                        # Now, average left and right counts, and get three prime value.
                        count_site_in_gene = (value_left_ta + value_right_ta) / 2
                        three_prime = (index_right_ta + 1) / gene_length

                        gene_sites_stats.append(
                            {"l": index_left_ta, "r": index_right_ta, "l_value": value_left_ta, "r_value": value_right_ta, '3-primeness': three_prime,
                             'count': count_site_in_gene})

                        total_counts_of_sites_in_gene += count_site_in_gene

                    # CPM math
                    total_counts_of_sites_in_gene_normalized = (total_counts_of_sites_in_gene * 10e6) / total_read_counts_in_sample
                    # Disruption in gene math
                    mean_three_prime = numpy.mean(numpy.array([i['3-primeness'] for i in gene_sites_stats]))
                    # Calculate the proportion of insertions that fell below the outlier threshold, important for knowing if insertions are effective
                    below_threshold = round(sum([i["count"] for i in gene_sites_stats if i['3-primeness'] <= outliers_threshold])/total_counts_of_sites_in_gene * 100, 2)

                # Didnt find any TA site :( life is sad
                # Lets be sad and take the mean of whole split, whatever
                else:
                    pass
                mean_coverage += numpy.mean(split)

            # Check that counts are not too small... Think about dead bacteria that dont actually contribute to fitness.
            if mean_coverage < 3:
                mean_coverage = 0

        return {'gene_callers_id': gene_callers_id,
                'sample_name': sample_name,
                'gene_coverage_values_per_nt': gene_coverage_values_per_nt,
                'mean_coverage': float(mean_coverage),
                'insertions': total_counts_of_sites_in_gene,
                'insertions_normalized': total_counts_of_sites_in_gene_normalized,
                'mean_disruption': mean_three_prime,
                'below_disruption': below_threshold}


    def get_gene_level_coverage_stats(self, split_name, contigs_db, min_cov_for_detection=0, outliers_threshold=1.5,
                                      zeros_are_outliers=False, mode=None, gene_caller_ids_of_interest=set([])):
        """THE function that returns the gene level coverage stats.

        Please note that the function returns a tuple of two items.

        Returns
        =======
        gene_coverage_stats : dict
            A nested dictionary that for each gene caller id, and sample name, contains
            a dictionary that looks like this:
                >>> (...): {'gene_callers_id': 0,
                            'sample_name': 'USA_0114C_007D',
                            'mean_coverage': 0.050724637681159424,
                            'detection': 0.050724637681159424,
                            'non_outlier_mean_coverage': 0.050724637681159424,
                            'non_outlier_coverage_std': 0.2194348395612568,
                            'gene_coverage_values_per_nt': array([1, 1, 1, ..., 0, 0, 0], dtype=uint16),
                            'non_outlier_positions': array([ True,  True,  True, ...,  True,  True,  True])
                            }, { (...)
                >>>

        failed_gene_caller_ids : set
            A set of gene caller ids for which anvi'o failed to recover coverage
            information for any reason.
        """


        failed_gene_caller_ids = set([])

        # sanity check
        if not isinstance(gene_caller_ids_of_interest, set):
            raise ConfigError("`gene_caller_ids_of_interest` must be of type `set`")

        # recover split coverage values from the auxiliary data file
        if split_name not in self.split_coverage_values_per_nt_dict:
            if not self.auxiliary_profile_data_available:
                raise ConfigError("You are trying to recover gene coverage stats dict for a single split, but (1) "
                                  "the split is not described in split coverage values per nucleotide dicts, and (2) "
                                  "you don't seem to have access to the auxiliary data file :/")

            self.split_coverage_values_per_nt_dict[split_name] = self.split_coverage_values.get(split_name)

        split_coverage = self.split_coverage_values_per_nt_dict[split_name]

        # identify entry ids for genes in `split_name`
        genes_in_splits_entries = contigs_db.split_name_to_genes_in_splits_entry_ids[split_name]

        # we have to go back, Kate :(
        if not genes_in_splits_entries:
            return {}, failed_gene_caller_ids

        output = {}

        # we will go through each gene entry in the split
        for genes_in_splits_entry in genes_in_splits_entries:
            e = contigs_db.genes_in_splits[genes_in_splits_entry]
            gene_callers_id, gene_start, gene_stop = e['gene_callers_id'], e['start_in_split'], e['stop_in_split']
            gene_length = gene_stop - gene_start

            # if the user requested to work only with a set of genes, check whether we are working with
            # one of those. see https://github.com/merenlab/anvio/issues/865 for details.
            if len(gene_caller_ids_of_interest) and gene_callers_id not in gene_caller_ids_of_interest:
                continue

            if gene_length <= 0:
                raise ConfigError("What? :( How! The gene with the caller id '%d' has a length of %d :/ We are done "
                                  "here!" % (gene_callers_id, gene_length))

            output[gene_callers_id] = dict([(sample_name, dict([('mean_coverage', 0), ('gene_detection', 0)])) for sample_name in self.p_meta['samples']])

            # the magic happens here:
            for sample_name in self.p_meta['samples']:
                if self.inseq_stats:
                    gene_coverage_stats = self.get_gene_level_coverage_stats_entry_for_inseq(gene_callers_id=gene_callers_id,
                                                                                             split_coverage=split_coverage,
                                                                                             sample_name=sample_name,
                                                                                             gene_start=gene_start,
                                                                                             gene_stop=gene_stop,
                                                                                             gene_length=gene_length,
                                                                                             outliers_threshold=outliers_threshold)
                else:
                    gene_coverage_stats = self.get_gene_level_coverage_stats_entry_for_default(gene_callers_id=gene_callers_id,
                                                                                               split_coverage=split_coverage,
                                                                                               sample_name=sample_name,
                                                                                               gene_start=gene_start,
                                                                                               gene_stop=gene_stop,
                                                                                               gene_length=gene_length,
                                                                                               outliers_threshold=outliers_threshold)

                if not gene_coverage_stats:
                    # if we got a None here, it means there was something wrong with this gene caller.
                    failed_gene_caller_ids.add(gene_callers_id)

                output[gene_callers_id][sample_name] = gene_coverage_stats

        # in case there are genes that failed, we need to clean up the output before
        # we return it
        if len(failed_gene_caller_ids):
            for gene_callers_id in failed_gene_caller_ids:
                output.pop(gene_callers_id)

        # now good to go
        return output, failed_gene_caller_ids


    def get_blank_variability_dict(self):
        """Returns an empty variability dictionary to be filled elsewhere"""
        d = {}

        for sample_name in self.p_meta['samples']:
            d[sample_name] = {'variability': {0: {}, 1: {}, 2: {}, 3: {}}, 'competing_nucleotides': {}}

        return d

    def get_blank_indels_dict(self):
        """Returns an empty indels dictionary to be filled elsewhere"""
        d = {}

        for sample_name in self.p_meta['samples']:
            d[sample_name] = {'indels': {}}

        return d


    def get_variability_information_for_split(self, split_name, skip_outlier_SNVs=False, return_raw_results=False):
        if not split_name in self.split_names:
            raise ConfigError("get_variability_information_for_split: The split name '%s' does not seem to be "
                               "represented in this profile database. Are you sure you are looking for it "
                               "in the right database?" % split_name)

        self.progress.new('Recovering variability information for split', discard_previous_if_exists=True)
        self.progress.update('...')

        profile_db = ProfileDatabase(self.profile_db_path)
        split_variability_information = list(profile_db.db.get_some_rows_from_table_as_dict(t.variable_nts_table_name, '''split_name = "%s"''' % split_name, error_if_no_data=False).values())
        profile_db.disconnect()

        if return_raw_results:
            return split_variability_information

        # they want pretty stuff...
        d = self.get_blank_variability_dict()

        for e in split_variability_information:
            frequencies = utils.get_variabile_item_frequencies(e, engine='NT')
            e['n2n1ratio'], e['consensus'], e['departure_from_consensus'] = utils.get_consensus_and_departure_data(frequencies)

            if skip_outlier_SNVs and e['cov_outlier_in_contig']:
                continue

            d[e['sample_id']]['variability'][e['base_pos_in_codon']][e['pos']] = e['departure_from_reference']
            d[e['sample_id']]['competing_nucleotides'][e['pos']] = e

        self.progress.end()

        return d

    def get_indels_information_for_split(self, split_name):
        if not split_name in self.split_names:
            raise ConfigError("get_indels_information_for_split: The split name '%s' does not seem to be "
                               "represented in this profile database. Are you sure you are looking for it "
                               "in the right database?" % split_name)

        self.progress.new('Recovering indels information for split', discard_previous_if_exists=True)
        self.progress.update('...')

        profile_db = ProfileDatabase(self.profile_db_path)
        split_indels_information = list(profile_db.db.get_some_rows_from_table_as_dict(t.indels_table_name, '''split_name = "%s"''' % split_name, error_if_no_data=False).values())
        profile_db.disconnect()

        d = self.get_blank_indels_dict()

        for i, e in enumerate(split_indels_information):
            d[e['sample_id']]['indels'][i] = e

        self.progress.end()

        return d


    def init_items_additional_data(self):
        items_additional_data = TableForItemAdditionalData(self.args)
        self.items_additional_data_keys, self.items_additional_data_dict = items_additional_data.get()


    def get_split_coverages_dict(self, use_Q2Q3_coverages=False, splits_mode=False, report_contigs=False):
        """ Returns an item coverages dictionary.

            Parameters:
                - use_Q2Q3_coverages: by default this function will return mean coverage instead of Q2Q3
                        normalized mean coverages.
                - splits_mode: this is important to understand. by default the function will return splits
                        coverage of which are going to be matching to contigs from which they're coming
                        from. splits_mode makes it so that the coverage of each split is for itself.
                - report_contigs: when this parameter is True, then anvi'o reports a dictionary where each
                        item is a contig and not a split.

            Returns:
                - A dictionary where each key is a split name or a contig name and each value is a dictionary
                        of samples and their coverages.
        """

        if self.p_meta['blank']:
            raise ConfigError("The anvi'o profile db %s seems to be a blank profile database. Blank "
                              "profiles do not have any coverage values of any sorts, so whatever you "
                              "were trying to do with this database will not work :/" % (self.profile_db_path))

        if splits_mode and report_contigs:
            self.progress.reset()
            raise ConfigError("--splits-mode and --report-contigs flags are incompatible. Pick one.")

        coverage_data_of_interest = 'mean_coverage_Q2Q3' if use_Q2Q3_coverages else 'mean_coverage'

        table_name = coverage_data_of_interest + '_' + ('splits' if splits_mode else 'contigs')

        profile_db = ProfileDatabase(self.profile_db_path)
        split_coverages_dict, _ = profile_db.db.get_view_data(table_name)
        profile_db.disconnect()

        if report_contigs:
            # if we are here it means the user is asking for coverages for contigs, not splits. easy peasy.
            contigs_db = ContigsDatabase(self.contigs_db_path)
            split_parents = contigs_db.db.get_table_as_dict(t.splits_info_table_name, columns_of_interest=['contig', 'parent'])
            contigs_db.disconnect()

            contig_coverages_dict = {}

            for split_name in split_coverages_dict:
                contig_name = split_parents[split_name]['parent']

                if contig_name in contig_coverages_dict:
                    continue

                contig_coverages_dict[contig_name] = split_coverages_dict[split_name]

            return contig_coverages_dict
        else:
            return split_coverages_dict


    def init_collection_profile(self, collection_name, calculate_Q2Q3_carefully=False):
        profile_db = ProfileDatabase(self.profile_db_path, quiet=True)

        # we only have a self.collections instance if the profile super has been inherited by summary super class.
        # the initialization of a collection profile should only be done through that module anyway. so we are
        # being cruel here, and sending the programmer back.
        if not hasattr(self, 'collections'):
            raise ConfigError("You are lost :/ You can only call `init_collection_profile` through an instance of "
                              "the `SummarizerSuperClass`. Go back and come another way.")

        # get trimmed collection and bins_info dictionaries
        collection, bins_info, self.split_names_in_profile_db_but_not_binned \
                    = self.collections.get_trimmed_dicts(collection_name, self.split_names)

        for bin_id in collection:
            self.collection_profile[bin_id] = {}

        table_names = [] if self.p_meta['blank'] else constants.essential_data_fields_for_anvio_profiles

        samples_template = dict([(s, []) for s in self.p_meta['samples']])

        if calculate_Q2Q3_carefully:
            self.run.warning("The anvi'o sumarizer class is instructed (hopefully by you) to calculate Q2Q3 mean "
                             "coverages carefully. This means, depending on the size of your dataset and the number "
                             "of contigs in your bins this step can take much much longer than usual, since anvi'o "
                             "will have to do a lot of sorting of very large arrays. But then you will get the best "
                             "mean coverage values for your populations (so brace yourself).",
                             header="💀 THINGS WILL TAKE LONGER 💀")

        self.progress.new(f"Collection profile for '{collection_name}'")
        for table_name in table_names:
            # if SNVs are not profiled, skip the `variability` table
            if table_name == 'variability' and not self.p_meta['SNVs_profiled']:
                continue

            self.progress.update(f"Populating view '{table_name}'")
            table_data, _ = profile_db.db.get_view_data(f'{table_name}_splits')

            for bin_id in collection:
                if calculate_Q2Q3_carefully and table_name == 'mean_coverage_Q2Q3':
                    self.collection_profile[bin_id][table_name] = {}
                    # we need to do something specific here.
                    for sample_name in samples_template:
                        nucleotide_level_coverage_values = numpy.array([])
                        for split_name in collection[bin_id]:
                            nucleotide_level_coverage_values = numpy.append(nucleotide_level_coverage_values, self.split_coverage_values.get(split_name)[sample_name])
                        stats = utils.CoverageStats(nucleotide_level_coverage_values)
                        self.collection_profile[bin_id][table_name][sample_name] = stats.mean_Q2Q3
                else:
                    # populate averages per bin
                    averages = copy.deepcopy(samples_template)

                    # These weights are used to properly account for differences in split lengths.
                    # Consider table_name == 'mean_coverage', for a bin with 2 splits. Without
                    # weighting, if one split is length 100 with coverage 100 and the other is length
                    # 900 wth coverage 500, the mean_coverage for this bin is (100 + 500)/2 = 300. But
                    # more accurately, mean_coverage of this bin is 100*[100/1000] + 500*[900/1000] =
                    # 460
                    weights = []

                    for split_name in collection[bin_id]:
                        if split_name not in table_data:
                            continue

                        weights.append(self.splits_basic_info[split_name]['length'])

                        for sample_name in samples_template:
                            averages[sample_name].append(table_data[split_name][sample_name])

                    # finalize averages per bin:
                    for sample_name in samples_template:
                        # weights is automatically normalized in numpy.average such that sum(weights) == 1
                        averages[sample_name] = numpy.average([a or 0 for a in averages[sample_name]], weights=weights)

                    self.collection_profile[bin_id][table_name] = averages

        # generating precent recruitment of each bin plus __splits_not_binned__ in each sample:
        coverage_table_data, _ = profile_db.db.get_view_data('mean_coverage_splits')

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


    def load_views(self, split_names_of_interest=None, omit_parent_column=False):
        profile_db = ProfileDatabase(self.profile_db_path)

        views_table = profile_db.db.get_table_as_dict(t.views_table_name)

        # if SNVs are not profiled, we should not have a view for `variability`. See the issue #1845.
        # At this point we could remove this comment and the next three lines of code since the variability
        # will never be in `views_table` if SNVs were not profiled, but I'm not deleting it anyway because
        # I guess I'm slowly becoming a code hoarder :(
        if not self.p_meta['SNVs_profiled']:
            if 'variability' in views_table:
                views_table.pop('variability')

        self.progress.new('Loading views%s' % (' for %d items' % len(split_names_of_interest) if split_names_of_interest else ''))
        for view in views_table:
            self.progress.update('for %s' % view)
            table_name = views_table[view]['target_table']

            data, header = profile_db.db.get_view_data(table_name,
                                                       split_names_of_interest=split_names_of_interest,
                                                       splits_basic_info=(None if omit_parent_column else self.splits_basic_info))

            self.views[view] = {'table_name': table_name,
                                'header': header,
                                'dict': data}

        self.progress.end()
        profile_db.disconnect()


class DatabasesMetaclass(ProfileSuperclass, ContigsSuperclass, object):
    """Essential data to load for a given run"""
    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        filesnpaths.is_file_exists(args.contigs_db)
        filesnpaths.is_file_exists(args.profile_db)

        utils.is_profile_db_and_contigs_db_compatible(args.profile_db, args.contigs_db)

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
        self.db_type = 'profile'

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if not os.path.exists(self.db_path):
            return

        self.meta = dbi(self.db_path, expecting=self.db_type).get_self_table()

        for key in ['min_contig_length', 'SNVs_profiled', 'SCVs_profiled', 'INDELs_profiled',
                    'merged', 'blank', 'items_ordered', 'report_variability_full', 'num_contigs',
                    'min_coverage_for_variability', 'max_contig_length', 'num_splits',
                    'total_length', 'skip_edges_for_variant_profiling']:
            try:
                self.meta[key] = int(self.meta[key])
            except:
                pass

        self.samples = set([s.strip() for s in self.meta['samples'].split(',')])
        for key in ['min_percent_identity', 'min_indel_fraction']:
            try:
                self.meta[key] = float(self.meta[key])
            except:
                pass


        # open the database
        self.db = db.DB(self.db_path, anvio.__profile__version__)

        self.run.info('Profile database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
        self.run.info('Samples', self.meta['samples'], quiet=self.quiet)


    def touch(self):
        """Creates an empty profile database on disk, and sets `self.db` to access to it.

        At some point self.db.disconnect() must be called to complete the creation of the new db."""

        is_db_ok_to_create(self.db_path, 'profile')

        self.db = db.DB(self.db_path, anvio.__profile__version__, new_database=True)

        # creating empty default tables
        self.db.create_table(t.item_additional_data_table_name, t.item_additional_data_table_structure, t.item_additional_data_table_types)
        self.db.create_table(t.item_orders_table_name, t.item_orders_table_structure, t.item_orders_table_types)
        self.db.create_table(t.layer_additional_data_table_name, t.layer_additional_data_table_structure, t.layer_additional_data_table_types)
        self.db.create_table(t.layer_orders_table_name, t.layer_orders_table_structure, t.layer_orders_table_types)
        self.db.create_table(t.variable_nts_table_name, t.variable_nts_table_structure, t.variable_nts_table_types)
        self.db.create_table(t.variable_codons_table_name, t.variable_codons_table_structure, t.variable_codons_table_types)
        self.db.create_table(t.indels_table_name, t.indels_table_structure, t.indels_table_types)
        self.db.create_table(t.views_table_name, t.views_table_structure, t.views_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)
        self.db.create_table(t.protein_abundances_table_name, t.protein_abundances_table_structure, t.protein_abundances_table_types)
        self.db.create_table(t.metabolite_abundances_table_name, t.metabolite_abundances_table_structure, t.metabolite_abundances_table_types)

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


class GenesDatabase:
    """To create an empty genes database and/or access one.

       This is an optional database that is always linked to a profile database, collection
       name, and a bin id. The purpose of it is to store coverage an detection statistics of
       genes within a bin across samples desribed within a profile database.

       Besides the performance gain during the initialization of gene coverage values, this
       database serves as a way to store collections and states that are specific to the
       '--gene-mode' operations"""

    def __init__(self, db_path, run=run, progress=progress, quiet=True):
        self.db = None
        self.db_path = db_path
        self.db_type = 'genes'

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if not os.path.exists(self.db_path):
            return

        self.meta = dbi(self.db_path, expecting=self.db_type).get_self_table()

        for key in ['min_cov_for_detection', 'zeros_are_outliers', 'gene_level_coverages_stored', 'items_ordered']:
            try:
                self.meta[key] = int(self.meta[key])
            except:
                pass

        for key in ['outliers_threshold']:
            try:
                self.meta[key] = float(self.meta[key])
            except:
                pass
        self.db = db.DB(self.db_path, anvio.__genes__version__)

        self.run.info('Genes database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)


    def touch(self):
        """Creates an empty genes database on disk, and sets `self.db` to access to it."""

        is_db_ok_to_create(self.db_path, 'genes')

        filesnpaths.gen_output_directory(os.path.dirname(self.db_path))

        self.db = db.DB(self.db_path, anvio.__genes__version__, new_database=True)

        # creating empty default tables
        self.db.create_table(t.item_additional_data_table_name, t.item_additional_data_table_structure, t.item_additional_data_table_types)
        self.db.create_table(t.item_orders_table_name, t.item_orders_table_structure, t.item_orders_table_types)
        self.db.create_table(t.layer_additional_data_table_name, t.layer_additional_data_table_structure, t.layer_additional_data_table_types)
        self.db.create_table(t.layer_orders_table_name, t.layer_orders_table_structure, t.layer_orders_table_types)
        self.db.create_table(t.gene_level_coverage_stats_table_name, t.gene_level_coverage_stats_table_structure, t.gene_level_coverage_stats_table_types)
        self.db.create_table(t.gene_level_inseq_stats_table_name, t.gene_level_inseq_stats_table_structure, t.gene_level_inseq_stats_table_types)
        self.db.create_table(t.collections_info_table_name, t.collections_info_table_structure, t.collections_info_table_types)
        self.db.create_table(t.collections_bins_info_table_name, t.collections_bins_info_table_structure, t.collections_bins_info_table_types)
        self.db.create_table(t.collections_splits_table_name, t.collections_splits_table_structure, t.collections_splits_table_types)
        self.db.create_table(t.collections_contigs_table_name, t.collections_contigs_table_structure, t.collections_contigs_table_types)
        self.db.create_table(t.states_table_name, t.states_table_structure, t.states_table_types)

        return self.db


    def create(self, meta_values={}):
        self.touch()

        self.db.set_meta_value('db_type', 'genes')

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value('gene_level_coverages_stored', False)
        self.db.set_meta_value('items_ordered', False)

        self.disconnect()

        self.run.info('Genes database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


class PanDatabase:
    """To create an empty pan database, and/or access to one."""
    def __init__(self, db_path, run=run, progress=progress, quiet=True):
        self.db = None
        self.db_path = db_path
        self.db_type = 'pan'

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.init()


    def init(self):
        if not os.path.exists(self.db_path):
            return

        self.meta = dbi(self.db_path, expecting=self.db_type).get_self_table()
        for key in ['num_genomes', 'gene_cluster_min_occurrence', 'use_ncbi_blast', 'exclude_partial_gene_calls',
                    'num_gene_clusters', 'num_genes_in_gene_clusters', 'gene_alignments_computed', 'items_ordered',
                    'reaction_network_ko_annotations_hash', 'reaction_network_kegg_database_release',
                    'reaction_network_modelseed_database_sha']:
            try:
                self.meta[key] = int(self.meta[key])
            except:
                pass

        for key in ['min_percent_identity', 'minbit', 'mcl_inflation',
                    'reaction_network_consensus_threshold', 'reaction_network_discard_ties']:
            try:
                self.meta[key] = float(self.meta[key])
            except:
                pass

        self.internal_genomes = [s.strip() for s in self.meta['internal_genome_names'].split(',')]
        self.external_genomes = [s.strip() for s in self.meta['external_genome_names'].split(',')]
        self.genomes = self.internal_genomes + self.external_genomes

        # open the database
        self.db = db.DB(self.db_path, anvio.__pan__version__)

        self.run.info('Pan database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
        self.run.info('Genomes', '%d found' % len(self.genomes), quiet=self.quiet)


    def touch(self):
        is_db_ok_to_create(self.db_path, self.db_type)

        self.db = db.DB(self.db_path, anvio.__pan__version__, new_database=True)

        # creating empty default tables for pan specific operations:
        self.db.create_table(t.pan_gene_clusters_table_name, t.pan_gene_clusters_table_structure, t.pan_gene_clusters_table_types)
        self.db.create_table(t.pan_reaction_network_reactions_table_name, t.pan_reaction_network_reactions_table_structure, t.pan_reaction_network_reactions_table_types)
        self.db.create_table(t.pan_reaction_network_metabolites_table_name, t.pan_reaction_network_metabolites_table_structure, t.pan_reaction_network_metabolites_table_types)
        self.db.create_table(t.pan_reaction_network_kegg_table_name, t.pan_reaction_network_kegg_table_structure, t.pan_reaction_network_kegg_table_types)

        # creating empty default tables for standard anvi'o pan dbs
        self.db.create_table(t.item_additional_data_table_name, t.item_additional_data_table_structure, t.item_additional_data_table_types)
        self.db.create_table(t.item_orders_table_name, t.item_orders_table_structure, t.item_orders_table_types)
        self.db.create_table(t.layer_additional_data_table_name, t.layer_additional_data_table_structure, t.layer_additional_data_table_types)
        self.db.create_table(t.layer_orders_table_name, t.layer_orders_table_structure, t.layer_orders_table_types)
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

        # know thyself
        self.db.set_meta_value('db_type', 'pan')

        self.disconnect()

        self.run.info('Pan database', 'A new database, %s, has been created.' % (self.db_path), quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


class ContigsDatabase:
    """To create an empty contigs database and/or access one."""

    def __init__(self, db_path, run=run, progress=progress, quiet=True, skip_init=False):
        self.db = None
        self.db_path = db_path
        self.db_type = 'contigs'

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.meta = {}

        if not skip_init:
            self.init()


    def init(self):
        if not os.path.exists(self.db_path):
            return

        self.meta = dbi(self.db_path, expecting=self.db_type).get_self_table()

        try:
            for key in ['split_length', 'kmer_size', 'total_length', 'num_splits', 'num_contigs',
                        'genes_are_called', 'splits_consider_gene_calls', 'scg_taxonomy_was_run',
                        'trna_taxonomy_was_run', 'external_gene_calls', 'external_gene_amino_acid_seqs', 'skip_predict_frame']:
                self.meta[key] = int(self.meta[key])
        except KeyError:
            raise ConfigError("Oh no :( There is a contigs database here at '%s', but it seems to be broken :( It is very "
                              "likely that the process that was trying to create this database failed, and left behind "
                              "this unfinished thingy (if you would like to picture its state you should imagine the baby "
                              "Voldemort at King's Cross). Well, anvi'o believes it is best if you make it go away with "
                              "fire, and try whatever you were trying before you got this error one more time with a "
                              "proper contigs database. End of sad news. Bye now." % self.db_path)

        # open database
        self.db = db.DB(self.db_path, anvio.__contigs__version__)

        self.meta['gene_callers'] = self.db.get_frequencies_of_values_from_a_column(t.genes_in_contigs_table_name, 'source')[::-1]
        self.meta['gene_function_sources'] = [s.strip() for s in self.meta['gene_function_sources'].split(',')] if self.meta['gene_function_sources'] else None

        # set a project name for the contigs database without any funny
        # characters to make sure it can be used programmatically later.
        self.meta['project_name_str'] = self.meta['project_name'].strip().translate({ord(c): "_" for c in "\"'!@#$%^&*()[]{};:,./<>?\|`~-=_+ "}).replace('__', '_').strip('_') \
                                if self.meta['project_name'] else '___'.join(['UNKNOWN', self.meta['contigs_db_hash']])

        if 'creation_date' not in self.meta:
            raise ConfigError("The contigs database ('%s') seems to be corrupted :/ This happens if the process that "
                               "that generates the database ends prematurely. Most probably, you will need to generate "
                               "the contigs database from scratch. Sorry!" % (self.db_path))

        self.run.info('Contigs database', 'An existing database, %s, has been initiated.' % self.db_path, quiet=self.quiet)
        self.run.info('Number of contigs', self.meta['num_contigs'], quiet=self.quiet)
        self.run.info('Number of splits', self.meta['num_splits'], quiet=self.quiet)
        self.run.info('Total number of nucleotides', self.meta['total_length'], quiet=self.quiet)
        self.run.info('Split length', self.meta['split_length'], quiet=self.quiet)


    def __del__(self):
        if self.db:
            self.db.disconnect()


    def get_date(self):
        return time.time()


    def get_hash(self):
        return 'hash' + str('%08x' % random.randrange(16**8))


    def touch(self, db_variant='unknown'):
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
        self.db.create_table(t.reaction_network_reactions_table_name, t.reaction_network_reactions_table_structure, t.reaction_network_reactions_table_types)
        self.db.create_table(t.reaction_network_metabolites_table_name, t.reaction_network_metabolites_table_structure, t.reaction_network_metabolites_table_types)
        self.db.create_table(t.reaction_network_kegg_table_name, t.reaction_network_kegg_table_structure, t.reaction_network_kegg_table_types)
        self.db.create_table(t.gene_amino_acid_sequences_table_name, t.gene_amino_acid_sequences_table_structure, t.gene_amino_acid_sequences_table_types)
        self.db.create_table(t.splits_info_table_name, t.splits_info_table_structure, t.splits_info_table_types)
        self.db.create_table(t.contigs_info_table_name, t.contigs_info_table_structure, t.contigs_info_table_types)
        self.db.create_table(t.nt_position_info_table_name, t.nt_position_info_table_structure, t.nt_position_info_table_types)
        self.db.create_table(t.scg_taxonomy_table_name, t.scg_taxonomy_table_structure, t.scg_taxonomy_table_types)
        self.db.create_table(t.trna_taxonomy_table_name, t.trna_taxonomy_table_structure, t.trna_taxonomy_table_types)
        self.db.create_table(t.nucleotide_additional_data_table_name, t.nucleotide_additional_data_table_structure, t.nucleotide_additional_data_table_types)
        self.db.create_table(t.amino_acid_additional_data_table_name, t.amino_acid_additional_data_table_structure, t.amino_acid_additional_data_table_types)

        if db_variant == 'trnaseq':
            self.db.create_table(t.trna_seed_feature_table_name, t.trna_seed_feature_table_structure, t.trna_seed_feature_table_types)

        return self.db


    def remove_gene_calls(self, gene_caller_ids_to_remove):
        if not self.meta['genes_are_called']:
            raise ConfigError("Well .. anvi'o can't remove gene calls from this contigs database since "
                              "genes were never called for it at the first place :(")

        if not isinstance(gene_caller_ids_to_remove, list):
            raise ConfigError("Gene caller ids to remove must be of type `list`.")

        gene_calls_in_db = set(self.db.get_single_column_from_table(t.genes_in_contigs_table_name,
                                                                    "gene_callers_id"))

        missing_gene_calls = [g for g in gene_caller_ids_to_remove if g not in gene_calls_in_db]
        if len(missing_gene_calls):
            raise ConfigError("Your contigs database has %d genes, but it's missing %d of %d gene calls "
                              "you want to remove from it :/ This doesn't make sense. Here is one of those "
                              "gene calls that were not in your database: %d" % \
                                    (len(gene_calls_in_db), len(missing_gene_calls), len(gene_caller_ids_to_remove), gene_caller_ids_to_remove[-1]))

        self.run.warning('%d gene calls %d is being removed from your contigs '
                         'database' % (len(gene_calls_in_db), len(gene_caller_ids_to_remove)))

        # tables from which the gene calls  to remove gene calls from:
        # TODO: Add functionality to remove data from reaction/metabolite tables associated with GCIDs.
        tables_dict = {
                    t.gene_function_calls_table_name: ('gene_callers_id', gene_caller_ids_to_remove),
                    t.gene_amino_acid_sequences_table_name: ('gene_callers_id', gene_caller_ids_to_remove),
                    t.genes_in_contigs_table_name: ('gene_callers_id', gene_caller_ids_to_remove),
                    t.genes_in_splits_table_name: ('gene_callers_id', gene_caller_ids_to_remove),
                    t.genes_taxonomy_table_name: ('gene_callers_id', gene_caller_ids_to_remove),
                    t.hmm_hits_table_name: ('gene_callers_id', gene_caller_ids_to_remove),
                }

        self.remove_data_from_db(tables_dict)


    def list_gene_caller_sources(self):
        self.run.warning(None, header="AVAILABLE GENE CALLERS", lc="green")

        gene_caller_sources = self.meta['gene_callers']

        if not len(gene_caller_sources):
            self.run.info_single("This contigs db does not have any gene calls :/")
        else:
            for gene_caller_source, num_genes in gene_caller_sources:
                self.run.info_single("'%s' (%s gene calls)" % (gene_caller_source, pp(num_genes)),
                                     nl_after = 1 if gene_caller_source == gene_caller_sources[-1][0] else 0)


    def list_function_sources(self):
        self.run.warning(None, header="AVAILABLE FUNCTIONAL ANNOTATION SOURCES", lc="green")

        gene_function_sources = sorted(self.meta['gene_function_sources']) if self.meta['gene_function_sources'] else []

        if not len(gene_function_sources):
            self.run.info_single('No functional annotations found in this contigs database :/', nl_after=1, mc='red')
        else:
            for source in gene_function_sources:
                num_annotations = self.db.get_row_counts_from_table(t.gene_function_calls_table_name, where_clause=f'source="{source}"')
                self.run.info_single(f'{source} ({pp(num_annotations)} annotations)', nl_after = 1 if source == gene_function_sources[-1] else 0)


    def list_available_hmm_sources(self):
        self.run.warning(None, header="AVAILABLE HMM SOURCES", lc="green")

        hmm_sources_dict = self.db.get_table_as_dict(t.hmm_hits_info_table_name)
        hmm_sources = sorted(list(hmm_sources_dict.keys()))

        if not len(hmm_sources_dict):
            self.run.info_single("This contigs db does not have HMMs :/")
        else:
            for source in hmm_sources:
                num_models = len(hmm_sources_dict[source]['genes'].split(','))
                num_hits = self.db.get_row_counts_from_table(t.hmm_hits_table_name, where_clause=f'source="{source}"')
                self.run.info_single(f"'{source}' ({P('model', num_models)} with {P('hit', num_hits)})",
                                     nl_after = 1 if source == hmm_sources[-1] else 0)


    def remove_data_from_db(self, tables_dict):
        """This is quite an experimental function to clean up tables in contigs databases. Use with caution.

           The expected tables dict should follow this structure:

           >>> tables = {
                            table_name: ('column_name', values_to_remove_list),
                            table_name: ('column_name', values_to_remove_list),
                            (...)
                        }
        """

        table_names_in_db = self.db.get_table_names()

        missing_table_names = [table_name for table_name in tables_dict if table_name not in table_names_in_db]
        if len(missing_table_names):
            raise ConfigError("Your tables dict contains table names that are not in the database. Like "
                              "these ones here: '%s'. WHO SENT YOU HERE? WHAT ARE YOU UP TO?" % (', '.join(missing_table_names)))

        for table_name in tables_dict:
            column = tables_dict[table_name][0]
            variables_to_remove = tables_dict[table_name][1]
            where_clause = '%s IN (%s)' % (column, ','.join(['{}'.format(x) for x in variables_to_remove]))

            self.db.remove_some_rows_from_table(table_name, where_clause)


    def create(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        contigs_fasta = A('contigs_fasta')
        db_variant = A('db_variant') or 'unknown'
        project_name = A('project_name')
        description_file_path = A('description')
        split_length = A('split_length')
        kmer_size = A('kmer_size')
        skip_gene_calling = A('skip_gene_calling')
        external_gene_calls_file_path = A('external_gene_calls')
        skip_mindful_splitting = A('skip_mindful_splitting')
        ignore_internal_stop_codons = A('ignore_internal_stop_codons')
        skip_predict_frame= A('skip_predict_frame')
        prodigal_translation_table = A('prodigal_translation_table')
        prodigal_single_mode = A('prodigal_single_mode')

        if external_gene_calls_file_path:
            filesnpaths.is_proper_external_gene_calls_file(external_gene_calls_file_path)

        if external_gene_calls_file_path and skip_gene_calling:
            raise ConfigError("You provided a file for external gene calls, and used requested gene calling to be "
                              "skipped. Please make up your mind.")

        if (external_gene_calls_file_path or skip_gene_calling) and prodigal_translation_table:
            raise ConfigError("You asked anvi'o to %s, yet you set a specific translation table for prodigal. These "
                              "parameters do not make much sense and anvi'o is kindly asking you to make up your "
                              "mind." % ('skip gene calling' if skip_gene_calling else 'use external gene calls'))

        if skip_predict_frame and not external_gene_calls_file_path:
            raise ConfigError("The flag `--skip-predict-frame` is only relevant if you are providing external gene calls. ")

        filesnpaths.is_file_fasta_formatted(contigs_fasta)
        contigs_fasta = os.path.abspath(contigs_fasta)

        # let's see if the user has provided extenral gene calls file with amino
        # acid sequences:
        external_gene_calls_include_amino_acid_sequences = external_gene_calls_file_path and 'aa_sequence' in utils.get_columns_of_TAB_delim_file(external_gene_calls_file_path)

        # let the user see what's up
        self.run.info('Input FASTA file', contigs_fasta)

        if not project_name:
            project_name = '.'.join(os.path.basename(os.path.abspath(contigs_fasta)).split('.')[:-1])
            project_name = project_name.translate({ord(c): "_" for c in "\"'!@#$%^&*()[]{};:,./<>?\|`~-=_+ "}).replace('__', '_')

            if project_name:
                self.run.warning(f"You are generating a new anvi'o contigs database, but you are not specifying a "
                                 f"project name for it. FINE. Anvi'o, in desperation, will use the input file name "
                                 f"to set the project name for this contigs database (i.e., '{project_name}'). If you are not "
                                 f"happy with that, feel free to kill and restart this process. If you are not happy "
                                 f"with this name, but you don't like killing things either, maybe next time you "
                                 f"should either name your FASTA files better, or use the `--project-name` parameter "
                                 f"to set your desired name.", "Anvi'o made things up for you", lc="green")
            else:
                raise ConfigError("Sorry, you must provide a project name for your contigs database :/ Anvi'o tried "
                                  "to make up one, but failed.")

        self.run.info('Name', project_name, mc='green')
        self.run.info('Description', os.path.abspath(description_file_path) if description_file_path else 'No description is given', mc='green')

        if description_file_path:
            filesnpaths.is_file_plain_text(description_file_path)
            description = open(os.path.abspath(description_file_path), 'r').read()
        else:
            description = ''

        # go through the FASTA file to make sure there are no surprises with deflines, sequence
        # lengths, or sequence characters.
        self.progress.new('Checking deflines and contig lengths')
        self.progress.update('tick tock ...')
        fasta = u.SequenceSource(contigs_fasta)

        # we only A, C, T, G, N, a, c, t, g, n
        character_regex = re.compile(r'^[ACGTNactgn]+$')

        while next(fasta):
            if not utils.check_contig_names(fasta.id, dont_raise=True):
                self.progress.end()
                raise ConfigError("At least one of the deflines in your FASTA File does not comply with the 'simple deflines' "
                                  "requirement of anvi'o. You can either use the script `anvi-script-reformat-fasta` to take "
                                  "care of this issue, or read this section in the tutorial to understand the reason behind "
                                  "this requirement (anvi'o is very upset for making you do this): %s" % \
                                       ('http://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file'))

            if len(fasta.seq) < kmer_size:
                self.progress.end()
                raise ConfigError("At least one of the contigs in your input FASTA '%s' is shorter than the k-mer size. The k "
                                  "is %d, and your contig is like %d :/ Anvi'o will not judge you for whatever you are doing "
                                  "with such short contigs, but the length of each contig must be at least as long as your `k` for "
                                  "k-mer analyis. You can use the script `anvi-script-reformat-fasta` to get rid of very short "
                                  "contigs if you like." % (contigs_fasta, kmer_size, len(fasta.seq)))

            if not bool(character_regex.search(fasta.seq)):
                raise ConfigError(f"Tough. {fasta.id} contains characters that are not any of A, C, T, G, N, a, c, t, g, n. To "
                                  f"save you from later headaches, anvi'o will not make a database. We recommend you identify "
                                  f"what's up with your sequences. If you want to continue, first reformat the FASTA file "
                                  f"(you could use anvi-script-reformat-fasta with the parameter --seq-type NT, and that would fix the problem)")

            try:
                int(fasta.id)
                is_int = True
            except:
                is_int = False

            if is_int:
                self.progress.end()
                raise ConfigError("At least one of the contigs in your FASTA file (well, this one to be precise: '%s') looks like "
                                  "a number. For reasons we can't really justify, anvi'o does not like those numeric names, and hereby "
                                  "asks you to make sure every contig name contains at least one alphanumeric character :/ Meanwhile we, "
                                  "the anvi'o developers, are both surprised by and thankful for your endless patience with such eccentric "
                                  "requests. You the real MVP." % fasta.id)

        fasta.close()
        self.progress.end()

        all_ids_in_FASTA = utils.get_all_ids_from_fasta(contigs_fasta)
        total_number_of_contigs = len(all_ids_in_FASTA)
        if total_number_of_contigs != len(set(all_ids_in_FASTA)):
            raise ConfigError("Every contig in the input FASTA file must have a unique ID. You know...")

        if split_length is None:
            raise ConfigError("Creating a new contigs database requires split length information to be "
                               "provided. But the ContigsDatabase class was called to create one without this "
                               "bit of information. Not cool.")

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
            raise ConfigError("We like our k-mer sizes between 2 and 8, sorry! (but then you can always change the "
                               "source code if you are not happy to be told what you can't do, let us know how it goes!).")

        if split_length < kmer_size + 1:
            raise ConfigError("Split size must be at least your k-mer size +1 (so in your case it can't be anything less "
                              "than %d)." % (kmer_size + 1))

        if skip_gene_calling:
            skip_mindful_splitting = True

        # create a blank contigs database on disk, and set the self.db
        self.touch(db_variant)

        # know thyself
        self.db.set_meta_value('db_type', 'contigs')
        self.db.set_meta_value('db_variant', db_variant)
        self.db.set_meta_value('project_name', project_name)
        self.db.set_meta_value('description', description)

        # this will be the unique information that will be passed downstream whenever this db is used:
        contigs_db_hash = self.get_hash()
        self.db.set_meta_value('contigs_db_hash', contigs_db_hash)

        # set split length variable in the meta table
        self.db.set_meta_value('split_length', split_length)

        # first things first: do the gene calling on contigs. this part is important. we are doing the
        # gene calling first, so we understand wher genes start and end. this information will guide the
        # arrangement of the breakpoint of splits
        genes_in_contigs_dict = {}
        contig_name_to_gene_start_stops = {}
        if not skip_gene_calling:
            # temporarily disconnect to perform gene calls
            self.db.disconnect()

            gene_calls_tables = TablesForGeneCalls(self.db_path, contigs_fasta, args=args, run=self.run, progress=self.progress, debug=anvio.DEBUG)

            # if the user provided a file for external gene calls, use it. otherwise do the gene calling yourself.
            if external_gene_calls_file_path:
                try:
                    gene_calls_tables.use_external_gene_calls_to_populate_genes_in_contigs_table(
                        input_file_path=external_gene_calls_file_path,
                        ignore_internal_stop_codons=ignore_internal_stop_codons,
                        skip_predict_frame=skip_predict_frame,
                    )
                except ConfigError as e:
                    if os.path.exists(self.db_path):
                        os.remove(self.db_path)
                    raise ConfigError(e.clear_text())
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

        # print some information for the user
        self.run.warning(None, header="CONTIGS DB CREATE REPORT", lc="cyan")
        self.run.info('Split Length', pp(split_length))
        self.run.info('K-mer size', kmer_size)
        self.run.info('Skip gene calling?', skip_gene_calling)
        self.run.info('External gene calls provided?', True if external_gene_calls_file_path else False, mc='green')

        if external_gene_calls_file_path:
            self.run.info('External gene calls file have AA sequences?', external_gene_calls_include_amino_acid_sequences, mc='green')
            self.run.info('Proper frames will be predicted?', (not skip_predict_frame), mc='green')
        else:
            self.run.info('Is prodigal run in single mode?', ('YES' if prodigal_single_mode else 'NO'), mc='green')

        self.run.info('Ignoring internal stop codons?', ignore_internal_stop_codons)
        self.run.info('Splitting pays attention to gene calls?', (not skip_mindful_splitting))

        # here we will process each item in the contigs fasta file.
        fasta = u.SequenceSource(contigs_fasta)
        db_entries_contig_sequences = []

        contigs_kmer_table = KMerTablesForContigsAndSplits('kmer_contigs', k=kmer_size)
        splits_kmer_table = KMerTablesForContigsAndSplits('kmer_splits', k=kmer_size)
        nt_positions_table = TableForNtPositions()
        contigs_info_table = TableForContigsInfo(split_length)
        splits_info_table = TableForSplitsInfo()

        recovered_split_lengths = []

        # THE INFAMOUS GEN CONTGS DB LOOP (because it is so costly, we call it South Loop)
        self.progress.new('The South Loop', progress_total_items=total_number_of_contigs)
        fasta.reset()
        while next(fasta):
            self.progress.increment()

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
                nt_position_info_list = self.compress_nt_position_info(contig_name, contig_length, genes_in_contig, genes_in_contigs_dict)
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
        self.db.set_meta_value('gene_level_taxonomy_source', None)
        self.db.set_meta_value('gene_function_sources', None)
        self.db.set_meta_value('genes_are_called', (not skip_gene_calling))
        self.db.set_meta_value('external_gene_calls', True if external_gene_calls_file_path else False)
        self.db.set_meta_value('external_gene_amino_acid_seqs', True if external_gene_calls_include_amino_acid_sequences else False)
        self.db.set_meta_value('skip_predict_frame', True if skip_predict_frame else False)
        self.db.set_meta_value('splits_consider_gene_calls', (not skip_mindful_splitting))
        self.db.set_meta_value('scg_taxonomy_was_run', False)
        self.db.set_meta_value('scg_taxonomy_database_version', None)
        self.db.set_meta_value('trna_taxonomy_was_run', False)
        self.db.set_meta_value('trna_taxonomy_database_version', None)
        self.db.set_meta_value('reaction_network_ko_annotations_hash', None)
        self.db.set_meta_value('reaction_network_kegg_database_release', None)
        self.db.set_meta_value('reaction_network_modelseed_database_sha', None)
        self.db.set_meta_value('reaction_network_consensus_threshold', None)
        self.db.set_meta_value('reaction_network_discard_ties', None)
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
        self.run.info("Average split length (what anvi'o gave back)", (int(round(numpy.mean(recovered_split_lengths)))) \
                                                                        if recovered_split_lengths \
                                                                            else "(Anvi'o did not create any splits)", quiet=self.quiet)


    def compress_nt_position_info(self, contig_name, contig_length, genes_in_contig, genes_in_contigs_dict):
        """Compress info regarding each nucleotide position in a given contig into a small int

        Every nucleotide position is represented by four bits depending on whether they occur in a
        coding open reading frame, and which base they correspond to in a codon.

             0000
             ||||
             ||| \
             |||   Third codon pos?
             || \
             ||   Second codon pos?
             | \
             |   First codon pos?
              \
                Whether the position is in an noncoding gene (call_type = NONCODING or UNKNOWN)

        8: int('1000', 2); nt position is in an noncoding gene
        4: int('0100', 2); nt position is in a coding gene, and is at the 1st position in the codon
        2: int('0010', 2); nt position is in a coding gene, and is at the 2nd position in the codon
        1: int('0001', 2); nt position is in a coding gene, and is at the 3rd position in the codon
        0: int('0000', 2); nt position not in a gene

        Notes
        =====
        - This code could be much faster by populating and returning a numpy array rather than a
          list
        - Upon testing (2020/05/19), this function is about 8% of runtime for
          anvi-gen-contigs-database (ignoring ORF prediction) so vectorization of this function is
          probably not worth it at this point in time.
        """

        # first we create a list of zeros for each position of the contig
        nt_position_info_list = [0] * contig_length

        # to keep track of things that didn't work
        gene_caller_ids_that_failed = []

        coding = constants.gene_call_types['CODING']

        for gene_unique_id, start, stop in genes_in_contig:
            gene_call = genes_in_contigs_dict[gene_unique_id]

            if (start > contig_length) or (stop > contig_length):
                if not anvio.DEBUG:
                    os.remove(self.db_path)
                raise ConfigError(f"Something is wrong with your external gene calls file. It seems "
                                  f"that the gene with gene callers id {gene_unique_id}, on contig "
                                  f"{contig_name}, has positions that go beyond the length of the contig. "
                                  f"Specifically, the length of the contig is {contig_length}, but the "
                                  f"gene starts at position {start} and goes to position {stop}. We've "
                                  f"removed the partially-created contigs database for you (but you can "
                                  f"see it if you re-run your command with the `--debug` flag).")

            if gene_call['call_type'] != coding:
                for nt_position in range(start, stop):
                    nt_position_info_list[nt_position] = 8
                continue

            # if the gene stop is identical to the contig length, we have to carefully assess
            # this situation. if we simply say,
            #
            #   >>> if stop == contig_length:
            #   >>>     continue
            #
            # then contigs that are solely composed of genes (i.e., gene = contig) ends up being
            # treated as intergenic regions. but if we skip this step, or say something like
            # this,
            #
            #   >>> if stop == contig_length + 1:
            #   >>>     continue
            #
            # then anvi'o explodes and dies a fiery death as explained at,
            #
            #    https://github.com/merenlab/anvio/issues/1943
            #
            # SO MUCH HISTORY HERE, BUT WE PAY OUR RISPEKS TO THOSE WHO SHARE REPRODUCIBLE TEST CASES.
            if stop == contig_length: # checking fo #1661
                if stop - start == contig_length: # realizing that it is #1943
                    pass # thanking the gene and sending it along.
                else: # finding out that it actually is #1661
                    continue # NEXT

            try:
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
            except IndexError:
                # FIXME: Please see https://github.com/merenlab/anvio/issues/2020 to see why this was necessary.
                #        Which also shows it is important to address this differently to make sure we better
                #        support eukaryotic organisms:
                gene_caller_ids_that_failed.append(gene_unique_id)

                # As a solution here we simply mark those nucleotide positions as of unknown nature,
                for nt_position in range(start, stop):
                    nt_position_info_list[nt_position] = 8

                # and move on to the next gene on the contig
                continue

        if len(gene_caller_ids_that_failed):
            progress.reset()
            run.warning(f"There were some problmes while anvi'o was trying to identify which nucleotides occur in which codon "
                        f"positions. These problems occurred in the following gene calls: {', '.join([str(g) for g in gene_caller_ids_that_failed])}. "
                        f"Please read this message for more information: https://github.com/merenlab/anvio/issues/2020#issuecomment-1341028673.",
                        header=f"SNAFU WITH CONTIG `{contig_name}`")

        return nt_position_info_list


    def disconnect(self):
        self.db.disconnect()
        self.db = None


class TRNASeqDatabase:
    """Used to create and/or access a tRNA-seq database"""

    def __init__(self, db_path, run=terminal.Run(), progress=terminal.Progress(), quiet=True):
        self.db = None
        self.db_path = db_path
        self.db_type = 'trnaseq'

        self.run = run
        self.progress = progress
        self.quiet = quiet

        self.meta_int_keys = [] # metadata to be stored as an int
        self.meta_float_keys = [] # metadata to be stored as a float

        self.meta = {}

        self.init()


    def init(self):
        if not os.path.exists(self.db_path):
            self.db = None
            return

        self.meta = dbi(self.db_path, expecting=self.db_type).get_self_table()

        self.db = db.DB(self.db_path, anvio.__trnaseq__version__)

        for key in self.meta_int_keys:
            try:
                self.meta[key] = int(self.meta[key])
            except:
                pass

        for key in self.meta_float_keys:
            try:
                self.meta[key] = float(self.meta[key])
            except:
                pass

        self.run.info("tRNA-seq database", f"An existing database, {self.db_path}, has been initiated.", quiet=self.quiet)


    def touch(self):
        is_db_ok_to_create(self.db_path, self.db_type)

        self.db = db.DB(self.db_path, anvio.__trnaseq__version__, new_database=True)

        self.db.create_table(t.trnaseq_sequences_table_name, t.trnaseq_sequences_table_structure, t.trnaseq_sequences_table_types)
        self.db.create_table(t.trnaseq_feature_table_name, t.trnaseq_feature_table_structure, t.trnaseq_feature_table_types)
        self.db.create_table(t.trnaseq_unconserved_table_name, t.trnaseq_unconserved_table_structure, t.trnaseq_unconserved_table_types)
        self.db.create_table(t.trnaseq_unpaired_table_name, t.trnaseq_unpaired_table_structure, t.trnaseq_unpaired_table_types)
        self.db.create_table(t.trnaseq_trimmed_table_name, t.trnaseq_trimmed_table_structure, t.trnaseq_trimmed_table_types)
        self.db.create_table(t.trnaseq_normalized_table_name, t.trnaseq_normalized_table_structure, t.trnaseq_normalized_table_types)
        self.db.create_table(t.trnaseq_modified_table_name, t.trnaseq_modified_table_structure, t.trnaseq_modified_table_types)

        return self.db


    def create(self, meta_values={}):
        self.touch()

        for key in meta_values:
            self.db.set_meta_value(key, meta_values[key])

        self.db.set_meta_value('creation_date', time.time())
        self.db.set_meta_value(self.db_type + '_db_hash', 'hash' + str('%08x' % random.randrange(16**8)))

        # know thyself
        self.db.set_meta_value('db_type', self.db_type)

        self.disconnect()

        self.run.info("%s database" % self.db_type, "A new database, %s, has been created." % self.db_path, quiet=self.quiet)


    def disconnect(self):
        self.db.disconnect()


####################################################################################################
#
#     TABLES
#
####################################################################################################


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
            raise ConfigError("Some contig names you declared do not seem to be present in the contigs "
                               "database :(")

        for contig_name in contigs_of_interest:
            self.counts_dict[contig_name] = self.get_AA_counts_dict(contig_names=set([contig_name]))['counts']


    def __AA_counts_for_genes(self):
        filesnpaths.is_file_exists(self.genes_of_interest_file_path)

        try:
            genes_of_interest = [int(line.strip()) for line in open(self.genes_of_interest_file_path).readlines()]
        except:
            raise ConfigError("Gene call ids in your genes of interest file does not resemble anvi'o gene "
                               "call ids (I tried to int them, and it didn't work!)")

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


####################################################################################################
#
#     HELPER FUNCTIONS
#
####################################################################################################

def is_db_ok_to_create(db_path, db_type):
    if not db_path.lower().endswith('.db'):
        raise ConfigError("Please make sure the file name for your new %s db has a '.db' extension. Anvi'o developers "
                           "apologize for imposing their views on how anvi'o databases should be named, and are "
                           "humbled by your cooperation." % db_type)

    filesnpaths.is_output_file_writable(db_path, ok_if_exists=False)


def get_auxiliary_data_path_for_profile_db(profile_db_path):
    return os.path.join(os.path.dirname(profile_db_path), 'AUXILIARY-DATA.db')


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
    description = open(os.path.abspath(description_file_path), 'r').read()

    update_description_in_db(anvio_db_path, description, run=run)


def update_description_in_db(anvio_db_path, description, run=run):
    """Updates the description in an anvi'o database"""

    if not isinstance(description, str):
        raise ConfigError("Description parameter must be of type `string`.")

    db_type = utils.get_db_type(anvio_db_path)

    anvio_db = db.DB(anvio_db_path, None, ignore_version=True)
    anvio_db.remove_meta_key_value_pair('description')
    anvio_db.set_meta_value('description', description)
    anvio_db.disconnect()

    run.info_single("The anvi'o %s database has just been updated with a description that contains %d words "
                    "and %d characters." % (db_type, len(description.split()), len(description)))


def do_hierarchical_clustering_of_items(anvio_db_path, clustering_configs, split_names=[], database_paths={}, input_directory=None, default_clustering_config=None, \
                                distance=constants.distance_metric_default, linkage=constants.linkage_method_default, run=run, progress=progress):
    """This is just an orphan function that computes hierarchical clustering w results
       and calls the `add_items_order_to_db` function with correct input.

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

        add_items_order_to_db(anvio_db_path=anvio_db_path,
                              order_name=config_name,
                              order_data=newick,
                              distance=distance,
                              linkage=linkage,
                              make_default=config_name == default_clustering_config,
                              run=run)


def add_items_order_to_db(anvio_db_path, order_name, order_data, order_data_type_newick=True, distance=None,
                          linkage=None, make_default=False, additional_data=None, dont_overwrite=False,
                          check_names_consistency=False, run=run):
    """Adds a new items order into an appropriate anvi'o db

       Here is a FIXME for future, smarter generations. This function should go away,
       and what its doing should be handled by a new items_order class in tables/miscdata.
    """

    if order_data_type_newick and (not distance or not linkage):
        raise ConfigError("You are trying to add a newick-formatted clustering dendrogram to the database without providing "
                          "distance and linkage data that generated this dendrogram :/")

    if not order_data_type_newick and (distance or linkage):
        raise ConfigError("Distance and linkage variables are only relevant if you are trying to add a newick-formatted "
                          "clustering dendrogram. But your function call suggests you are not.")

    # let's learn who we are dealing with:
    db_type = utils.get_db_type(anvio_db_path)

    # replace clustering id with a text that contains distance and linkage information
    if order_data_type_newick:
        order_name = ':'.join([order_name, distance, linkage])
    else:
        order_name = ':'.join([order_name, 'NA', 'NA'])

    # check names consistency if the user asked for it
    if check_names_consistency:
        if order_data_type_newick:
            names_in_data = sorted(utils.get_names_order_from_newick_tree(order_data))
        else:
            names_in_data = sorted([n.strip() for n in order_data.split(',')])

        names_in_db = sorted(utils.get_all_item_names_from_the_database(anvio_db_path))

        if not len(names_in_db):
            raise ConfigError(f"Your {db_type} database does not have any item names stored, but whoever called this "
                              f"function asked for a check between the item names in the items order and item names in "
                              f"the database. Well. It will not happen. A proper way to deal with it is to set the "
                              f"variable `check_names_consistency` to `False` when calling this function. This message "
                              f"is not to be seen by a user, so if you are a user, it means someone screwed up something "
                              f"somewhere, and those of us at the anvi'o headquarters for snafu handling sincerely apologize "
                              f"sincerely on their behalf :(")

        names_in_db_not_in_data = set(names_in_db) - set(names_in_data)
        if names_in_db_not_in_data:
            raise ConfigError("Ehem. There is something wrong with the incoming items order data here :/ Basically,\
                               the names found in your input data do not match to the item names found in the\
                               database. For example, this item \"%s\" is in your database, but not in your input data\
                               " % next(iter(names_in_db_not_in_data)))

        names_in_data_not_in_db = set(names_in_data) - set(names_in_db)
        if names_in_data_not_in_db:
            raise ConfigError("Ehem. There is something wrong with the incoming items order data here :/ Basically,\
                               the names found in your input data do not match to the item names found in the\
                               database. For example, this item \"%s\" is in your input data, but not in your database\
                               " % next(iter(names_in_data_not_in_db)))

    # additional data is JSON formatted entry
    # for now it will only contain collapsed node information.
    # in future we may extend this column to include other annotations
    if not additional_data:
        additional_data = json.dumps({})
    else:
        additional_data = json.dumps(additional_data)

    anvio_db = DBClassFactory().get_db_object(anvio_db_path)

    if t.item_orders_table_name not in anvio_db.db.get_table_names():
        raise ConfigError("You can't add a new items order into this %s database (%s). You know why? Becasue it doesn't "
                          "have a table for 'item_order' :(" % (db_type, anvio_db_path))

    try:
        available_item_orders = anvio_db.db.get_meta_value('available_item_orders').split(',')
    except:
        available_item_orders = []

    if order_name in available_item_orders:
        if dont_overwrite:
            raise ConfigError("The order name '%s' is already in the database, and you are not allowed to overwrite that. "
                              "Probably it is time for you to come up with a new name?" % (order_name))
        else:
            run.warning('Clustering for "%s" is already in the database. It will be replaced with the new content.' % (order_name))

        anvio_db.db._exec('''DELETE FROM %s where name = "%s"''' % (t.item_orders_table_name, order_name))
    else:
        available_item_orders.append(order_name)

    anvio_db.db._exec('''INSERT INTO %s VALUES (?,?,?,?)''' % t.item_orders_table_name, tuple([order_name, 'newick' if order_data_type_newick else 'basic', order_data, additional_data]))

    anvio_db.db.set_meta_value('available_item_orders', ','.join(available_item_orders))
    anvio_db.db.set_meta_value('items_ordered', True)

    try:
        anvio_db.db.get_meta_value('default_item_order')
        default_item_order_is_set = True
    except:
        default_item_order_is_set = False

    if make_default or not default_item_order_is_set:
        anvio_db.db.set_meta_value('default_item_order', order_name)

    anvio_db.disconnect()

    run.info('New items order', '"%s" (type %s) has been added to the database...' % (order_name, 'newick' if order_data_type_newick else 'basic'))


def get_item_orders_from_db(anvio_db_path):
    anvio_db = DBClassFactory().get_db_object(anvio_db_path)

    utils.is_pan_or_profile_db(anvio_db_path, genes_db_is_also_accepted=True)

    if not anvio_db.meta['items_ordered']:
        return ([], {})

    if not anvio_db.meta['available_item_orders']:
        # this means the database thinks that the items are ordered, but in fact
        # there are no item orders anwyhere to be found. this means that the
        # clustering failed somehow, and we're dealing with a database that is not
        # quite accurate. but since we can't fix the problem that led the database
        # to find itself in this state, we can do the next best thing and avoid
        # going any further and return an empty set as if items are NOT ordered
        # (special thanks to Rose Kantor for helping us to diagnose this):
        return ([], {})

    available_item_orders = sorted([s.strip() for s in anvio_db.meta['available_item_orders'].split(',')])
    item_orders_dict = anvio_db.db.get_table_as_dict(t.item_orders_table_name)

    for item_order in item_orders_dict:
        if item_orders_dict[item_order]['type'] == 'basic':
            try:
                item_orders_dict[item_order]['data'] = item_orders_dict[item_order]['data'].split(',')
            except:
                raise ConfigError("Something is wrong with the basic order `%s` in this %s database :(" % (item_order, utils.get_db_type(anvio_db_path)))

    return (available_item_orders, item_orders_dict)


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
        if 'mean_coverage:euclidean:ward' in item_orders_dict:
            default_item_order = 'mean_coverage:euclidean:ward'
        else:
            default_item_order = list(item_orders_dict.keys())[0]
            run.warning('`get_default_item_order_name` function is concerned, because nothing in the item_orders '
                        'dict matched to the desired order class "%s". So the order literally set to "%s" '
                        '(a class of "%s") randomly as the default order. Good luck :/' % (default_item_order_requested,
                                                                                     default_item_order,
                                                                                     default_item_order.split(':')[0]))
        return default_item_order
    elif len(matching_item_order_names) == 1:
        return matching_item_order_names[0]
    else:
        default_item_order = matching_item_order_names[0]
        run.warning('`get_default_item_order_name` function is concerned, because there were multiple entries '
                    'in the item_orders dict matched to the desired default order class "%s". So it set '
                    'the first of all %d matching item_orders, which happened to be the "%s", as the '
                    'default. We hope that will not screw up your mojo :/' % (default_item_order_requested,
                                                                              len(matching_item_order_names),
                                                                              default_item_order))
        return default_item_order
