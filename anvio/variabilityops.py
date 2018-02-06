# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of single nucleotide variation"""


from collections import Counter

import os
import sys
import copy
import random
import numpy as np

from scipy.stats import entropy

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.variability as variability
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = ['Alon Shaiber']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


pp = terminal.pretty_print
progress = terminal.Progress()
run = terminal.Run(width=62)


class VariabilitySuper(object):
    def __init__(self, args={}, p=progress, r=run):
        self.args = args

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        if args.engine not in variability_engines:
            raise ConfigError("You are doing something wrong :/ Focus '%s' does not correspond to an available engine." % args.engine)
        self.data = {}

        self.splits_of_interest = set([])
        self.samples_of_interest = set([])

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.bin_id = A('bin_id', null)
        self.collection_name = A('collection_name', null)
        self.splits_of_interest_path = A('splits_of_interest', null)
        self.min_departure_from_reference = A('min_departure_from_reference', float) or 0
        self.max_departure_from_reference = A('max_departure_from_reference', float) or 1
        self.min_departure_from_consensus = A('min_departure_from_consensus', float) or 0
        self.max_departure_from_consensus = A('max_departure_from_consensus', float) or 1
        self.min_occurrence = A('min_occurrence', int) or 1
        self.num_positions_from_each_split = A('num_positions_from_each_split', int) or 0
        self.min_scatter = A('min_scatter', int) or 0
        self.min_coverage_in_each_sample = A('min_coverage_in_each_sample', int) or 0
        self.profile_db_path = A('profile_db', null)
        self.contigs_db_path = A('contigs_db', null)
        self.quince_mode = A('quince_mode', bool)
        self.skip_comprehensive_variability_scores = A('skip_comprehensive_variability_scores', bool) or False
        self.output_file_path = A('output_file', null)
        self.samples_of_interest_path = A('samples_of_interest', null)
        self.genes_of_interest_path = A('genes_of_interest', null)
        self.include_contig_names_in_output = A('include_contig_names', null)
        self.include_split_names_in_output = A('include_split_names', null)
        self.gene_caller_id = A('gene_caller_id', null)

        self.substitution_scoring_matrices = None
        self.merged_split_coverage_values = None
        self.unique_pos_identifier = 0
        self.split_name_position_dict = {}
        self.unique_pos_id_to_entry_id = {}
        self.contig_sequences = None
        self.input_file_path = None

        if self.engine not in variability_engines:
            raise ConfigError("The superclass is inherited with an unknown engine. Anvi'o needs an adult :(")

        self.comprehensive_stats_headers = []
        self.comprehensive_variability_scores_computed = False

        # Initialize the contigs super
        filesnpaths.is_file_exists(self.contigs_db_path)
        dbops.ContigsSuperclass.__init__(self, self.args, r=self.run, p=self.progress)
        self.init_contig_sequences()


    def init_commons(self):
        self.progress.new('Init')

        self.progress.update('Checking the output file path ..')
        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)

        self.progress.update('Checking the samples of interest ..')
        if self.samples_of_interest_path:
            filesnpaths.is_file_tab_delimited(self.samples_of_interest_path, expected_number_of_fields=1)
            self.samples_of_interest = set([s.strip() for s in open(self.samples_of_interest_path).readlines()])
        else:
            self.samples_of_interest = set([])

        self.progress.update('Any genes of interest?')
        if self.genes_of_interest_path and self.gene_caller_id:
            self.progress.end()
            raise ConfigError("You can't provide a gene caller id from the command line, and a list of gene caller ids\
                               as a file at the same time, obviously.")

        if self.gene_caller_id:
            try:
                self.gene_caller_id = int(self.gene_caller_id)
            except:
                raise ConfigError("Anvi'o does not like your gene caller id '%s'..." % str(self.gene_caller_id))

            self.genes_of_interest = set([self.gene_caller_id])
        elif self.genes_of_interest_path:
            filesnpaths.is_file_tab_delimited(self.genes_of_interest_path, expected_number_of_fields=1)

            try:
                self.genes_of_interest = set([int(s.strip()) for s in open(self.genes_of_interest_path).readlines()])
            except ValueError:
                self.progress.end()
                raise ConfigError("Well. Anvi'o was working on your genes of interest .. and ... those gene IDs did not\
                                   look like anvi'o gene caller ids :/ Anvi'o is now sad.")
        else:
            self.genes_of_interest = set([])

        self.progress.update('Setting up genes of interest data ..')
        if self.genes_of_interest:
            # check for genes that do not appear in the contigs database
            bad_gene_caller_ids = [g for g in self.genes_of_interest if g not in self.gene_callers_id_to_split_name_dict]
            if bad_gene_caller_ids:
                self.progress.end()
                raise ConfigError("The gene caller id you provided is not known to this contigs database. You have only 2 lives\
                                   left. 2 more mistakes, and anvi'o will automatically uninstall itself. Yes, seriously :(")

            # if we know gene names, we know split names. set split names straight:
            self.splits_of_interest = list(set([self.gene_callers_id_to_split_name_dict[g] for g in self.genes_of_interest]))

        self.progress.update('Making sure you are not playing games ..')
        if self.engine not in ['NT', 'AA']:
            raise ConfigError("Anvi'o doesn't know what to do with a engine on '%s' yet :/" % self.engine)

        # set items of interest while you are at it.
        self.items = constants.amino_acids if self.engine == 'AA' else constants.nucleotides

        self.progress.update('Making sure our databases are here ..')
        if not self.profile_db_path:
            raise ConfigError('You need to provide a profile database.')

        if not self.contigs_db_path:
            raise ConfigError('You need to provide a contigs database.')

        self.progress.update('Making sure our databases are compatible ..')
        dbops.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        if self.min_coverage_in_each_sample and not self.quince_mode:
            self.progress.end()
            raise ConfigError("When you sepecify a coverage value through --min-coverage-in-each-sample, you must also\
                                use --quince-mode flag, since the former parameter needs to know the coverage values in all\
                                samples even if variation is reported for only one sample among otheres. This is the only way\
                                to figure out whether variation is not reported for other samples due to low or zero coverage,\
                                or there was no variation to report despite the high coverage. Anvi'o could turn --quince-mode\
                                flat automatically for you, but then it is much better if you have full control and understaning\
                                of what is going on.")

        if self.quince_mode:
            self.progress.update('Accessing auxiliary data file ...')
            auxiliary_data_file_path = os.path.join(os.path.dirname(self.profile_db_path), 'AUXILIARY-DATA.h5')
            if not os.path.exists(auxiliary_data_file_path):
                raise ConfigError("Anvi'o needs the auxiliary data file to run this program with '--quince-mode' flag.\
                                    However it wasn't found at '%s' :/" % auxiliary_data_file_path)
            self.merged_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(auxiliary_data_file_path, None, ignore_hash=True)

        self.progress.update('Attempting to get our splits of interest sorted ..')
        if self.collection_name:
            # the user wants to go with the collection id path. fine. we will get our split names from
            # the profile database.
            if not self.bin_id:
                self.progress.end()
                raise ConfigError('When you declare a collection id, you must also declare a bin name\
                                    (from which the split names of interest will be acquired)')
            if self.collection_name and self.splits_of_interest_path:
                self.progress.end()
                raise ConfigError("You declared a collection id and one or more bin names so anvi'o can find out\
                                    splits of interest, but you also have specified informaiton for split names?\
                                    This is confusing. You should choose one way or another :/")

            self.splits_of_interest = ccollections.GetSplitNamesInBins(self.args).get_split_names_only()
        else:
            # OK. no collection id. we will go oldschool. we whope to find what we are looking for in
            # self.splits_of_interst_path  at this point (which may have been filled through the command
            # line client), or in self.splits_of_interest (which may have been filled in by another program)
            if not self.splits_of_interest:
                if not self.splits_of_interest_path:
                    self.progress.end()
                    raise ConfigError('You did not declare a source for split names. You either should give me\
                                        a file with split names you are interested in, or a collection id and\
                                        bin name so I can learn split names from the profile database.')
                filesnpaths.is_file_exists(self.splits_of_interest_path)
                self.splits_of_interest = set([c.strip().replace('\r', '') for c in open(self.splits_of_interest_path).readlines()])

        self.input_file_path = '/' + '/'.join(os.path.abspath(self.profile_db_path).split('/')[:-1])

        self.progress.update('Reading the data ...')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        self.sample_ids = sorted(list(profile_db.samples)) # we set this now, but we will overwrite it with args.samples_of_interest if necessary

        if not profile_db.meta['SNVs_profiled']:
            self.progress.end()
            raise ConfigError("Well well well. It seems SNVs were not characterized for this profile database.\
                                Sorry, there is nothing to report here!")

        # populate substitution scoring matrices
        self.progress.end()
        import anvio.data.SSMs as SSMs
        self.substitution_scoring_matrices = SSMs.get(self.engine)
        self.progress.new('Init')

        ##################### LOAD ENGINE-SPECIFIC DATA #####################
        # data is one of them, since they will be read from different tables.
        # another one is the substitution scoring matrices.
        if self.engine == 'NT':
            self.data = profile_db.db.get_table_as_dict(t.variable_nts_table_name)
            self.progress.end()

        elif self.engine == 'AA':
            if not profile_db.meta['AA_frequencies_profiled']:
                raise ConfigError("It seems AA frequencies were not characterized for this profile database.\
                                    There is nothing to report here for AAs!")

            self.data = profile_db.db.get_table_as_dict(t.variable_aas_table_name)
            self.progress.end()

            # append split_name information
            for e in list(self.data.values()):
                e['split_name'] = self.gene_callers_id_to_split_name_dict[e['corresponding_gene_call']]

        # we're done here. bye.
        profile_db.disconnect()


    def filter(self, filter_name, test_func):
        self.progress.new('Filtering based on "%s"' % filter_name)

        entry_ids_to_remove, counter = set([]), 0

        for entry_id in self.data:
            if counter % 1000 == 0:
                self.progress.update('identifying entries to remove :: %s' % pp(counter))

            counter += 1

            if test_func(self.data[entry_id]):
                entry_ids_to_remove.add(entry_id)
                continue

        self.progress.end()

        self.remove_entries_from_data(entry_ids_to_remove, reason=filter_name)


    def check_if_data_is_empty(self):
        if not len(self.data):
            self.progress.end()
            self.run.info_single('Nothing left in the variability data to work with. Quitting :/', 'red', 1, 1)
            sys.exit()


    def apply_advanced_filters(self):
        if self.min_departure_from_consensus:
            self.run.info('Min departure from consensus', self.min_departure_from_consensus)
            self.filter('max departure from consensus', lambda x: x['departure_from_consensus'] < self.min_departure_from_consensus)

        if self.max_departure_from_consensus < 1:
            self.run.info('Max departure from consensus', self.max_departure_from_consensus)
            self.filter('max departure from consensus', lambda x: x['departure_from_consensus'] > self.max_departure_from_consensus)

        self.check_if_data_is_empty()


    def apply_preliminary_filters(self):
        self.run.info('Variability data', '%s entries in %s splits across %s samples'\
                % (pp(len(self.data)), pp(len(self.splits_basic_info)), pp(len(self.sample_ids))))

        self.run.info('Samples in the profile db', ', '.join(sorted(self.sample_ids)))
        if self.samples_of_interest:
            samples_missing_from_db = [sample for sample in self.samples_of_interest if sample not in self.sample_ids]

            if len(samples_missing_from_db):
                raise ConfigError('One or more samples you are interested in seem to be missing from\
                                    the profile database: %s' % ', '.join(samples_missing_from_db))

            self.run.info('Samples of interest', ', '.join(sorted(list(self.samples_of_interest))))
            self.sample_ids = sorted(list(self.samples_of_interest))
            self.filter('samples of interest', lambda x: x['sample_id'] not in self.samples_of_interest)

        if self.genes_of_interest:
            self.run.info('Num genes of interest', pp(len(self.genes_of_interest)))
            self.filter('genes of interest', lambda x: x['corresponding_gene_call'] not in self.genes_of_interest)

        if self.splits_of_interest:
            self.run.info('Num splits of interest', pp(len(self.splits_of_interest)))
            self.filter('splits of interest', lambda x: x['split_name'] not in self.splits_of_interest)

        # let's report the number of positions reported in each sample before filtering any futrher:
        num_positions_each_sample = Counter([v['sample_id'] for v in list(self.data.values())])
        self.run.info('Total number of variable positions in samples', '; '.join(['%s: %s' % (s, num_positions_each_sample[s]) for s in sorted(self.sample_ids)]))

        if self.min_departure_from_reference:
            self.run.info('Min departure from reference', self.min_departure_from_reference)
            self.filter('min departure from reference', lambda x: x['departure_from_reference'] < self.min_departure_from_reference)

        if self.max_departure_from_reference < 1:
            self.run.info('Max departure from reference', self.max_departure_from_reference)
            self.filter('max departure from reference', lambda x: x['departure_from_reference'] > self.max_departure_from_reference)

        # this is the first time we are going through each data entry. before things get thorny down below,
        # we will use this opportunity to add gene_length as well as unique_pos_identifier_str for each entry.
        # a more granular design may be necessary for exposure later on.
        for entry_id in self.data:
            v = self.data[entry_id]
            if self.engine == 'NT':
                v['unique_pos_identifier_str'] = '_'.join([v['split_name'], str(v['pos'])])
            if self.engine == 'AA':
                v['unique_pos_identifier_str'] = '_'.join([v['split_name'], str(v['corresponding_gene_call']), str(v['codon_order_in_gene'])])

            v['gene_length'] = self.get_gene_length(v['corresponding_gene_call'])

        if self.min_occurrence == 1:
            return

        if self.min_occurrence > 1:
            self.run.info('Min occurrence requested', self.min_occurrence)

        self.progress.new('Processing positions further')

        self.progress.update('counting occurrences of each position across samples ...')
        unique_pos_identifier_str_occurrences = {}
        for entry_id in self.data:
            v = self.data[entry_id]
            if v['unique_pos_identifier_str'] in unique_pos_identifier_str_occurrences:
                unique_pos_identifier_str_occurrences[v['unique_pos_identifier_str']] += 1
            else:
                unique_pos_identifier_str_occurrences[v['unique_pos_identifier_str']] = 1

        self.progress.update('identifying entries that occurr in less than %d samples ...' % (self.min_occurrence))
        entry_ids_to_remove = set([])
        for entry_id in self.data:
            v = self.data[entry_id]
            if not unique_pos_identifier_str_occurrences[v['unique_pos_identifier_str']] >= self.min_occurrence:
                entry_ids_to_remove.add(entry_id)

        self.progress.end()

        self.remove_entries_from_data(entry_ids_to_remove, reason="preliminary filters")


    def set_unique_pos_identification_numbers(self):
        self.progress.new('Further processing')
        self.progress.update('re-setting unique identifiers to track split/position pairs across samples')

        for entry_id in self.data:
            v = self.data[entry_id]
            v['unique_pos_identifier'] = self.get_unique_pos_identification_number(v['unique_pos_identifier_str'])
            v['contig_name'] = self.splits_basic_info[v['split_name']]['parent']

        self.progress.end()


    def get_unique_pos_identification_number(self, unique_pos_identifier_str):
        if unique_pos_identifier_str in self.split_name_position_dict:
            return self.split_name_position_dict[unique_pos_identifier_str]
        else:
            self.split_name_position_dict[unique_pos_identifier_str] = self.unique_pos_identifier
            self.unique_pos_identifier += 1
            return self.split_name_position_dict[unique_pos_identifier_str]


    def gen_unique_pos_identifier_to_entry_id_dict(self):
        self.progress.new('Generating the `unique pos` -> `entry id` dict')
        self.progress.update('...')

        self.unique_pos_id_to_entry_id = {}

        for entry_id in self.data:
            v = self.data[entry_id]
            u = v['unique_pos_identifier_str']
            if u in self.unique_pos_id_to_entry_id:
                self.unique_pos_id_to_entry_id[u].add(entry_id)
            else:
                self.unique_pos_id_to_entry_id[u] = set([entry_id])

        self.progress.end()


    def insert_additional_fields(self, keys=[]):
        if not len(keys):
            keys = list(self.data.keys())

        for key in keys:
            e = self.data[key]
            item_frequencies = utils.get_variabile_item_frequencies(e, self.engine)
            e['n2n1ratio'], e['consensus'], e['departure_from_consensus'] = utils.get_consensus_and_departure_data(item_frequencies)

            # this is where we will make use of the substitution scoring matrices framework.
            if self.engine == 'NT':
                competing_items = list(e['competing_nts'])
            elif self.engine == 'AA':
                competing_items = variability.get_competing_items(e['reference'], item_frequencies)

                if not competing_items:
                    # this is a position that did have variation with SNVs, but all of them turned out to be synonymous.
                    # we will just list this one as itself.
                    competing_items = [item_frequencies[0][0]] * 2

                e['competing_aas'] = ''.join(competing_items)
            else:
                raise ConfigError("You cray :( Ain't nobody got an engine for %s." % self.engine)

            for m in self.substitution_scoring_matrices:
                try:
                    e[m] = self.substitution_scoring_matrices[m][competing_items[0]][competing_items[1]]
                except KeyError:
                    e[m] = None


    def filter_based_on_scattering_factor(self):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from eachother."""

        if self.min_scatter == 0:
            return

        num_samples = len(self.sample_ids)
        if self.min_scatter > num_samples / 2:
            raise ConfigError('Minimum scatter (%d) can not be more than half of the number of samples\
                                (%d) :/' % (self.min_scatter, num_samples))

        self.run.info('Min scatter', self.min_scatter)

        # we need the unique pos_id to entry id dict filled for this function:
        self.gen_unique_pos_identifier_to_entry_id_dict()

        self.progress.new('Examining scatter')
        self.progress.update('...')

        entry_ids_to_remove = set([])

        for unique_pos_id in self.unique_pos_id_to_entry_id:
            entry_ids = self.unique_pos_id_to_entry_id[unique_pos_id]

            # find how many samples it occurs:
            num_occurrence = len(entry_ids)

            # if the number of samples it occurs more
            scatter = num_occurrence if num_occurrence < num_samples - num_occurrence else num_samples - num_occurrence
            if scatter < self.min_scatter:
                entry_ids_to_remove.update(entry_ids)

        self.progress.end()

        self.remove_entries_from_data(entry_ids_to_remove, reason="minimum scatter")


    def remove_entries_from_data(self, entry_ids_to_remove=set([]), reason="unknown reason"):
        """Safely remove entries from self.data"""

        if not entry_ids_to_remove:
            return

        self.progress.new('Data removal for "%s"' % reason)
        self.progress.update('...')

        num_entries_before = len(self.data)

        # when data is removed, one of the dictionaries to also update is `self.unique_pos_id_to_entry_id`,
        # however, not at all stages of the process this dictionary is present. hence, we need to determine
        # whether we need to work with it early on:
        unique_pos_id_to_entry_id_needs_updating = 'unique_pos_identifier_str' in self.data[next(iter(entry_ids_to_remove))]

        self.progress.update('removing %s entries from data...' % pp(len(entry_ids_to_remove)))
        unique_pos_ids_to_remove = set([])
        for entry_id in entry_ids_to_remove:
            if unique_pos_id_to_entry_id_needs_updating:
                unique_pos_ids_to_remove.add(self.data[entry_id]['unique_pos_identifier_str'])
            self.data.pop(entry_id)

        if self.unique_pos_id_to_entry_id:
            self.progress.update('removing %s unique positions...' % pp(len(unique_pos_ids_to_remove)))
            for unique_pos_id in unique_pos_ids_to_remove:
                self.unique_pos_id_to_entry_id.pop(unique_pos_id)

        num_entries_after = len(self.data)

        self.progress.end()

        self.run.info('Remaining entries after "%s"' % (reason),
                      '%s (%s was removed)' % (pp(num_entries_after),
                                               pp(num_entries_before - num_entries_after)),
                      mc='green')

        self.check_if_data_is_empty()


    def filter_based_on_minimum_coverage_in_each_sample(self):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from eachother."""

        if self.min_coverage_in_each_sample < 1:
            return

        self.run.info('Min coverage in all samples', '%dX' % self.min_coverage_in_each_sample)

        # we need to make sure we have an up-to-date dictionary for unque position to entry id conversion:
        self.gen_unique_pos_identifier_to_entry_id_dict()

        self.progress.new('Examining coverage of each variable position in each sample')
        self.progress.update('...')

        entry_ids_to_remove = set([])

        for unique_pos_id in self.unique_pos_id_to_entry_id:
            entry_ids = self.unique_pos_id_to_entry_id[unique_pos_id]

            min_coverage_in_a_sample = min([self.data[entry_id]['coverage'] for entry_id in entry_ids])

            if min_coverage_in_a_sample < self.min_coverage_in_each_sample:
                entry_ids_to_remove.update(entry_ids)

        self.progress.end()

        self.remove_entries_from_data(entry_ids_to_remove, "min cov in all samples")


    def filter_based_on_num_positions_from_each_split(self):
        if self.num_positions_from_each_split > 0:
            self.run.info('Num positions to keep from each split', self.num_positions_from_each_split)
        else:
            self.run.info('Num positions to keep from each split', '(all positions)')
            return

        self.progress.new('Filtering based on -n')

        self.progress.update('Generating splits and positions tuples ...')
        splits_and_positions = set([(v['split_name'], v['unique_pos_identifier']) for v in list(self.data.values())])
        unique_positions_to_remove = set([])

        self.progress.update('Generating positions in splits dictionary ...')
        positions_in_splits = {}
        for split_name, position in splits_and_positions:
            if split_name not in positions_in_splits:
                positions_in_splits[split_name] = set([])

            positions_in_splits[split_name].add(position)

        self.progress.update('Randomly subsampling from splits with num positions > %d' % self.num_positions_from_each_split)
        for split_name in positions_in_splits:
            if self.num_positions_from_each_split and len(positions_in_splits[split_name]) > self.num_positions_from_each_split:
                positions_to_keep = set(random.sample(positions_in_splits[split_name], self.num_positions_from_each_split))
                for pos in (positions_in_splits[split_name] - positions_to_keep):
                    unique_positions_to_remove.add(pos)

        self.progress.update('Identifying entry ids to remove ...')
        entry_ids_to_remove = set([entry_id for entry_id in self.data if self.data[entry_id]['unique_pos_identifier'] in unique_positions_to_remove])

        self.progress.end()

        self.remove_entries_from_data(entry_ids_to_remove, reason="max num positions from each split")


    def compute_comprehensive_variability_scores(self):
        """
            Comprehensive stats - we compute scores that take into consideration the entire vector of variability
            and not only the two most competing items (and thus it is comprehensive).
            Currently the scores that are included are: site-entropy, site-Kullback-Leibler divergence (both a 
            raw score and a normalized score (see below)), and weighted substitution scores (i.e. BLOSUM).
            
            site-entropy - the entropy of the items in a single site (in a sample).

            Kullback-Leibler divergence raw - the Kullback-Leibler divergence of the frequencies in a sample
            compared to the raw frequencies of the sum of occurances in the same site accross samples.

            Kullback-Leibler divergence normalized - the Kullback-Leibler divergence of the frequencies in a sample
            compared to the frequencies of the sum of normalized occurances in the same site accross samples. Where
            the normalization is such that the occernce of items is normalized to sum to one in each sample. This method
            eliminates the effect of coverage on the score. The disadvantage of this method is that if there is a sample with
            low coverage then any noise (like a single sequencing error) could have a major effect. It is recommended 
            to use this score in combination with the --min-coverage-in-each-sample.
            
            Weighted substitution scores - the weights per substitution score is weighted by the joint frequency of the items
            i.e. sum(S_{i,j}*pi*pj) where i does not equal j (i.e. the substitution of an item with itself is not considered)
        """

        if self.skip_comprehensive_variability_scores:
            self.run.warning("Anvi'o will skip comprehensive variability score computations.")
            return

        # FIXME: In the future we should separate this function into two:
        #   1. scores that require quince mode - these are scores that are computed according 
        #       to the occurence accross samples, and hence they require quince mode
        #       for example: Kullbak-Leibler divergence
        #   2. scores that dont require quince mode - these are calculated per sample (e.g. entropy, BLOSUM)
        # For now we just require quince mode for all scores
        if not self.quince_mode:
            self.run.warning("Comprehensive variability score computations can only be done with `--quince-mode`")
            return

        unique_positions_and_frequencies_dict = self.get_unique_positions_and_frequencies_dict()

        self.progress.new('Comprehensive stats')
        self.progress.update('...')

        comprehensive_stats = {}
        self.comprehensive_stats_headers = [m + '_weighted' for m in self.substitution_scoring_matrices] + \
                                           ['entropy', 'kullback_leibler_divergence_raw', 'kullback_leibler_divergence_normalized']

        # converting the substitution matrices from dict of dicts to numpy matrix (to allow vector operations later on)
        # we keep track of the indices of the items for the following reasons:
        #    1. in case the substitution matrix has only a sub-set of the items
        #    2. to make sure the order of items is compatible between this matrix and the items frequency vectors (and thus allow vector operations - see below)
        self.progress.update('Initializing numpy formatted substitution matrices...')
        item_indices_for_substitution_scoring_matrices = {}
        numpy_matrices_for_substitution_scoring_matrices = {}
        for m in self.substitution_scoring_matrices:
            item_indices_for_substitution_scoring_matrices[m] = [self.items.index(item) for item in self.substitution_scoring_matrices[m].keys()]

            items = list(self.substitution_scoring_matrices[m].keys())
            num_items = len(items)
            numpy_matrices_for_substitution_scoring_matrices[m] = np.zeros([num_items, num_items])
            for i in range(num_items):
                for j in range(num_items):
                    if i == j:
                        # We set the substitution score of an item with itself to zero. This way we only consider substitutions to other items
                        # (and dont consider substitution of an item with itself)
                        numpy_matrices_for_substitution_scoring_matrices[m][i][j] = 0
                    else:
                        numpy_matrices_for_substitution_scoring_matrices[m][i][j] = self.substitution_scoring_matrices[m][items[i]][items[j]]

        self.progress.update('Running comprehensive stats across all unique positions X samples')
        unique_positions = list(unique_positions_and_frequencies_dict.keys())
        num_unique_positions = len(unique_positions)
        for unique_pos_index in range(num_unique_positions):
            unique_pos = unique_positions[unique_pos_index]

            if unique_pos_index % 100 == 0:
                self.progress.update('Running comprehensive stats: %s of %s ...' % (pp(unique_pos_index + 1), pp(num_unique_positions)))

            comprehensive_stats = {}

            # first create a vector of frequencies for all samples
            list_of_sample_frequencies = []
            for sample_id in self.sample_ids:
                list_of_sample_frequencies.append([unique_positions_and_frequencies_dict[unique_pos][sample_id][f] for f in self.items])
                comprehensive_stats[sample_id] = {}

            # create a numpy array for all samples, and sum sample frequencies
            list_of_sample_frequencies = np.array(list_of_sample_frequencies)
            sum_sample_frequencies = sum(list_of_sample_frequencies)

            # normalizing frequencies (so they add to 1) (notice that by deviding by a column vector we divide each row of the matrix by the sum that
            # corresponds to that row)
            list_of_sample_frequencies_normalized = np.divide(np.array(list_of_sample_frequencies), np.sum(list_of_sample_frequencies, axis=1)[:, np.newaxis])
            sum_sample_frequencies_normalized = sum(list_of_sample_frequencies_normalized)

            for i in range(0, len(self.sample_ids)):
                sample_id = self.sample_ids[i]

                # compute entropy
                comprehensive_stats[sample_id]['entropy'] = entropy(list_of_sample_frequencies[i])

                # compute Kullback-Leibler divergence for raw counts
                kullback_leibler_divergence_raw = entropy(list_of_sample_frequencies[i], sum_sample_frequencies)
                # compute Kullback-Leibler divergence for normalized counts (normalized to sum to 1 in each sample)
                kullback_leibler_divergence_normalized = entropy(list_of_sample_frequencies_normalized[i], sum_sample_frequencies_normalized)

                comprehensive_stats[sample_id]['kullback_leibler_divergence_raw'] = kullback_leibler_divergence_raw
                comprehensive_stats[sample_id]['kullback_leibler_divergence_normalized'] = kullback_leibler_divergence_normalized

                # computing weighted substitution score
                # we compute a weighted score for each existing substitution matrix
                if not comprehensive_stats[sample_id]['entropy']:
                    # if the entropy is zero it means there are no substitutions in this sample (i.e. this row is here due to quince mode)
                    # so no need to calculate substitution score
                    for m in self.substitution_scoring_matrices:
                        comprehensive_stats[sample_id][m + '_weighted'] = None
                else:
                    for m in self.substitution_scoring_matrices:
                        # here we subsample AND reorder our sample frequencies based on the items that
                        # appear in the substitution matrix `m`. see the code above with the for loop
                        # to remember how they are set.
                        S = list_of_sample_frequencies[i][item_indices_for_substitution_scoring_matrices[m]]
                        S = S / sum(S)

                        if np.count_nonzero(S) > 1:
                            # because the substitution matrix might hold only a subset of the items,
                            # we could have positions in which the entropy is greater than zero, but 
                            # only due to items that are not in the substitution matrix (for example stop codon is not in BLOSUM)
                            # hence to avoid devision by zero we make sure there is more than one nonzero frequency

                            # a normalization score is needed since we dont consider a substitution of an item with itself.
                            # hence, the sum of frequencies wouldn't sum to 1, and so to make sure they sum to 1
                            # we multiply by this normalization factor. We want these to sum to 1 because otherwise these
                            # are not valid weights.
                            normalization_factor = 1 / (1 - np.sum(np.square(S)))
                            # element-wise multiplication of the substitution matrix with the outer-product of the frequency vector with itself
                            # this means that each substitution score between two items is weighted by the product of the weight of the items
                            A = normalization_factor * np.sum(np.multiply(numpy_matrices_for_substitution_scoring_matrices[m], np.outer(S, S)))
                            comprehensive_stats[sample_id][m + '_weighted'] = A
                        else:
                            comprehensive_stats[sample_id][m + '_weighted'] = None

            # update self.data with comprehensive stats
            for entry_id in self.unique_pos_id_to_entry_id[unique_pos]:
                e = self.data[entry_id]
                s = e['sample_id']

                for h in self.comprehensive_stats_headers:
                    e[h] = comprehensive_stats[s][h]

        self.progress.end()

        self.comprehensive_variability_scores_computed = True


    def get_unique_positions_and_frequencies_dict(self):
        """From the self.data object, creates a dict that contains item frequencies for
           each sample for each unique position identifier."""

        unique_positions_and_frequencies_dict = {}

        template = dict.fromkeys(self.items, 0)

        if not self.unique_pos_id_to_entry_id:
            self.gen_unique_pos_identifier_to_entry_id_dict()

        self.progress.new('The unique positions and frequencies dict')
        self.progress.update('generating ..')

        for entry_ids in self.unique_pos_id_to_entry_id.values():
            unique_pos_identifier = self.data[list(entry_ids)[0]]['unique_pos_identifier_str']
            unique_positions_and_frequencies_dict[unique_pos_identifier] = {}

            for entry_id in entry_ids:
                v = self.data[entry_id]
                unique_positions_and_frequencies_dict[unique_pos_identifier][v['sample_id']] = copy.deepcopy(template)

                for item in self.items:
                    if v[item]:
                        unique_positions_and_frequencies_dict[unique_pos_identifier][v['sample_id']][item] = v[item]

        self.progress.end()

        return unique_positions_and_frequencies_dict


    def process(self):
        self.init_commons()

        self.apply_preliminary_filters()

        self.set_unique_pos_identification_numbers() # which allows us to track every unique position across samples

        self.filter_based_on_scattering_factor()

        self.filter_based_on_num_positions_from_each_split()

        self.insert_additional_fields()

        self.apply_advanced_filters()

        if self.quince_mode: # will be very costly...
            self.recover_base_frequencies_for_all_samples()

        self.filter_based_on_minimum_coverage_in_each_sample()

        self.compute_comprehensive_variability_scores()


    def get_gene_length(self, gene_callers_id):
            if gene_callers_id in self.gene_lengths:
                return self.gene_lengths[gene_callers_id]
            else:
                return -1


    def get_unique_pos_identifier_to_corresponding_gene_id(self):
        self.progress.update('populating a dict to track corresponding gene ids for each unique position')
        unique_pos_identifier_to_corresponding_gene_id = {}
        for entry in list(self.data.values()):
            unique_pos_identifier_to_corresponding_gene_id[entry['unique_pos_identifier']] = entry['corresponding_gene_call']

        return unique_pos_identifier_to_corresponding_gene_id


    def get_unique_pos_identifier_to_codon_order_in_gene(self):
        self.progress.update('populating a dict to track codon order in genes for each unique position')
        unique_pos_identifier_to_codon_order_in_gene = {}
        for entry in list(self.data.values()):
            unique_pos_identifier_to_codon_order_in_gene[entry['unique_pos_identifier']] = entry['codon_order_in_gene']

        return unique_pos_identifier_to_codon_order_in_gene


    def report(self):
        self.progress.new('Reporting variability data')

        if self.engine == 'NT':
            new_structure = [t.variable_nts_table_structure[0]] + \
                             ['unique_pos_identifier', 'gene_length'] + \
                             [x for x in t.variable_nts_table_structure[1:] if x != 'split_name'] + \
                             list(self.substitution_scoring_matrices.keys()) + \
                             ['consensus', 'departure_from_consensus', 'n2n1ratio'] + \
                             self.comprehensive_stats_headers
        elif self.engine == 'AA':
            new_structure = [t.variable_nts_table_structure[0]] + \
                             ['unique_pos_identifier', 'gene_length'] + \
                             [x for x in t.variable_aas_table_structure[1:] if x != 'split_name'] + \
                             list(self.substitution_scoring_matrices.keys()) + \
                             ['competing_aas', 'consensus', 'departure_from_consensus', 'n2n1ratio'] + \
                             self.comprehensive_stats_headers


        if self.include_contig_names_in_output:
            new_structure.append('contig_name')

        if self.include_split_names_in_output:
            new_structure.append('split_name')

        self.progress.update('exporting variable positions table as a TAB-delimited file ...')

        utils.store_dict_as_TAB_delimited_file(self.data, self.args.output_file, new_structure)
        self.progress.end()

        self.run.info('Num entries reported', pp(len(self.data)))
        self.run.info('Output File', self.output_file_path)
        self.run.info('Num %s positions reported' % self.engine, pp(len(set([e['unique_pos_identifier'] for e in list(self.data.values())]))))


class VariableNtPositionsEngine(dbops.ContigsSuperclass, VariabilitySuper):
    """This is the main class to make sense and report variability for a given set of splits,
       or a bin in a collection, across multiple or all samples. The user can scrutinize the
       nature of the variable positions to be reported dramatically given the ecology and/or
       other biologically-relevant considerations, or purely technical limitations such as
       minimum coverage of a given nucleotide position or the ratio of the competing nts at a
       given position. The default entry to this class is the `anvi-gen-variability-profile`
       program."""

    def __init__(self, args={}, p=progress, r=run):
        self.run = r
        self.progress = p

        self.engine = 'NT'

        # Init Meta
        VariabilitySuper.__init__(self, args=args, r=self.run, p=self.progress)


    def recover_base_frequencies_for_all_samples(self):
        """this function populates variable_nts_table dict with entries from samples that have no
           variation at nucleotide positions reported in the table"""

        self.progress.new('Recovering NT data')

        samples_wanted = self.samples_of_interest if self.samples_of_interest else self.sample_ids
        splits_wanted = self.splits_of_interest if self.splits_of_interest else set(self.splits_basic_info.keys())
        next_available_entry_id = max(self.data.keys()) + 1

        unique_pos_identifier_to_corresponding_gene_id = self.get_unique_pos_identifier_to_corresponding_gene_id()

        unique_pos_identifier_to_codon_order_in_gene = self.get_unique_pos_identifier_to_codon_order_in_gene()

        self.progress.update('creating a dict to track missing base frequencies for each sample / split / pos')
        split_pos_to_unique_pos_identifier = {}
        splits_to_consider_dict = {}
        for split_name in splits_wanted:
            splits_to_consider_dict[split_name] = {}
            split_pos_to_unique_pos_identifier[split_name] = {}

        self.progress.update('populating the dict to track missing base frequencies for each sample / split / pos')
        for entry_id in self.data:
            v = self.data[entry_id]
            p = v['pos']
            d = splits_to_consider_dict[v['split_name']]
            u = split_pos_to_unique_pos_identifier[v['split_name']]

            if p in d:
                d[p].remove(v['sample_id'])
            else:
                d[p] = copy.deepcopy(samples_wanted)
                d[p].remove(v['sample_id'])

            if p not in u:
                u[p] = v['unique_pos_identifier']

        split_names_to_consider = list(splits_to_consider_dict.keys())
        num_splits = len(split_names_to_consider)
        for split_index in range(num_splits):
            split = split_names_to_consider[split_index]
            self.progress.update('Accessing split covs, updating variable pos dict (%s of %s)' % (pp(split_index + 1), pp(num_splits)))

            split_coverage_across_samples = self.merged_split_coverage_values.get(split)

            split_info = self.splits_basic_info[split]
            contig_name_name = split_info['parent']

            for pos in splits_to_consider_dict[split]:
                unique_pos_identifier = split_pos_to_unique_pos_identifier[split][pos]
                contig_name_seq = self.contig_sequences[contig_name_name]['sequence']
                pos_in_contig = split_info['start'] + pos
                base_at_pos = contig_name_seq[pos_in_contig]
                corresponding_gene_call = unique_pos_identifier_to_corresponding_gene_id[unique_pos_identifier]
                gene_length = self.get_gene_length(corresponding_gene_call)
                codon_order_in_gene = unique_pos_identifier_to_codon_order_in_gene[unique_pos_identifier]

                in_partial_gene_call, in_complete_gene_call, base_pos_in_codon = self.get_nt_position_info(contig_name_name, pos_in_contig)

                for sample in splits_to_consider_dict[split][pos]:
                    self.data[next_available_entry_id] = {'contig_name': contig_name_name,
                                                          'departure_from_reference': 0,
                                                          'reference': base_at_pos,
                                                          'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0,
                                                          'pos': pos,
                                                          'pos_in_contig': pos_in_contig,
                                                          'in_partial_gene_call': in_partial_gene_call,
                                                          'in_complete_gene_call': in_complete_gene_call,
                                                          'base_pos_in_codon': base_pos_in_codon,
                                                          'coverage': split_coverage_across_samples[sample][pos],
                                                          'sample_id': sample,
                                                          'cov_outlier_in_split': 0,
                                                          'cov_outlier_in_contig': 0,
                                                          'competing_nts': base_at_pos + base_at_pos,
                                                          'unique_pos_identifier': unique_pos_identifier,
                                                          'unique_pos_identifier_str': '%s_%d' % (split, pos),
                                                          'corresponding_gene_call': corresponding_gene_call,
                                                          'gene_length': gene_length,
                                                          'codon_order_in_gene': codon_order_in_gene,
                                                          'split_name': split}
                    self.data[next_available_entry_id][base_at_pos] = split_coverage_across_samples[sample][pos]
                    self.insert_additional_fields([next_available_entry_id])
                    next_available_entry_id += 1

        self.progress.end()


class VariableAAPositionsEngine(dbops.ContigsSuperclass, VariabilitySuper):
    def __init__(self, args={}, p=progress, r=run):
        self.run = r
        self.progress = p

        self.engine = 'AA'

        # Init Meta
        VariabilitySuper.__init__(self, args=args, r=self.run, p=self.progress)


    def recover_base_frequencies_for_all_samples(self):
        self.progress.new('Recovering AA data')

        samples_wanted = self.samples_of_interest if self.samples_of_interest else self.sample_ids
        splits_wanted = self.splits_of_interest if self.splits_of_interest else set(self.splits_basic_info.keys())
        next_available_entry_id = max(self.data.keys()) + 1

        unique_pos_identifier_str_to_consenus_codon = {}
        unique_pos_identifier_str_to_unique_pos_identifier = {}
        for e in list(self.data.values()):
            upi = e['unique_pos_identifier_str']
            unique_pos_identifier_str_to_consenus_codon[upi] = e['reference']
            unique_pos_identifier_str_to_unique_pos_identifier[upi] = e['unique_pos_identifier']

        self.progress.update('creating a dict to track missing AA frequencies for each sample / split / pos')

        splits_to_consider_dict = {}
        for split_name in splits_wanted:
            splits_to_consider_dict[split_name] = {}

        self.progress.update('populating the dict to track missing AA frequencies for each sample / split / pos')
        for entry_id in self.data:
            v = self.data[entry_id]
            gene_codon_key = '%d_%d' % (v['corresponding_gene_call'], v['codon_order_in_gene'])
            d = splits_to_consider_dict[v['split_name']]

            if gene_codon_key in d:
                d[gene_codon_key].remove(v['sample_id'])
            else:
                d[gene_codon_key] = copy.deepcopy(samples_wanted)
                d[gene_codon_key].remove(v['sample_id'])

        split_names_to_consider = list(splits_to_consider_dict.keys())
        num_splits = len(split_names_to_consider)
        for split_index in range(num_splits):
            split_name = split_names_to_consider[split_index]
            self.progress.update('Accessing split covs, updating variable pos dict (%s of %s)' % (pp(split_index + 1), pp(num_splits)))

            split_coverage_across_samples = self.merged_split_coverage_values.get(split_name)

            split_info = self.splits_basic_info[split_name]
            contig_name = split_info['parent']

            for gene_codon_key in splits_to_consider_dict[split_name]:
                corresponding_gene_call, codon_order_in_gene = [int(k) for k in gene_codon_key.split('_')]
                gene_length = self.get_gene_length(corresponding_gene_call)

                for sample_name in splits_to_consider_dict[split_name][gene_codon_key]:
                    unique_pos_identifier_str = '_'.join([split_name, str(corresponding_gene_call), str(codon_order_in_gene)])
                    reference_codon = unique_pos_identifier_str_to_consenus_codon[unique_pos_identifier_str]

                    self.data[next_available_entry_id] = {'unique_pos_identifier_str': unique_pos_identifier_str,
                                                          'unique_pos_identifier': unique_pos_identifier_str_to_unique_pos_identifier[unique_pos_identifier_str],
                                                          'sample_id': sample_name,
                                                          'split_name': split_name,
                                                          'contig_name': contig_name,
                                                          'corresponding_gene_call': corresponding_gene_call,
                                                          'gene_length': gene_length,
                                                          'codon_order_in_gene': codon_order_in_gene,
                                                          'departure_from_reference': 0,
                                                          'coverage': None,
                                                          'reference': reference_codon}

                    # DEALING WITH COVERAGE ##################################################################
                    # some very cool but expensive shit is going on here, let me break it down for poor souls of the future.
                    # what we want to do is to learn the coverage of this codon in the sample. all we have is the corresponding
                    # gene call id, and the order of this codon in the gene. so here how it goes:
                    #
                    # learn the gene call
                    gene_call = self.genes_in_contigs_dict[corresponding_gene_call]

                    # the following dict converts codon orders into nt positions in contig for a geven gene call
                    codon_order_to_nt_positions_in_contig = utils.get_codon_order_to_nt_positions_dict(gene_call)

                    # so the nucleotide positions for this codon in the contig is the following:
                    nt_positions_for_codon_in_contig = codon_order_to_nt_positions_in_contig[codon_order_in_gene]

                    # but we need to convert those positions to the context of this split. so here is the start pos:
                    split_start = self.splits_basic_info[split_name]['start']

                    # here we map nt positions from the contig context to split context using the start position
                    nt_positions_for_codon_in_split = [p - split_start for p in nt_positions_for_codon_in_contig]

                    # we acquire coverages that match to these positions
                    coverages = split_coverage_across_samples[sample_name][nt_positions_for_codon_in_split]
                    coverage = int(round(sum(coverages) / 3))

                    # and finally update the data table
                    self.data[next_available_entry_id]['coverage'] = coverage

                    # DEALING WITH AAs ##################################################################
                    # here we need to put all the codons into the data table for this sample
                    for codon in set(constants.codon_to_AA.values()):
                        self.data[next_available_entry_id][codon] = 0

                    # and finally update the frequency of the reference codon with the coverage (WHICH IS VERY BAD,
                    # WE HAVE NO CLUE WHAT IS THE ACTUAL COVERAGE OF TRIPLICATE LINKMERS):
                    self.data[next_available_entry_id][reference_codon] = coverage

                    # insert additional fields for this newly added data point
                    self.insert_additional_fields([next_available_entry_id])

                    next_available_entry_id += 1

        self.progress.end()


class ConsensusSequences(VariableNtPositionsEngine, VariableAAPositionsEngine):
    def __init__(self, args={}, p=progress, r=run):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.engine = A('engine', null)
        self.tab_delimited_output = A('tab_delimited', null)

        if not self.engine:
            raise ConfigError("You somehow managed to call the ConsensusSequences class with an args object that does not\
                               contain an engine variable. Not appropriate.")

        if self.engine not in variability_engines:
            raise ConfigError("Anvi'o does not know how to make sene of the variability engine '%s :/" % self.engine)

        if self.engine == 'NT':
            VariableNtPositionsEngine.__init__(self, args=args, r=self.run, p=self.progress)
        elif self.engine == 'AA':
            VariableAAPositionsEngine.__init__(self, args=args, r=self.run, p=self.progress)

        self.sequence_variants_in_samples_dict = {}


    def populate_seqeunce_variants_in_samples_dict(self):
        """Populates the main dictionary that keeps track of variants for each sample."""

        # no data no play.
        if not self.data:
            raise ConfigError("ConsensusSequences class is upset because it doesn't have any data. There can be two reasons\
                               to this. One, anvi'o variability engines reported nothing (in which case you should have gotten\
                               an error much earler). Two, you are a programmer and failed to call the 'process()' on your\
                               instance from this class. Do you see how the second option is much more likely? :/")

        # learn about the gene seqeunces per gene call.
        gene_sequences = {}
        for gene_callers_id in self.genes_of_interest:
            _, d = self.get_sequences_for_gene_callers_ids([gene_callers_id])
            gene_sequences[gene_callers_id] = d[gene_callers_id]['sequence'].lower()

        # here we populate a dictionary with all the right items but witout any real data.
        sample_names = set([e['sample_id'] for e in self.data.values()])
        for sample_name in sample_names:
            self.sequence_variants_in_samples_dict[sample_name] = {}

            for gene_callers_id in self.genes_of_interest:
                self.sequence_variants_in_samples_dict[sample_name][gene_callers_id] = {'sequence_as_list': list(gene_sequences[gene_callers_id]),
                                                                                        'num_changes': 0, 'gene_callers_id': gene_callers_id,
                                                                                        'in_pos_1': 0, 'in_pos_2': 0, 'in_pos_3': 0}

        # here we will go through every single variant in our data, and correct replace
        # some items in sequences for each sample based on variability infomration.
        self.progress.new('Populating sequence variants in samples data')
        self.progress.update('processing %d variants ...' % len(self.data))
        for entry in list(self.data.values()):
            sample_name = entry['sample_id']
            gene_callers_id = entry['corresponding_gene_call']

            # the dict item we will be playing with
            d = self.sequence_variants_in_samples_dict[sample_name][gene_callers_id]

            reference = entry['reference']
            consensus = entry['consensus']

            if reference != consensus:
                codon_order = entry['codon_order_in_gene']
                base_pos_in_codon = entry['base_pos_in_codon']

                nt_position_to_update = ((codon_order * 3) + base_pos_in_codon) - 1

                # update the entry.
                d['sequence_as_list'][nt_position_to_update] = consensus
                d['num_changes'] += 1
                d['in_pos_%d' % base_pos_in_codon] += 1

        self.progress.end()


    def get_formatted_consensus_sequence_entry(self, key, sample_name, gene_callers_id):
        """Gets a sample naem and gene callers id, returns a well-formatted dict for sequence
           entry using the `self.sequence_variants_in_samples_dict`.

           `key` must be unique string identifier."""

        F = self.sequence_variants_in_samples_dict[sample_name][gene_callers_id]
        return {'key': key,
                'sample_name':  sample_name,
                'gene_caller_id': gene_callers_id,
                'num_changes': F['num_changes'],
                'in_pos_1': F['in_pos_1'],
                'in_pos_2': F['in_pos_2'],
                'in_pos_3': F['in_pos_3'],
                'sequence': ''.join(F['sequence_as_list'])}


    def report(self):
        if not self.sequence_variants_in_samples_dict:
            self.populate_seqeunce_variants_in_samples_dict()

        self.progress.new('Generating the report')
        self.progress.update('...')

        output_d = {}
        counter = 1
        for sample_name in self.sequence_variants_in_samples_dict:
            for gene_callers_id in self.sequence_variants_in_samples_dict[sample_name]:
                d = self.get_formatted_consensus_sequence_entry('e%.7d' % (counter), sample_name, gene_callers_id)
                counter += 1

                if self.tab_delimited_output:
                    output_d[d['key']] = d
                else:
                    key = '|'.join([d['key'],
                                    'sample_name:%s' % d['sample_name'],
                                    'gene_caller_id:%d' % d['gene_caller_id'],
                                    'num_changes:%s' % d['num_changes'],
                                    'in_pos_1:%d' % d['in_pos_1'],
                                    'in_pos_2:%d' % d['in_pos_2'],
                                    'in_pos_3:%d' % d['in_pos_3']])
                    output_d[key] = d['sequence']

        if self.tab_delimited_output:
            utils.store_dict_as_TAB_delimited_file(output_d, self.output_file_path, headers=['entry_id', 'sample_name', 'gene_caller_id', 'num_changes', 'in_pos_1', 'in_pos_2', 'in_pos_3', 'sequence'])
        else:
            utils.store_dict_as_FASTA_file(output_d, self.output_file_path)

        self.progress.end()

        self.run.info('Num genes reported', pp(len(self.genes_of_interest)))
        self.run.info('Num sequences reported', pp(len(self.sequence_variants_in_samples_dict)))
        self.run.info('Output File', self.output_file_path)


class VariabilityNetwork:
    def __init__(self, args={}, p=progress, r=run):
        self.args = args

        self.run = r
        self.progress = p

        self.samples = None
        self.samples_information_dict = None
        self.data = None

        A = lambda x, t: t(args.__dict__[x]) if x in args.__dict__ else None
        null = lambda x: x
        self.input_file_path = A('input_file', null)
        self.samples_information_path = A('samples_information', null)
        self.max_num_unique_positions = A('max_num_unique_positions', int)
        self.output_file_path = A('output_file', null)

        filesnpaths.is_output_file_writable(self.output_file_path)

        if self.samples_information_path:
            filesnpaths.is_file_tab_delimited(self.samples_information_path)
            self.samples_information_dict = utils.get_TAB_delimited_file_as_dictionary(self.samples_information_path)
            num_attributes = len(list(self.samples_information_dict.values())[0])

            self.run.info('samples_information', '%d attributes read for %d samples' % (num_attributes, len(self.samples_information_dict)))

        if self.input_file_path:
            filesnpaths.is_file_tab_delimited(self.input_file_path)
            self.progress.new('Reading the input file')
            self.progress.update('...')
            self.data = utils.get_TAB_delimited_file_as_dictionary(self.input_file_path)
            self.progress.end()

            self.run.info('input_file', '%d entries read' % len(self.data))


    def generate(self):
        if not self.data:
            raise ConfigError("There is nothing to report. Either the input file you provided was empty, or you\
                                haven't filled in the variable positions data into the class.")

        if self.max_num_unique_positions < 0:
            raise ConfigError("Max number of unique positions cannot be less than 0.. Obviously :/")

        self.samples = sorted(list(set([e['sample_id'] for e in list(self.data.values())])))
        self.run.info('samples', '%d found: %s.' % (len(self.samples), ', '.join(self.samples)))

        if self.samples_information_dict:
            samples_missing_in_information_dict = [s for s in self.samples if s not in self.samples_information_dict]
            if len(samples_missing_in_information_dict):
                raise ConfigError("The sample names you provided in the samples information data is not a subset of\
                                    sample names found in the variable positions data :/ Essentially, every sample name\
                                    appears in the variability data must be present in the samples information data,\
                                    however, you are missing these ones from your samples information: %s."\
                                                % (', '.join(samples_missing_in_information_dict)))

        self.unique_variable_nt_positions = set([e['unique_pos_identifier'] for e in list(self.data.values())])
        self.run.info('unique_variable_nt_positions', '%d found.' % (len(self.unique_variable_nt_positions)))

        if self.max_num_unique_positions and len(self.unique_variable_nt_positions) > self.max_num_unique_positions:
            self.unique_variable_nt_positions = set(random.sample(self.unique_variable_nt_positions, self.max_num_unique_positions))
            self.run.info('unique_variable_nt_positions', 'Unique positions are subsampled to %d' % self.max_num_unique_positions, mc='red')

        self.progress.new('Samples dict')
        self.progress.update('Creating an empty one ...')
        samples_dict = {}
        for sample_name in self.samples:
            samples_dict[sample_name] = {}
            for unique_variable_position in self.unique_variable_nt_positions:
                samples_dict[sample_name][unique_variable_position] = 0

        self.progress.update('Updating the dictionary with data')
        for entry in list(self.data.values()):
            sample_id = entry['sample_id']
            pos = entry['unique_pos_identifier']
            frequency = entry['departure_from_reference']

            samples_dict[sample_id][pos] = float(frequency)


        self.progress.update('Generating the network file')
        utils.gen_gexf_network_file(sorted(list(self.unique_variable_nt_positions)), samples_dict, self.output_file_path, sample_mapping_dict=self.samples_information_dict)
        self.progress.end()

        self.run.info('network_description', self.output_file_path)


variability_engines = {'NT': VariableNtPositionsEngine,
                       'AA': VariableAAPositionsEngine}
