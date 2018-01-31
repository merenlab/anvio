# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of single nucleotide variation"""


import os
import sys
import copy
import random
import numpy as np
import pandas as pd

from scipy.stats import entropy

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = ['Alon Shaiber']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


pd.options.display.max_columns=100
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

        if self.gene_caller_id is not None:
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

        # Set items of interest while you are at it. Ensure self.items is sorted alphabetically.
        # This is required for resolving ties in coverage alphabetically, which is described in the
        # docstring of self.insert_additional_fields.
        self.items = constants.amino_acids if self.engine == 'AA' else list(constants.nucleotides)
        self.items = sorted(self.items)

        self.progress.update('Making sure our databases are here ..')
        if not self.profile_db_path:
            raise ConfigError('You need to provide a profile database.')

        if not self.contigs_db_path:
            raise ConfigError('You need to provide a contigs database.')

        self.progress.update('Making sure our databases are compatible ..')
        utils.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

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
            auxiliary_data_file_path = dbops.get_auxiliary_data_path_for_profile_db(self.profile_db_path)
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
            self.data = profile_db.db.get_table_as_dataframe(t.variable_nts_table_name, table_structure=t.variable_nts_table_structure)
            self.progress.end()

        elif self.engine == 'AA':
            if not profile_db.meta['AA_frequencies_profiled']:
                raise ConfigError("It seems AA frequencies were not characterized for this profile database.\
                                    There is nothing to report here for AAs!")
            self.data = profile_db.db.get_table_as_dataframe(t.variable_aas_table_name)
            self.progress.end()

            # append split_name information
            self.data["split_name"] = self.data["corresponding_gene_call"].apply(lambda x: self.gene_callers_id_to_split_name_dict[x])

        # we're done here. bye.
        profile_db.disconnect()


    def check_if_data_is_empty(self):
        if not len(self.data):
            self.progress.end()
            self.run.info_single('Nothing left in the variability data to work with. Quitting :/', 'red', 1, 1)
            sys.exit()


    def apply_advanced_filters(self):
        if self.min_departure_from_consensus:
            self.run.info('Min departure from consensus', self.min_departure_from_consensus)
            self.progress.new('Filtering based on min departure from consensus')
            entries_before = len(self.data.index)
            self.data = self.data[self.data['departure_from_consensus'] >= self.min_departure_from_consensus]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="min departure from consensus")

        if self.max_departure_from_consensus < 1:
            self.run.info('Max departure from consensus', self.max_departure_from_consensus)
            self.progress.new('Filtering based on max departure from consensus')
            entries_before = len(self.data.index)
            self.data = self.data[self.data['departure_from_consensus'] <= self.max_departure_from_consensus]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="max departure from consensus")

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
            self.progress.new('Filtering based on samples of interest')
            entries_before = len(self.data.index)
            self.data = self.data.loc[self.data["sample_id"].isin(self.samples_of_interest)]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="samples of interest")

        if self.genes_of_interest:
            self.run.info('Num genes of interest', pp(len(self.genes_of_interest)))
            self.progress.new('Filtering based on genes of interest')
            entries_before = len(self.data.index)
            self.data = self.data.loc[self.data["corresponding_gene_call"].isin(self.genes_of_interest)]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="genes of interest")

        if self.splits_of_interest:
            self.run.info('Num splits of interest', pp(len(self.splits_of_interest)))
            self.progress.new('Filtering based on splits of interest')
            entries_before = len(self.data.index)
            self.data = self.data.loc[self.data["split_name"].isin(self.splits_of_interest)]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="splits of interest")

        # let's report the number of positions reported in each sample before filtering any further:
        num_positions_each_sample = dict(self.data.sample_id.value_counts())
        self.run.info('Total number of variable positions in samples', '; '.join(['%s: %s' % (s, num_positions_each_sample.get(s, 0)) for s in sorted(self.sample_ids)]))

        if self.min_departure_from_reference:
            self.run.info('Min departure from reference', self.min_departure_from_reference)
            self.progress.new('Filtering based on min departure from reference')
            entries_before = len(self.data.index)
            self.data = self.data.loc[self.data["departure_from_reference"] >= self.min_departure_from_reference]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="min departure from reference")

        if self.max_departure_from_reference < 1:
            self.run.info('Max departure from reference', self.max_departure_from_reference)
            self.progress.new('Filtering based on max departure from reference')
            entries_before = len(self.data.index)
            self.data = self.data.loc[self.data["departure_from_reference"] <= self.max_departure_from_reference]
            entries_after = len(self.data.index)
            self.progress.end()
            self.report_change_in_entry_number(entries_before, entries_after, reason="max departure from reference")

        if self.engine == 'NT':
            self.data['unique_pos_identifier_str'] = self.data['split_name'] + "_" + self.data['pos'].astype(str)
        if self.engine == 'AA':
            self.data['unique_pos_identifier_str'] = self.data['split_name'] + "_" + self.data['corresponding_gene_call'].astype(str) + "_" + self.data['codon_order_in_gene'].astype(str)

        # this could go anywhere now
        self.data['gene_length'] = self.data['corresponding_gene_call'].apply(self.get_gene_length)

        if self.min_occurrence == 1:
            return

        if self.min_occurrence > 1:
            self.run.info('Min occurrence requested', self.min_occurrence)

        self.progress.new('Processing positions further')

        self.progress.update('counting occurrences of each position across samples ...')

        occurrences = self.data["unique_pos_identifier_str"].value_counts()
        entries_before = len(self.data.index)
        self.data = self.data[self.data["unique_pos_identifier_str"].isin(occurrences[occurrences >= self.min_occurrence].index)]
        entries_after = len(self.data.index)
        self.progress.end()
        self.report_change_in_entry_number(entries_before, entries_after, reason="min occurrence")


    def set_unique_pos_identification_numbers(self):
        self.progress.new('Further processing')
        self.progress.update('re-setting unique identifiers to track split/position pairs across samples')

        self.data['unique_pos_identifier'] = self.data['unique_pos_identifier_str'].apply(self.get_unique_pos_identification_number)
        self.data['contig_name'] = self.data['split_name'].apply(lambda split: self.splits_basic_info[split]['parent'])

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


    def insert_additional_fields(self, entry_ids=[]):
        """
        This calculates the following columns: consensus, n2n1 ratio, competing_aas,
        departure_from_consensus, and all substitution scoring matrices.

        NOTE For defining the "consensus" column and the "competing_aas" column (or related column
        for different --engine values), it is important to make it explict how we resolve ties. If
        you finish reading this comment and still do not understand, then I have failed you. It's
        an issue because suppose Ala and Trp are tied for sharing the most reads. Is the consensus
        Ala or Trp? Is competing_aas AlaTrp or TrpAla? In a separate example, if Ser is most
        common and Gly and Glu are tied for second, should Gly or Glu be a part of competing_aas
        and should it be GlxSer or SerGlx? There are three rules that define our conventions:

            1. Competing_aas ALWAYS appear in alphabetical order. Even if Cys is most common, and
               Ala is second most commond, competing_aas = AlaCys.  
            2. Ties are always resolved alphabetically. If there is a 3-way tie for second between
               His, Met, and Thr, the item including in competing_aas will be His.
            3. If the coverage of the second-most common item is 0, the most common is paired with
               itself.

        NOTE if the coverage is 0 (which is stupid, but it can exist if both
        --min-coverage-in-each-sample is 0 and --quince-mode is active), departure_from_consensus
        and n2n1ratio are set to 0.

        NOTE To get a first glance at variation, we display competing_nts during inspect mode in the
        interactive interface, which means they are stored. I (Evan) don't know how the
        competing_nts are computed during this process, but after looking at 3 datasets (Infant Gut,
        E. Faecalis; TARA Oceans, HIMB083; and Mushroom Spring, Synechococcus), there were no "N"
        counts (--engine NT), which makes me think that consensus and competing_nts are calculated
        without consideration of "N" values. In contrast, the variability table is not prejudiced to
        N, and treats it like any other item in self.items. This means it can be the consensus
        value, or be one of the items in competing_nts, etc.
        """

        # First, we just make sure that whatever operations we have performed on self.data, we have
        # not altered the types of coverage values from int. I am learning that with pandas
        # DataFrames, it is WORTH being explicit with types, even at the cost of redundancy.
        self.data.loc[:, self.items] = self.data.loc[:, self.items].astype(int)

        if not len(entry_ids):
            entry_ids = list(self.data.index)

        # index entries with and without non-zero coverage (introduced by --quince-mode)
        coverage_zero = self.data.index.isin(entry_ids) & (self.data["coverage"] == 0)
        coverage_nonzero = self.data.index.isin(entry_ids) & (self.data["coverage"] > 0)

        # rank the items of each entry from 1 - 21 (for --engine AA) based on item coverage.
        # method="first" ensures alphabetic ordering in the case of ties. Convert the rank DataFrame
        # into a numpy array, and find the order of indices that sort each entry's items based on
        # their ranks. type(ranks) = pd.DataFrame, type(item_index_order) = numpy array
        ranks = self.data.loc[entry_ids, self.items].rank(ascending=False, axis=1, method="first").astype(int)
        item_index_order = np.argsort(ranks.values, axis=1)

        # the first and second most common items, according to the 2nd convention in the docstring,
        # are now defined for each entry in two pandas Series.
        items_first_and_second  = self.data.loc[entry_ids, self.items].columns[item_index_order[:,:2]]

        # we also calculate the coverage values for the first and second most common items
        sorted_coverage = np.sort(self.data.loc[entry_ids, self.items].values, axis=1)
        coverages_first =  pd.Series(sorted_coverage[:,-1], index=self.data.loc[entry_ids].index)
        coverages_second = pd.Series(sorted_coverage[:,-2], index=self.data.loc[entry_ids].index)

        # define consensus as the first most common item (hence `[:,0]`)
        self.data.loc[entry_ids, "consensus"] = items_first_and_second[:,0]

        # if the coverage is zero, departure_from_consensus = 0
        self.data.loc[coverage_zero, "departure_from_consensus"] = 0
        self.data.loc[coverage_nonzero, "departure_from_consensus"] = \
                    (self.data.loc[coverage_nonzero, "coverage"] - coverages_first) / self.data.loc[coverage_nonzero, "coverage"]

        # if the coverage is zero, n2n1ratio = 0
        self.data.loc[coverage_zero, "n2n1ratio"] = 0
        self.data.loc[coverage_nonzero, "n2n1ratio"] = coverages_second[coverage_nonzero] / coverages_first[coverage_nonzero]

        # we name the competing_items column differently depending on the engine used
        competing_items = "competing_aas" if self.engine=="AA" else "competing_nts" if self.engine=="NT" else "competing_codons"

        # This step computes the "competing_items" column. First, convert items_first_and_second
        # [which has shape (len(entry_ids), 2)] into a numpy array, then sort them alphabetically
        # (in order to comply with convention 1 in docstring). Then, sum along the sorted axis.
        # Since the elements are of type str, the sum operator concatenates the strings together.
        # Finally, if the coverage of the most common item is equal to the total coverage, we pair
        # the most common item with itself.
        self.data.loc[entry_ids, competing_items] = np.sum(np.sort(items_first_and_second.values, axis=1), axis=1) # V/\
        self.data.loc[self.data.index.isin(entry_ids) & (self.data["coverage"] == self.data[self.items].max(axis=1)), competing_items] = self.data["consensus"]*2

        # Loop through each SSM, filling each corresponding column entry by entry using the `apply`
        # operator. Instead of using self.substitution_scoring_matrices[m], we speed things up by
        # converting to a dictionary we call `substitution_scoring_matrix` that takes as its input
        # an entry from competing_items instead of having to pull from the first AND second most
        # common item, and then indexing a triple-nested dictionary.
        for m in self.substitution_scoring_matrices:
            substitution_scoring_matrix = utils.convert_SSM_to_single_accession(self.substitution_scoring_matrices[m])
            self.data.loc[entry_ids, m] = self.data.loc[entry_ids, competing_items].apply(lambda x: substitution_scoring_matrix.get(x, None))


    def filter_based_on_scattering_factor(self):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from eachother."""

        if self.min_scatter == 0:
            return

        # FIXME THIS FUNCTION DOES NOT DO WHAT IT IS SUPPOSED TO.
        raise ConfigError("Woops. The function that handles --min-scatter doesn't do \
                           what we thought it did. This will be fixed soon. Sorry for \
                           the inconvenience.")

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


    def report_change_in_entry_number(self, num_before, num_after, reason="unknown reason", added=False):
        """
        Reports how many entries were removed during a filtering step. If added=True, then it
        reports the number gained (quince mode)
        """
        changed = "removed" if not added else "added"

        self.run.info('Entries after "%s"' % (reason),
                         '%s (%s were %s)' % (pp(num_after),
                                              pp(abs(num_before - num_after)),
                                              changed),
                      mc='green')

        self.check_if_data_is_empty()


    def filter_based_on_minimum_coverage_in_each_sample(self):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from eachother."""

        if self.min_coverage_in_each_sample < 1:
            return

        self.run.info('Min coverage in all samples', '%dX' % self.min_coverage_in_each_sample)

        self.progress.new('Examining coverage of each variable position in each sample')
        self.progress.update('...')

        num_entries_before = len(self.data.index)

        # get the minimum coverage for each unique_pos_identifier
        min_cov_each_pos = self.data.groupby("unique_pos_identifier")["coverage"].min()

        # get a list of all unique_pos_identifiers with min cov > self.min_coverage_in_each_sample].index
        pos_to_keep = list(min_cov_each_pos[min_cov_each_pos >= self.min_coverage_in_each_sample].index)

        # only include entries with unique_pos_identifier in pos_to_keep
        self.data = self.data[self.data["unique_pos_identifier"].isin(pos_to_keep)]

        num_entries_after = len(self.data.index)

        self.progress.end()

        self.report_change_in_entry_number(num_entries_before, num_entries_after, reason="min cov for all samples")


    def filter_based_on_num_positions_from_each_split(self):
        if self.num_positions_from_each_split > 0:
            self.run.info('Num positions to keep from each split', self.num_positions_from_each_split)
        else:
            self.run.info('Num positions to keep from each split', '(all positions)')
            return

        self.progress.new('Filtering based on -n')

        self.progress.update('Randomly subsampling from splits with num positions > %d' % self.num_positions_from_each_split)

        subsample_func = lambda x: pd.Series(x.unique()) if len(x.unique()) <= self.num_positions_from_each_split else\
                                   pd.Series(np.random.choice(x.unique(), size=self.num_positions_from_each_split, replace=False))
        unique_positions_to_keep = self.data.groupby('split_name')['unique_pos_identifier'].apply(subsample_func)

        entries_before = len(self.data.index)
        self.data = self.data[self.data['unique_pos_identifier'].isin(unique_positions_to_keep)]
        entries_after = len(self.data.index)

        self.progress.end()
        self.report_change_in_entry_number(entries_before, entries_after, reason="num positions each split")


    def compute_comprehensive_variability_scores(self):
        """
            Comprehensive stats are defined as scores that take into consideration the entire vector of variability and
            not only the two most competing items (and thus it is comprehensive).  Currently the scores that are
            included are: site-entropy, site-Kullback-Leibler divergence (both a raw score and a normalized score (see
            below)), and weighted substitution scores (i.e. BLOSUM).
            
            site-entropy 
            ============
            The entropy of the items at a single position in a single sample. This means that each entry of the
            variability table receives its own site-entropy score (in the table we just called it "entropy"). If a site
            has no coverage, which can happen if --quince-mode is enabled, the value of entropy is -np.inf, something we
            maintain until exporting the table as a tab-delimited file, at which point we recast them to something
            reasonable.

            Kullback-Leibler divergence raw
            ===============================
            the Kullback-Leibler divergence of the frequencies in a sample compared to the raw frequencies of the sum of
            occurrences in the same site accross samples.

            Kullback-Leibler divergence normalized
            ======================================
            The Kullback-Leibler divergence of the frequencies in a sample compared to the frequencies of the sum of
            normalized occurances in the same site accross samples. Where the normalization is such that the occernce of
            items is normalized to sum to one in each sample. This method eliminates the effect of coverage on the
            score. The disadvantage of this method is that if there is a sample with low coverage then any noise (like a
            single sequencing error) could have a major effect. It is recommended to use this score in combination with
            the --min-coverage-in-each-sample.
            
            Weighted substitution scores
            ============================
            The weights per substitution score is weighted by the joint frequency of the items i.e. sum(S_{i,j}*pi*pj)
            where i does not equal j (i.e. the substitution of an item with itself is not considered)
        """

        if self.skip_comprehensive_variability_scores:
            self.run.warning("Anvi'o will skip comprehensive variability score computations.")
            return

        # suppress division by zero runtime warning
        np.seterr(invalid='ignore')

        if not self.quince_mode:
            self.run.warning("Only some comprehensive variability score computations can only be done without `--quince-mode`")

        self.progress.new("Comprehensive stats")
        self.progress.update("Those that don't require --quince-mode")

        self.comprehensive_stats_headers = [m + '_weighted' for m in self.substitution_scoring_matrices] + \
                                           ['entropy', 'kullback_leibler_divergence_raw', 'kullback_leibler_divergence_normalized']

        # Pandas is fun, but numpy is fast. Here we convert the coverage table information from the DataFrame to a
        # numpy array. The transpose is required because scipy.stats entropy function calculates along an
        # unspecifiable axis that we must conform to. But before any of this is done we order the entries according
        # to unique_pos_identifier (and for a given unique_pos_identifier, entries are ordered alphabetically by
        # sample_id). The reason for this is aesthetic but also required for vectorized operations that occur after
        # self.progress.update("Those that do require --quince-mode")
        self.data = self.data.sort_values(by=["unique_pos_identifier","sample_id"])
        coverage_table = self.data[self.items].T.astype(int).as_matrix()

        # Now we compute the entropy, which is defined at a per position, per sample basis. There is a reason we
        # pass coverage_table instead of a normalized table. If we pass a normalized table scipy.stat.entropy complains
        # that there is division by zero for entries introduced by --quince-mode that have 0 coverage.  By passing
        # coverage_table instead, scipy.stats.entropy does the normalization itself ensures that such entries return
        # -inf instead of raising an error.
        self.data["entropy"] = entropy(coverage_table)

        # Next we worry ourselves with weighted substitution scores. We convert each SSM to a numpy array and store
        # each of them in a dictionary with keys equal to the SSM names, e.g. BLOSUM90. The index of numpy arrays
        # are organized alphabetically. For example, for an amino acid substitution matrix, the substitution
        # Ala->Ala is indexed by [0,0], whereas Val->Val is indexed by [-1,-1]. We decided it makes sense to set the
        # substitution score of an item with itself to zero. This way we only consider substitutions to other items
        # (and don't consider substitution of an item with itself)
        numpy_substitution_scoring_matrices = {SMM: np.array([[matrix[i][j] if j != i else 0 for j in sorted(matrix[i])] for i in sorted(matrix)]) for SMM, matrix in self.substitution_scoring_matrices.items()}

        # We loop through each of the substitution scoring matrices available. It's possible the items in the
        # substitution matrices aren't a complete set of the items we have coverage data for. For example, BLOSUM90
        # doesn't have substitution scores for the STP codon (yet we have coverage data for STP). Available items
        # can vary from SSM to SSM. To deal with this, we subset the coverage_table to only include available items
        # for a given SSM. We then transform this subset of the coverage data to a frequency table that's
        # potentially unique to the SSM. To avoid calculating frequency tables unnecessarily, the previous SSM's
        # items are stored in `previous_indices` so the table is only recalculated if the current SSM items vary
        # from the previous SSM items.
        previous_indices = None
        for SMM, matrix in numpy_substitution_scoring_matrices.items():

            # initialize the entropy column in self.data with np.nans
            self.data[SMM + "weighted"] = np.nan

            # We find which items in self.items are in the SSM and record the index at which they appear in
            # self.items. By definition this is also the index ordering of coverage_table. For example, the index of
            # Ala is 0 so the coverage data for Ala is coverage_table[0,:]. Hence, the subset of coverage_table for
            # the given SSM is coverage_table[indices, :].
            indices = sorted([self.items.index(item) for item in self.substitution_scoring_matrices[SMM]])
            if indices != previous_indices:

                # To get the frequency table we divide each column by the sum of the rows (the array called
                # total_coverage). But we must be careful, because it's possible for some entries of total_coverage
                # to be 0, which can occur if both a) --quince-mode is enabled and b) --min-coverage-in-each-sample
                # is 0 (the default). If this is the case we set the total_coverage entry to -1 so that the
                # frequency for each item in the entry becomes 0/-1 = 0, instead of producing a NaN due to division
                # by zero. 
                total_coverage = np.sum(coverage_table[indices,:], axis=0)
                freq_table = np.divide(coverage_table[indices,:], np.where(total_coverage==0, -1, total_coverage))

                # While in this loop, we define a Boolean array with length equal to the number of entries in
                # freq_table. It is True if the substitution score should be calculated and False otherwise. What
                # can cause an entry not to be calculated? First of all, if there is no coverage data then there is
                # no substitution to report. Secondly, since we don't consider the substitution of an item with
                # itself (Remember? We set all diagonals in each SSM to 0), we don't report scores for entries that
                # only have 1 item with non-zero coverage. Both of these conditions are encapsulated nicely with the
                # np.count_nonzero function. We then subset the freq_table to only include these entries but we hang
                # onto the Boolean array for when we put the scores into self.data
                keep = np.count_nonzero(freq_table, axis=0) > 1
                freq_table = freq_table[:, keep]

                # The last thing we do in here is calculate a normalization factor for the entries, since it is
                # unique to each freq_table. A normalization score is needed since we don't consider a substitution
                # of an item with itself. Hence, the sum of frequencies doesn't sum to 1, and so to make sure they
                # sum to 1 we multiply by this normalization factor. We want these to sum to 1 because otherwise
                # these are not valid weights.
                normalization = 1 / (1 - np.sum(np.square(freq_table), axis=0))

            # This is legitimately legendary status array broadcasting. I can't explain how it works but the
            # quantity sum(S_{i,j}*pi*pj) is being calculated for each entry, where pi is the frequency of the ith
            # item, pj is the frequency of the jth item, and S_{i,j} is the substitution matrix score for the i->j
            # substitution. If you remember, we already set i=j entries to zero so they contribute zero weight to
            # the score.
            self.data.loc[keep, SMM + "_weighted"] = normalization * np.sum(freq_table[np.newaxis, :].T * freq_table[:, np.newaxis].T * matrix, axis=(1,2))

        if not self.quince_mode:
            self.progress.end()
            self.comprehensive_variability_scores_computed = True
            return

        self.progress.update("Those that do require --quince-mode")

        # Due to --quince-mode every unique position identifier has the same number of entries (this screams
        # vectorization). We abuse this to make a 3 dimensional numpy array, coverage_by_pos. The zeroth axis
        # indexes the items, the first indexes the unique_pos_identifier, and the second indexes the sample number.
        # This reshaping operation works only because at the start of this function we ordered entries by
        # unique_pos_identifier.
        numpy_dimension = (len(self.items), self.data["unique_pos_identifier"].nunique(), len(self.sample_ids))
        coverage_by_pos = coverage_table.reshape(numpy_dimension)

        # We also define a normalized version of coverage_by_pos so that each entries defines the frequency
        # (probability) of a certain item occuring, rather than the raw counts. This is used for the normalized
        # Kullback-Leibler divergence. If the entry is all zeros, the frequencies for the entry are all defined to
        # be 0. NOTE If this is not done, any position where one or more of the samples has 0 coverage yields a
        # normalized kullback-leibler value of inf. I think it makes more sense this way.
        counts_per_pos_per_sample = np.sum(coverage_by_pos, axis=0)
        #freq_by_pos = np.divide(coverage_by_pos, np.where(counts_per_pos_per_sample==0, -1, counts_per_pos_per_sample))
        freq_by_pos = np.divide(coverage_by_pos, counts_per_pos_per_sample)

        # The entropy from scipy.stats operates on axis = 0, so the returned array dimension is flattened in the
        # zeroth dimension compared to the input array. For example, if the input is shape (X,Y,Z), the output is
        # shape (Y,Z). That means when we pass coverage_by_pos, an entropy array is returned with the zeroth axis
        # indexed by unique_pos_identifiers and the first axis indexed by sample_ids. The entropy function insists
        # that our reference distributions (the normalized and unnormalized mean frequency counts across samples,
        # see docstring for more details) must ALSO have the shape (X,Y,Z). The issue with this is that by
        # calculating the mean over samples we collapse the Z dimension so the shape of our reference distribution
        # array is (X,Y). We solve this issue by stacking Z identical reference distribution arrays to create a
        # pseudo-second axis so the final shape is (X,Y,Z). coverage_summed_over_samples is the reference
        # distribution array for the raw Kullback-Leibler divergence and freq_summed_over_samples is for the
        # normalized Kullback-Leibler divergence.
        coverage_summed_over_samples = np.repeat(np.sum(coverage_by_pos, axis=2)[:,:,np.newaxis], len(self.sample_ids), axis=2)
        freq_summed_over_samples     = np.repeat(np.sum(freq_by_pos,     axis=2)[:,:,np.newaxis], len(self.sample_ids), axis=2)

        # As mentioned in last comment, a 2D array is returned by entropy, which we flatten into 1D.  The flattened
        # array is equal to the length of entries in self.data. Furthermore, the order of entropy values is the same
        # order as the entries they correspond to, so all we do is assign the arrays to new columns in self.data and
        # we're done.
        self.data['kullback_leibler_divergence_raw'] = entropy(coverage_by_pos, coverage_summed_over_samples).flatten()
        self.data['kullback_leibler_divergence_normalized'] = entropy(freq_by_pos, freq_summed_over_samples).flatten()

        # Vape Nation V/\
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

        # key = unique_pos_identifier, val = corresponding_gene_call
        return self.data[["unique_pos_identifier","corresponding_gene_call"]].\
               drop_duplicates().set_index("unique_pos_identifier").to_dict()["corresponding_gene_call"]


    def get_unique_pos_identifier_to_codon_order_in_gene(self):
        self.progress.update('populating a dict to track codon order in genes for each unique position')
        # key = unique_pos_identifier, val = codon_order_in_gene
        return self.data[["unique_pos_identifier","codon_order_in_gene"]].\
               drop_duplicates().set_index("unique_pos_identifier").to_dict()["codon_order_in_gene"]


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

        # Update entry_id with sequential numbers based on the final ordering of the data:
        self.data.reset_index(drop=True, inplace=True)
        self.data["entry_id"] = self.data.index

        self.progress.update('exporting variable positions table as a TAB-delimited file ...')
        utils.store_dataframe_as_TAB_delimited_file(self.data, self.args.output_file, columns=new_structure)
        self.progress.end()

        self.run.info('Num entries reported', pp(len(self.data.index)))
        self.run.info('Output File', self.output_file_path)
        self.run.info('Num %s positions reported' % self.engine, self.data["unique_pos_identifier"].nunique())


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
        next_available_entry_id = self.data["entry_id"].max() + 1

        unique_pos_identifier_to_corresponding_gene_id = self.get_unique_pos_identifier_to_corresponding_gene_id()

        unique_pos_identifier_to_codon_order_in_gene = self.get_unique_pos_identifier_to_codon_order_in_gene()
        self.progress.update('creating a dict to track missing base frequencies for each sample / split / pos')
        split_pos_to_unique_pos_identifier = {}
        splits_to_consider_dict = {}
        for split_name in splits_wanted:
            splits_to_consider_dict[split_name] = {}
            split_pos_to_unique_pos_identifier[split_name] = {}

        self.progress.update('populating the dict to track missing base frequencies for each sample / split / pos')
        for entry_id, v in self.data.iterrows():
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
        new_entries = {}
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
                    new_entries[next_available_entry_id] = {'entry_id': next_available_entry_id,
                                                            'contig_name': contig_name_name,
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
                    new_entries[next_available_entry_id][base_at_pos] = split_coverage_across_samples[sample][pos]
                    next_available_entry_id += 1

        # convert to pandas DataFrame (its much faster to build and convert a dictionary than to build DataFrame row by row)
        new_entries = pd.DataFrame(new_entries).T
        new_entries.set_index("entry_id", drop=False, inplace=True)

        # concatenate new columns to self.data
        self.data = pd.concat([self.data, new_entries])

        # fill in additional fields for new entries
        self.insert_additional_fields(list(new_entries.index))
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
        next_available_entry_id = self.data["entry_id"].max() + 1

        unique_pos_identifier_str_to_consenus_codon = {}
        unique_pos_identifier_str_to_unique_pos_identifier = {}
        for _, e in self.data.iterrows():
            upi = e['unique_pos_identifier_str']
            unique_pos_identifier_str_to_consenus_codon[upi] = e['reference']
            unique_pos_identifier_str_to_unique_pos_identifier[upi] = e['unique_pos_identifier']

        self.progress.update('creating a dict to track missing AA frequencies for each sample / split / pos')

        splits_to_consider_dict = {}
        for split_name in splits_wanted:
            splits_to_consider_dict[split_name] = {}

        self.progress.update('populating the dict to track missing AA frequencies for each sample / split / pos')
        for entry_id, v in self.data.iterrows():
            gene_codon_key = '%d_%d' % (v['corresponding_gene_call'], v['codon_order_in_gene'])
            d = splits_to_consider_dict[v['split_name']]

            if gene_codon_key in d:
                d[gene_codon_key].remove(v['sample_id'])
            else:
                d[gene_codon_key] = copy.deepcopy(samples_wanted)
                d[gene_codon_key].remove(v['sample_id'])

        split_names_to_consider = list(splits_to_consider_dict.keys())
        num_splits = len(split_names_to_consider)
        new_entries = {}
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

                    new_entries[next_available_entry_id] = {'entry_id': next_available_entry_id,
                                                            'unique_pos_identifier_str': unique_pos_identifier_str,
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
                    new_entries[next_available_entry_id]['coverage'] = coverage

                    # DEALING WITH AAs ##################################################################
                    # here we need to put all the codons into the data table for this sample
                    for codon in set(constants.codon_to_AA.values()):
                        new_entries[next_available_entry_id][codon] = 0

                    # and finally update the frequency of the reference codon with the coverage (WHICH IS VERY BAD,
                    # WE HAVE NO CLUE WHAT IS THE ACTUAL COVERAGE OF TRIPLICATE LINKMERS):
                    new_entries[next_available_entry_id][reference_codon] = coverage

                    next_available_entry_id += 1

        # convert to pandas DataFrame (its faster to build and convert a dictionary than to build
        # DataFrame row by row).
        new_entries = pd.DataFrame(new_entries).T

        # before concatenating the new entries, store the self.data column order. Also, check that
        # no columns exist in new_entries but not in self.data. This is unacceptable, and could have
        # happened if code for new_entries was changed or if the workflow in process() is
        # significantly reworked.
        column_order = self.data.columns.tolist()
        if len([x for x in new_entries.columns.tolist() if x not in self.data.columns.tolist()]):
            raise ValueError("Columns found in new_entries exist that aren't in self.data.")

        # concatenate new columns to self.data
        entries_before = len(self.data.index)
        self.data = pd.concat([self.data, new_entries])
        new_entries.set_index("entry_id", drop=False, inplace=True)
        self.data = self.data[column_order]
        entries_after = len(self.data.index)

        # fill in additional fields for new entries. insert_additional_fields takes a list of
        # entry_ids to consider for self.data, which here is provided from new_entries (what I'm
        # saying is new_entries is not passed, only the entry_id's in new_entries
        self.insert_additional_fields(list(new_entries["entry_id"]))

        self.progress.end()

        self.report_change_in_entry_number(entries_before, entries_after, reason="quince mode", added=True)


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
        if not len(self.data):
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
        sample_names = set(self.data['sample_id'])
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
        for idx, entry in self.data.iterrows():
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
