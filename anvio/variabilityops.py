# -*- coding: utf-8

"""Classes to make sense of single nucleotide variation"""

from __future__ import division
from collections import Counter

import os
import sys
import copy
import random

import anvio
import anvio.tables as t
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


pp = terminal.pretty_print
progress = terminal.Progress()
run = terminal.Run(width = 60)


class VariablePositionsEngine:
    """This is the main class to make sense and report variability for a given set of splits,
       or a bin in a collection, across multiple or all samples. The user can scrutinize the
       nature of the variable positions to be reported dramatically given the ecology and/or
       other biologically-relevant considerations, or purely technical limitations such as
       minimum coverage of a given nucleotide position or the ratio of the competing nts at a
       given position. The default entry to this class is the `anvi-gen-variability-profile`
       program."""

    def __init__(self, args = {}, p=progress, r=run):
        self.args = args
        self.splits_of_interest = set([])
        self.samples_of_interest = set([])

        self.run = r
        self.progress = p

        A = lambda x, t: t(args.__dict__[x]) if args.__dict__.has_key(x) else None
        null = lambda x: x
        self.bin_id = A('bin_id', null)
        self.collection_id = A('collection_id', null)
        self.splits_of_interest_path = A('splits_of_interest', null)
        self.min_ratio = A('min_ratio', float)
        self.min_occurrence = A('min_occurrence', int)
        self.num_positions_from_each_split = A('num_positions_from_each_split', int)
        self.min_scatter = A('min_scatter', int)
        self.min_coverage_in_each_sample = A('min_coverage_in_each_sample', int)
        self.profile_db_path = A('profile_db', null)
        self.contigs_db_path = A('contigs_db', null)
        self.quince_mode = A('quince_mode', bool)
        self.output_file_path = A('output_file', null)
        self.samples_of_interest_path = A('samples_of_interest', null)

        self.variable_positions_table = {} 
        self.merged_split_coverage_values = None
        self.unique_pos_identifier = 0
        self.split_name_position_dict = {}
        self.unique_pos_id_to_entry_id = {}
        self.contig_sequences = None
        self.input_file_path = None


    def init(self):
        self.progress.new('Init')

        self.progress.update('Checking the output file path ..')
        if self.output_file_path:
            filesnpaths.is_output_file_writable(self.output_file_path)

        self.progress.update('Checking the samples of interest ..')
        if self.samples_of_interest_path:
            filesnpaths.is_file_exists(self.samples_of_interest_path)
            self.samples_of_interest = set([s.strip() for s in open(self.samples_of_interest_path).readlines()])
        else:
            self.samples_of_interest = set([])

        self.progress.update('Making sure our databases are here, and they are compatible ..')
        if not self.profile_db_path:
            raise ConfigError, 'You need to provide a profile database.'

        if not self.contigs_db_path:
            raise ConfigError, 'You need to provide a contigs database.'

        if self.min_coverage_in_each_sample and not self.quince_mode:
            self.progress.end()
            raise ConfigError, "When you sepecify a coverage value through --min-coverage-in-each-sample, you must also\
                                use --quince-mode flag, since the former parameter needs to know the coverage values in all\
                                samples even if variation is reported for only one sample among otheres. This is the only way\
                                to figure out whether variation is not reported for other samples due to low or zero coverage,\
                                or there was no variation to report despite the high coverage. Anvi'o could turn --quince-mode\
                                flat automatically for you, but then it is much better if you have full control and understaning\
                                of what is going on."

        dbops.is_profile_db_and_contigs_db_compatible(self.profile_db_path, self.contigs_db_path)

        self.progress.update('Attempting to get our splits of interest sorted ..')
        if self.collection_id:
            # the user wants to go with the collection id path. fine. we will get our split names from
            # the profile database.
            if not self.bin_id:
                self.progress.end()
                raise ConfigError, 'When you declare a collection id, you must also declare a bin name\
                                    (from which the split names of interest will be acquired)'
            if self.splits_of_interest or self.splits_of_interest_path:
                self.progress.end()
                raise ConfigError, "You declared a collection id and one or more bin names so anvi'o can find out\
                                    splits of interest, but you also have specified informaiton for split names?\
                                    This is confusing. You should choose one way or another :/"

            self.splits_of_interest = ccollections.GetSplitNamesInBins(self.args).get_split_names_only()
        else:
            # OK. no collection id. we will go oldschool. we whope to find what we are looking for in
            # self.splits_of_interst_path  at this point (which may have been filled through the command
            # line client), or in self.splits_of_interest (which may have been filled in by another program)
            if not self.splits_of_interest:
                if not self.splits_of_interest_path:
                    self.progress.end()
                    raise ConfigError, 'You did not declare a source for split names. You either should give me\
                                        a file with split names you are interested in, or a collection id and\
                                        bin name so I can learn split names from the profile database.'
                filesnpaths.is_file_exists(self.splits_of_interest_path)
                self.splits_of_interest = set([c.strip().replace('\r', '') for c in open(self.splits_of_interest_path).readlines()])

        self.input_file_path = '/' + '/'.join(os.path.abspath(self.profile_db_path).split('/')[:-1])

        self.progress.update('Reading variable positions table ...')
        profile_db = dbops.ProfileDatabase(self.profile_db_path)
        self.sample_ids = profile_db.samples # we set this now, but we will overwrite it with args.samples_of_interest if necessary
        self.variable_positions_table = profile_db.db.get_table_as_dict(t.variable_positions_table_name)
        db_hash = profile_db.meta['contigs_db_hash']
        profile_db.disconnect()

        if self.quince_mode:
            self.progress.update('Accessing auxiliary data file ...')
            auxiliary_data_file_path = os.path.join(os.path.dirname(self.profile_db_path), 'AUXILIARY-DATA.h5')
            if not os.path.exists(auxiliary_data_file_path):
                raise ConfigError, "Anvi'o needs the aixiliary data file to run this program with '--quince-more' flag.\
                                    However it wasn't found at '%s' :/" % auxiliary_data_file_path
            self.merged_split_coverage_values = auxiliarydataops.AuxiliaryDataForSplitCoverages(auxiliary_data_file_path, db_hash)

        self.progress.update('Reading splits info table ...')
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        self.splits_info_table = contigs_db.db.get_table_as_dict(t.splits_info_table_name)
        self.num_splits_in_db = len(self.splits_info_table)
        if self.quince_mode:
            self.progress.update('Reading contig sequences table ...')
            self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
        contigs_db.disconnect()

        self.progress.end()

        self.process_variable_positions_table()

        self.set_unique_pos_identification_numbers() # which allows us to track every unique position across samples

        self.filter_based_on_scattering_factor()
        
        self.filter_based_on_num_positions_from_each_split()

        if self.quince_mode: # will be very costly...
            self.recover_base_frequencies_for_all_samples() 

        self.filter_based_on_minimum_coverage_in_each_sample()


    def filter(self, filter_name, test_func):
        self.progress.new('Filtering based on "%s"' % filter_name)
        num_entries_before_filter = len(self.variable_positions_table)

        entry_ids_to_remove, counter = set([]), 0

        for entry_id in self.variable_positions_table:
            if counter % 1000 == 0:
                self.progress.update('identifying entries to remove :: %s' % pp(counter))

            counter += 1

            if test_func(self.variable_positions_table[entry_id]):
                entry_ids_to_remove.add(entry_id)
                continue

        self.progress.update('removing %s entries from table ...' % pp(len(entry_ids_to_remove)))
        for entry_id in entry_ids_to_remove:
            self.variable_positions_table.pop(entry_id)

        num_entries_after_filter = len(self.variable_positions_table)

        self.progress.end()

        self.run.info('Remaining entries after "%s" filter' % filter_name,
                      '%s (filter removed %s entries)' % (pp(num_entries_after_filter),
                                                          pp(num_entries_before_filter - num_entries_after_filter)),
                      mc = 'green')

        self.check_if_variable_table_is_empty()


    def check_if_variable_table_is_empty(self):
        if not len(self.variable_positions_table):
            self.progress.end()
            self.run.info_single('No variable positions left to work with. Quitting.', 'red', 1, 1)
            sys.exit()


    def process_variable_positions_table(self):
        self.run.info('Variability table', '%s entries in %s splits across %s samples'\
                % (pp(len(self.variable_positions_table)), pp(self.num_splits_in_db), pp(len(self.sample_ids))))

        self.run.info('Samples in the profile db', ', '.join(sorted(self.sample_ids)))
        if self.samples_of_interest:
            samples_missing_from_db = [sample for sample in self.samples_of_interest if sample not in self.sample_ids]

            if len(samples_missing_from_db):
                raise ConfigError, 'One or more samples you are interested in seem to be missing from\
                                    the profile database: %s' % ', '.join(samples_missing_from_db)

            self.run.info('Samples of interest', ', '.join(sorted(list(self.samples_of_interest))))
            self.sample_ids = self.samples_of_interest
            self.filter('samples of interest', lambda x: x['sample_id'] not in self.samples_of_interest)

        if self.splits_of_interest:
            self.run.info('Num splits of interest', pp(len(self.splits_of_interest)))
            self.filter('splits of interest', lambda x: x['split_name'] not in self.splits_of_interest)

        # let's report the number of positions reported in each sample before filtering any futrher:
        num_positions_each_sample = Counter([v['sample_id'] for v in self.variable_positions_table.values()])
        self.run.info('Total number of variable positions in samples', '; '.join(['%s: %s' % (s, num_positions_each_sample[s]) for s in sorted(self.sample_ids)]))

        if self.min_ratio:
            self.run.info('Min departure from consensus ratio', self.min_ratio)
            self.filter('n2/n1', lambda x: x['departure_from_consensus'] < self.min_ratio)

        for entry_id in self.variable_positions_table:
            v = self.variable_positions_table[entry_id]
            v['unique_position_id'] = '_'.join([v['split_name'], str(v['pos'])])

        if self.min_occurrence == 1:
            return

        if self.min_occurrence > 1:
            self.run.info('Min occurrence requested', self.min_occurrence)

        self.progress.new('Processing positions further')

        self.progress.update('counting occurrences of each position across samples ...')
        unique_position_id_occurrences = {}
        for entry_id in self.variable_positions_table:
            v = self.variable_positions_table[entry_id]
            if unique_position_id_occurrences.has_key(v['unique_position_id']):
                unique_position_id_occurrences[v['unique_position_id']] += 1
            else:
                unique_position_id_occurrences[v['unique_position_id']] = 1

        self.progress.update('identifying entries that occurr in less than %d samples ...' % (self.min_occurrence))
        entry_ids_to_remove = set([])
        for entry_id in self.variable_positions_table:
            v = self.variable_positions_table[entry_id]
            if not unique_position_id_occurrences[v['unique_position_id']] >= self.min_occurrence:
                entry_ids_to_remove.add(entry_id)

        self.progress.update('removing %s entries from table ...' % pp(len(entry_ids_to_remove)))
        for entry_id in entry_ids_to_remove:
            self.variable_positions_table.pop(entry_id)

        self.progress.end()

        self.check_if_variable_table_is_empty()


    def get_unique_pos_identification_number(self, unique_position_id):
        if unique_position_id in self.split_name_position_dict:
            return self.split_name_position_dict[unique_position_id]
        else:
            self.split_name_position_dict[unique_position_id] = self.unique_pos_identifier
            self.unique_pos_identifier += 1
            return self.split_name_position_dict[unique_position_id]


    def gen_unique_pos_identifier_to_entry_id_dict(self):
        self.progress.new('Generating the `unique pos` -> `entry id` dict')
        self.progress.update('...')

        self.unique_pos_id_to_entry_id = {}

        for entry_id in self.variable_positions_table:
            v = self.variable_positions_table[entry_id]
            u = v['unique_position_id']
            if u in self.unique_pos_id_to_entry_id:
                self.unique_pos_id_to_entry_id[u].add(entry_id)
            else:
                self.unique_pos_id_to_entry_id[u] = set([entry_id])

        self.progress.end()


    def set_unique_pos_identification_numbers(self):
        self.progress.new('Further processing')
        self.progress.update('re-setting unique identifiers to track split/position pairs across samples')

        for entry_id in self.variable_positions_table:
            v = self.variable_positions_table[entry_id]
            v['unique_pos_identifier'] = self.get_unique_pos_identification_number(v['unique_position_id'])
            v['parent'] = self.splits_info_table[v['split_name']]['parent']

        self.progress.end()


    def filter_based_on_scattering_factor(self):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from eachother."""

        if self.min_scatter == 0:
            return

        num_samples = len(self.sample_ids)
        if self.min_scatter > num_samples / 2:
            raise ConfigError, 'Minimum scatter (%d) can not be more than half of the number of samples\
                                (%d) :/' % (self.min_scatter, num_samples)

        self.run.info('Min scatter', self.min_scatter)

        num_entries_before_filter = len(self.variable_positions_table)

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

        self.progress.update('removing %s entries from table ...' % pp(len(entry_ids_to_remove)))
        for entry_id in entry_ids_to_remove:
            self.variable_positions_table.pop(entry_id)

        num_entries_after_filter = len(self.variable_positions_table)

        self.progress.end()

        self.run.info('Remaining entries after "minimum scatter" filter',
                      '%s (filter removed %s entries)' % (pp(num_entries_after_filter),
                                                          pp(num_entries_before_filter - num_entries_after_filter)),
                      mc = 'green')

        self.check_if_variable_table_is_empty()


    def filter_based_on_minimum_coverage_in_each_sample(self):
        """To remove any unique entry from the variable positions table that describes a variable position
           and yet is not helpful to distinguish samples from eachother."""

        if self.min_coverage_in_each_sample < 1:
            return

        self.run.info('Min coverage in all samples', '%dX' % self.min_coverage_in_each_sample)

        num_entries_before_filter = len(self.variable_positions_table)

        # we need to make sure we have an up-to-date dictionary for unque position to entry id conversion:
        self.gen_unique_pos_identifier_to_entry_id_dict()

        self.progress.new('Examining coverage of each variable position in each sample')
        self.progress.update('...')

        entry_ids_to_remove = set([])

        for unique_pos_id in self.unique_pos_id_to_entry_id:
            entry_ids = self.unique_pos_id_to_entry_id[unique_pos_id]

            min_coverage_in_a_sample = min([self.variable_positions_table[entry_id]['coverage'] for entry_id in entry_ids])

            if min_coverage_in_a_sample < self.min_coverage_in_each_sample:
                entry_ids_to_remove.update(entry_ids)

        self.progress.update('removing %s entries from table ...' % pp(len(entry_ids_to_remove)))
        for entry_id in entry_ids_to_remove:
            self.variable_positions_table.pop(entry_id)

        num_entries_after_filter = len(self.variable_positions_table)

        self.progress.end()

        self.run.info('Remaining entries after "minimum cov in all samples" filter',
                      '%s (filter removed %s entries)' % (pp(num_entries_after_filter),
                                                          pp(num_entries_before_filter - num_entries_after_filter)),
                      mc = 'green')

        self.check_if_variable_table_is_empty()


    def filter_based_on_num_positions_from_each_split(self):
        self.run.info('Num positions to keep from each split (-n)', self.num_positions_from_each_split)

        num_entries_before_filter = len(self.variable_positions_table)

        self.progress.new('Filtering based on -n')

        self.progress.update('Generating splits and positions tuples ...')
        splits_and_positions = set([(v['split_name'], v['unique_pos_identifier']) for v in self.variable_positions_table.values()])
        unique_positions_to_remove = set([])

        self.progress.update('Generating positions in splits dictionary ...')
        positions_in_splits = {}
        for split_name, position in splits_and_positions:
            if not positions_in_splits.has_key(split_name):
                positions_in_splits[split_name] = set([])

            positions_in_splits[split_name].add(position)

        self.progress.update('Randomly subsampling from splits with num positions > %d' % self.num_positions_from_each_split)
        for split_name in positions_in_splits:
            if self.num_positions_from_each_split and len(positions_in_splits[split_name]) > self.num_positions_from_each_split:
                positions_to_keep = set(random.sample(positions_in_splits[split_name], self.num_positions_from_each_split))
                for pos in (positions_in_splits[split_name] - positions_to_keep):
                    unique_positions_to_remove.add(pos)

        self.progress.update('Identifying entry ids to remove ...')
        entry_ids_to_remove = set([entry_id for entry_id in self.variable_positions_table if self.variable_positions_table[entry_id]['unique_pos_identifier'] in unique_positions_to_remove])

        self.progress.update('Removing %d positions ...' % len(unique_positions_to_remove))
        for entry_id in entry_ids_to_remove:
            self.variable_positions_table.pop(entry_id)

        num_entries_after_filter = len(self.variable_positions_table)

        self.progress.end()

        self.run.info('Remaining entries after "-n" filter',
                      '%s (filter removed %s entries)' % (pp(num_entries_after_filter),
                                                          pp(num_entries_before_filter - num_entries_after_filter)),
                      mc = 'green')

        self.check_if_variable_table_is_empty()


    def recover_base_frequencies_for_all_samples(self):
        """this function populates variable_positions_table dict with entries from samples that have no
           variation at nucleotide positions reported in the table"""

        self.progress.new('Recovering base frequencies for all')

        samples_wanted = self.samples_of_interest if self.samples_of_interest else self.sample_ids
        splits_wanted = self.splits_of_interest if self.splits_of_interest else set(self.splits_info_table.keys())
        next_available_entry_id = max(self.variable_positions_table.keys()) + 1

        self.progress.update('creating a dicts to track missing base frequencies for each sample / split / pos')
        split_pos_to_unique_pos_identifier = {}
        splits_to_consider = {}
        for split_name in splits_wanted:
            splits_to_consider[split_name] = {}
            split_pos_to_unique_pos_identifier[split_name] = {}

        self.progress.update('populating the dict to track missing base frequencies for each sample / split / pos')
        for entry_id in self.variable_positions_table:
            v = self.variable_positions_table[entry_id]
            p = v['pos']
            d = splits_to_consider[v['split_name']]
            u = split_pos_to_unique_pos_identifier[v['split_name']]

            if d.has_key(p):
                d[p].remove(v['sample_id'])
            else:
                d[p] = copy.deepcopy(samples_wanted)
                d[p].remove(v['sample_id'])

            if not u.has_key(p):
                u[p] = v['unique_pos_identifier']

        counter = 0
        for split in splits_to_consider:
            counter += 1
            self.progress.update('accessing split coverages and updating variable positions dict :: %s' % pp(counter))

            split_coverage_across_samples = self.merged_split_coverage_values.get(split)

            split_info = self.splits_info_table[split]

            for pos in splits_to_consider[split]:
                parent_seq = self.contig_sequences[split_info['parent']]['sequence']
                pos_in_contig = split_info['start'] + pos
                base_at_pos = parent_seq[pos_in_contig]
                for sample in splits_to_consider[split][pos]:
                    self.variable_positions_table[next_available_entry_id] = {'parent': split_info['parent'],
                                                                              'departure_from_consensus': 0,
                                                                              'consensus': base_at_pos,
                                                                              'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0,
                                                                              'pos': pos,
                                                                              'pos_in_contig': pos_in_contig, 
                                                                              'coverage': split_coverage_across_samples[sample][pos],
                                                                              'sample_id': sample,
                                                                              'competing_nts': base_at_pos + base_at_pos,
                                                                              'unique_pos_identifier': split_pos_to_unique_pos_identifier[split][pos],
                                                                              'unique_position_id': '%s_%d' % (split, pos),
                                                                              'split_name': split}
                    self.variable_positions_table[next_available_entry_id][base_at_pos] = split_coverage_across_samples[sample][pos]
                    next_available_entry_id += 1

        self.progress.end()

    def report(self):
        self.progress.new('Reporting')

        new_structure = [t.variable_positions_table_structure[0]] + ['unique_pos_identifier'] + t.variable_positions_table_structure[1:] + ['parent']

        self.progress.update('exporting variable positions table as a TAB-delimited file ...')

        utils.store_dict_as_TAB_delimited_file(self.variable_positions_table, self.args.output_file, new_structure)
        self.progress.end()

        self.run.info('Num entries reported', pp(len(self.variable_positions_table)))
        self.run.info('Output File', self.args.output_file) 
        self.run.info('Num nt positions reported', pp(len(set([e['unique_pos_identifier'] for e in self.variable_positions_table.values()]))))


class VariabilityNetwork:
    def __init__(self, args = {}, p=progress, r=run):
        self.args = args

        self.run = r
        self.progress = p

        self.samples = None
        self.samples_information_dict = None
        self.variable_positions_table = None

        A = lambda x, t: t(args.__dict__[x]) if args.__dict__.has_key(x) else None
        null = lambda x: x
        self.input_file_path = A('input_file', null)
        self.samples_information_path = A('samples_information', null)
        self.max_num_unique_positions = A('max_num_unique_positions', int)
        self.output_file_path = A('output_file', null)

        filesnpaths.is_output_file_writable(self.output_file_path)

        if self.samples_information_path:
            filesnpaths.is_file_tab_delimited(self.samples_information_path)
            self.samples_information_dict = utils.get_TAB_delimited_file_as_dictionary(self.samples_information_path)
            num_attributes = len(self.samples_information_dict.values()[0])

            self.run.info('samples_information', '%d attributes read for %d samples' % (num_attributes, len(self.samples_information_dict)))

        if self.input_file_path:
            filesnpaths.is_file_tab_delimited(self.input_file_path)
            self.progress.new('Reading the input file')
            self.progress.update('...')
            self.variable_positions_table = utils.get_TAB_delimited_file_as_dictionary(self.input_file_path)
            self.progress.end()

            self.run.info('input_file', '%d entries read' % len(self.variable_positions_table))


    def generate(self):
        if not self.variable_positions_table:
            raise ConfigError, "There is nothing to report. Either the input file you provided was empty, or you\
                                haven't filled in the variable positions data into the class."

        if self.max_num_unique_positions < 0:
            raise ConfigError, "Max number of unique positions cannot be less than 0.. Obviously :/"

        self.samples = sorted(list(set([e['sample_id'] for e in self.variable_positions_table.values()])))
        self.run.info('samples', '%d found: %s.' % (len(self.samples), ', '.join(self.samples)))

        if self.samples_information_dict:
            samples_missing_in_information_dict = [s for s in self.samples if s not in self.samples_information_dict]
            if len(samples_missing_in_information_dict):
                raise ConfigError, "The sample names you provided in the samples information data is not a subset of\
                                    sample names found in the variable positions data :/ Essentially, every sample name\
                                    appears in the variability data must be present in the samples information data,\
                                    however, you are missing these ones from your samples information: %s."\
                                                % (', '.join(samples_missing_in_information_dict))

        self.unique_variable_positions = set([e['unique_pos_identifier'] for e in self.variable_positions_table.values()])
        self.run.info('unique_variable_positions', '%d found.' % (len(self.unique_variable_positions)))

        if self.max_num_unique_positions and len(self.unique_variable_positions) > self.max_num_unique_positions:
            self.unique_variable_positions = set(random.sample(self.unique_variable_positions, self.max_num_unique_positions))
            self.run.info('unique_variable_positions', 'Unique positions are subsampled to %d' % self.max_num_unique_positions, mc = 'red')

        self.progress.new('Samples dict')
        self.progress.update('Creating an empty one ...')
        samples_dict = {}
        for sample_name in self.samples:
            samples_dict[sample_name] = {}
            for unique_variable_position in self.unique_variable_positions:
                samples_dict[sample_name][unique_variable_position] = 0

        self.progress.update('Updating the dictionary with data')
        for entry in self.variable_positions_table.values():
            sample_id = entry['sample_id']
            pos = entry['unique_pos_identifier']
            frequency = entry['departure_from_consensus']

            samples_dict[sample_id][pos] = float(frequency)


        self.progress.update('Generating the network file')
        utils.gen_gexf_network_file(sorted(list(self.unique_variable_positions)), samples_dict, self.output_file_path, sample_mapping_dict = self.samples_information_dict)
        self.progress.end()

        self.run.info('network_description', self.output_file_path)
