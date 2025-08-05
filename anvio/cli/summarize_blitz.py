#!/usr/bin/env python
# -*- coding: utf-8
"""Script to quickly summarize many single profiles."""

import sys
import numpy as np
from anvio.argparse import ArgumentParser

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.terminal as terminal
import anvio.constants as constants

from anvio.errors import ConfigError, FilesNPathsError
from anvio.utils.files import AppendableFile
from anvio.dbinfo import is_blank_profile, is_profile_db_and_contigs_db_compatible, is_profile_db_merged


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ivagljiva']
__requires__ = ['single-profile-db', 'contigs-db']
__provides__ = ['quick-summary']
__description__ = "FAST summary of many anvi'o single profile databases (without having to use the program anvi-merge)."
__resources__ = []


class RapidSummarizer:
    def __init__(self, args, skip_sanity_check=False, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_paths = A('profile_dbs')
        self.collection_name = A('collection_name')
        self.output_file_path = A('output_file') or self.collection_name + '-SUMMARY-BLITZ.txt'
        self.stats_to_summarize = A('stats_to_summarize')

        if self.stats_to_summarize:
            self.stats_to_summarize = self.stats_to_summarize.split(',')
        else:
            self.stats_to_summarize = ["detection", "mean_coverage_Q2Q3"]

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Num profile DBs', len(self.profile_db_paths))
        self.run.info('Collection', self.collection_name)
        self.run.info('Summarizing stats', ", ".join(self.stats_to_summarize))

        if not skip_sanity_check:
            self.sanity_check()


    def sanity_check(self):
        self.progress.new('Sanity checks')
        self.progress.update('Are profile and contigs dbs compatible?')
        [is_profile_db_and_contigs_db_compatible(p, self.contigs_db_path) for p in self.profile_db_paths]
        self.progress.update('Are profile dbs single?')
        if any([is_profile_db_merged(p) for p in self.profile_db_paths]):
            raise ConfigError("At least one of the profile dbs you provided is a merged profile. Unfortunately, "
                              "this program only works for single profiles, so get rid of the merged ones (please and "
                              "thank you, says anvi'o).")
        self.progress.update('Are profile dbs blank?')
        if any([is_blank_profile(p) for p in self.profile_db_paths]):
            raise ConfigError("At least one of the profile dbs you provided is blank, so it has no data for us "
                              "to use. Please make sure to give this program only non-blank profiles.")
        self.progress.end()

        # did user request stats that exist?
        for s in self.stats_to_summarize:
            if s not in constants.essential_data_fields_for_anvio_profiles:
                stat_str = ", ".join(constants.essential_data_fields_for_anvio_profiles)
                raise ConfigError(f"The statistic you requested, {s}, does not exist. Here are the options to "
                                  f"choose from: {stat_str}")


    def process(self):
        self.progress.new('Loading data from contigs db')
        contigs_db = db.DB(self.contigs_db_path, None, ignore_version=True)
        splits_basic_info = contigs_db.get_table_as_dict(t.splits_info_table_name)
        contigs_db.disconnect()

        self.progress.update("Loading collection info from profile db")
        profile_db = db.DB(self.profile_db_paths[0], None, ignore_version=True)
        collection = {s:b for s,b in profile_db.get_some_columns_from_table(t.collections_splits_table_name, 'split,bin_name',
                                                                            where_clause=f"""collection_name = '{self.collection_name}'""")}
        self.progress.end()

        if not collection:
            raise ConfigError(f"Someone did an oopsie :) The collection '{self.collection_name}' does not exist "
                              f"in the profile database at {self.profile_db_paths[0]}.")
        split_list = ','.join(["'%s'" % split_name for split_name in collection.keys()])
        splits_where_clause = f'''item IN ({split_list})'''

        # prepare for appending to output file
        output_file = AppendableFile(self.output_file_path, append_type=dict, fail_if_file_exists=True)
        header_list = ['unique_id', 'bin_name', 'sample'] + self.stats_to_summarize
        key = 0         # arbitrary int to serve as unique index in output

        # obtain statistics from all profile dbs
        self.progress.new("Summarizing statistics from profile dbs", progress_total_items=len(self.profile_db_paths))
        for profile_db_path in self.profile_db_paths:
            profile_db = db.DB(profile_db_path, None, ignore_version=True)
            sample_id = profile_db.get_meta_value('sample_id')
            self.progress.update("[%d of %d] Processing sample %s" % (self.progress.progress_current_item + 1, len(self.profile_db_paths), sample_id))

            split_length_dict = {}
            stats_per_bin = {} # 3-level dict containing weighted avg of statistics per bin per sample, keyed by bin and then by sample and then by statistic name

            for stat in self.stats_to_summarize:
                if stat == 'variability' and not profile_db.get_meta_value('SNVs_profiled'):
                    self.run.warning(f"SNVs are not profiled for sample {sample_id}, hence the variability statistic will not be summarized.")
                    continue

                stat_table_name = stat + "_splits"
                stat_table_splits = {i:v for i,v in profile_db.get_some_columns_from_table(stat_table_name, 'item,value', where_clause=splits_where_clause)}

                # compute statistic per bin, normalized by split length
                all_vals_per_bin = {} # list of stat values for each split in bin
                all_lens_per_bin = {} # list of split lengths (weights) for each split in bin

                for split,stat_val in stat_table_splits.items():
                    bin_for_split = collection[split]
                    if bin_for_split not in all_vals_per_bin:
                        all_vals_per_bin[bin_for_split] = []
                        all_lens_per_bin[bin_for_split] = []
                    all_vals_per_bin[bin_for_split].append(stat_val)

                    if split not in split_length_dict:
                        split_length_dict[split] = splits_basic_info[split]['length']
                    all_lens_per_bin[bin_for_split].append(split_length_dict[split])

                for bin_name in all_vals_per_bin.keys():
                    if bin_name not in stats_per_bin:
                        stats_per_bin[bin_name] = {}
                    if sample_id not in stats_per_bin[bin_name]:
                        stats_per_bin[bin_name][sample_id] = {}
                    stats_per_bin[bin_name][sample_id][stat] = np.average(all_vals_per_bin[bin_name], weights=all_lens_per_bin[bin_name])

            profile_db.disconnect()

            # print stats to output file
            output_dict = {}  # 2-level dictionary keyed by arbitrary integer, each entry is one line in the output file
            for bin_name, sample_dict in stats_per_bin.items():
                for sample, stat_dict in sample_dict.items():
                    output_dict[key] = {}
                    output_dict[key]['bin_name'] = bin_name
                    output_dict[key]['sample'] = sample
                    for stat, val in stat_dict.items():
                        output_dict[key][stat] = val
                    key += 1
            output_file.append(output_dict, key_header="unique_id", headers=header_list)

            self.progress.increment()
            self.progress.reset()

        self.progress.end()
        self.run.info("Output file", self.output_file_path)


def main():
    args = get_args()

    try:
        r = RapidSummarizer(args)
        r.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)


def get_args():
    parser = ArgumentParser(description=__description__)

    parser.add_argument('profile_dbs', metavar = 'SINGLE_PROFILE(S)', nargs='+',
                        help = "Anvo'o single profiles to summarize. All profiles should be associated with the same contigs db.")

    parser.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))
    parser.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name', {'required': True,
                                 'help': "Name of the collection you wish to summarize. The SAME collection will be summarized "
                                         "across all of your input profiles. This collection must be defined in at least the "
                                         "first profile in the argument list."}))
    parser.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    parser.add_argument(*anvio.A('stats-to-summarize'), **anvio.K('stats-to-summarize'))

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
