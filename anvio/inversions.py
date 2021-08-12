# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to characterize Florian's inversions"""

import os
import argparse
import numpy as np

import anvio
import anvio.tables as t
import anvio.dbinfo as dbi
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
P = terminal.pluralize
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)


class Inversions:
    def __init__(self, args, skip_sanity_check=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.profile_db_paths = A('profile_dbs')

        if not skip_sanity_check:
            self.sanity_check()

        # we will generate our splits info and contigs to splits dicts here.
        split_names = utils.get_all_item_names_from_the_database(self.profile_db_paths[0])
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.splits_basic_info = contigs_db.db.smart_get(t.splits_info_table_name, column='split', data=split_names)
        contigs_db.disconnect()

        # next, we will generate a dictionary to convert contig names to split names
        self.contig_name_to_split_names = {}
        for split_name in sorted(self.splits_basic_info.keys()):
            contig_name = self.splits_basic_info[split_name]['parent']

            if contig_name not in self.contig_name_to_split_names:
                self.contig_name_to_split_names[contig_name] = []

            self.contig_name_to_split_names[contig_name].append(split_name)

        # let's have a variable of convenience:
        self.contig_names = sorted(list(self.contig_name_to_split_names.keys()))


    def process_db(self, profile_db_path):
        self.progress.append(" recovering coverages")

        profile_db = dbops.ProfileSuperclass(argparse.Namespace(profile_db=profile_db_path, contigs_db=self.contigs_db_path), r=run_quiet, p=progress_quiet)
        sample_id = profile_db.p_meta['sample_id']

        # FIXME: this will need to have more reasonable defaults
        #        that are also parameterized
        min_cov = 10
        min_stretch_length = 50
        min_distance_between_independent_stretches = 2000

        coverage_stretches_in_contigs = {}
        for contig_name in self.contig_names:
            contig_coverage = np.array([])

            split_names = self.contig_name_to_split_names[contig_name]

            for i in range(len(split_names)):
                split_name = split_names[i]
                split_coverages = auxiliarydataops.AuxiliaryDataForSplitCoverages(profile_db.auxiliary_data_path, profile_db.p_meta['contigs_db_hash']).get(split_name)
                contig_coverage = np.concatenate((contig_coverage, split_coverages[sample_id]), axis=None)

            # now we know the `contig_coverage`. it is time to break it into stretches
            # of 'high coverage' regions (as in coverage > `min_cov`), and store that
            # information into the dictionary `coverage_stretches_in_contigs`
            coverage_stretches_in_contigs[contig_name] = []

            # to find regions of high coverage, we first need to 'pad' our array to ensure it always
            # starts and ends with 'low coverage'.
            regions_of_contig_covered_enough = np.hstack([[False], contig_coverage >= min_cov, [False]])

            regions_of_contig_covered_enough_diff = np.diff(regions_of_contig_covered_enough.astype(int))
            cov_stretch_start_positions = np.where(regions_of_contig_covered_enough_diff == 1)[0]
            cov_stretch_end_positions = np.where(regions_of_contig_covered_enough_diff == -1)[0]

            # at this stage, `cov_stretch_start_positions` and `cov_stretch_end_positions` contain pairs of
            # positions that match to the begining and end of stretches. we will remove those that are too
            # short to be considered, and store the start/end positions for the remaining stretches of
            # high coverage into the dictionary `coverage_stretches_in_contigs`
            for i in range(0, len(cov_stretch_start_positions)):
                cov_stretch_start, cov_stretch_end = cov_stretch_start_positions[i], cov_stretch_end_positions[i]

                if (cov_stretch_end - cov_stretch_start) >= min_stretch_length:
                    coverage_stretches_in_contigs[contig_name].append((cov_stretch_start, cov_stretch_end),)

            # now it is time to merge those stretches of coverage if they are close to one another to avoid
            # over-splitting areas of coverage due to short regions with low-coverage in the middle like this,
            # where we wish to identify A and B together in a single stretch:
            #
            #                A         B
            #
            #                -         -
            #               ---        --
            #              -----      -----
            #             --------   --------
            #           -----------------------
            # -----------------------------------------------
            coverage_stretches_in_contigs[contig_name] = utils.merge_stretches(coverage_stretches_in_contigs[contig_name],
                                                                               min_distance_between_independent_stretches=min_distance_between_independent_stretches)


            # DO SOMETHING WITH `coverage_stretches_in_contigs` :)


    def process(self):
        self.progress.new("Processing coverages", progress_total_items=len(self.profile_db_paths))
        self.progress.update('...')

        for profile_db_path in self.profile_db_paths:
            profile_name = dbi.DBInfo(profile_db_path).get_self_table()['sample_id']
            self.progress.update(f"{profile_name} ...", increment=True)
            self.process_db(profile_db_path)

        self.progress.end()


    def sanity_check(self):
        utils.is_contigs_db(self.contigs_db_path)

        if len(set([os.path.abspath(p) for p in self.profile_db_paths])) != len(self.profile_db_paths):
            raise ConfigError("We don't like redundant profile databases. Every database should appear "
                              "only once in the list of databases to be processed. RULES, yo.")

        [utils.is_profile_db_and_contigs_db_compatible(p, self.contigs_db_path) for p in self.profile_db_paths]

        bad_profile_dbs = [p for p in self.profile_db_paths if dbi.DBInfo(p).variant != 'inversions']
        if len(bad_profile_dbs):
            if len(bad_profile_dbs) == len(self.profile_db_paths):
                if len(bad_profile_dbs) == 1:
                    summary = "The only profile database you have here is also the one with the wrong variant :( CONGRATULATIONS."
                else:
                    summary = f"But none of the {P('database', len(self.profile_db_paths))} here have the right variant :("
            else:
                summary = (f"Of the total of {P('database', len(self.profile_db_paths))} "
                           f"you are working with, {len(bad_profile_dbs)} {P('does not', len(bad_profile_dbs), alt='do not')} have "
                           f"the right variant :/")

            raise ConfigError(f"To report inversions, you must provide this program with one or more profile "
                              f"databases that are of variant 'inversions'. {summary} PRO TIP: a given profile database is "
                              f"of variant type 'inversions' only if it was generated by `anvi-profile` with the additional "
                              f"flag `--fetch-filter inversions`. PLUS, you can always use the program `anvi-db-info` to "
                              f"learn about the variant of any anvi'o database. BYE NOW.")


    def report(self):
        raise ConfigError("Nothing to report yet :(")
