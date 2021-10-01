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
from anvio.sequencefeatures import Palindromes


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
        self.min_palindrome_length = A('min_palindrome_length') or 10
        self.max_num_mismatches = A('max_num_mismatches') or 0
        self.verbose = A('verbose')

        if not skip_sanity_check:
            self.sanity_check()

        # we will generate our splits info and contigs to splits dicts here.
        split_names = utils.get_all_item_names_from_the_database(self.profile_db_paths[0])
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.splits_basic_info = contigs_db.db.smart_get(t.splits_info_table_name, column='split', data=split_names)
        self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
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
        """Function that does everything"""

        profile_name = dbi.DBInfo(profile_db_path).get_self_table()['sample_id']

        self.progress.new(f"Processing '{profile_name}'")

        ################################################################################
        self.progress.update("Recovering the coverage data")
        ################################################################################

        profile_db = dbops.ProfileSuperclass(argparse.Namespace(profile_db=profile_db_path, contigs_db=self.contigs_db_path), r=run_quiet, p=progress_quiet)
        sample_id = profile_db.p_meta['sample_id']

        # FIXME: this will need to have more reasonable defaults
        #        that are also parameterized
        min_cov = 10
        min_stretch_length = 50
        min_distance_between_independent_stretches = 2000
        num_nts_to_pad_stretches = 100

        ################################################################################
        self.progress.update("Computing coverage stretches")
        ################################################################################
        # populate coverage stretches in contigs based on coverage data in this
        # particular profile_db. we will then go through each stretch to find
        # those that include palindromic sequences
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

            # we also know the contig length here, so let's keep that in mind:
            contig_length = len(contig_coverage)

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
            # extend start and stop positions of merged stretches to ENSURE we are not
            # missing important information because bioinformatics.
            coverage_stretches_in_contigs[contig_name] = [(0 if (e[0] - num_nts_to_pad_stretches < 0) else e[0] - num_nts_to_pad_stretches,
                                                           contig_length if (e[1] + num_nts_to_pad_stretches) > contig_length else e[1] + num_nts_to_pad_stretches) \
                                                                for e in coverage_stretches_in_contigs[contig_name]]

        ################################################################################
        self.progress.update("Identifying palindromes within stretches")
        ################################################################################
        # time to go through each stretch and look for palindromes
        # first, we will set up the Palindromes class
        _args = argparse.Namespace(min_palindrome_length=self.min_palindrome_length, max_num_mismatches=self.max_num_mismatches)
        P = Palindromes(_args,
                        run=run_quiet,
                        progress=progress_quiet)
        P.verbose = False

        # now we can go through all the stretches to look for palindromes
        for contig_name in coverage_stretches_in_contigs:
            contig_sequence = self.contig_sequences[contig_name]['sequence']
            for start, stop in coverage_stretches_in_contigs[contig_name]:
                stretch_sequence = contig_sequence[start:stop]
                sequence_name = f"{contig_name}_{start}_{stop}"

                P.find(stretch_sequence, sequence_name=sequence_name, display_palindromes=False)

                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.warning(None, header=f"Palindromes in {sequence_name}")
                    self.run.info_single(f"Sequence {stretch_sequence}", cut_after=0)
                    [p.display() for p in P.palindromes[sequence_name]]

                if not len(P.palindromes[sequence_name]):
                    continue

                for inversion_candidate in P.palindromes[sequence_name]:
                    pass

        self.progress.end()


    def process(self):
        for profile_db_path in self.profile_db_paths:
            self.process_db(profile_db_path)


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
