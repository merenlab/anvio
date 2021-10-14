# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to characterize Florian's inversions"""

import argparse
import numpy as np

import anvio
import anvio.tables as t
import anvio.dbinfo as dbi
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.bamops as bamops
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
PL = terminal.pluralize
run_quiet = terminal.Run(verbose=False)
progress_quiet = terminal.Progress(verbose=False)


class Inversions:
    def __init__(self, args, skip_sanity_check=False, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bams_and_profiles_file_path = A('bams_and_profiles')

        if not self.bams_and_profiles_file_path:
            raise ConfigError("Sorry, you can't get an instance of this class without a `--bams-and-profiles` arguemnt.")

        # get these filled in immediately
        self.contigs_db_path, self.profile_db_bam_file_pairs = utils.get_bams_and_profiles_txt_as_data(self.bams_and_profiles_file_path)
        self.profile_db_paths = [e['profile_db_path'] for e in self.profile_db_bam_file_pairs.values()]

        # palindrome search parameters
        self.min_palindrome_length = A('min_palindrome_length') or 10
        self.max_num_mismatches = A('max_num_mismatches') or 0
        self.min_distance_palindrome = A('min-distance') or 50

        # be talkative or not
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


    def process_db(self, entry_name, profile_db_path, bam_file_path):
        """Function that does everything.

        `entry_name` is the entry name in bams and profiles file.
        """

        self.progress.new(f"Processing '{entry_name}'")

        ################################################################################
        self.progress.update("Recovering the coverage data")
        ################################################################################

        profile_db = dbops.ProfileSuperclass(argparse.Namespace(profile_db=profile_db_path, contigs_db=self.contigs_db_path), r=run_quiet, p=progress_quiet)
        sample_id = profile_db.p_meta['sample_id']

        # here we open our bam file with an inversions fetch filter.
        # we will access to it later when it is time to get the FWD/FWD and
        # REV/REV reads.
        bam_file = bamops.BAMFileObject(bam_file_path, 'rb')
        bam_file.fetch_filter = 'inversions'

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
        contig_coverages = {}
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

            contig_coverages[contig_name] = contig_coverage

        ################################################################################
        self.progress.update("Getting ready to process stretches")
        ################################################################################
        # time to go through each stretch and look for palindromes
        # first, we will set up the Palindromes class
        _args = argparse.Namespace(min_palindrome_length=self.min_palindrome_length,
                                   max_num_mismatches=self.max_num_mismatches,
                                   min_distance=self.min_distance_palindrome)

        P = Palindromes(_args,
                        run=run_quiet,
                        progress=progress_quiet)
        P.verbose = False

        # now we can go through all the stretches to look for palindromes. this is a LOOOOOONG loop.
        # down below, we will got through each contig name, find stretches of good coverage of FWD/FWD
        # and REV/REV reads (since their coverage values are stored in the profile db of 'inversions'
        # type), find palindromes in those sequences that match to those coverage stretches, build some
        # constructs, and then go through every FWD/FWD and REV/REV read from the BAM file to see if
        # our constructs occur in any of them, which is the only 100% proof of an active inversion.
        for contig_name in coverage_stretches_in_contigs:
            contig_sequence = self.contig_sequences[contig_name]['sequence']
            for start, stop in coverage_stretches_in_contigs[contig_name]:
                stretch_sequence_coverage = contig_coverages[contig_name][start:stop]
                stretch_sequence = contig_sequence[start:stop]
                sequence_name = f"{contig_name}_{start}_{stop}"

                # before we go any further, let's print out the sequence in consideration
                # for the user if they used `--verbose`
                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.warning(None, header=f"Palindromes in {sequence_name}", lc='yellow', nl_before=3)
                    self.run.info_single(f"Sequence {stretch_sequence}", cut_after=0)
                    self.run.info_single("Coverage:", nl_before=1, nl_after=1)
                    self.plot_coverage(f"{sequence_name}", stretch_sequence_coverage)

                ################################################################################
                self.progress.update(f"{contig_name}: looking for palindromes")
                ################################################################################
                P.find(stretch_sequence, sequence_name=sequence_name, display_palindromes=False)

                if not len(P.palindromes[sequence_name]):
                    # there is no palindrome in this one
                    if anvio.DEBUG or self.verbose:
                        self.progress.reset()
                        self.run.info_single("No palindromes in this one :/", mc="red")
                    continue
                else:
                    if anvio.DEBUG or self.verbose:
                        self.progress.reset()
                        self.run.info_single(f"The sequence has {PL('palindrome', len(P.palindromes[sequence_name]))}:", mc="green")

                ################################################################################
                self.progress.update(f"{contig_name}: building constructs")
                ################################################################################
                # this is important. here we are getting ready to test each our inversion candidate
                # by reconstructing Florian's imaginary sequences. in the next step we will see if
                # any of these sequences are in any of the FWD/FWD or REV/REV reads
                inversion_candidates = []
                for inversion_candidate in P.palindromes[sequence_name]:
                    region_A_start = inversion_candidate.first_start - 6
                    region_A_end = inversion_candidate.first_start
                    region_A = stretch_sequence[region_A_start:region_A_end]

                    region_B_start = inversion_candidate.first_end
                    region_B_end = inversion_candidate.first_end + 6
                    region_B = stretch_sequence[region_B_start:region_B_end]

                    region_C_start = inversion_candidate.second_start - 6
                    region_C_end = inversion_candidate.second_start
                    region_C = stretch_sequence[region_C_start:region_C_end]

                    region_D_start = inversion_candidate.second_end
                    region_D_end = inversion_candidate.second_end + 6
                    region_D = stretch_sequence[region_D_start:region_D_end]

                    construct_v1_left = region_A + inversion_candidate.first_sequence + utils.rev_comp(region_C)
                    construct_v1_right = utils.rev_comp(region_B) + utils.rev_comp(inversion_candidate.second_sequence) + region_D

                    construct_v2_left = region_A + inversion_candidate.second_sequence + utils.rev_comp(region_C)
                    construct_v2_right = utils.rev_comp(region_B) + utils.rev_comp(inversion_candidate.first_sequence) + region_D

                    # update the palindrome instance with its constructs
                    inversion_candidate.v1_left = construct_v1_left
                    inversion_candidate.v1_right = construct_v1_right
                    inversion_candidate.v2_left = construct_v2_left
                    inversion_candidate.v2_right = construct_v2_right

                    if anvio.DEBUG or self.verbose:
                        self.progress.reset()
                        inversion_candidate.display()
                        self.run.info("Construct v1 left", construct_v1_left, mc="cyan")
                        self.run.info("Construct v1 right", construct_v1_right, mc="cyan")
                        self.run.info("Construct v2 left", construct_v2_left, mc="cyan")
                        self.run.info("Construct v2 right", construct_v2_right, mc="cyan")

                    inversion_candidates.append(inversion_candidate)

                # here we have, for a given `contig_name` and `start` and `stop` positions of a stretch in it, we have our inversion candidates,
                ################################################################################
                self.progress.update(f"{contig_name}[{start}:{stop}]: testing constructs")
                ################################################################################
                true_inversions = []
                for read in bam_file.fetch_only(contig_name, start=start, end=stop):
                    hit = False
                    for inversion_candidate in inversion_candidates:
                        if inversion_candidate.v1_left in read.query_sequence:
                            hit = True
                        elif inversion_candidate.v1_right in read.query_sequence:
                            hit = True
                        elif inversion_candidate.v2_left in read.query_sequence:
                            hit = True
                        elif inversion_candidate.v2_right in read.query_sequence:
                            hit = True
                    if hit:
                        true_inversions.append(inversion_candidate)
                        break

                if anvio.DEBUG or self.verbose:
                    if not len(true_inversions):
                        self.progress.reset()
                        self.run.info_single(f"No true inversions in this one: none of the REV/REV or FWD/FWD reads "
                                             f"had any of the constructs in {PL('inversion candidate', len(inversion_candidates))}.",
                                             mc="red", nl_before=1)
                    else:
                        self.progress.reset()
                        self.run.info_single(f"The following {PL('inversion candidate', len(inversion_candidates))} "
                                             f"had one or more perfect matches to their constructs in REV/REV or "
                                             f"FWD/FWD reads from the BAM file:", mc="green", nl_before=1)
                        for true_inversion in true_inversions:
                            true_inversion.display()

        self.progress.end()


    def plot_coverage(self, sequence_name, coverage, num_bins=100):
        try:
            import plotext as plt
        except:
            self.run.warning("You don't have the `plotext` library to plot data :/ You can "
                             "install it by running `pip install plotext` in your anvi'o "
                             "environment.", header="NO PLOT FOR YOU :(")

        plt.clp()
        plt.title(f"{sequence_name}")
        plt.xlabel("Position")
        plt.ylabel("Coverage")
        plt.plot(coverage, fillx = True)
        plt.plotsize(self.progress.terminal_width, 25)
        plt.canvas_color("cloud")
        plt.axes_color("cloud")
        plt.ticks_color("iron")

        try:
            plt.show()
        except OSError:
            self.run.warning("Redirecting things into files and working with funny TTYs confuse "
                             "the plotting services. Is ok tho.", header="NO PLOT FOR YOU :(")
            pass


    def process(self):
        for entry_name in self.profile_db_bam_file_pairs:
            bam_file_path = self.profile_db_bam_file_pairs[entry_name]['bam_file_path']
            profile_db_path = self.profile_db_bam_file_pairs[entry_name]['profile_db_path']
            self.process_db(entry_name, profile_db_path, bam_file_path)


    def sanity_check(self):
        bad_profile_dbs = [p for p in self.profile_db_paths if dbi.DBInfo(self.profile_db_paths[0]).get_self_table()['fetch_filter'] != 'inversions']
        if len(bad_profile_dbs):
            if len(bad_profile_dbs) == len(self.profile_db_paths):
                if len(bad_profile_dbs) == 1:
                    summary = "The only profile database you have here is also the one with the wrong 'fetch filter' :( CONGRATULATIONS."
                else:
                    summary = f"But none of the {PL('database', len(self.profile_db_paths))} here have the right 'fetch filter' :("
            else:
                summary = (f"Of the total of {PL('database', len(self.profile_db_paths))} "
                           f"you are working with, {len(bad_profile_dbs)} {PL('does not', len(bad_profile_dbs), alt='do not')} have "
                           f"the right 'fetch filter' :/")

            raise ConfigError(f"To report inversions, you must provide this program with one or more profile "
                              f"databases that have the fetch filter 'inversions'. {summary} PRO TIP: a given profile database "
                              f"will have a fetch filter 'inversions' only if it was generated by `anvi-profile` with the additional "
                              f"flag `--fetch-filter inversions`. PLUS, you can always use the program `anvi-db-info` to "
                              f"learn about whether a fetch filter was used when anvi'o profiled the data. BYE NOW.")


    def report(self):
        raise ConfigError("Nothing to report yet :(")
