# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to characterize Florian's inversions"""

import os
import copy
import argparse
import numpy as np
from collections import OrderedDict

import anvio
import anvio.tables as t
import anvio.dbinfo as dbi
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.bamops as bamops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
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

        # the primary data structure that will be filled by this
        # class once everything is initialized and the user calls
        # the member function `self.process`:
        self.inversions = {}

        # the purpose of this is to keep track of every stretch where
        # REV/REV or FWD/FWD reads indicated some activity. this is
        # regardless of whether we found a true inversion or not
        self.stretches_considered = {}

        # the purpose of this is to report all the consensus inversions
        self.consensus_inversions = []

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bams_and_profiles_file_path = A('bams_and_profiles')
        self.output_directory = A('output_dir') or 'INVERSIONS-OUTPUT'

        if not self.bams_and_profiles_file_path:
            raise ConfigError("Sorry, you can't get an instance of this class without a `--bams-and-profiles` argument.")

        # get these filled in immediately
        self.contigs_db_path, self.profile_db_bam_file_pairs = utils.get_bams_and_profiles_txt_as_data(self.bams_and_profiles_file_path)
        self.profile_db_paths = [e['profile_db_path'] for e in self.profile_db_bam_file_pairs.values()]

        # params to identify regions of interest. if you are studying the code, don't forget to read
        # the information stored in the help menu of the program about these parameters
        self.min_coverage_to_define_stretches = A('min_coverage_to_define_stretches') or 10
        self.min_stretch_length = A('min_stretch_length') or 50
        self.min_distance_between_independent_stretches = A('min_distance_between_independent_stretches') or 2000
        self.num_nts_to_pad_a_stretch = A('num_nts_to_pad-a_stretch') or 100

        # palindrome search parameters
        self.palindrome_search_algorithm = A('palindrome_search_algorithm')
        self.min_palindrome_length = A('min_palindrome_length') or 10
        self.max_num_mismatches = A('max_num_mismatches') or 0
        self.min_distance_palindrome = A('min-distance') or 50
        self.min_mismatch_distance_to_first_base = A('min_mismatch_distance_to_first_base') or 1

        # parameters to survey inversions
        self.process_only_inverted_reads = A('process_only_inverted_reads')
        self.check_all_palindromes = A('check_all_palindromes')

        # be talkative or not
        self.verbose = A('verbose')

        # debugging mode:
        self.only_report_from = A('only_report_from')

        if self.only_report_from:
            self.verbose = True

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

        Once this function is done, it will populate the `self.inversions` dictionary
        with all the inversions.
        """

        self.progress.new(f"Processing '{entry_name}'")

        self.inversions[entry_name] = []

        ################################################################################
        self.progress.update("Recovering the coverage data")
        ################################################################################

        profile_db = dbops.ProfileSuperclass(argparse.Namespace(profile_db=profile_db_path, contigs_db=self.contigs_db_path), r=run_quiet, p=progress_quiet)
        sample_id = profile_db.p_meta['sample_id']

        # here we open our bam file with an inversions fetch filter.
        # we will access to it later when it is time to get the FWD/FWD and
        # REV/REV reads.
        bam_file = bamops.BAMFileObject(bam_file_path, 'rb')

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
            # of 'high coverage' regions (as in coverage > `self.min_coverage_to_define_stretches`), and store that
            # information into the dictionary `coverage_stretches_in_contigs`
            coverage_stretches_in_contigs[contig_name] = []

            # we also know the contig length here, so let's keep that in mind:
            contig_length = len(contig_coverage)

            # to find regions of high coverage, we first need to 'pad' our array to ensure it always
            # starts and ends with 'low coverage'.
            regions_of_contig_covered_enough = np.hstack([[False], contig_coverage >= self.min_coverage_to_define_stretches, [False]])

            regions_of_contig_covered_enough_diff = np.diff(regions_of_contig_covered_enough.astype(int))
            cov_stretch_start_positions = np.where(regions_of_contig_covered_enough_diff == 1)[0]
            cov_stretch_end_positions = np.where(regions_of_contig_covered_enough_diff == -1)[0]

            # at this stage, `cov_stretch_start_positions` and `cov_stretch_end_positions` contain pairs of
            # positions that match to the begining and end of stretches. we will remove those that are too
            # short to be considered, and store the start/end positions for the remaining stretches of
            # high coverage into the dictionary `coverage_stretches_in_contigs`
            for i in range(0, len(cov_stretch_start_positions)):
                cov_stretch_start, cov_stretch_end = cov_stretch_start_positions[i], cov_stretch_end_positions[i]

                if (cov_stretch_end - cov_stretch_start) >= self.min_stretch_length:
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
                                                                               min_distance_between_independent_stretches=self.min_distance_between_independent_stretches)
            # extend start and stop positions of merged stretches to ENSURE we are not
            # missing important information because bioinformatics.
            coverage_stretches_in_contigs[contig_name] = [(0 if (e[0] - self.num_nts_to_pad_a_stretch< 0) else e[0] - self.num_nts_to_pad_a_stretch,
                                                           contig_length if (e[1] + self.num_nts_to_pad_a_stretch) > contig_length else e[1] + self.num_nts_to_pad_a_stretch) \
                                                                for e in coverage_stretches_in_contigs[contig_name]]

            contig_coverages[contig_name] = contig_coverage

        ################################################################################
        self.progress.update("Getting ready to process stretches")
        ################################################################################
        # time to go through each stretch and look for palindromes
        # first, we will set up the Palindromes class
        _args = argparse.Namespace(min_palindrome_length=self.min_palindrome_length,
                                   max_num_mismatches=self.max_num_mismatches,
                                   min_distance=self.min_distance_palindrome,
                                   palindrome_search_algorithm=self.palindrome_search_algorithm,
                                   min_mismatch_distance_to_first_base=self.min_mismatch_distance_to_first_base)

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

                # if the user wants to learn about only a single sequence, we only
                # focus on that one and prematurely go to the next stretch unless
                # there is a match
                if self.only_report_from and sequence_name != self.only_report_from:
                    continue

                # before we go any further, let's print out the sequence in consideration
                # for the user if they used `--verbose`
                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.warning(None, header=f"Palindromes in {sequence_name}", lc='yellow', nl_before=3)
                    self.run.info_single(f"Sequence {stretch_sequence}", cut_after=0)
                    self.run.info_single(f"Coverage in {sample_id}:", nl_before=1, nl_after=1)
                    self.plot_coverage(f"{sequence_name}", stretch_sequence_coverage)

                # make a record of the stretch that is about to be considered for having
                # palindromes and later true inversions
                self.stretches_considered[f"{entry_name}_{sequence_name}"] = {'sequence_name': sequence_name,
                                                                              'sample_name': entry_name,
                                                                              'contig_name': contig_name,
                                                                              'start_stop': f"{start}-{stop}",
                                                                              'max_coverage': int(max(stretch_sequence_coverage)),
                                                                              'num_palindromes_found': 0,
                                                                              'true_inversions_found': False}

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
                    num_palindromes_found = len(P.palindromes[sequence_name])
                    self.stretches_considered[f"{entry_name}_{sequence_name}"]['num_palindromes_found'] = num_palindromes_found

                    if anvio.DEBUG or self.verbose:
                        self.progress.reset()
                        self.run.info_single(f"The sequence has {PL('palindrome', num_palindromes_found)}:", mc="green")

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

                # here we have, for a given `contig_name`, `start` and `stop` positions of a stretch in it,
                # we have our inversion candidates,
                ################################################################################
                self.progress.update(f"{contig_name}[{start}:{stop}]: true inv testing w/constructs")
                ################################################################################
                true_inversions = self.get_true_inversions_in_stretch(inversion_candidates, bam_file, contig_name, start, stop)

                if anvio.DEBUG or self.verbose:
                    if true_inversions:
                        self.progress.reset()
                        self.run.info_single(f"Of the {PL('inversion candidate', len(inversion_candidates))} above, "
                                             f"anvi'o found the following inversion(s) to have at least one perfect match "
                                             f"to their constructs in the BAM file:", mc="green", nl_before=1)

                        for true_inversion in true_inversions:
                            true_inversion.display()
                    else:
                        self.progress.reset()
                        self.run.info_single(f"No true inversions in this one: none of the REV/REV or FWD/FWD reads "
                                             f"had any of the constructs in {PL('inversion candidate', len(inversion_candidates))}.",
                                             mc="red", nl_before=1)

                # if we are here, it means we took care of one region.
                # it is time to update our global data dictionaries with
                # everything we know about this inversion, and the stretch
                # it belongs:
                if len(true_inversions):
                    for true_inversion in true_inversions:
                        self.stretches_considered[f"{entry_name}_{sequence_name}"]['true_inversions_found'] = True

                        d = OrderedDict({'entry_id': true_inversion.sequence_name,
                                         'sample_name': entry_name,
                                         'contig_name': contig_name,
                                         'first_seq': true_inversion.first_sequence,
                                         'midline': true_inversion.midline,
                                         'second_seq': utils.rev_comp(true_inversion.second_sequence),
                                         'first_start': true_inversion.first_start + start,
                                         'first_end': true_inversion.first_end + start,
                                         'second_start': true_inversion.second_start + start,
                                         'second_end': true_inversion.second_end + start,
                                         'num_mismatches': true_inversion.num_mismatches,
                                         'num_gaps': true_inversion.num_gaps,
                                         'length': true_inversion.length,
                                         'distance': true_inversion.distance})

                        self.inversions[entry_name].append(d)

        self.progress.end()


    def get_true_inversions_in_stretch(self, inversion_candidates, bam_file, contig_name, start, stop):
        """Survey a bunch of palindromes with 'constructs' to find true/active inversions.
    
        A true inversion is one with evidence from short reads that it has activity.
        """
    
        true_inversions_in_stretch = []
    
        if self.process_only_inverted_reads:
            bam_file.fetch_filter = 'inversions'
        else:
            bam_file.fetch_filter = None

        if anvio.DEBUG or self.verbose:
            self.progress.reset()
            self.run.info_single(f"Testing if any of the palindromes resolve to true inversions:", mc="green", nl_before=1, nl_after=1)

        # fetching reads first
        reads = [r.query_sequence for r in bam_file.fetch_only(contig_name, start=start, end=stop)]

        total_num_inversions = len(inversion_candidates)
        current_inversion = 0
        for inversion_candidate in inversion_candidates:
            current_inversion += 1
            num_reads_considered = 0
            match = False
            evidence = ''
            for read in reads:
                num_reads_considered += 1
                if inversion_candidate.v1_left in read:
                    match = True
                    evidence = 'v1_left'
                    break
                elif inversion_candidate.v1_right in read:
                    match = True
                    evidence = 'v1_right'
                    break
                elif inversion_candidate.v2_left in read:
                    evidence = 'v2_left'
                    match = True
                    break
                elif inversion_candidate.v2_right in read:
                    evidence = 'v2_right'
                    match = True
                    break

            if match:
                # we found an inversion candidate that has at least one confirmed
                # construct. We add this one into the list of true inversions:
                true_inversions_in_stretch.append(inversion_candidate)
    
                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.info_single(f"👍 Candidate {current_inversion} of {total_num_inversions}: confirmed by {evidence} "
                                         f"after {num_reads_considered} reads.", mc="yellow", level=2)

                if not self.check_all_palindromes:
                    # if the user is not interested in testing of every single palindrome
                    # found in the stretch of interest to see whether there may be more
                    # inversion candidates, we return the current list which includes only
                    # one confirmed inversion
                    return true_inversions_in_stretch
            else:
                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.info_single(f"👎 Candidate {current_inversion} of {total_num_inversions}: no confirmation "
                                         f"after processing {num_reads_considered} reads.", mc="red", level=2)
    
        # we can simply go with `return true_inversions_in_stretch` here, but the rest of the
        # lines are for posterity:
    
        num_true_inversions_in_stretch = len(true_inversions_in_stretch)
    
        if self.check_all_palindromes:
            if num_true_inversions_in_stretch == 1:
                # if we are here, we have gone through each read and inversion candidate
                # identified in this strecth, and our `true_inversions_in_stretch` contains
                # only one inversion. which means this loop could have taken a fraction of
                # what it did to search for additional inversions when there was none.
                return true_inversions_in_stretch
            elif num_true_inversions_in_stretch > 1:
                # if we are here, it means we actually did find more than one active inversion
                # in a single stretch, and `self.check_all_palindromes` worked well to discover
                # more
                return true_inversions_in_stretch
            else:
                # if we are here, it means we did everything we can and there was nothing
                return []
        else:
            # the user elected not to check all palindromes in the stretch, and if we are here
            # there was no match for any of the constructs for any of the inversoin candidates
            # for this strecth (because if there was one, it would have been already returned)
            return []


    def compute_consensus_inversions(self):
        """Compute a final, consensus list of unique inversions.

        By identifying redundant inversions and reporting only one
        of them, this function reports a final list of inversions
        identifications of which are informed by invdividual samples
        and their coverages in them, but will continue their lives
        as strong and independent inversions.
        """

        # to do this, we need to get all the start positions for all
        # inversions across all samples, cluster them if their start
        # positions are too close to one another to be an different
        # inversion, and choose a single representative for each
        # cluster. first, we get a simpler version of all entries:
        all_entries = []

        self.progress.new('Computing consensus inversions')
        self.progress.update('...')

        for sample_name in self.inversions:
            for inv in self.inversions[sample_name]:
                all_entries.append((sample_name, inv['contig_name'], inv['first_start'], inv['length']), )

        # now it is time to identify clusters. the following state
        # machine does that:
        self.progress.update('Running the state machine')
        clusters = []
        while 1:
            if not len(all_entries):
                break

            entry = all_entries.pop(0)
            cluster = [entry]
            sample, contig_name, start, length = entry
            matching_entries = []

            for i in range(0, len(all_entries)):
                n_sample, n_contig_name, n_start, n_length = all_entries[i]
                if n_contig_name == contig_name and n_start > start - 5 and n_start < start + 5:
                    matching_entries.append(i)

            # add all matching entries
            for i in sorted(matching_entries, reverse=True):
                cluster.append(all_entries.pop(i))

            clusters.append(cluster)

        # now we know our clusters. time to collect all the entires
        # that match to them to have a final list of consensus
        # inversions that occur in at least one sample.
        self.progress.update('Selecting representatives')
        for cluster in clusters:
            num_samples = len(cluster)
            sample_names = ','.join(sorted([x[0] for x in cluster]))
            sample, contig_name, start, length = sorted(cluster, key=lambda x: x[3], reverse=True)[0]

            consensus_found = False
            for sample_name in self.inversions:
                if consensus_found:
                    break

                for inv in self.inversions[sample_name]:
                    if inv['contig_name'] == contig_name and inv['first_start'] == start:
                        consensus_inversion = copy.deepcopy(inv)
                        consensus_inversion['num_samples'] = num_samples
                        consensus_inversion['sample_names'] = sample_names
                        self.consensus_inversions.append(consensus_inversion)
                        consensus_found = True
                        break

        self.progress.end()


    def plot_coverage(self, sequence_name, coverage, num_bins=100):
        try:
            import plotext as plt
        except:
            self.run.warning("You don't have the `plotext` library to plot data :/ You can "
                             "install it by running `pip install plotext` in your anvi'o "
                             "environment.", header="NO PLOT FOR YOU :(")
            return

        try:
            plt.clp()
            plt.title(f"{sequence_name}")
            plt.xlabel("Position")
            plt.ylabel("Coverage")
            plt.plot(coverage, fillx = True)
            plt.plotsize(self.progress.terminal_width, 25)
            plt.canvas_color("cloud")
            plt.axes_color("cloud")
            plt.ticks_color("iron")
        except Exception as e:
            self.run.warning(f"Something bad happen when anvi'o atempted to plot the coverage data :/ "
                             f"Your program will continue running, but here is the error message from "
                             f"the library: \"{e}\".", nl_after=0)
            return

        try:
            plt.show()
        except OSError:
            self.run.warning("Redirecting things into files and working with funny TTYs confuse "
                             "the plotting services. Is ok tho.", header="NO PLOT FOR YOU :(")
            return


    def process(self):
        """Do the processing of everything"""

        w = self.run.width
        self.run.width += 15
        self.run.info("[General] Input BAMs and profiles file", os.path.abspath(self.bams_and_profiles_file_path))
        self.run.info("[General] Number of samples to process", len(self.profile_db_bam_file_pairs))
        self.run.info("[General] Be talkative (--verbose)?", "True" if self.verbose else "False", nl_after=1)

        self.run.info("[Defining stretches] Min FF/RR coverage to qualify", self.min_coverage_to_define_stretches)
        self.run.info("[Defining stretches] Min length", self.min_stretch_length)
        self.run.info("[Defining stretches] Min dist between independent stretches", self.min_distance_between_independent_stretches)
        self.run.info("[Defining stretches] Num nts to pad a stretch", self.num_nts_to_pad_a_stretch, nl_after=1)

        self.run.info('[Finding palindromes] Algorithm', self.palindrome_search_algorithm or "[will be dynamically determined based on sequence length]", mc="red")
        self.run.info("[Finding palindromes] Min palindrome length", self.min_palindrome_length)
        self.run.info("[Finding palindromes] Max num mismatches", self.max_num_mismatches)
        self.run.info("[Finding palindromes] Min mismatch distance to first base", self.min_mismatch_distance_to_first_base)
        self.run.info("[Finding palindromes] Min distance betwen each sequence", self.min_distance_palindrome, nl_after=1)

        self.run.info("[Confirming inversions] Check all palindromes in a stretch?",  "True" if self.check_all_palindromes else "False")
        self.run.info("[Confirming inversions] Process only inverted reads?",  "True" if self.process_only_inverted_reads else "False", nl_after=1)

        if self.only_report_from:
            self.run.info("[Debug] Anvi'o will only report data for:",  self.only_report_from, mc="red", nl_after=1)

        self.run.width = w

        for entry_name in self.profile_db_bam_file_pairs:
            bam_file_path = self.profile_db_bam_file_pairs[entry_name]['bam_file_path']
            profile_db_path = self.profile_db_bam_file_pairs[entry_name]['profile_db_path']

            # populate `self.inversions` with inversions associated with `entry_name`
            self.process_db(entry_name, profile_db_path, bam_file_path)

        # here we know every single inversion. The same inversion site might be
        # found as true inversion in multiple samples when bams and profiles file
        # includes more than one sample. To avoid redundancy, we wish to report
        # a concensus file that describes only unique inversions.
        self.compute_consensus_inversions()


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

        filesnpaths.check_output_directory(self.output_directory)
        filesnpaths.gen_output_directory(self.output_directory)


    def report(self):
        """Reporting per-sample as well as consensus inversions, along with other reporting files"""

        headers = []
        # just learning headers here. NO JUDGING.
        for entry_name in self.inversions:
            if len(self.inversions[entry_name]):
                headers = list(self.inversions[entry_name][0].keys())
                break

        self.run.warning(None, header="REPORTING OUTPUTS", lc="green")

        # report inversions per sample
        for entry_name in self.inversions:
            if len(self.inversions[entry_name]):
                output_path = os.path.join(self.output_directory, f'INVERSIONS-IN-{entry_name}.txt')
                with open(output_path, 'w') as output:
                    output.write('\t'.join(headers) + '\n')
                    for v in self.inversions[entry_name]:
                        output.write('\t'.join([f"{v[k]}" for k in headers]) + '\n')
                self.run.info(f'Inversions in {entry_name}', os.path.abspath(output_path), mc='green')
            else:
                self.run.info(f'Inversions in {entry_name}', 'No true inversions in this one :/', mc='red')

        # report consensus inversions
        output_path = os.path.join(self.output_directory, 'INVERSIONS-CONSENSUS.txt')
        headers = ['contig_name', 'first_seq', 'midline', 'second_seq', 'first_start', 'first_end', 'second_start', 'second_end', 'num_mismatches', 'num_gaps', 'length', 'distance', 'num_samples', 'sample_names']
        with open(output_path, 'w') as output:
            output.write('\t'.join(headers) + '\n')
            for v in self.consensus_inversions:
                output.write('\t'.join([f"{v[k]}" for k in headers]) + '\n')
        self.run.info('Reporting file for consensus inversions', os.path.abspath(output_path), mc='green')

        # report all stretches
        output_path = os.path.join(self.output_directory, 'ALL-STRETCHES-CONSIDERED.txt')
        headers = ['entry_id', 'sequence_name', 'sample_name', 'contig_name', 'start_stop', 'max_coverage', 'num_palindromes_found', 'true_inversions_found']
        utils.store_dict_as_TAB_delimited_file(self.stretches_considered, output_path, headers=headers)
        self.run.info('Reporting file on all stretches considered', os.path.abspath(output_path), nl_before=1, nl_after=1)
