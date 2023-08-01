# -*- coding: utf-8
# pylint: disable=line-too-long
"""A module to characterize Florian's inversions"""

import os
import copy
import argparse
import numpy as np
from collections import OrderedDict, Counter

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing



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
from anvio.sequencefeatures import Palindromes, PrimerSearch


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

        # in which we will store the genomic context that surrounds
        # consensus inversions for downstream fun
        self.genomic_context_surrounding_consensus_inversions = {}

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.bams_and_profiles_file_path = A('bams_and_profiles')
        self.output_directory = A('output_dir') or 'INVERSIONS-OUTPUT'

        if not self.bams_and_profiles_file_path:
            raise ConfigError("Sorry, you can't get an instance of this class without a `--bams-and-profiles` argument.")

        # get these filled in immediately
        self.contigs_db_path, self.profile_db_bam_file_pairs = utils.get_bams_and_profiles_txt_as_data(self.bams_and_profiles_file_path)
        self.profile_db_paths = [e['profile_db_path'] for e in self.profile_db_bam_file_pairs.values()]
        self.raw_r1_r2_reads_are_present = all([('r1' in v) and ('r1' in v) for v in self.profile_db_bam_file_pairs.values()])

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

        # the purpose of this variable is to determine how much of the palindrome for a
        # given inversion should be used to build a primer to search for short reads. the
        # longer it is, the more specific the primers will be when surveying ALL FASTQ READS
        # in a (meta)genome sequencing dataset. but if it is too long, then there will
        # not be enough reads to keep depending on the minimum short read length of the
        # sequencing.
        self.oligo_primer_base_length = A('oligo_primer_base_length') or 12

        # this variable tells us how long is the oligonucleotide we should be focusing on
        # to measure proportions of activity
        self.oligo_length = 6

        # compute inversion activity across samples?
        self.skip_compute_inversion_activity = A('skip_compute_inversion_activity') or False
        self.pre_computed_inversions_path = A('pre_computed_inversions')

        # stop inversion activity computation early for testing?
        self.end_primer_search_after_x_hits = A('end_primer_search_after_x_hits')

        # skip learning about the genomic context that surrounds inversions?
        self.skip_recovering_genomic_context = A('skip_recovering_genomic_context')
        self.gene_caller_to_consider_in_context = A('gene_caller') or 'prodigal'
        self.num_genes_to_consider_in_context = A('num_genes_to_consider_in_context') or 3

        # be talkative or not
        self.verbose = A('verbose')

        # performance
        self.num_threads = int(A('num_threads')) if A('num_threads') else 1

        # focus mode:
        self.target_contig = A('target_contig')
        self.target_region_start = A('target_region_start')
        self.target_region_end = A('target_region_end')

        # these are the keys we are interested in finding in input files offered to reconstruct
        # inverson profiles via the --pre-computed-inversions flag. NOTE that these keys are not ALL
        # keys that are used to build inversion profiles in the code, but the minimum set that
        # co-occur both sample-specific and consensus inversion reports. this way, the user can
        # attempt to characterize the activity of inversions found in a single sample if they wish:
        self.essential_keys_to_describe_inversions = [('contig_name', str), ('first_seq', str), ('midline', str), ('second_seq', str),
                                                      ('first_start', int), ('first_end', int), ('first_oligo_primer', str),
                                                      ('first_oligo_reference', str), ('second_start', int), ('second_end', int),
                                                      ('second_oligo_primer', str), ('second_oligo_reference', str), ('num_mismatches', int),
                                                      ('num_gaps', int), ('length', int), ('distance', int)]

        if self.target_contig:
            self.verbose = True

        if not skip_sanity_check:
            self.sanity_check()

        # we will generate our splits info, contigs to splits dicts, and check a few things to learn more about the
        # contigs db.
        split_names = utils.get_all_item_names_from_the_database(self.profile_db_paths[0])
        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)
        self.splits_basic_info = contigs_db.db.smart_get(t.splits_info_table_name, column='split', data=split_names)
        self.contig_sequences = contigs_db.db.get_table_as_dict(t.contig_sequences_table_name)
        self.genes_are_called_in_contigs_db = contigs_db.meta['genes_are_called']
        self.genes_annotated_with_functions_in_contigs_db = contigs_db.meta['gene_function_sources'] is not None and len(contigs_db.meta['gene_function_sources']) > 0
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

        if self.target_contig:
            if self.target_contig not in self.contig_names:
                raise ConfigError(f"You asked anvi'o to focus on a single contig, named {self.target_contig}, "
                                  f"but there doesn't seem to be a contig in this database with that name :/")
            else:
                self.contig_names = [self.target_contig]


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
        auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(profile_db.auxiliary_data_path, profile_db.p_meta['contigs_db_hash'])
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
                split_coverages = auxiliary_db.get(split_name)
                contig_coverage = np.concatenate((contig_coverage, split_coverages[sample_id]), axis=None)

            # if the user is asking us to focus only a particular stretch in the contig
            # we are going to grant their wish by setting the values in `contig_coverage`
            # to zero that are outside of those regions
            if self.target_region_start or self.target_region_end:
                if self.target_region_start:
                    contig_coverage[0:self.target_region_start] = 0.0
                if self.target_region_end:
                    contig_coverage[self.target_region_end:] = 0.0

            # now we know the `contig_coverage`
            contig_coverages[contig_name] = contig_coverage

            # and all we want to do now via the downstream code in this loop is to break it into stretches
            # of 'high coverage' regions of FWD/FWD or REV/REV reads (as in coverage > `self.min_coverage_to_define_stretches`),
            # and store that  information into the dictionary `coverage_stretches_in_contigs` for further processing.
            # so here is the blank entry that will be filled soon:
            coverage_stretches_in_contigs[contig_name] = []

            # but we don't want to go through any of this if the contig has no coverage at all. so here we will test that first.
            # speicial thanks goes to Andrea Watson who identified this edge case in https://github.com/merenlab/anvio/issues/1970
            if not max(contig_coverage) > 0:
                continue

            # if we are here, we're good to go. let's keep the contig lenght in a separate variable:
            contig_length = len(contig_coverage)

            # to find regions of high coverage, we first need to 'pad' our array to ensure it always
            # starts and ends with 'low coverage'.
            regions_of_contig_covered_enough = np.hstack([[False], contig_coverage >= self.min_coverage_to_define_stretches, [False]])

            # but if there aren't any regions covered enough we want to stop:
            if not regions_of_contig_covered_enough.any():
                continue

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

            # and if there are no coverage stretches long enough, stop:
            if not coverage_stretches_in_contigs[contig_name]:
                continue

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

                    if (anvio.DEBUG or self.verbose) and not anvio.QUIET:
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

                if (anvio.DEBUG or self.verbose) and not anvio.QUIET:
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

                # if there are no true inversions, go back.
                if not len(true_inversions):
                    continue

                # here we have, for a given `contig_name`, `start` and `stop` positions of a stretch in it,
                # and one or more ture inversions. the next step will require us to count them in the original
                # FASTQ files to see their ratio across samples, but for that we need primers to get oligos
                ################################################################################
                self.progress.update(f"{contig_name}[{start}:{stop}]: oligo determination")
                ################################################################################
                self.update_inversions_with_primer_sequences(true_inversions, contig_sequence, start)


                # if we are here, it means we took care of one region.
                # it is time to update our global data dictionaries with
                # everything we know about this inversion, and the stretch
                # it belongs:
                if len(true_inversions):
                    for inv in true_inversions:
                        self.stretches_considered[f"{entry_name}_{sequence_name}"]['invs_found'] = True

                        d = OrderedDict({'entry_id': inv.sequence_name,
                                         'sample_name': entry_name,
                                         'contig_name': contig_name,
                                         'first_seq': inv.first_sequence,
                                         'midline': inv.midline,
                                         'second_seq': utils.rev_comp(inv.second_sequence),
                                         'first_start': inv.first_start + start,
                                         'first_end': inv.first_end + start,
                                         'first_oligo_primer': inv.first_oligo_primer,
                                         'first_oligo_reference': inv.first_oligo_reference,
                                         'second_start': inv.second_start + start,
                                         'second_end': inv.second_end + start,
                                         'second_oligo_primer': inv.second_oligo_primer,
                                         'second_oligo_reference': inv.second_oligo_reference,
                                         'num_mismatches': inv.num_mismatches,
                                         'num_gaps': inv.num_gaps,
                                         'length': inv.length,
                                         'distance': inv.distance})

                        self.inversions[entry_name].append(d)

        self.progress.end()

        self.run.info(f"[Inversions found] In sample {entry_name}", f"{len(self.inversions[entry_name])}", lc="yellow")


    def update_inversions_with_primer_sequences(self, inversions, contig_sequence, start):
        """Takes a set of true inversions and updates them with oligo primers.

        What happens here is this. For true inversions, we want to go back to original
        raw reads (not only mapped reads) and find out the proportion of the inversion
        states for each true inversion across all samples. For this, we need a primer
        that guides our search through short reads.

        Parameters
        ==========
        inversions : list
            List of inversion objects.
        contig_sequence : str
            The contig sequence in which the inversion is found.
        start : int
            The start position of the coverage stretch in which the inversion is found.
            This value is important to normalize relative start stop positions to the
            contig sequence so upstream and downstream genomic context can be recovered.

        Returns
        =======
        None
            But it updates the inversion data structure with `first_oligo_primer` and `second_oligo_primer`.
        """

        for inv in inversions:
            # we will first identify upstream and downstream genomic regions that are before the first palindrome
            # in the inversion and after the second palindrome in the inversion
            first_genomic_region_start = inv.first_start + start - self.oligo_primer_base_length
            first_genomic_region_end = inv.first_start + start
            first_genomic_region = contig_sequence[first_genomic_region_start:first_genomic_region_end]

            second_genomic_region_start = inv.second_end + start
            second_genomic_region_end = inv.second_end + start + self.oligo_primer_base_length + 1
            second_genomic_region = contig_sequence[second_genomic_region_start:second_genomic_region_end]

            # then, we will replace nucleotides found on the contig with `.` character where
            # there are mismatches in the palindrome sequence
            first_with_mismatches = ''.join(['.' if inv.midline[i] == 'x' else inv.first_sequence[i] for i in range(0, inv.length)])
            second_with_mismatches = ''.join(['.' if inv.midline[i] == 'x' else inv.second_sequence[i] for i in range(0, inv.length)])

            # here we first update the inversion object with primers
            inv.first_oligo_primer = first_genomic_region + first_with_mismatches
            inv.second_oligo_primer = utils.rev_comp(second_genomic_region) + second_with_mismatches

            # and finally we update the inversion object with reference oligos
            # found in the original contig
            inv.first_oligo_reference = contig_sequence[first_genomic_region_end + inv.length:first_genomic_region_end + inv.length + self.oligo_length]
            inv.second_oligo_reference = utils.rev_comp(contig_sequence[second_genomic_region_start - inv.length - self.oligo_length:second_genomic_region_start - inv.length])

        return


    def test_inversion_candidates_using_short_reads(self, inversion_candidates, reads):
        """Takes a set of inversion candidates and short reads to look for true inversions"""

        true_inversions = []

        total_num_inversions = len(inversion_candidates)
        current_inversion = 0
        for inversion_candidate in inversion_candidates:
            current_inversion += 1
            num_reads_considered = 0
            match = False
            evidence = ''

            # we are going to be using these to variables to ensure both left and right
            # construct are supported by short reads
            evidence_right = None
            evidence_left = None

            if inversion_candidate.num_mismatches == 0:
                # Here we take advantage of construct symmetry. If the inversion candidate has no
                # mismatches, we only needs to test for the v1 constructs.
                for read in reads:
                    num_reads_considered += 1
                    if not evidence_left:
                        if inversion_candidate.v1_left in read:
                            evidence += 'v1_left and '
                            evidence_left = True

                            if evidence_right:
                                match = True
                                evidence_right = False
                                evidence_left = False
                                break
                    elif not evidence_right:
                        if inversion_candidate.v1_right in read:
                            evidence += 'v1_right'
                            evidence_right = True

                            if evidence_left:
                                match = True
                                evidence_right = False
                                evidence_left = False
                                break
            else:
                # Unfortunately, the inversion candidate has some mismatches, which requires testing
                # for v1 _and_ v2 constructs.
                for read in reads:
                    num_reads_considered += 1
                    if not evidence_left:
                        if inversion_candidate.v1_left in read:
                            evidence += 'v1_left and '
                            evidence_left = True

                            if evidence_right:
                                match = True
                                evidence_right = False
                                evidence_left = False
                                break
                        elif inversion_candidate.v2_left in read:
                            evidence += 'v2_left and'
                            evidence_left = True

                            if evidence_right:
                                match = True
                                evidence_right = False
                                evidence_left = False
                                break
                    if not evidence_right:
                        if inversion_candidate.v1_right in read:
                            evidence += 'v1_right'
                            evidence_right = True

                            if evidence_left:
                                match = True
                                evidence_right = False
                                evidence_left = False
                                break
                        elif inversion_candidate.v2_right in read:
                            evidence += 'v2_right'
                            evidence_right = True

                            if evidence_left:
                                match = True
                                evidence_right = False
                                evidence_left = False
                                break

            if match:
                # we found an inversion candidate that has at least one confirmed
                # construct. We add this one into the list of true inversions:
                true_inversions.append(inversion_candidate)

                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.info_single(f"👍 Candidate {current_inversion} of {total_num_inversions}: confirmed by {evidence} "
                                         f"after {num_reads_considered} reads.", mc="yellow", level=2)

                if not self.check_all_palindromes:
                    # if the user is not interested in testing of every single palindrome
                    # found in the stretch of interest to see whether there may be more
                    # inversion candidates, we return the current list which includes only
                    # one confirmed inversion
                    return true_inversions
            else:
                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.info_single(f"👎 Candidate {current_inversion} of {total_num_inversions}: no confirmation "
                                         f"after processing {num_reads_considered} reads.", mc="red", level=2)

        return true_inversions


    def get_true_inversions_in_stretch(self, inversion_candidates, bam_file, contig_name, start, stop):
        """Survey a bunch of palindromes with 'constructs' to find true/active inversions.

        A true inversion is one with evidence from short reads that it has activity.
        """

        # here we will do something sneaky. if the user wants us to test inversions only using inverted reads,
        # we will abide and move on. as inverted reads will be a fraction of all reads that cover a given
        # region in the BAM file, it is not as costly as going through all of the regular reads. that said,
        # inverted reads will often miss active inversions, and going through the costly route of using all
        # reads will be necessary to be sure. the following if/else block covers that with the consideration
        # of one more key parameter: `self.check_all_palindromes`. the code is intentionally redundant below
        # to ensure maximum readability and to make it easier to follow this relatively complex logic.
        if self.process_only_inverted_reads:
            # if we are here, we will only use the inverted reads and move on. this is the simple case.
            bam_file.fetch_filter = 'inversions'
            reads = [r.query_sequence for r in bam_file.fetch_only(contig_name, start=start, end=stop)]

            if anvio.DEBUG or self.verbose:
                self.progress.reset()
                self.run.info_single(f"Testing {len(inversion_candidates)} palindromes with {len(reads)} REV/REV and FWD/FWD reads to find "
                                     f"true inversions:", mc="green", nl_before=1, nl_after=1)

            true_inversions_in_stretch = self.test_inversion_candidates_using_short_reads(inversion_candidates, reads)
        else:
            # if we are here, we want to use the inverted reads first, and if we find nothing, we want to
            # try all reads. BUT THERE IS ONE MORE CONSIDERATION: if we have multiple palindromes, and if
            # the user asked for `self.check_all_palindromes`, it means we may have a match to one of them
            # in inverted reads, but even if we could have found more if we were to use all reads, we would
            # not be checking them since `if not len(true_inversions_in_stretch)` would no longer be True.
            # Which means, if the user is asking for `self.check_all_palindromes`, we can't do the sneaky
            # thing we wished to do and must work with all reads from the get go :/ first, we will take
            # care of the default case where the user did not ask to check all palindromes:
            if not self.check_all_palindromes:
                bam_file.fetch_filter = 'inversions'
                reads = [r.query_sequence for r in bam_file.fetch_only(contig_name, start=start, end=stop)]
                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.info_single(f"Testing {len(inversion_candidates)} palindromes with {len(reads)} REV/REV and FWD/FWD reads to find "
                                         f"true inversions:", mc="green", nl_before=1, nl_after=1)

                true_inversions_in_stretch = self.test_inversion_candidates_using_short_reads(inversion_candidates, reads)

                if not len(true_inversions_in_stretch):
                    # if we are here, it means inverted reads did not match any of the constructs, and next
                    # we will try all reads.
                    bam_file.fetch_filter = None
                    reads = [r.query_sequence for r in bam_file.fetch_only(contig_name, start=start, end=stop)]

                    if anvio.DEBUG or self.verbose:
                        self.progress.reset()
                        self.run.info_single(f"REV/REV and FWD/FWD reads didn't work :( Now testing the same {len(inversion_candidates)} palindromes "
                                             f"with {len(reads)} regular reads to look for true inversions:", mc="green", nl_before=1, nl_after=1)

                    true_inversions_in_stretch = self.test_inversion_candidates_using_short_reads(inversion_candidates, reads)
            else:
                # so the user wants us to test test all palindromes. this is the case where we shouldn't use
                # inverted reads to speed things up when possible, and simply use all reads.
                bam_file.fetch_filter = None
                reads = [r.query_sequence for r in bam_file.fetch_only(contig_name, start=start, end=stop)]

                if anvio.DEBUG or self.verbose:
                    self.progress.reset()
                    self.run.info_single(f"Testing {len(inversion_candidates)} palindromes with {len(reads)} regular reads to find "
                                         f"true inversions:", mc="green", nl_before=1, nl_after=1)

                true_inversions_in_stretch = self.test_inversion_candidates_using_short_reads(inversion_candidates, reads)

        return true_inversions_in_stretch



    def recover_genomic_context_surrounding_inversions(self):
        """Learn about what surrounds the consensus inversion sites"""

        # we are not wanted
        if self.skip_recovering_genomic_context:
            return

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet)

        # are there genes?
        if not contigs_db.meta['genes_are_called']:
            self.run.warning("There are no gene calls in your contigs database, therefore there is context to "
                             "learn about :/ Your reports will not include a file to study the genomic context "
                             "that surrounds consensus inversions.")

            contigs_db.disconnect()

            return

        # are there functions?
        function_sources_found = contigs_db.meta['gene_function_sources'] or []
        if not len(function_sources_found):
            self.run.warning("There are no functions for genes in your contigs database :/ Your reports on the "
                             "genomic context that surrounds consensus inversions will not have any functions "
                             "for gnes. PITY.")


        self.progress.new('Recovering genomic context surrounding inversions', progress_total_items=len(self.consensus_inversions))
        self.progress.update('...')

        # now we will go through each consensus inversion to populate `self.genomic_context_surrounding_consensus_inversions`
        # with gene calls and functions
        gene_calls_per_contig = {}
        inversions_with_no_gene_calls_around = set([])
        for entry in self.consensus_inversions:
            inversion_id = entry['inversion_id']

            self.progress.update(f"{inversion_id}", increment=True)

            contig_name = entry['contig_name']
            first_start = entry['first_start']
            second_end = entry['second_end']

            # lazy recovery of gene calls in contigs of relevance. should be useful for
            # metagenomes, but useless for isolate genomes with a single contig.
            if contig_name not in gene_calls_per_contig:
                where_clause = f'''contig="{contig_name}" and source="{self.gene_caller_to_consider_in_context}"'''
                gene_calls_per_contig[contig_name] = contigs_db.db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, where_clause=where_clause, error_if_no_data=False)

            gene_calls_in_contig = gene_calls_per_contig[contig_name]

            if not len(gene_calls_in_contig):
                # well, we don't have anything to work with here. let's keep this rebel in mind and
                # move on to the next inversion
                inversions_with_no_gene_calls_around.add(inversion_id)
                continue

            # here we will find out what is hte closes genes to the beginning of the first palindrome
            # inversion and end of the second palindrome
            min_distance_to_first_start, min_distance_to_second_end = float('inf'), float('inf')
            closest_gene_call_to_first_start, closest_gene_call_to_second_end = None, None
            for gene_callers_id, gene_call in gene_calls_in_contig.items():
                if abs(gene_call['start'] - first_start) < min_distance_to_first_start:
                    closest_gene_call_to_first_start = gene_callers_id
                    min_distance_to_first_start = abs(gene_call['start'] - first_start)

                if abs(gene_call['start'] - second_end) < min_distance_to_second_end:
                    closest_gene_call_to_second_end = gene_callers_id
                    min_distance_to_second_end = abs(gene_call['start'] - second_end)

            # now we can recover gene calls of interest for our inversions:
            _range = range(closest_gene_call_to_first_start - self.num_genes_to_consider_in_context,
                           closest_gene_call_to_second_end + self.num_genes_to_consider_in_context)
            gene_caller_ids_of_interest = [c for c in _range if c in gene_calls_in_contig]

            # if there are funtion sources, let's recover them for our genes of interest
            if function_sources_found:
                where_clause = '''gene_callers_id IN (%s)''' % (', '.join([f"{str(g)}" for g in gene_caller_ids_of_interest]))
                hits = list(contigs_db.db.get_some_rows_from_table_as_dict(t.gene_function_calls_table_name, where_clause=where_clause, error_if_no_data=False).values())
            else:
                # so none of these genes have any functions? WELL FINE.
                hits = None

            # we are now ready
            c = []
            for gene_callers_id in gene_caller_ids_of_interest:
                gene_call = gene_calls_in_contig[gene_callers_id]
                gene_call['gene_callers_id'] = gene_callers_id

                # if there are any functions at all, add that to the dict
                if hits:
                    gene_call['functions'] = [h for h in hits if h['gene_callers_id'] == gene_callers_id]

                c.append(gene_call)

            # done! `c` now goes to live its best life as a part of the main class
            self.genomic_context_surrounding_consensus_inversions[inversion_id] = copy.deepcopy(c)

        contigs_db.disconnect()

        self.progress.end()

        self.run.info(f"[Genomic Context] Searched for {PL('inversion', len(self.consensus_inversions))}",
                      f"Recovered for {len(self.genomic_context_surrounding_consensus_inversions)}",
                      nl_before=1,
                      lc="yellow")

        if len(inversions_with_no_gene_calls_around):
            self.run.warning(f"There were one or more inversions that did not have any valid "
                             f"{self.gene_caller_to_consider_in_context} gene calls around them. This may happen if an inversion "
                             f"is occurring in a contig that happens to have no gene calls (either due to its too short, "
                             f"or because you need a Nobel prize). So results from these weird inversions will not appear "
                             f"in your final reports. Here is the list in case you would like to track them down: "
                             f"{', '.join(inversions_with_no_gene_calls_around)}.")

        if not len(self.genomic_context_surrounding_consensus_inversions):
            self.run.warning(f"Even though anvi'o went through all {PL('inversion', len(self.consensus_inversions))} "
                             f"it was unable to recover any genomic context for any of them. So your final reports will "
                             f"not include any insights into the surrounding genomic context of inversions (but otherwise "
                             f"you will be fine).")


    def compute_consensus_inversions(self):
        """Compute a final, consensus list of unique inversions.

        By identifying redundant inversions and reporting only one
        of them, this function reports a final list of inversions
        identifications of which are informed by invdividual samples
        and their coverages in them, but will continue their lives
        as strong and independent inversions.
        """

        total_num_inversions = sum([len(l) for l in self.inversions.values()])
        if total_num_inversions == 0:
            raise ConfigError("You called a function to compute consensus inversions with "
                              "zero inversions. Even though everyone knows that no inversions "
                              "lead to no consensus inversions, and no consensus inversions "
                              "lead to no 🍰")

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
        inversion_counter = 0
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
                        inversion_counter += 1
                        consensus_inversion = OrderedDict({'inversion_id': f'INV_{inversion_counter:04}' })
                        consensus_inversion.update(copy.deepcopy(inv))
                        consensus_inversion['num_samples'] = num_samples
                        consensus_inversion['sample_names'] = sample_names
                        self.consensus_inversions.append(consensus_inversion)
                        consensus_found = True
                        break

        self.progress.end()

        self.run.info(f"[Consensus Inversions] Across {PL('sample', len(self.inversions))}", f"{len(self.consensus_inversions)}", nl_before=1, lc="yellow")


    @staticmethod
    def compute_inversion_activity_for_sample(input_queue, output_queue, samples_dict, primers_dict, oligo_length=6, end_primer_search_after_x_hits=None, run=run_quiet, progress=progress_quiet):
        """Go back to the raw metagenomic reads to compute activity of inversions for a single sample.

        Returns
        =======
        sample_counts : list of tuples
            A list that summarizes frequencies of all oligos for each primer within the sample. Each
            tuple in the list holds information for a single oligo and follows the order,

                >>> (sample_name, inversion_id, oligo_primer, oligo, frequency, relative_abundance)

            where,

                - sample_name: sample name as written in the first column of bams-and-profiles-txt
                - inversion_id: the simplified name (i.e, INV_0001) that matches entries in consensus output file
                - oligo_primer: frequencies reported whether for the first or second palindrome
                - oligo: the actual oligonucleotides of `self.oligo_length` bases
                - frequency: the frequency of oligo in sample
                - relative_abundance: witin-sample relative abundance of the frequency
        """

        while True:
            sample_name = input_queue.get(True)

            # the `samples_dict` knows all samples, `sample_name` knows the sample we are
            # interested in in this thread. we will subsample the `samples_dict` first
            # because otherwise PrimerSearch will search primers for every sample in
            # `samples_dict`
            samples_dict_for_sample = {sample_name: samples_dict[sample_name]}

            # setup the args object
            args = argparse.Namespace(samples_dict=samples_dict_for_sample,
                                      primers_dict=primers_dict,
                                      min_remainder_length=oligo_length,
                                      only_keep_remainder=True)

            # if the user is testing:
            if end_primer_search_after_x_hits:
                args.stop_after = end_primer_search_after_x_hits

            s = PrimerSearch(args, run=run, progress=progress)
            sample_dict, primer_hits = s.process_sample(sample_name)

            # we now have results for a single sample. prepare for return.
            sample_counts = []
            for primer_name in primers_dict:
                oligos = s.get_sequences(primer_name, primer_hits, target='remainders')
                num_oligos = len(oligos)
                oligos_frequency_dict = Counter(oligos)
                reads_found = False

                for oligo, frequency in oligos_frequency_dict.items():
                    sample_counts.append((sample_name, primer_name, oligo, oligo == primers_dict[primer_name]['oligo_reference'], frequency, frequency / num_oligos))
                    if frequency:
                        reads_found = True

                # if the reference oligo has no frequency but reads were found for other oligo
                # then add reference oligo with a frequency of 0
                oligo_reference = primers_dict[primer_name]['oligo_reference']
                if oligo_reference not in oligos_frequency_dict and reads_found:
                    sample_counts.append((sample_name, primer_name, oligo_reference, 'True', 0, 0))

            output_queue.put(sample_counts)


    def populate_consensus_inversions_from_input_file(self):
        """Get the consensus inversions from a previously generated output file"""

        inversions_dict = utils.get_TAB_delimited_file_as_dictionary(self.pre_computed_inversions_path)

        self.consensus_inversions = []

        for inversion_id in sorted(list(inversions_dict.keys())):
            entry = OrderedDict({'inversion_id': inversion_id})
            for tpl in self.essential_keys_to_describe_inversions:
                entry[tpl[0]] = tpl[1](inversions_dict[inversion_id][tpl[0]])
            self.consensus_inversions.append(entry)


    def compute_inversion_activity(self):
        """Go back to the raw metagenomic reads to compute activity of inversions"""

        if self.skip_compute_inversion_activity or not self.raw_r1_r2_reads_are_present:
            return

        if not len(self.consensus_inversions):
            self.run.info_single("Compute inversion activity function is speaking: There are no consensus inversions to "
                                 "compute in-sample activity :/", mc="red")

        sample_names = list(self.profile_db_bam_file_pairs.keys())
        num_samples = len(sample_names)

        # let the user know what is going on
        msg = (f"Now anvi'o will compute in-sample activity of consensus {PL('inversion', len(self.consensus_inversions))} "
               f"across {PL('sample', num_samples)}. Brace yourself and please note that this can "
               f"take a very long time since for each sample, anvi'o will go through each short read to search for two "
               f"sequences per inversion. IF IT COMES TO A POINT WHERE you (or your job on your HPC) can't continue running "
               f"it, this process can be killed without any loss of data from the previous steps, as your primary output "
               f"files must have already been reported. You can always skip this step and search for individual primers "
               f"listed in the consensus output file using the program `anvi-search-primers` with the parameter "
               f"`--min-remainder-length 6` and the flag `--only-report-remainders` to explore inversion activity "
               f"manually")
        self.run.warning(None, header="PERFORMANCE NOTE", lc="yellow")
        if num_samples > self.num_threads:
            self.run.info_single(f"You have {PL('sample', num_samples)} but {PL('thread', self.num_threads)}. Not all samples will be processed "
                                 f"in parallel. Just FYI. {msg}.", level=0, nl_after=1)
        elif self.num_threads > num_samples:
            self.run.info_single(f"You have {PL('sample', num_samples)} but {PL('thread', self.num_threads)}. Since only samples are run in "
                                 f"parallel, the additional {PL('thread', self.num_threads - num_samples)} you have there is not really "
                                 f"useful for anything. So anvi'o will set your number of threads to {num_samples}. So. {msg}.", level=0, nl_after=1)
            self.num_threads = num_samples
        else:
            self.run.info_single(f"{msg}.", level=0, nl_after=1)

        # here we will need to reconstruct a samples_dict and primers_dict to pass to the class
        # `PrimerSearch`. for this we first need to generate a list of primers. for each
        # consensus inversion, which at this point are described in a list like this,
        #
        #    >>> [OrderedDict([('inversion_id', 'INV_0001'),
        #    >>>               ('first_oligo_primer', 'GAGCAAAGATCATGTTTCAAAA.ACGTTC'),
        #    >>>               ('second_oligo_primer', 'TGTTTCAAAA.ACGTTCCCATTGTGTATTC'),
        #    >>>               (...)
        #    >>>               ('num_samples', 1),
        #    >>>               ('sample_names', 'xP3')]),
        #    >>>  OrderedDict([('inversion_id', 'INV_0002'),
        #    >>>               ('first_oligo_primer', 'GGCAGAAATGCCAAGT.CTATCAGAACTT'),
        #    >>>               ('second_oligo_primer', 'AAGT.CTATCAGAACTTAGAGTAGAGCACT'),
        #    >>>               (...)
        #    >>>               ('num_samples', 2),
        #    >>>               ('sample_names', 'xP3,xP2')]),
        #    >>>  OrderedDict([('inversion_id', 'INV_0003'),
        #    >>>               ('first_oligo_primer', 'AGATGTTTCAAAAAACGTTCGTCTATTTGAAC'),
        #    >>>               ('second_oligo_primer', 'AAACGTTCGTCTATTTGAACAAAAACACGTTTT'),
        #    >>>               (...)
        #    >>>               ('num_samples', 1),
        #    >>>               ('sample_names', 'xP3')]),
        #    >>>  (...)]
        #
        # there will have to be two primers, one for `first_oligo_primer` and one for `second_oligo_primer`.
        # so our data structure for primers will include two entries for one consensus inversions that
        # will need to be decomposed later.
        primers_dict = {}
        for entry in self.consensus_inversions:
            inversion_id = entry['inversion_id']
            primers_dict[inversion_id + '-first_oligo_primer'] = {'primer_sequence': entry['first_oligo_primer'],
                                                                  'oligo_reference': entry['first_oligo_reference']}
            primers_dict[inversion_id + '-second_oligo_primer'] = {'primer_sequence': entry['second_oligo_primer'],
                                                                   'oligo_reference': entry['second_oligo_reference']}

        ##########################################################################################
        # MULTITHREADING
        ##########################################################################################

        # setup the input/output queues
        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        # put all the sample names in our input queue
        for sample_name in sample_names:
            input_queue.put(sample_name)

        # engage the proletariat, our hard-working wage-earner class
        workers = []
        for i in range(self.num_threads):
            worker = multiprocessing.Process(target=Inversions.compute_inversion_activity_for_sample,
                                             args=(input_queue,
                                                   output_queue,
                                                   self.profile_db_bam_file_pairs,
                                                   primers_dict,
                                                   self.oligo_length,
                                                   self.end_primer_search_after_x_hits),
                                             kwargs=({'progress': self.progress if self.num_threads == 1 else progress_quiet}))
            workers.append(worker)
            worker.start()


        # these if blocks for progress is an ugly hack, but they are serving a very useful purpose.
        # if the user is working with a single thread, we would like them to see th eactual PrimerSearch
        # activity that comes from `compute_inversion_activity_for_sample`. But when that function is run
        # in multiple threads, we don't want any progress report from that guy, and report an overall
        # progress from this main function. So if there is a single thread, we pass our own progress
        # object to `compute_inversion_activity_for_sample`, which passes it to `PrimerSearch`, which
        # tells us what it is up to on our terminal. If we have multiple threads, we pass a `progress_quiet`
        # to `compute_inversion_activity_for_sample`, which omits `PrimerSearch` progress messages, and
        # instead here we give user an overall idea about how their process is going. Thus, we need to
        # keep track of thread numbers here :/
        if self.num_threads > 1:
            self.progress.new('Inversion activity', progress_total_items=num_samples)
            self.progress.update(f"Processing {PL('sample', num_samples)} and {PL('primer', len(primers_dict))} in {PL('thread', self.num_threads)}.")

        oligo_frequencies = []
        num_samples_processed = 0
        while num_samples_processed < num_samples:
            try:
                oligo_frequencies_for_one_sample = output_queue.get()
                if oligo_frequencies_for_one_sample:
                    oligo_frequencies.extend(oligo_frequencies_for_one_sample)

                num_samples_processed += 1
                self.progress.increment(increment_to=num_samples_processed)
                if self.num_threads > 1:
                    if num_samples_processed < num_samples:
                        self.progress.update(f"Samples processed: {num_samples_processed} of {num_samples}. Still working ...")
                    else:
                        self.progress.update("All done!")
            except KeyboardInterrupt:
                self.run.info_single("Recieved SIGINT, terminating all processes... Don't believe anything you see "
                                     "below this and sanitize all the output files with fire.", nl_before=1, nl_after=1)
                break

        if self.num_threads > 1:
            self.progress.end()

        # always double-tap?
        for worker in workers:
            worker.terminate()

        ##########################################################################################
        # /MULTITHREADING
        ##########################################################################################

        # we must be ready with all the oligo frequencies by now. it is time to report this
        output_path = os.path.join(self.output_directory, 'INVERSION-ACTIVITY.txt')
        headers = ['sample', 'inversion_id', 'oligo_primer', 'oligo', 'reference', 'frequency_count', 'relative_abundance']
        with open(output_path, 'w') as output:
            output.write('\t'.join(headers) + '\n')
            for e in oligo_frequencies:
                inversion_id, oligo_primer = '-'.join(e[1].split('-')[:-1]), e[1].split('-')[-1]
                output.write(f"{e[0]}\t{inversion_id}\t{oligo_primer}\t{e[2]}\t{e[3]!s}\t{e[4]}\t{e[5]:.3f}\n")

        self.run.info('Long-format reporting file for inversion activity', output_path, mc='green')


    def plot_coverage(self, sequence_name, coverage, num_bins=100):
        if anvio.QUIET:
            return

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
        self.run.info("[General] Input BAMs and profiles file", self.bams_and_profiles_file_path)
        self.run.info("[General] Number of samples to process", len(self.profile_db_bam_file_pairs))
        self.run.info("[General] R1/R2 for raw reads present?", "True" if self.raw_r1_r2_reads_are_present else "False")
        self.run.info("[General] Be talkative (--verbose)?", "True" if self.verbose else "False", nl_after=1)

        # do we have a previously computed list of consensus inversions to focus for inversion
        # activity calculations?
        if self.pre_computed_inversions_path:
            self.run.warning("Anvi'o is taking a shortcut to calculate inversion activity using the inversions you have "
                             "provided in the 'consensus inversions' file. There is nothing for you to be concerned about "
                             "-- except the fact that some very fancy coding is at play here and catastrophic failures "
                             "are not a remote possibility. Still, anvi'o prints this message in green for positive "
                             "vibes 🥲", header="🌈 PRE-COMPUTED INVERSIONS FOUND 🌈", lc="green")

            self.populate_consensus_inversions_from_input_file()
            self.compute_inversion_activity()

            # WE'RE DOnE HERE. DoNe.
            return

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

        if not self.skip_recovering_genomic_context:
            self.run.info("[Genomic context] Recover and report genomic context?",  "True", mc="green")
            self.run.info("[Genomic context] Gene caller to use", self.gene_caller_to_consider_in_context)
            self.run.info("[Genomic context] Number of genes to consider",  self.num_genes_to_consider_in_context, nl_after=1)
        else:
            self.run.info("[Genomic context] Recover and report genomic context?",  "False", mc="red", nl_after=1)

        # are we to compute inversion activity by going through raw reads?
        inversion_activity_will_be_computed = self.raw_r1_r2_reads_are_present and not self.skip_compute_inversion_activity
        self.run.info("[Inversion activity] Compute inversion activity?",  "True" if inversion_activity_will_be_computed else "False", mc=("green" if inversion_activity_will_be_computed else "red"))
        if not inversion_activity_will_be_computed:
            if not self.raw_r1_r2_reads_are_present:
                self.run.info("[Inversion activity] Not computing because",  "`bams-and-profiles-txt` does not contain raw reads for r1/r2", nl_after=1)
            elif self.skip_compute_inversion_activity:
                self.run.info("[Inversion activity] Not computing because",  "The user asked not to do it :/", nl_after=1)
            else:
                self.run.info("[Inversion activity] Not computing because",  "Anvi'o has no idea what it is doing", nl_after=1)
        else:
            self.run.info("[Inversion activity] Number of threads", self.num_threads, mc=("green" if self.num_threads > 1 else "red"))
            if self.end_primer_search_after_x_hits:
                self.run.info("[Inversion activity] Oligo primer base length", self.oligo_primer_base_length)
                self.run.info("[Inversion activity Debug] Num hits to end primer search",  self.end_primer_search_after_x_hits, mc="red", nl_after=1)
            else:
                self.run.info("[Inversion activity] Oligo primer base length", self.oligo_primer_base_length, nl_after=1)

        if self.target_contig:
            self.run.info("[Targetly] Only focus on the contig:",  self.target_contig, mc="red")
        if self.target_region_start:
            self.run.info("[Targetly] START of the region of interest:",  self.target_region_start, mc="red")
        if self.target_region_end:
            self.run.info("[Targetly] END of the region of interest:",  self.target_region_end, mc="red")

        self.run.width = w

        for entry_name in self.profile_db_bam_file_pairs:
            bam_file_path = self.profile_db_bam_file_pairs[entry_name]['bam_file_path']
            profile_db_path = self.profile_db_bam_file_pairs[entry_name]['profile_db_path']

            # populate `self.inversions` with inversions associated with `entry_name`
            self.process_db(entry_name, profile_db_path, bam_file_path)

        # this is time to check `self.inversions`. if there is nothing in it, then there is no
        # need to continue with the analysis, reporting, etc.
        total_num_inversions = sum([len(l) for l in self.inversions.values()])
        if total_num_inversions == 0:
            samples_note = f"across your {self.profile_db_bam_file_pairs} samples" \
                                if len(self.profile_db_bam_file_pairs) > 1 else "in your single sample"
            self.run.warning(f"The code came all the way down here with zero inversions {samples_note}. "
                             f"Zero. Zip. Which means either there truly are no inversions in your data, "
                             f"or anvi'o couldn't find them given the search criteria you have defined "
                             f"(or criteria defined de facto by the default parameters of the program). "
                             f"Either ways, this is a goodbye.", header="NO INVERSIONS FOUND")
            return

        # here we know every single inversion. The same inversion site might be
        # found as true inversion in multiple samples when bams and profiles file
        # includes more than one sample. To avoid redundancy, we wish to report
        # a concensus file that describes only unique inversions.
        self.compute_consensus_inversions()

        # here we will try to generate some insights into the genomic context that
        # surround inversions, so this data can be reported along with other files.
        self.recover_genomic_context_surrounding_inversions()

        # we want to report now because the very last step can take a long time, and
        # if the user kills the process, we don't want them to go home empty handed.
        self.report()

        # here we have our consensus inversions in `self.consensus_inversions`. It is time
        # to go back to the raw reads and compute their activity IF r1/r2 files are provided
        # AND the user didn't ask this step to be skipped.
        self.compute_inversion_activity()


    def sanity_check(self):
        """Basic checks for a smooth operation"""

        if self.pre_computed_inversions_path:
            if self.skip_compute_inversion_activity:
                raise ConfigError("You can't provide consensus inversions to calculate inversion activity, and then ask "
                                  "anvi'o to skip calculating inversion activity :/")

            if not self.raw_r1_r2_reads_are_present:
                raise ConfigError("You asked anvi'o to calculate inversion activity across samples, but your bams-and-profiles-txt "
                                  "does not include raw R1/R2 reads :(")

        if not self.skip_recovering_genomic_context:
            if not dbops.ContigsDatabase(self.contigs_db_path, run=run_quiet, progress=progress_quiet).meta['genes_are_called']:
                raise ConfigError("Your parameter setup asks anvi'o to recover genomic context of active inversions at the end "
                                  "but the contigs database does not have any genes called :/ For the sake of being explicit "
                                  "about it anvi'o requests you to use the flag `--skip-recovering-genomic-context` so it is clear "
                                  "to everyone that this step is meant to be skipped.")

            if 'RNA' in self.gene_caller_to_consider_in_context:
                raise ConfigError(f"Anvi'o truly hates to do this, but '{self.gene_caller_to_consider_in_context}' is really not a "
                                  f"very good source for gene calls to consider to understand genomic context for inversions :/ Please "
                                  f"write to anvi'o developers if you think they're wrong.")

            gene_callers = [g[0] for g in dbops.ContigsDatabase(self.contigs_db_path).meta['gene_callers'] if 'RNA' not in g[0]]
            if self.gene_caller_to_consider_in_context not in gene_callers:
                raise ConfigError(f"The gene caller '{self.gene_caller_to_consider_in_context}' is not one of the gene callers "
                                  f"your contigs database knows about :/ Please use the `--gene-caller` parameter to select one "
                                  f"that fits .. such as {PL('this one', len(gene_callers), alt='one of these')}: {', '.join(gene_callers)}.")

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

        if self.num_threads < 0:
            raise ConfigError("{self.num_threads} for number of threads? You must be joking, Mr. Feynman.")

        if self.num_genes_to_consider_in_context < 1:
            raise ConfigError("The number of genes to consider when recovering genomic context around active inversions can't be less than "
                              "one :/ If you need to consider less than 1 gene around your inversions, you should use the flag "
                              "`--skip-recovering-genomic-context` instead of confusing anvi'o :(")

        if self.target_region_start and self.target_region_end:
            if self.target_region_end <= self.target_region_start:
                raise ConfigError(f"The end position of the target region you wish to focus on must be larger "
                                  f"than the start poisition of the same region .. for obvious reasons. But "
                                  f"{pp(self.target_region_end)} is not larger than {pp(self.target_region_start)}.")

        if (self.target_region_start and self.target_region_start < 0) or (self.target_region_end and self.target_region_end < 0):
            raise ConfigError("Anvi'o kindly advices you to check yourself before you wreck yourself. BUT, it leaves "
                              "HOW you should check yourself completely up to you as a mystery for you to solve.")


    def report(self):
        """Reporting per-sample as well as consensus inversions, along with other reporting files"""

        headers = []
        # just learning headers here. NO JUDGING.
        for entry_name in self.inversions:
            if len(self.inversions[entry_name]):
                headers = list(self.inversions[entry_name][0].keys())
                break

        self.run.warning(None, header="REPORTING OUTPUTS", lc="green")

        ################################################################################################
        # Per sample inversions
        ################################################################################################
        for entry_name in self.inversions:
            if len(self.inversions[entry_name]):
                output_path = os.path.join(self.output_directory, f'INVERSIONS-IN-{entry_name}.txt')
                with open(output_path, 'w') as output:
                    output.write('\t'.join(headers) + '\n')
                    for v in self.inversions[entry_name]:
                        output.write('\t'.join([f"{v[k]}" for k in headers]) + '\n')
                self.run.info(f'Inversions in {entry_name}', output_path, mc='green')
            else:
                self.run.info(f'Inversions in {entry_name}', 'No true inversions in this one :/', mc='red')

        headers = ['inversion_id', 'contig_name', 'first_seq', 'midline', 'second_seq', 'first_start',
                   'first_end', 'second_start', 'second_end', 'num_mismatches', 'num_gaps', 'length',
                   'distance', 'num_samples', 'sample_names', 'first_oligo_primer',  'first_oligo_reference',
                   'second_oligo_primer', 'second_oligo_reference']

        ################################################################################################
        # Consensus inversions
        ################################################################################################
        output_path = os.path.join(self.output_directory, 'INVERSIONS-CONSENSUS.txt')
        with open(output_path, 'w') as output:
            output.write('\t'.join(headers) + '\n')
            for v in self.consensus_inversions:
                output.write('\t'.join([f"{v[k]}" for k in headers]) + '\n')
        self.run.info('Reporting file for consensus inversions', output_path, mc='green', nl_before=1)

        ################################################################################################
        # All stretches considered
        ################################################################################################
        output_path = os.path.join(self.output_directory, 'ALL-STRETCHES-CONSIDERED.txt')
        headers = ['entry_id', 'sequence_name', 'sample_name', 'contig_name', 'start_stop', 'max_coverage', 'num_palindromes_found', 'true_inversions_found']
        utils.store_dict_as_TAB_delimited_file(self.stretches_considered, output_path, headers=headers)
        self.run.info('Reporting file on all stretches considered', output_path, nl_before=1, nl_after=1)

        ################################################################################################
        # Genomic context
        ################################################################################################
        self.report_genomic_context_surrounding_inversions()


    def report_genomic_context_surrounding_inversions(self):
        """Reports two long-format output files for genes and functions around inversion"""

        if self.skip_recovering_genomic_context:
            return

        if not len(self.genomic_context_surrounding_consensus_inversions):
            return

        # we are in business
        genes_output_path = os.path.join(self.output_directory, 'INVERSIONS-CONSENSUS-SURROUNDING-GENES.txt')
        functions_output_path = os.path.join(self.output_directory, 'INVERSIONS-CONSENSUS-SURROUNDING-FUNCTIONS.txt')

        genes_output_headers = ["gene_callers_id", "start", "stop", "direction", "partial", "call_type", "source", "version", "contig"]
        functions_output_headers = ["gene_callers_id", "source", 'accession', 'function']

        with open(genes_output_path, 'w') as genes_output, open(functions_output_path, 'w') as functions_output:
            genes_output.write("inversion_id\tentry_type\t%s\n" % '\t'.join(genes_output_headers))
            functions_output.write("inversion_id\t%s\n" % '\t'.join(functions_output_headers))

            for v in self.consensus_inversions:

                # this is the inversion id we are working with here:
                inversion_id = v['inversion_id']

                # but we will first create fake gene call entries for inversions:
                d = dict([(h, '') for h in genes_output_headers])

                # fill in non-empty data for the first palindrome in inversion and insert it:
                d['contig'] = v['contig_name']
                d['start'] = v['first_start']
                d['stop'] = v['first_end']
                genes_output.write(f"{inversion_id}\tFIRST_PALINDROME\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                # fill in non-empty data for the second palindrome in inversion and insert it:
                d['contig'] = v['contig_name']
                d['start'] = v['second_start']
                d['stop'] = v['second_end']
                genes_output.write(f"{inversion_id}\tSECOND_PALINDROME\t%s\n" % '\t'.join([f"{d[h]}" for h in genes_output_headers]))

                if inversion_id not in self.genomic_context_surrounding_consensus_inversions:
                    # we don't have anything else to report on this one :/
                    continue

                for gene_call in self.genomic_context_surrounding_consensus_inversions[inversion_id]:
                    genes_output.write(f"{inversion_id}\tGENE\t%s\n" % '\t'.join([f"{gene_call[h]}" for h in genes_output_headers]))

                    if 'functions' in gene_call:
                        for hit in gene_call['functions']:
                            functions_output.write(f"{inversion_id}\t{hit['gene_callers_id']}\t{hit['source']}\t{hit['accession'].split('!!!')[0]}\t{hit['function'].split('!!!')[0]}\n")
                    else:
                        functions_output.write(f"{inversion_id}\t{gene_call['gene_callers_id']}\t\t\t\n")

        self.run.info('Reporting file on gene context', genes_output_path)
        self.run.info('Reporting file on functional context', functions_output_path, nl_after=1)
