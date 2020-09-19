# -*- coding: utf-8
# pylint: disable=line-too-long
""" Classes for tRNA-seq dataset operations. anvi-trnaseq is the default client using this module. """

import gc
import os
import shutil
import itertools
import numpy as np
import pickle as pkl
import multiprocessing as mp

from bisect import bisect
from collections import OrderedDict, deque
from itertools import combinations, product

import anvio
import anvio.utils as utils
import anvio.tables as tables
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.trnaidentifier as trnaidentifier

from anvio.errors import ConfigError
from anvio.dbops import TRNASeqDatabase
from anvio.agglomeration import Agglomerator
from anvio.sequence import Aligner, Dereplicator
from anvio.constants import THREEPRIME_VARIANTS, TRNA_FEATURE_NAMES


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


unambiguous_nucs = ('A', 'C', 'G', 'T')
nuc_int_dict = {nuc: i for i, nuc in enumerate(unambiguous_nucs, start=1)}
int_nuc_dict = {i: nuc for i, nuc in enumerate(unambiguous_nucs, start=1)}


class UniqueSeq:
    __slots__ = ('seq_string',
                 'representative_name',
                 'has_complete_feature_set',
                 'input_count',
                 'identification_method',
                 'acceptor_length',
                 'extra_fiveprime_length')

    def __init__(self,
                 seq_string,
                 representative_name,
                 input_count,
                 identification_method=None,
                 acceptor_length=None,
                 has_complete_feature_set=None,
                 extra_fiveprime_length=None):
        """A dereplicated tRNA-seq read, with information from tRNA feature profiling"""

        self.seq_string = seq_string
        self.representative_name = representative_name
        self.input_count = input_count
        # If dealing with tRNA, identification method 0 = profiled, 1 = mapped
        self.identification_method = identification_method
        self.acceptor_length = acceptor_length
        self.has_complete_feature_set = has_complete_feature_set
        self.extra_fiveprime_length = extra_fiveprime_length


class TrimmedSeq:
    __slots__ = ('seq_string',
                 'unique_seqs',
                 'has_complete_feature_set',
                 'input_count',
                 'unique_with_extra_fiveprime_count',
                 'input_with_extra_fiveprime_count',
                 'representative_name',
                 'input_acceptor_variant_count_dict',
                 'identification_method',
                 'normalized_seq_count')

    def __init__(self, seq_string, unique_seqs, skip_init=False):
        """A tRNA sequence with bases trimmed 5' of the acceptor stem and 3' of the discriminator"""

        self.seq_string = seq_string
        self.unique_seqs = unique_seqs # list of UniqueSeq objects
        self.has_complete_feature_set = unique_seqs[0].has_complete_feature_set
        self.normalized_seq_count = 0

        if skip_init:
            self.input_count = None
            self.unique_with_extra_fiveprime_count = None
            self.input_with_extra_fiveprime_count = None
            self.representative_name = None
            self.input_acceptor_variant_count_dict = None
            self.identification_method = None
        else:
            self.init()


    def init(self):
        """Set attributes representative of a final set of input `UniqueSeq` objects"""

        self.input_count = sum([unique_seq.input_count for unique_seq in self.unique_seqs])

        self.unique_with_extra_fiveprime_count = sum([1 if unique_seq.extra_fiveprime_length else 0
                                                      for unique_seq in self.unique_seqs])
        self.input_with_extra_fiveprime_count = sum([unique_seq.input_count if unique_seq.extra_fiveprime_length else 0
                                                     for unique_seq in self.unique_seqs])

        # The representative name is chosen as follows:
        # 1. Most abundant full-length tRNA (no extra 5' bases), ignoring acceptor sequence
        # 2. Most abundant longer-than-full-length tRNA
        # 3. Most abundant fragmentary tRNA
        # Sort such that the first sequence is the most abundant longest and the last is the least abundant shortest.
        unique_seqs = sorted(self.unique_seqs,
                             key=lambda unique_seq: (-unique_seq.extra_fiveprime_length, -unique_seq.input_count))

        if unique_seqs[0].extra_fiveprime_length > 0:
            # If there is also a unique sequence that was ultimately trimmed down
            # to the same sequence as the sequence with extra 5' bases, it must be a full-length sequence.
            if unique_seqs[-1].extra_fiveprime_length == 0:
                # Sort such that the last sequence is the most abundant shortest.
                representative_name = sorted(unique_seqs,
                                             key=lambda unique_seq: (-unique_seq.extra_fiveprime_length,
                                                                     unique_seq.input_count))[-1].representative_name
            else:
                representative_name = unique_seqs[0].representative_name
        else:
            # ALL unique sequences are EITHER full-length OR a fragment.
            representative_name = unique_seqs[0].representative_name

        self.representative_name = representative_name

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for unique_seq in self.unique_seqs:
            if unique_seq.acceptor_length: # unique_seq need not have an acceptor
                acceptor_seq_string = unique_seq.seq_string[-unique_seq.acceptor_length: ]
                input_acceptor_variant_count_dict[acceptor_seq_string] += unique_seq.input_count
        self.input_acceptor_variant_count_dict = input_acceptor_variant_count_dict

        identification_methods = set(unique_seq.identification_method for unique_seq in self.unique_seqs)
        if len(identification_methods) == 1:
            self.identification_method = identification_methods.pop()
        else:
            raise ConfigError("A TrimmedSeq should not be made from UniqueSeq objects "
                              "with different identification methods. "
                              "Trimmed tRNA sequences will EITHER be formed from "
                              "\"profiled\" tRNA sequences or \"mapped\" tRNA sequences, "
                              "because they are of different lengths and are fragments from different parts of the tRNA.")


class NormalizedSeq:
    __slots__ = ('trimmed_seqs',
                 'representative_name',
                 'seq_string',
                 'has_complete_feature_set',
                 'start_positions',
                 'end_positions',
                 'input_count',
                 'input_with_extra_fiveprime_count',
                 'input_acceptor_variant_count_dict',
                 'trimmed_seqs_mapped_without_extra_fiveprime_count',
                 'input_seqs_mapped_without_extra_fiveprime_count',
                 'trimmed_seqs_mapped_with_extra_fiveprime_count',
                 'input_seqs_mapped_with_extra_fiveprime_count',
                 'coverages',
                 'unique_coverages',
                 'modified_seq')

    def __init__(self, trimmed_seqs, start_positions=None, end_positions=None, skip_init=False):
        """A longer tRNA sequence consolidated from shorter tRNA fragments"""

        self.trimmed_seqs = trimmed_seqs # list of TrimmedSeq objects
        for trimmed_seq in trimmed_seqs:
            trimmed_seq.normalized_seq_count += 1
        self.representative_name = trimmed_seqs[0].representative_name
        self.seq_string = trimmed_seqs[0].seq_string
        self.has_complete_feature_set = trimmed_seqs[0].has_complete_feature_set
        if start_positions and end_positions:
            self.start_positions = start_positions
            self.end_positions = end_positions
        elif (not start_positions) and (not end_positions):
            # Trimmed seqs were dereplicated from the 3' end of the normalized sequence.
            normalized_seq_length = len(self.seq_string)
            self.start_positions = [normalized_seq_length - len(trimmed_seq.seq_string) for trimmed_seq in self.trimmed_seqs]
            self.end_positions = [normalized_seq_length] * len(trimmed_seqs)
        else:
            self.start_positions = None
            self.end_positions = None

        # It is useful to know which ModifiedSeq, if any, encompasses this NormalizedSeq.
        self.modified_seq = None

        if skip_init:
            self.input_count = None
            self.input_with_extra_fiveprime_count = None
            self.input_acceptor_variant_count_dict = None
            self.trimmed_seqs_mapped_without_extra_fiveprime_count = None
            self.input_seqs_mapped_without_extra_fiveprime_count = None
            self.trimmed_seqs_mapped_with_extra_fiveprime_count = None
            self.input_seqs_mapped_with_extra_fiveprime_count = None
            self.coverages = None
            self.unique_coverages = None
        else:
            self.init()


    def init(self):
        """Set the attributes representative of a finalized list of `TrimmedSeq` objects"""

        self.input_count = sum([trimmed_seq.input_count for trimmed_seq in self.trimmed_seqs])

        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for trimmed_seq in self.trimmed_seqs:
            for acceptor_seq_string, input_count in trimmed_seq.input_acceptor_variant_count_dict.items():
                if input_count > 0:
                    input_acceptor_variant_count_dict[acceptor_seq_string] += input_count
        self.input_acceptor_variant_count_dict = input_acceptor_variant_count_dict

        input_with_extra_fiveprime_count = 0
        trimmed_seqs_mapped_without_extra_fiveprime_count = 0
        input_seqs_mapped_without_extra_fiveprime_count = 0
        trimmed_seqs_mapped_with_extra_fiveprime_count = 0
        input_seqs_mapped_with_extra_fiveprime_count = 0
        coverages = np.zeros(len(self.seq_string), dtype=int)
        unique_coverages = np.zeros(len(self.seq_string), dtype=int)
        for trimmed_seq, start_position, end_position in zip(self.trimmed_seqs, self.start_positions, self.end_positions):
            if trimmed_seq.identification_method == 1: # 1 => mapped
                # TrimmedSeqs are comprised of EITHER profiled (0) OR mapped (1) UniqueSeqs.
                if trimmed_seq.unique_with_extra_fiveprime_count == 0:
                    trimmed_seqs_mapped_without_extra_fiveprime_count += 1
                    input_seqs_mapped_without_extra_fiveprime_count += trimmed_seq.input_count
                else:
                    input_with_extra_fiveprime_count += trimmed_seq.input_with_extra_fiveprime_count
                    trimmed_seqs_mapped_with_extra_fiveprime_count += 1
                    input_seqs_mapped_with_extra_fiveprime_count += trimmed_seq.input_count
            else:
                input_with_extra_fiveprime_count += trimmed_seq.input_with_extra_fiveprime_count

            coverages[start_position: end_position] += trimmed_seq.input_count
            if trimmed_seq.normalized_seq_count == 1:
                unique_coverages[start_position: end_position] += trimmed_seq.input_count
        self.trimmed_seqs_mapped_without_extra_fiveprime_count = trimmed_seqs_mapped_without_extra_fiveprime_count
        self.input_seqs_mapped_without_extra_fiveprime_count = input_seqs_mapped_without_extra_fiveprime_count
        self.trimmed_seqs_mapped_with_extra_fiveprime_count = trimmed_seqs_mapped_with_extra_fiveprime_count
        self.input_seqs_mapped_with_extra_fiveprime_count = input_seqs_mapped_with_extra_fiveprime_count
        self.coverages = coverages
        self.unique_coverages = unique_coverages


class ModifiedSeq:
    __slots__ = ('normalized_seqs_with_substitutions',
                 'all_normalized_seqs',
                 'representative_name',
                 'substitution_indices',
                 'normalized_seqs_with_deletions',
                 'deletion_indices',
                 'coverages',
                 'unique_coverages',
                 'deletion_coverages',
                 'deletion_unique_coverages',
                 'substitution_coverages_array',
                 'substitution_unique_coverages_array',
                 'input_count',
                 'input_with_extra_fiveprime_count',
                 'input_acceptor_variant_count_dict',
                 'input_seqs_mapped_without_extra_fiveprime_count',
                 'input_seqs_mapped_with_extra_fiveprime_count',
                 'consensus_seq_string')

    def __init__(self, normalized_seqs_with_substitutions, substitution_indices, init_substitutions=True):
        """A tRNA sequence with sites of predicted modification-induced substitutions and deletions

        Parameters
        ==========
        normalized_seqs : list
            NormalizedSeq objects representing sequences with distinct modification-induced substitutions
            Sequences with modification-induced deletions must be added later.
            The first sequence in the list should be longest or tied for longest.

        substitution_indices : list
            Indices of modification-induced substitutions
            These indices are in the reference frame of the modified sequence/first sequence from the `normalized_seqs` argument.

        init_substitutions : bool, True
            Triggers the analysis of added normalized sequences, which should contain substitutions but not deletions
            Normalized sequences with deletions should be added
            after sequences with substitutions have been added and the `init_substitutions` method run.
            Deletions occur at and around substitution sites and are thus identified after substitutions.
        """

        self.normalized_seqs_with_substitutions = normalized_seqs_with_substitutions
        self.all_normalized_seqs = []
        for normalized_seq in self.normalized_seqs_with_substitutions:
            self.all_normalized_seqs.append(normalized_seq)
            normalized_seq.modified_seq = self
        self.representative_name = self.all_normalized_seqs[0].representative_name
        self.substitution_indices = substitution_indices


        if init_substitutions:
            self.init_substitutions()
        else:
            self.normalized_seqs_with_deletions = None
            self.deletion_indices = None
            self.coverages = None
            self.unique_coverages = None
            self.deletion_coverages = None
            self.deletion_unique_coverages = None
            self.substitution_coverages_array = None
            self.substitution_unique_coverages_array = None
            self.input_count = None
            self.input_with_extra_fiveprime_count = None
            self.input_acceptor_variant_count_dict = None
            self.input_seqs_mapped_without_extra_fiveprime_count = None
            self.input_seqs_mapped_with_extra_fiveprime_count = None
            self.consensus_seq_string = None


    def init_substitutions(self):
        """This method analyzes normalized sequences with substitutions (and without deletions)

        This method should be called after adding normalized sequences with substitutions
        and before adding normalized sequences with deletions (using `add_normalized_seq_with_deletion`).
        The full complement of substitutions is used in identifying deletions,
        as deletions occur at and around substitution sites.
        """

        # Allow normalized sequences with deletions to be added
        # by changing attributes initialized with None to empty lists.
        self.normalized_seqs_with_deletions = []
        self.deletion_indices = []
        # Deletion coverages are stored in a list rather than numpy array
        # to facilitate the insertion of new deletion positions as sequences with deletions are added.
        self.deletion_coverages = []
        self.deletion_unique_coverages = []

        modified_seq_length = len(self.all_normalized_seqs[0].seq_string)
        coverages = np.zeros(modified_seq_length, dtype=int)
        unique_coverages = np.zeros(modified_seq_length, dtype=int)

        substitution_coverages_array = np.zeros((len(self.substitution_indices), len(nuc_int_dict)), dtype=int)
        substitution_unique_coverages_array = np.zeros((len(self.substitution_indices), len(nuc_int_dict)), dtype=int)

        # Set coverage and input sequence count attributes.
        processed_trimmed_seq_names = []
        input_count = 0
        input_with_extra_fiveprime_count = 0
        input_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        input_seqs_mapped_without_extra_fiveprime_count = 0
        input_seqs_mapped_with_extra_fiveprime_count = 0
        for normalized_seq in self.all_normalized_seqs:
            normalized_seq_length = len(normalized_seq.seq_string)
            normalized_seq_start_in_modified_seq = modified_seq_length - normalized_seq_length
            for (trimmed_seq,
                 trimmed_seq_start_in_normalized_seq,
                 trimmed_seq_stop_in_normalized_seq) in zip(normalized_seq.trimmed_seqs,
                                                            normalized_seq.start_positions,
                                                            normalized_seq.end_positions):
                if trimmed_seq.representative_name in processed_trimmed_seq_names:
                    continue
                processed_trimmed_seq_names.append(trimmed_seq.representative_name)

                trimmed_seq_input_count = trimmed_seq.input_count

                input_count += trimmed_seq_input_count

                if trimmed_seq.identification_method == 1: # 1 => mapped
                    # TrimmedSeqs are comprised of EITHER profiled (0) OR mapped (1) UniqueSeqs.
                    if trimmed_seq.unique_with_extra_fiveprime_count == 0:
                        input_seqs_mapped_without_extra_fiveprime_count += trimmed_seq_input_count
                    else:
                        input_with_extra_fiveprime_count += trimmed_seq.input_with_extra_fiveprime_count
                        input_seqs_mapped_with_extra_fiveprime_count += trimmed_seq_input_count
                else:
                    input_with_extra_fiveprime_count += trimmed_seq.input_with_extra_fiveprime_count

                for acceptor_variant_seq, acceptor_variant_input_count in trimmed_seq.input_acceptor_variant_count_dict.items():
                    input_acceptor_variant_count_dict[acceptor_variant_seq] += acceptor_variant_input_count

                # Check if the trimmed sequence only occurs in the present normalized sequence or in others as well.
                if trimmed_seq.normalized_seq_count == 1:
                    is_trimmed_seq_unique_to_normalized_seq = True
                else:
                    is_trimmed_seq_unique_to_normalized_seq = False

                # Augment the coverage of nucleotides by the number of input sequences (reads) in the trimmed sequence.
                trimmed_seq_start_in_modified_seq = normalized_seq_start_in_modified_seq + trimmed_seq_start_in_normalized_seq
                trimmed_seq_stop_in_modified_seq = trimmed_seq_start_in_modified_seq + len(trimmed_seq.seq_string)
                coverages[trimmed_seq_start_in_modified_seq: trimmed_seq_stop_in_modified_seq] += trimmed_seq_input_count
                if is_trimmed_seq_unique_to_normalized_seq:
                    unique_coverages[trimmed_seq_start_in_modified_seq: trimmed_seq_stop_in_modified_seq] += trimmed_seq_input_count

                # Augment the coverage of the specific nucleotides from the trimmed sequence located at substitution sites.
                for n, substitution_index in enumerate(self.substitution_indices):
                    # Check that the position of the substitution is covered by the present trimmed sequence.
                    if not (trimmed_seq_start_in_modified_seq <= substitution_index < trimmed_seq_stop_in_modified_seq):
                        continue
                    nuc_int = nuc_int_dict[normalized_seq.seq_string[substitution_index - normalized_seq_start_in_modified_seq]]
                    substitution_number_in_modified_seq = self.substitution_indices.index(substitution_index)
                    substitution_coverages_array[n, nuc_int - 1] += trimmed_seq_input_count
                    if is_trimmed_seq_unique_to_normalized_seq:
                        substitution_unique_coverages_array[n, nuc_int - 1] += trimmed_seq_input_count
        self.coverages = coverages
        self.unique_coverages = unique_coverages
        self.substitution_coverages_array = substitution_coverages_array
        self.substitution_unique_coverages_array = substitution_unique_coverages_array
        self.input_count = input_count
        self.input_with_extra_fiveprime_count = input_with_extra_fiveprime_count
        self.input_acceptor_variant_count_dict = input_acceptor_variant_count_dict
        self.input_seqs_mapped_without_extra_fiveprime_count = input_seqs_mapped_without_extra_fiveprime_count
        self.input_seqs_mapped_with_extra_fiveprime_count = input_seqs_mapped_with_extra_fiveprime_count


    def add_normalized_seq_with_deletion(self,
                                         normalized_seq_with_deletion,
                                         normalized_seq_start_in_modified_seq,
                                         deletion_indices_in_modified_seq,
                                         deletion_indices_in_normalized_seq,
                                         substitution_indices_in_modified_seq,
                                         substitution_indices_in_normalized_seq):
        """Add a normalized sequence with deletions to the modified sequence.

        This method should only be called after `init_substitutions` has been run.

        This method ensures that trimmed sequences comprising the normalized sequence
        do not contribute coverage data to the modified sequence
        when trimmed sequences are already represented in normalized sequences already added to the modified sequence.

        For a trimmed sequence to cover a deletion, it cannot end in the deletion.
        In other words, it must have nucleotides remaining on either side of the deletion site.

        Deletion coverage is equivalent to nucleotide coverage,
        the count of input sequences (reads) containing the deletion.
        As with nucleotides, both total and unique coverage is calculated for deletions,
        since input sequences may map to multiple normalized and modified sequences.

        Parameters
        ==========
        normalized_seq_with_deletion : NormalizedSeq object
            One or more nucleotides in the modified sequence are deleted in this sequence.

        normalized_seq_start_in_modified_seq : int
            The index at which the normalized sequence starts in the modified sequence
            The normalized sequence cannot start with a deletion,
            as it doesn't make sense for the 5' end of a sequence to be identified as a deletion
            when it could instead be a truncated tRNA fragment.

        deletion_indices_in_modified_seq : list
            The list of positions in the modified sequence that are deleted in the normalized sequence
            These are the positions at which nucleotides are deleted in the normalized sequence.
            This list must be the same length as `deletion_indices_in_normalized_seq`, as the deletions must correspond.

        deletion_indices_in_normalized_seq : list
            The list of indices marking deletions in the normalized sequence
            These are the positions of the remaining nucleotides on the 5' side of deletions.
            If there are multiple adjacent deletions,
            there should be consecutive identical entries in this list for the same index 5' of the deletions.
            This list must be the same length as `deletion_indices_in_modified_seq`, as the deletions must correspond.

        substitution_indices_in_modified_seq : list
            The list of indices with nucleotide substitutions in the modified sequence
            This list must be the same length as `substitution_indices_in_normalized_seq`, as the deletions must correspond.

        substitution_indices_in_normalized_seq : list
            The list of indices with nucleotide substitutions in the normalized sequence
            This list must be the same length as `substitution_indices_in_modified_seq`, as the deletions must correspond.
        """

        assert len(deletion_indices_in_modified_seq) == len(deletion_indices_in_normalized_seq)
        assert len(substitution_indices_in_modified_seq) == len(substitution_indices_in_normalized_seq)

        # Get a nonredundant set of trimmed sequences already contained in the modified sequence.
        processed_trimmed_seq_names = set()
        for normalized_seq in self.all_normalized_seqs:
            for trimmed_seq in normalized_seq.trimmed_seqs:
                processed_trimmed_seq_names.add(trimmed_seq.representative_name)

        self.normalized_seqs_with_deletions.append(normalized_seq_with_deletion)
        self.all_normalized_seqs.append(normalized_seq_with_deletion)
        normalized_seq.modified_seq = self

        input_count = 0
        input_with_extra_fiveprime_count = 0
        input_seqs_mapped_without_extra_fiveprime_count = 0
        input_seqs_mapped_with_extra_fiveprime_count = 0
        for (trimmed_seq,
             trimmed_seq_start_in_normalized_seq,
             trimmed_seq_stop_in_normalized_seq) in zip(normalized_seq_with_deletion.trimmed_seqs,
                                                        normalized_seq_with_deletion.start_positions,
                                                        normalized_seq_with_deletion.end_positions):
            if trimmed_seq.representative_name in processed_trimmed_seq_names:
                continue

            trimmed_seq_input_count = trimmed_seq.input_count

            input_count += trimmed_seq_input_count

            if trimmed_seq.identification_method == 1: # 1 => mapped
                # TrimmedSeqs are comprised of EITHER profiled (0) OR mapped (1) UniqueSeqs.
                if trimmed_seq.unique_with_extra_fiveprime_count == 0:
                    input_seqs_mapped_without_extra_fiveprime_count += trimmed_seq_input_count
                else:
                    input_with_extra_fiveprime_count += trimmed_seq.input_with_extra_fiveprime_count
                    input_seqs_mapped_with_extra_fiveprime_count += trimmed_seq_input_count
            else:
                input_with_extra_fiveprime_count += trimmed_seq.input_with_extra_fiveprime_count

            for acceptor_variant_seq, acceptor_variant_input_count in trimmed_seq.input_acceptor_variant_count_dict:
                self.input_acceptor_variant_count_dict[acceptor_variant_seq] += acceptor_variant_input_count

            # Check if the trimmed sequence only occurs in the present normalized sequence or in others as well.
            if trimmed_seq.normalized_seq_count == 1:
                is_trimmed_seq_unique_to_normalized_seq = True
            else:
                is_trimmed_seq_unique_to_normalized_seq = False

            # Find the deletions covered by the present trimmed sequence.
            deletion_indices_in_normalized_seq_covered_by_trimmed_seq = []
            deletion_indices_in_modified_seq_covered_by_trimmed_seq = []
            for trimmed_seq_index, normalized_seq_index in enumerate(range(trimmed_seq_start_in_normalized_seq + 1, trimmed_seq_stop_in_normalized_seq)):
                try:
                    deletion_number_in_normalized_seq = deletion_indices_in_normalized_seq.index(normalized_seq_index)
                except ValueError:
                    continue
                if normalized_seq_index + 1 == trimmed_seq_stop_in_normalized_seq:
                    # The trimmed sequence only covers the 5' but not the 3' nucleotide flanking the deletion site,
                    # so the sequence cannot be said to cover the deletion.
                    continue
                deletion_indices_in_normalized_seq_covered_by_trimmed_seq.append(normalized_seq_index)
                deletion_indices_in_modified_seq_covered_by_trimmed_seq.append(deletion_indices_in_modified_seq[deletion_number_in_normalized_seq])

            # Augment the coverage of the deletions in the trimmed sequence.
            for deletion_index_in_modified_seq in deletion_indices_in_modified_seq_covered_by_trimmed_seq:
                # Add a record of this deletion in the modified sequence if it is encountered for the first time.
                try:
                    deletion_number_in_modified_seq = self.deletion_indices.index(deletion_index_in_modified_seq)
                except ValueError:
                    deletion_number_in_modified_seq = bisect(self.deletion_indices, deletion_index_in_modified_seq)
                    self.deletion_indices.insert(deletion_number_in_modified_seq, deletion_index_in_modified_seq)
                    self.deletion_coverages.insert(deletion_number_in_modified_seq, 0)
                    self.deletion_unique_coverages.insert(deletion_number_in_modified_seq, 0)

                self.deletion_coverages[deletion_number_in_modified_seq] += trimmed_seq_input_count
                if is_trimmed_seq_unique_to_normalized_seq:
                    self.deletion_unique_coverages[deletion_number_in_modified_seq] += trimmed_seq_input_count

            # Augment the coverage of undeleted nucleotides by the number of input sequences (reads) in the trimmed sequence.
            undeleted_indices = np.setdiff1d(np.array(range(len(trimmed_seq.seq_string)
                                                            + len(deletion_indices_in_modified_seq_covered_by_trimmed_seq)))
                                             + trimmed_seq_start_in_normalized_seq + normalized_seq_start_in_modified_seq,
                                             np.array(deletion_indices_in_modified_seq_covered_by_trimmed_seq))
            self.coverages[undeleted_indices] += trimmed_seq_input_count
            if is_trimmed_seq_unique_to_normalized_seq:
                self.unique_coverages[undeleted_indices] += trimmed_seq_input_count

            # Augment the coverage of the specific nucleotides from the trimmed sequence located at substitution sites.
            for substitution_number_in_normalized_seq, substitution_index_in_normalized_seq in enumerate(substitution_indices_in_modified_seq):
                substitution_index_in_modified_seq = substitution_indices_in_modified_seq[substitution_number_in_normalized_seq]

                # Check that the position of the substitution in the normalized sequence is covered by the present trimmed sequence.
                if not (trimmed_seq_start_in_normalized_seq <= substitution_index_in_normalized_seq < trimmed_seq_stop_in_normalized_seq):
                    continue

                nuc_int = nuc_int_dict[normalized_seq.seq_string[substitution_index_in_normalized_seq]]
                substitution_number_in_modified_seq = self.substitution_indices.index(substitution_index_in_modified_seq)
                self.substitution_coverages_array[substitution_number_in_modified_seq, nuc_int] += trimmed_seq_input_count
                if is_trimmed_seq_unique_to_normalized_seq:
                    self.substitution_unique_coverages_array[substitution_number_in_modified_seq, nuc_int] += trimmed_seq_input_count
        self.input_count += input_count
        self.input_with_extra_fiveprime_count = input_with_extra_fiveprime_count
        self.input_seqs_mapped_without_extra_fiveprime_count = input_seqs_mapped_without_extra_fiveprime_count
        self.input_seqs_mapped_with_extra_fiveprime_count = input_seqs_mapped_with_extra_fiveprime_count


    def set_consensus_seq_string(self):
        """Set a consensus sequence from the nucleotides with the highest coverage at each position.

        There is no point calling this method until all normalized sequences with substitutions and deletions have been added.
        """

        consensus_seq_string = self.all_normalized_seqs[0].seq_string
        for substitution_index, nuc_coverages in zip(self.substitution_indices, self.substitution_coverages_array):
            max_index = nuc_coverages.argmax()
            nuc_int = max_index + 1
            consensus_seq_string = (consensus_seq_string[: substitution_index]
                                    + int_nuc_dict[nuc_int]
                                    + consensus_seq_string[substitution_index + 1: ])
        self.consensus_seq_string = consensus_seq_string


class TRNASeqDataset:
    # Column headers for supplementary tables written to text files
    UNIQUED_NONTRNA_HEADER = ["representative_name", "input_count", "sequence"]
    TRIMMED_ENDS_HEADER = ["representative_name", "unique_name", "fiveprime_sequence", "threeprime_sequence", "input_seq_count"]
    NORMALIZED_FRAGMENTS_HEADER = ["representative_name", "trimmed_name", "start", "end"]


    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        """Class for processing a tRNA-seq dataset"""

        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # Argument group 1A: MANDATORY
        self.input_fasta_path = A('fasta_file')
        self.project_name = A('project_name')
        self.output_dir = os.path.abspath(A('output_dir')) if A('output_dir') else None

        # Argument group 1B: MUNDANE
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.description_file_path = os.path.abspath(A('description')) if A('description') else None

        # Argument group 1C: PERFORMANCE
        self.num_threads = A('num_threads')
        self.skip_fasta_check = A('skip_fasta_check')
        self.write_checkpoints = A('write_checkpoints')
        self.load_checkpoint = A('load_checkpoint')
        self.write_buffer_size = A('write_buffer_size')
        self.alignment_target_chunk_size = A('alignment_target_chunk_size')
        self.fragment_mapping_query_chunk_length = A('fragment_mapping_query_chunk_length')

        # Argument group 1D: ADVANCED
        self.feature_param_path = os.path.abspath(A('feature_param_file')) if A('feature_param_file') else None
        self.min_trna_fragment_size = A('min_trna_fragment_size')
        self.agglomeration_max_mismatch_freq = A('agglomeration_max_mismatch_freq')
        self.min_modification_count = A('min_modification_count')
        self.min_modification_fraction = A('min_modification_fraction')
        self.max_deletion_size = A('max_deletion_size')

        # Argument group 1E: MINUTIAE
        self.alignment_progress_interval = A('alignment_progress_interval')
        self.agglomeration_progress_interval = A('agglomeration_progress_interval')

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using --fasta-file or -f.")
        if not self.project_name:
            raise ConfigError("Please set a project name using --project-name or -n.")
        if not self.output_dir:
            raise ConfigError("Please provide an output directory using --output-dir or -o.")

        self.trnaseq_db_path = os.path.join(self.output_dir, self.project_name + "-TRNASEQ.db")

        # Supplementary text file paths
        self.uniqued_nontrna_path = os.path.join(self.output_dir, self.project_name + "-UNIQUED_NONTRNA.txt")
        self.uniqued_trna_path = os.path.join(self.output_dir, self.project_name + "-UNIQUED_TRNA.txt")
        self.trimmed_ends_path = os.path.join(self.output_dir, self.project_name + "-TRIMMED_ENDS.txt")
        self.normalized_fragments_path = os.path.join(self.output_dir, self.project_name + "-NORMALIZED_FRAGMENTS.txt")

        # Intermediate pickle file paths
        self.profile_unique_trna_seqs_path = os.path.join(self.output_dir, "UNIQUE_TRNA_SEQS-PROFILE_CHECKPOINT.pkl")
        self.profile_unique_nontrna_seqs_path = os.path.join(self.output_dir, "UNIQUE_NONTRNA_SEQS-PROFILE_CHECKPOINT.pkl")
        self.threeprime_normalization_unique_trna_seqs_path = os.path.join(self.output_dir, "UNIQUE_TRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.threeprime_normalization_unique_nontrna_seqs_path = os.path.join(self.output_dir, "UNIQUE_NONTRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.threeprime_normalization_trimmed_trna_seqs_path = os.path.join(self.output_dir, "TRIMMED_TRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.threeprime_normalization_normalized_trna_seqs_path = os.path.join(self.output_dir, "NORMALIZED_TRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.fragment_mapping_unique_trna_seqs_path = os.path.join(self.output_dir, "UNIQUE_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
        self.fragment_mapping_unique_nontrna_seqs_path = os.path.join(self.output_dir, "UNIQUE_NONTRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
        self.fragment_mapping_trimmed_trna_seqs_path = os.path.join(self.output_dir, "TRIMMED_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
        self.fragment_mapping_normalized_trna_seqs_path = os.path.join(self.output_dir, "NORMALIZED_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")

        self.unique_nontrna_seqs = []
        self.unique_trna_seqs = []
        self.trimmed_trna_seqs = []
        self.normalized_trna_seqs = []
        self.modified_trna_seqs = []

        self.counts_of_normalized_seqs_containing_trimmed_seqs = [] # same length as self.trimmed_trna_seqs
        # "Multiplicity" of a trimmed tRNA sequence
        # = number of normalized sequences that the sequence is in * number of input sequences represented by the trimmed sequence
        self.multiplicities_of_trimmed_seqs_among_normalized_seqs = [] # same length as self.trimmed_trna_seqs
        self.average_multiplicities_of_normalized_seqs = [] # same length as self.normalized_trna_seqs
        self.average_multiplicities_of_modified_seqs = []


    def sanity_check(self):
        """Check user inputs before proceeding."""

        if os.path.exists(self.output_dir):
            if len(os.listdir(self.output_dir)) == 0:
                # `anvi-run-workflow` creates the output directory even before `anvi-trnaseq` is called.
                pass
            elif self.overwrite_output_destinations:
                if self.load_checkpoint:
                    raise ConfigError("You cannot use `--load-checkpoint` in conjunction with `--overwrite-output-destinations`. "
                                      "Starting at a checkpoint requires loading intermediate files "
                                      "written to the output directory in a previous `anvi-trnaseq` run, "
                                      "but this directory would be removed with `--overwrite-output-destinations`.")
                shutil.rmtree(self.output_dir)
            else:
                if not self.load_checkpoint:
                    raise ConfigError("The directory that was specified by --output-dir or -o, %s, already exists. "
                                      "Use the flag --overwrite-output-destinations to overwrite this directory." % self.output_dir)
        missing_intermediate_files = False
        if self.load_checkpoint == 'profile':
            if (not os.path.exists(self.profile_unique_trna_seqs_path)
                or not os.path.exists(self.profile_unique_nontrna_seqs_path)):
                missing_intermediate_files = True
        elif self.load_checkpoint == 'threeprime_normalization':
            if (not os.path.exists(self.threeprime_normalization_unique_trna_seqs_path)
                or not os.path.exists(self.threeprime_normalization_unique_nontrna_seqs_path)
                or not os.path.exists(self.threeprime_normalization_trimmed_trna_seqs_path)
                or not os.path.exists(self.threeprime_normalization_normalized_trna_seqs_path)):
                missing_intermediate_files = True
        elif self.load_checkpoint == 'fragment_mapping':
            if (not os.path.exists(self.fragment_mapping_unique_trna_seqs_path)
                or not os.path.exists(self.fragment_mapping_unique_nontrna_seqs_path)
                or not os.path.exists(self.fragment_mapping_trimmed_trna_seqs_path)
                or not os.path.exists(self.fragment_mapping_normalized_trna_seqs_path)):
                missing_intermediate_files = True
        else:
            if not os.path.exists(self.output_dir):
                # `anvi-run-workflow` creates the output directory even before `anvi-trnaseq` is called.
                os.mkdir(self.output_dir)
        if missing_intermediate_files:
            raise ConfigError("Intermediate files needed for running `anvi-trnaseq` with `--load-checkpoint %s` are missing. "
                              "You should probably run `anvi-trnaseq` from the beginning without `--load-checkpoint`. "
                              "To generate necessary intermediate files for future use of `--load-checkpoint`, use the flag `--write-checkpoints`."
                              % self.load_checkpoint)
        filesnpaths.is_output_dir_writable(self.output_dir)

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = self.description_file_path.read()
        else:
            self.description = None
        self.run.info("Description", self.description_file_path if self.description_file_path else "No description given")

        if not 1 < self.num_threads < mp.cpu_count():
            ConfigError("The number of threads to use must be a positive integer "
                        "less than or equal to %d. Try again!" % mp.cpu_count())

        self.run.info("Input FASTA file", self.input_fasta_path)

        # UNCOMMENT
        # if not self.load_checkpoint:
        #     utils.check_fasta_id_uniqueness(self.input_fasta_path)

        if not self.skip_fasta_check:
            self.progress.new("Checking input FASTA defline format")
            self.progress.update("...")

            utils.check_fasta_id_formatting(self.input_fasta_path)

            self.progress.end()

            self.run.info_single("FASTA deflines were found to be anvi'o-compliant", mc='green')


    def create_trnaseq_db(self):
        """Create an empty tRNA-seq database."""

        meta_values = {'project_name': self.project_name,
                       'description': self.description if self.description else '_No description is provided_'}
        TRNASeqDatabase(self.trnaseq_db_path, quiet=False).create(meta_values)


    def get_unique_input_seqs(self):
        """Dereplicate input reads."""

        self.progress.new("Finding unique input sequences")
        self.progress.update("...")

        fasta = fastalib.SequenceSource(self.input_fasta_path)
        names = []
        seqs = []
        input_seq_count = 0
        while next(fasta):
            names.append(fasta.id)
            seqs.append(fasta.seq)
            input_seq_count += 1
        fasta.close()
        self.input_seq_count = input_seq_count

        clusters = Dereplicator(names, seqs, progress=self.progress).full_length_dereplicate()

        unique_input_seqs = [UniqueSeq(cluster.representative_seq_string,
                                       cluster.member_names[0],
                                       len(cluster.member_names)) for cluster in clusters]

        self.progress.end()

        return unique_input_seqs


    def profile_trna(self, unique_input_seqs):
        """Profile tRNA features in input sequences.

        Appends UniqueSeq objects representing profiled tRNA sequences to `self.unique_trna_seqs`
        Appends leftover UniqueSeq objects representing unprofiled tRNA sequences to `self.unique_nontrna_seqs`

        Parameters
        ==========
        unique_input_seqs : list
            List of UniqueSeq objects
        """

        self.progress.new("Profiling tRNA features in input sequences")
        self.progress.update("...")

        processed_read_count = 0
        processed_seq_count = 0
        trna_read_count = 0
        unique_trna_count = 0
        trna_containing_anticodon_read_count = 0
        full_length_trna_read_count = 0
        trna_with_one_to_three_extra_fiveprime_bases_read_count = 0
        trna_with_more_than_three_extra_fiveprime_bases_read_count = 0
        trna_with_extrapolated_fiveprime_feature_read_count = 0
        trna_with_threeprime_cca_read_count = 0
        trna_with_threeprime_cc_read_count = 0
        trna_with_threeprime_c_read_count = 0
        trna_with_threeprime_ccan_ccann_read_count = 0

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        processes = [mp.Process(target=profile_worker, args=(input_queue, output_queue))
                     for _ in range(self.num_threads)]
        for p in processes:
            p.start()

        write_point_iterator = iter([self.write_buffer_size * (i + 1)
                                     for i in range(len(unique_input_seqs) // self.write_buffer_size)]
                                    + [len(unique_input_seqs), None])
        write_point = next(write_point_iterator)
        fetched_profile_count = 0
        unique_input_seqs_to_write_dict = {}
        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        for unique_input_seq in unique_input_seqs:
            input_queue.put((unique_input_seq.seq_string, unique_input_seq.representative_name))
            unique_input_seqs_to_write_dict[unique_input_seq.representative_name] = unique_input_seq
            processed_read_count += unique_input_seq.input_count
            processed_seq_count += 1

            if processed_seq_count == write_point:
                # Write a chunk of sequence results.
                write_point = next(write_point_iterator)

                # List of entries for each tRNA-seq database table
                trnaseq_sequences_table_entries = []
                trnaseq_info_table_entries = []
                trnaseq_features_table_entries = []
                trnaseq_unconserved_table_entries = []
                trnaseq_unpaired_table_entries = []

                while fetched_profile_count < processed_seq_count:
                    trna_profile = output_queue.get()
                    fetched_profile_count += 1
                    output_name = trna_profile.name
                    output_seq = trna_profile.input_seq

                    unique_seq = unique_input_seqs_to_write_dict.pop(output_name)

                    if not trna_profile.is_predicted_trna:
                        self.unique_nontrna_seqs.append(unique_seq)
                        continue

                    unique_seq.identification_method = 0 # 0 => profiled
                    unique_seq.acceptor_length = len(trna_profile.acceptor_variant_string)
                    unique_seq.has_complete_feature_set = trna_profile.has_complete_feature_set
                    unique_seq.extra_fiveprime_length = trna_profile.num_extra_fiveprime
                    self.unique_trna_seqs.append(unique_seq)

                    unique_trna_count += 1
                    output_seq_length = len(output_seq)

                    num_replicates = unique_seq.input_count
                    trna_read_count += num_replicates

                    # Recover nucleotides that did not fit expectation,
                    # either by not being the expected nucleotide or type of nucleotide
                    # or by not base pairing in a stem.
                    unconserved_info = trna_profile.get_unconserved_positions()
                    unpaired_info = trna_profile.get_unpaired_positions()

                    if trna_profile.anticodon_seq:
                        trna_containing_anticodon_read_count += num_replicates
                    if trna_profile.has_complete_feature_set:
                        full_length_trna_read_count += num_replicates
                    if trna_profile.num_in_extrapolated_fiveprime_feature > 0:
                        trna_with_extrapolated_fiveprime_feature_read_count += num_replicates

                    if trna_profile.num_extra_threeprime > 0:
                        trna_with_threeprime_ccan_ccann_read_count += num_replicates
                    elif trna_profile.acceptor_variant_string == 'CCA':
                        trna_with_threeprime_cca_read_count += num_replicates
                    elif trna_profile.acceptor_variant_string == 'CC':
                        trna_with_threeprime_cc_read_count += num_replicates
                    elif trna_profile.acceptor_variant_string == 'C':
                        trna_with_threeprime_c_read_count += num_replicates

                    if trna_profile.num_extra_fiveprime > 3:
                        trna_with_more_than_three_extra_fiveprime_bases_read_count += num_replicates
                    elif trna_profile.num_extra_fiveprime > 0:
                        trna_with_one_to_three_extra_fiveprime_bases_read_count += num_replicates

                    trnaseq_sequences_table_entries.append((output_name, num_replicates, output_seq))

                    # The alpha and beta regions of the D loop vary in length.
                    # Record their start and stop positions in the sequence if they were profiled.
                    # These are included in the info table rather than the feature table,
                    # because they are subfeatures of the D loop feature,
                    # and the positions of the features in the feature table are not overlapping.
                    alpha_start = trna_profile.alpha_start if trna_profile.alpha_start else '??'
                    alpha_stop = trna_profile.alpha_stop - 1 if trna_profile.alpha_stop else '??'
                    beta_start = trna_profile.beta_start if trna_profile.beta_start else '??'
                    beta_stop = trna_profile.beta_stop - 1 if trna_profile.beta_stop else '??'

                    trnaseq_info_table_entries.append(
                        (output_name,
                         trna_profile.has_complete_feature_set,
                         trna_profile.anticodon_seq,
                         trna_profile.anticodon_aa,
                         output_seq_length,
                         # Zero-based start index of identified tRNA features within the input sequence.
                         output_seq_length - len(trna_profile.profiled_seq),
                         # Stop index of features (real stop position, not Pythonic stop index for slicing).
                         output_seq_length - trna_profile.num_extra_threeprime - 1,
                         trna_profile.num_conserved,
                         trna_profile.num_unconserved,
                         trna_profile.num_paired,
                         trna_profile.num_unpaired,
                         trna_profile.num_in_extrapolated_fiveprime_feature,
                         trna_profile.num_extra_fiveprime,
                         trna_profile.num_extra_threeprime,
                         alpha_start,
                         alpha_stop,
                         beta_start,
                         beta_stop))

                    trnaseq_features_table_entries.append(
                        (output_name, )
                        # When tRNA features were not found at the 5' end of the read,
                        # their start and stop indices also were not found.
                        + tuple(['?' * 2 for _ in range((len(TRNA_FEATURE_NAMES) - len(trna_profile.features)))]) * 2
                        + tuple(itertools.chain(*zip(
                            [str(feature.start_index) if hasattr(feature, 'start_index')
                             else ','.join(map(str, feature.start_indices))
                             for feature in trna_profile.features],
                            # Convert Pythonic stop index for slicing to real stop position of feature.
                            [str(feature.stop_index - 1) if hasattr(feature, 'stop_index')
                             else ','.join(map(str, [stop_index - 1 for stop_index in feature.stop_indices]))
                             for feature in trna_profile.features]))))

                    for unconserved_tuple in unconserved_info:
                        trnaseq_unconserved_table_entries.append((output_name, ) + unconserved_tuple)

                    for unpaired_tuple in unpaired_info:
                        trnaseq_unpaired_table_entries.append((output_name, ) + unpaired_tuple)

                if len(trnaseq_sequences_table_entries) > 0:
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('sequences', ','.join('?' * len(tables.trnaseq_sequences_table_structure))),
                                             trnaseq_sequences_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('basic_info', ','.join('?' * len(tables.trnaseq_info_table_structure))),
                                             trnaseq_info_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('features', ','.join('?' * len(tables.trnaseq_features_table_structure))),
                                             trnaseq_features_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('unconserved_nucleotides', ','.join('?' * len(tables.trnaseq_unconserved_table_structure))),
                                             trnaseq_unconserved_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('unpaired_nucleotides', ','.join('?' * len(tables.trnaseq_unpaired_table_structure))),
                                             trnaseq_unpaired_table_entries)

                self.progress.update("%d of %d unique sequences have been profiled"
                                     % (fetched_profile_count, len(unique_input_seqs)))

        for p in processes:
            p.terminate()
            p.join()

        # Profiled seqs were added to the output queue as they were processed, so sort by name.
        self.unique_trna_seqs.sort(key=lambda unique_seq: unique_seq.representative_name)
        self.unique_nontrna_seqs.sort(key=lambda unique_seq: unique_seq.representative_name)

        trnaseq_db.db.set_meta_value('input_reads_processed', processed_read_count)
        trnaseq_db.db.set_meta_value('unique_input_seqs_processed', processed_seq_count)
        trnaseq_db.db.set_meta_value('trna_reads', trna_read_count)
        trnaseq_db.db.set_meta_value('unique_trna_seqs', unique_trna_count)
        trnaseq_db.db.set_meta_value('trna_reads_containing_anticodon', trna_containing_anticodon_read_count)
        trnaseq_db.db.set_meta_value('full_length_trna_reads', full_length_trna_read_count)
        trnaseq_db.db.set_meta_value('trna_with_one_to_three_extra_fiveprime_bases', trna_with_one_to_three_extra_fiveprime_bases_read_count)
        trnaseq_db.db.set_meta_value('trna_with_more_than_three_extra_fiveprime_bases', trna_with_more_than_three_extra_fiveprime_bases_read_count)
        trnaseq_db.db.set_meta_value('trna_reads_with_extrapolated_fiveprime_feature', trna_with_extrapolated_fiveprime_feature_read_count)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_cca', trna_with_threeprime_cca_read_count)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_cc', trna_with_threeprime_cc_read_count)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_c', trna_with_threeprime_c_read_count)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_ccan_ccann', trna_with_threeprime_ccan_ccann_read_count)
        trnaseq_db.disconnect()

        self.progress.end()

        self.run.info("Reads processed", processed_read_count)
        self.run.info("Unique sequences processed", processed_seq_count)
        self.run.info("Reads profiled as tRNA", trna_read_count)
        self.run.info("Unique profiled tRNA sequences", unique_trna_count)
        self.run.info("Profiled reads with anticodon", trna_containing_anticodon_read_count)
        self.run.info("Profiled reads spanning acceptor stem", full_length_trna_read_count)
        self.run.info("Profiled reads with 1-3 extra 5' bases", trna_with_one_to_three_extra_fiveprime_bases_read_count)
        self.run.info("Profiled reads with >3 extra 5' bases", trna_with_more_than_three_extra_fiveprime_bases_read_count)
        self.run.info("Profiled reads with extrapolated 5' feature", trna_with_extrapolated_fiveprime_feature_read_count)
        self.run.info("Profiled reads ending in 3'-CCA", trna_with_threeprime_cca_read_count)
        self.run.info("Profiled reads ending in 3'-CC", trna_with_threeprime_cc_read_count)
        self.run.info("Profiled reads ending in 3'-C", trna_with_threeprime_c_read_count)
        self.run.info("Profiled reads ending in 3'-CCAN/CCANN", trna_with_threeprime_ccan_ccann_read_count)


    def trim_ends(self, unique_trna_seqs):
        """Trim any nucleotides 5' of the acceptor stem and 3' of the discriminator.

        Appends TrimmedSeq objects formed from input UniqueSeq objects to `self.trimmed_trna_seqs`

        Parameters
        ==========
        unique_trna_seqs : list
            List of UniqueSeq objects
        """

        self.progress.new("Trimming the 3' and 5' ends of sequences")
        self.progress.update("...")

        representative_names = [unique_seq.representative_name for unique_seq in unique_trna_seqs]
        trimmed_seq_strings = [unique_seq.seq_string[unique_seq.extra_fiveprime_length: len(unique_seq.seq_string) - unique_seq.acceptor_length]
                               for unique_seq in unique_trna_seqs]

        clusters = Dereplicator(representative_names,
                                trimmed_seq_strings,
                                extras=unique_trna_seqs,
                                progress=self.progress).full_length_dereplicate()

        trimmed_seqs = [TrimmedSeq(cluster.representative_seq_string, cluster.member_extras) for cluster in clusters]

        self.trimmed_trna_seqs.extend(trimmed_seqs)

        self.trimmed_trna_seqs.sort(key=lambda trimmed_seq: trimmed_seq.representative_name)

        self.progress.end()


    def dereplicate_threeprime(self):
        """Dereplicate trimmed tRNA sequences from the 3' end of longer trimmed sequences.

        EXAMPLE:
        normalized tRNA (trimmed tRNA 1): TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        trimmed tRNA 2                  :                       AATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        trimmed tRNA 3                  :                                     GCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        """

        self.progress.new("Dereplicating trimmed tRNA sequences from the 3' end")
        self.progress.update("...")

        representative_names = [trimmed_seq.representative_name for trimmed_seq in self.trimmed_trna_seqs]
        # Reverse sequence orientation to dereplicate from the 3' end.
        reversed_seq_strings = [trimmed_seq.seq_string[::-1] for trimmed_seq in self.trimmed_trna_seqs]
        clusters = Dereplicator(representative_names,
                                reversed_seq_strings,
                                extras=self.trimmed_trna_seqs,
                                progress=self.progress).prefix_dereplicate()

        # Skip initialization of NormalizedSeq objects,
        # as additional TrimmedSeq constituents are added after mapping unprofiled fragments to profiled tRNA.
        self.normalized_trna_seqs = [NormalizedSeq(cluster.member_extras, skip_init=True) for cluster in clusters]

        self.progress.end()


    def map_fragments(self):
        """Map unprofiled tRNA fragments to longer profiled tRNA sequences.

        If the specified minimum fragment length is shorter than the minimum length of profiled tRNA fragments --
        which does not happen with the default settings, as both lengths are 25,
        but can happen when the user adjusts `--min-trna-fragment-size` downward --
        then mapped fragments may occur at the 3' end of a tRNA, but will not include 3' acceptor variants (CCA, CC, C, etc.),
        as these 3' extensions were trimmed off the mapping targets.

        EXAMPLE:
        normalized tRNA:                 (GT)TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        mapped tRNA 1 (extra 5' bases) :   T TCCGTGATAGTTTAATGGTCAGAATGG
        mapped tRNA 2 (interior)  :                 TAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGG
        """

        self.progress.new("Mapping unprofiled tRNA fragments to profiled tRNA")

        self.progress.update("Retrieving queries from unprofiled sequences")
        query_length_intervals = []
        max_query_length = max(map(len, [seq.seq_string for seq in self.unique_nontrna_seqs]))
        # By default, avoid sequences shorter than 25 nucleotides, the minimum length of a profiled 3' fragment of tRNA.
        interval_start_length = self.min_trna_fragment_size
        interval_stop_length = interval_start_length + self.fragment_mapping_query_chunk_length
        if interval_stop_length + self.fragment_mapping_query_chunk_length > max_query_length:
            # Group the last chunk in with the second to last chunk if it is less than full size.
            interval_stop_length = max_query_length + 1
        query_length_intervals.append((interval_start_length, interval_stop_length))
        query_names = []
        query_seqs = []
        query_name_chunks = [query_names]
        query_seq_chunks = [query_seqs]
        for nontrna_index, seq in sorted([t for t in zip(range(len(self.unique_nontrna_seqs)), self.unique_nontrna_seqs)
                                          if len(t[1].seq_string) >= self.min_trna_fragment_size],
                                         key=lambda t: len(t[1].seq_string)):
            if len(seq.seq_string) < interval_stop_length:
                query_names.append((seq.representative_name, nontrna_index))
                query_seqs.append(seq.seq_string)
            else:
                interval_start_length = interval_stop_length
                interval_stop_length = interval_stop_length + self.fragment_mapping_query_chunk_length
                if interval_stop_length + self.fragment_mapping_query_chunk_length > max_query_length:
                    interval_stop_length = max_query_length + 1
                query_length_intervals.append((interval_start_length, interval_stop_length))
                query_names = [(seq.representative_name, nontrna_index)]
                query_seqs = [seq.seq_string]
                query_name_chunks.append(query_names)
                query_seq_chunks.append(query_seqs)

        self.progress.update("Retrieving targets from profiled tRNA")
        # Leftover non-tRNA sequences are mapped to normalized tRNA sequences
        # with extra 5' bases added when present in underlying unique tRNA sequences.
        # Multiple targets for each normalized sequence are therefore produced for different 5' sequence extensions.
        target_names = []
        target_seqs = []
        for normalized_seq_index, normalized_seq in enumerate(self.normalized_trna_seqs):
            normalized_name = normalized_seq.representative_name
            normalized_seq_string = normalized_seq.seq_string
            # The longest trimmed sequence (the first in the list) is by design
            # the only one of the profiled trimmed sequences forming the normalized sequence that may have extra 5' bases.
            longest_trimmed_seq = normalized_seq.trimmed_seqs[0]
            if longest_trimmed_seq.unique_with_extra_fiveprime_count > 0:
                fiveprime_seq_string_set = set()
                for unique_seq in longest_trimmed_seq.unique_seqs:
                    if unique_seq.extra_fiveprime_length > 0:
                        fiveprime_seq_string_set.add(unique_seq.seq_string[: unique_seq.extra_fiveprime_length])

                # Avoid creating superfluous target sequences that are subsequences of other target sequences
                # due to a 5' extension of a normalized sequence being a subsequence of a longer 5' extension of the same normalized sequence.
                fiveprime_seq_strings = sorted(fiveprime_seq_string_set, key=lambda s: -len(s))
                fiveprime_seq_string_additions = [fiveprime_seq_strings[0]]
                for fiveprime_seq_string in fiveprime_seq_strings[1: ]:
                    fiveprime_seq_string_length = len(fiveprime_seq_string)
                    for fiveprime_seq_string_addition in fiveprime_seq_string_additions:
                        if fiveprime_seq_string == fiveprime_seq_string_addition[-fiveprime_seq_string_length: ]:
                            break
                    else:
                        fiveprime_seq_string_additions.append(fiveprime_seq_string)

                for fiveprime_index, fiveprime_seq_string in enumerate(fiveprime_seq_strings):
                    # Use an index to distinguish otherwise equivalent targets with different 5' extensions of the same length.
                    target_names.append((normalized_seq_index, len(fiveprime_seq_string), fiveprime_index))
                    target_seqs.append(fiveprime_seq_string + normalized_seq_string)
            else:
                target_names.append((normalized_seq_index, 0, 0)) # no extra 5' bases
                target_seqs.append(normalized_seq_string)

        self.progress.end()


        interval_index = 0
        nontrna_indices = []
        for query_names, query_seqs in zip(query_name_chunks, query_seq_chunks):
            self.progress.new("Mapping %d unprofiled sequences of length %d-%d to profiled tRNA"
                              % (len(query_names), query_length_intervals[interval_index][0], query_length_intervals[interval_index][1] - 1))

            aligned_query_dict, aligned_target_dict = Aligner(query_names, # aligned_target_dict is not used for anything
                                                              query_seqs,
                                                              target_names,
                                                              target_seqs,
                                                              num_threads=self.num_threads,
                                                              progress=self.progress).align(max_mismatch_freq=0,
                                                                                            target_chunk_size=self.alignment_target_chunk_size,
                                                                                            query_progress_interval=self.alignment_progress_interval)
            del aligned_target_dict
            gc.collect()

            self.progress.update("Processing alignments")

            for query_name, aligned_query in aligned_query_dict.items():
                if len(aligned_query.alignments) == 0:
                    continue

                nontrna_index = query_name[1]

                trimmed_seq = None

                for alignment in aligned_query.alignments:
                    reference_alignment_start = alignment.target_start
                    reference_alignment_end = alignment.target_start + alignment.alignment_length

                    normalized_seq_index, reference_fiveprime_length, _ = alignment.aligned_target.name # extra 5' index doesn't matter now

                    normalized_end_position = reference_alignment_end - reference_fiveprime_length
                    if normalized_end_position < 0:
                        # Ignore queries that align entirely to extra 5' bases.
                        continue

                    normalized_seq = self.normalized_trna_seqs[normalized_seq_index]

                    if not trimmed_seq:
                        nontrna_indices.append(nontrna_index)

                        unique_mapped_seq = self.unique_nontrna_seqs[nontrna_index]
                        unique_mapped_seq.identification_method = 1 # 1 => mapped
                        unique_mapped_seq.acceptor_length = 0
                        unique_mapped_seq.has_complete_feature_set = False

                        # Assume that 5' extensions are the same for the query regardless of the reference.
                        # This could be false in the unlikely cases of
                        # 1. tRNA profiling erroneously identifying the end of the acceptor stem
                        # or 2. the query mapping to different places at the end of the acceptor stem in different tRNAs.
                        if reference_fiveprime_length - reference_alignment_start > 0:
                            unique_mapped_seq.extra_fiveprime_length = reference_fiveprime_length - reference_alignment_start
                            normalized_start_position = 0
                        else:
                            unique_mapped_seq.extra_fiveprime_length = 0
                            normalized_start_position = reference_alignment_start - reference_fiveprime_length

                        trimmed_seq = TrimmedSeq(unique_mapped_seq.seq_string[unique_mapped_seq.extra_fiveprime_length: ],
                                                 [unique_mapped_seq])
                        self.trimmed_trna_seqs.append(trimmed_seq)

                        normalized_seq.trimmed_seqs.append(trimmed_seq)
                        trimmed_seq.normalized_seq_count += 1
                        normalized_seq.start_positions.append(normalized_start_position)
                        normalized_seq.end_positions.append(normalized_end_position)
                    else:
                        for prev_trimmed_seq in normalized_seq.trimmed_seqs[::-1]:
                            # Ensure that the trimmed sequence maps to the normalized sequence only once.
                            # Multiple targets can be created from the same normalized sequence for different 5' extensions.
                            if prev_trimmed_seq.identification_method == 0:
                                normalized_seq.trimmed_seqs.append(trimmed_seq)
                                trimmed_seq.normalized_seq_count += 1
                                normalized_seq.start_positions.append(normalized_start_position)
                                normalized_seq.end_positions.append(normalized_end_position)
                                break
                            if trimmed_seq.representative_name == prev_trimmed_seq.representative_name:
                                break
            interval_index += 1

            del aligned_query_dict
            gc.collect()

            self.progress.end()

        for normalized_seq in self.normalized_trna_seqs:
            normalized_seq.init()

        interior_mapped_count = 0
        fiveprime_mapped_count = 0
        for nontrna_index in sorted(nontrna_indices, reverse=True):
            unique_mapped_seq = self.unique_nontrna_seqs.pop(nontrna_index)

            if unique_mapped_seq.extra_fiveprime_length > 0:
                fiveprime_mapped_count += unique_mapped_seq.input_count
            else:
                interior_mapped_count += unique_mapped_seq.input_count

        self.run.info("Mapped reads without extra 5' tRNA bases", interior_mapped_count)
        self.run.info("Mapped reads with extra 5' tRNA bases", fiveprime_mapped_count)


    def find_modifications(self):
        self.progress.new("Finding modifications")

        # Cluster normalized tRNA sequences.
        # Clusters agglomerate sequences that differ from at least one other sequence in the cluster
        # by no more than 2 substitutions per 71 aligned positions (by default) in a gapless end-to-end alignment.
        agglomerator = Agglomerator([seq.representative_name for seq in self.normalized_trna_seqs],
                                    [seq.seq_string for seq in self.normalized_trna_seqs],
                                    num_threads=self.num_threads,
                                    progress=self.progress)
        # Provide a priority function for seeding clusters
        # that favors fully profiled tRNA over "longer" tRNA without a full set of profiled features.
        # Such incompletely profiled longer tRNA includes tRNA-tRNA chimeras,
        # and some of these have a long 5' section that is a long 3' fragment of tRNA,
        # which can cause other shorter normalized sequences to agglomerate by aligning to the 5' section of the chimera.
        full_length_trna_dict = {seq.representative_name: seq.has_complete_feature_set
                                 for seq in self.normalized_trna_seqs}
        agglomerator.agglomerate(max_mismatch_freq=self.agglomeration_max_mismatch_freq,
                                 priority_function=lambda aligned_reference: (-full_length_trna_dict[aligned_reference.name],
                                                                              -len(aligned_reference.seq_string),
                                                                              -len(aligned_reference.alignments),
                                                                              aligned_reference.name),
                                 alignment_target_chunk_size=self.alignment_target_chunk_size,
                                 alignment_progress_interval=self.alignment_progress_interval,
                                 agglomeration_progress_interval=self.agglomeration_progress_interval)

        agglomerated_aligned_reference_dict = agglomerator.agglomerated_aligned_reference_dict

        self.progress.update("Separating modification-induced substitutions from \"inter-strain\" variants")

        normalized_seq_dict = {seq.representative_name: seq for seq in self.normalized_trna_seqs}
        num_nuc_bins = len(unambiguous_nucs) + 1 # There is one bin (value 0) for end gaps in the alignment.
        names_of_normalized_seqs_assigned_to_modified_seqs = []
        for reference_name, aligned_reference in agglomerated_aligned_reference_dict.items():
            # A modification requires at least 3 different nucleotides to be detected,
            # and each normalized sequence differs by at least 1 nucleotide (substitution or gap),
            # so for a cluster to form a modified sequence, it must contain at least 3 normalized sequences.
            if len(aligned_reference.alignments) < 2:
                continue

            aligned_reference_length = len(aligned_reference.seq_string)

            valid_aligned_queries = []
            for alignment in aligned_reference.alignments:
                # Normalized tRNA sequences should only align at the 3' end.
                # Alignments to the interior of the sequence can occur when the reference is a tRNA-tRNA chimera.
                if aligned_reference_length != alignment.target_start + alignment.alignment_length:
                    continue

                query_name = alignment.aligned_query.name
                # The normalized sequence query may have agglomerated with another reference as well.
                # If the query formed a modified sequence,
                # it would form the same modified sequence when starting with this agglomeration.
                if query_name in names_of_normalized_seqs_assigned_to_modified_seqs:
                    continue

                valid_aligned_queries.append(normalized_seq_dict[query_name])

            # Confirm that at > 1 query passed the filters,
            # so at least 3 normalized sequences are still in the cluster.
            if len(valid_aligned_queries) < 2:
                continue

            sequence_array = np.zeros((len(aligned_reference.alignments) + 1, aligned_reference_length), dtype=int)
            # Rather than using the ASCII representation of each character,
            # which saves some time in converting the sequence string to a numpy array,
            # constrain the integer representation to the smallest possible range of integers
            # to speed up the bincount method used to determine the number of unique nucleotides at an alignment position.
            sequence_array[0, :] += [nuc_int_dict[nuc] for nuc in aligned_reference.seq_string]
            for i, aligned_query in enumerate(valid_aligned_queries, start=1):
                sequence_array[i, aligned_reference_length - len(aligned_query.seq_string): ] += [nuc_int_dict[nuc] for nuc in aligned_query.seq_string]

            normalized_seqs = np.array([normalized_seq_dict[reference_name]] + valid_aligned_queries)

            # Find positions in the alignment with nucleotide variability.
            alignment_position_unique_nuc_counts = (np.bincount((sequence_array + np.arange(aligned_reference_length, dtype=int) * num_nuc_bins).ravel(),
                                                                minlength=aligned_reference_length * num_nuc_bins).reshape(-1, num_nuc_bins)[:, 1:] != 0).sum(axis=1)
            three_four_nuc_alignment_positions = (alignment_position_unique_nuc_counts > 2).nonzero()[0]

            # Modification sites must have at least 3 nucleotides.
            if not three_four_nuc_alignment_positions.size:
                continue

            two_nuc_alignment_positions = (alignment_position_unique_nuc_counts == 2).nonzero()[0]
            clusters = deque(((sequence_array,
                               normalized_seqs,
                               three_four_nuc_alignment_positions), ))
            for alignment_position in two_nuc_alignment_positions:
                next_clusters = deque() # Make a new object with each iteration rather than clearing it.

                while clusters:
                    sequence_array, normalized_seqs, three_four_nuc_alignment_positions = clusters.pop()

                    # A modification requires at least 3 different nucleotides to be detected,
                    # and each normalized sequence differs by at least 1 nucleotide (substitution or gap),
                    # so for a cluster to form a modified sequence, it must contain at least 3 normalized sequences.
                    if normalized_seqs.size < 3:
                        continue

                    aligned_nucs = sequence_array[:, alignment_position]
                    nuc_counts = np.bincount(aligned_nucs, minlength=num_nuc_bins)[1: ]

                    if (nuc_counts != 0).sum() < 2:
                        # There are now < 2 nucleotides at the alignment position in the (derived) cluster under consideration.
                        # 2 different nucleotides are needed to distinguish single nucleotide variants.
                        next_clusters.appendleft((sequence_array,
                                                  normalized_seqs,
                                                  three_four_nuc_alignment_positions))
                        continue

                    # Add a new cluster for each nucleotide variant to the stack of clusters to process
                    # if the new cluster contains at least 3 sequences.
                    represented_nucs = nuc_counts.nonzero()[0] + 1
                    for nuc in represented_nucs:
                        split_cluster_seq_indices = (aligned_nucs == nuc).nonzero()[0]
                        if split_cluster_seq_indices.size > 2:
                            next_clusters.appendleft((sequence_array[split_cluster_seq_indices, :],
                                                      normalized_seqs[split_cluster_seq_indices],
                                                      three_four_nuc_alignment_positions))
                if next_clusters:
                    clusters = next_clusters
                else:
                    break
            if not clusters:
                continue

            # Check alignment positions previously found to have 3-4 nucleotides.
            # Further split (derived) clusters when positions now have 2 nucleotides.
            next_clusters = deque()
            while clusters:
                sequence_array, normalized_seqs, three_four_nuc_alignment_positions = clusters.pop()
                candidates_to_remove = []

                for i, alignment_position in enumerate(three_four_nuc_alignment_positions):
                    aligned_nucs = sequence_array[:, alignment_position]
                    nuc_counts = np.bincount(aligned_nucs, minlength=num_nuc_bins)[1: ]
                    # At least 3 different nucleotides are needed at a position to predict a modification.
                    represented_nucs = nuc_counts.nonzero()[0] + 1
                    if represented_nucs.size < 2:
                        candidates_to_remove.append(i)
                    elif represented_nucs.size == 2:
                        candidates_to_remove.append(i)
                        split_three_four_nuc_alignment_positions = np.delete(three_four_nuc_alignment_positions, candidates_to_remove)
                        for nuc in represented_nucs:
                            split_cluster_seq_indices = (aligned_nucs == nuc).nonzero()[0]
                            # At least 3 normalized sequences are needed to form a modified sequence.
                            if split_cluster_seq_indices.size > 2:
                                clusters.appendleft((sequence_array[split_cluster_seq_indices, :],
                                                     normalized_seqs[split_cluster_seq_indices],
                                                     split_three_four_nuc_alignment_positions))
                        # Reevaluate previous alignment positions in the split clusters.
                        break
                else:
                    # At least 1 position was discounted as no longer having 3-4 different nucleotides,
                    # but these positions had fewer than 2 nucleotides,
                    # and so did not cause the cluster to be split into new clusters.
                    # Therefore, do not cycle through the remaining positions again to find those with fewer than 3 nucleotides.
                    if candidates_to_remove:
                        next_clusters.appendleft((sequence_array,
                                                  normalized_seqs,
                                                  np.delete(three_four_nuc_alignment_positions, candidates_to_remove)))
                    else:
                        next_clusters.appendleft((sequence_array,
                                                  normalized_seqs,
                                                  three_four_nuc_alignment_positions))
            if next_clusters:
                clusters.clear()
                while next_clusters:
                    sequence_array, normalized_seqs, modification_positions = next_clusters.pop()

                    unique_coverage_array = np.zeros(sequence_array.shape, dtype=int)
                    array_width = sequence_array.shape[1]
                    added_and_rejected_trimmed_seq_names = []
                    for n, normalized_seq in enumerate(normalized_seqs):
                        normalized_seq_start_in_array = array_width - len(normalized_seq.seq_string)
                        for trimmed_seq, trimmed_seq_start_in_normalized_seq, trimmed_seq_stop_in_normalized_seq in zip(normalized_seq.trimmed_seqs, normalized_seq.start_positions, normalized_seq.end_positions):
                            if trimmed_seq.representative_name in added_and_rejected_trimmed_seq_names:
                                continue

                            if trimmed_seq.normalized_seq_count == 1:
                                unique_coverage_array[n, normalized_seq_start_in_array + trimmed_seq_start_in_normalized_seq: normalized_seq_start_in_array + trimmed_seq_stop_in_normalized_seq] += trimmed_seq.input_count
                                added_and_rejected_trimmed_seq_names.append(trimmed_seq.representative_name)
                                continue

                            num_normalized_seqs_in_cluster_containing_trimmed_seq = 1
                            for p, other_normalized_seq in enumerate(normalized_seqs):
                                if n == p:
                                    continue
                                for other_trimmed_seq in other_normalized_seq.trimmed_seqs:
                                    if trimmed_seq.representative_name == other_trimmed_seq.representative_name:
                                        num_normalized_seqs_in_cluster_containing_trimmed_seq += 1
                            if num_normalized_seqs_in_cluster_containing_trimmed_seq == trimmed_seq.normalized_seq_count:
                                unique_coverage_array[n, normalized_seq_start_in_array + trimmed_seq_start_in_normalized_seq: normalized_seq_start_in_array + trimmed_seq_stop_in_normalized_seq] += trimmed_seq.input_count
                                added_and_rejected_trimmed_seq_names.append(trimmed_seq.representative_name)
                            elif num_normalized_seqs_in_cluster_containing_trimmed_seq < trimmed_seq.normalized_seq_count:
                                added_and_rejected_trimmed_seq_names.append(trimmed_seq.representative_name)
                            else:
                                raise ConfigError("The number of normalized sequences containing the trimmed sequence was somehow miscalculated.")

                    clusters.appendleft((sequence_array,
                                         unique_coverage_array,
                                         normalized_seqs,
                                         modification_positions))
            else:
                continue

            # Check that nucleotide positions with at least 3 different nucleotides meet the criteria for modifications.
            while clusters:
                sequence_array, unique_coverage_array, normalized_seqs, modification_positions = clusters.pop()
                candidates_to_remove = []

                for i, alignment_position in enumerate(modification_positions):
                    if unique_coverage_array[:, alignment_position].sum() < self.min_modification_count:
                        # There is a minimum coverage threshold at the position required to predict a modification.
                        candidates_to_remove.append(i)
                        continue

                    aligned_nucs = sequence_array[:, alignment_position]
                    nuc_counts = np.bincount(aligned_nucs, minlength=num_nuc_bins)[1: ]
                    max_coverage = 0
                    total_coverage = 0
                    for nuc_int, nuc_count in enumerate(nuc_counts, start=1):
                        if nuc_count > 0:
                            nuc_coverage = unique_coverage_array[(aligned_nucs == nuc_int).nonzero()[0], alignment_position].sum()
                            if nuc_coverage > max_coverage:
                                max_coverage = nuc_coverage
                            total_coverage += nuc_coverage
                    if 1 - max_coverage / total_coverage < self.min_modification_fraction:
                        # There is a minimum fraction of minority nucleotide coverage required to predict a modification.
                        candidates_to_remove.append(i)
                        continue

                if not candidates_to_remove:
                    normalized_seqs = sorted([seq for seq in normalized_seqs], key=lambda seq: (-len(seq.seq_string), seq.representative_name))
                    start = aligned_reference_length - len(normalized_seqs[0].seq_string)
                    modification_positions -= start
                    modified_seq = ModifiedSeq(normalized_seqs, modification_positions.tolist())
                    for normalized_seq in normalized_seqs:
                        names_of_normalized_seqs_assigned_to_modified_seqs.append(normalized_seq.representative_name)
                    self.modified_trna_seqs.append(modified_seq)
                elif len(candidates_to_remove) < modification_positions.size:
                    normalized_seqs = sorted([seq for seq in normalized_seqs], key=lambda seq: (-len(seq.seq_string), seq.representative_name))
                    start = aligned_reference_length - len(normalized_seqs[0].seq_string)
                    modification_positions = np.delete(modification_positions, candidates_to_remove)
                    modification_positions -= start
                    modified_seq = ModifiedSeq(normalized_seqs, modification_positions.tolist())
                    for normalized_seq in normalized_seqs:
                        names_of_normalized_seqs_assigned_to_modified_seqs.append(normalized_seq.representative_name)
                    self.modified_trna_seqs.append(modified_seq)

        for i, modified_seq in enumerate(self.modified_trna_seqs[::-1]):
            modified_seq.set_consensus_seq_string()

        self.progress.end()


    def calc_normalization_stats(self):
        self.progress.new("Calculating normalized tRNA stats")
        self.progress.update("...")

        # Count the normalized sequences containing each trimmed sequence.
        normalized_count_dict = OrderedDict([(trimmed_seq.representative_name, 0) for trimmed_seq in self.trimmed_trna_seqs])
        for normalized_seq in self.normalized_trna_seqs:
            for trimmed_seq in normalized_seq.trimmed_seqs:
                normalized_count_dict[trimmed_seq.representative_name] += 1
        self.counts_of_normalized_seqs_containing_trimmed_seqs = [normalized_count for normalized_count in normalized_count_dict.values()]

        # Find the "multiplicity" of trimmed sequences among normalized sequences.
        multiplicity_dict = OrderedDict()
        for trimmed_seq, normalized_count_item in zip(self.trimmed_trna_seqs, normalized_count_dict.items()):
            trimmed_representative_name, normalized_count = normalized_count_item
            multiplicity_dict[trimmed_representative_name] = trimmed_seq.input_count * normalized_count
        self.multiplicities_of_trimmed_seqs_among_normalized_seqs = [multiplicity for multiplicity in multiplicity_dict.values()]

        # Find the "average multiplicity" of normalized and modified sequences.
        for normalized_seq in self.normalized_trna_seqs:
            multiplicity_sum = 0
            for trimmed_seq in normalized_seq.trimmed_seqs:
                multiplicity_sum += multiplicity_dict[trimmed_seq.representative_name]
            self.average_multiplicities_of_normalized_seqs.append(round(multiplicity_sum / normalized_seq.input_count, 1))

        for modified_seq in self.modified_trna_seqs:
            multiplicity_sum = 0
            for normalized_seq in modified_seq.all_normalized_seqs:
                for trimmed_seq in normalized_seq.trimmed_seqs:
                    multiplicity_sum += multiplicity_dict[trimmed_seq.representative_name]
            self.average_multiplicities_of_modified_seqs.append(round(multiplicity_sum / modified_seq.input_count, 1))

        self.progress.end()


    def write_trimmed_table(self):
        self.progress.new("Writing tRNA-seq database table of trimmed tRNA sequences")
        self.progress.update("...")

        trimmed_table_entries = []
        for trimmed_seq, normalized_seq_count in zip(self.trimmed_trna_seqs, self.counts_of_normalized_seqs_containing_trimmed_seqs):
            trimmed_table_entries.append(
                (trimmed_seq.representative_name,
                 len(trimmed_seq.unique_seqs),
                 trimmed_seq.input_count,
                 trimmed_seq.seq_string,
                 normalized_seq_count,
                 trimmed_seq.unique_with_extra_fiveprime_count,
                 trimmed_seq.input_with_extra_fiveprime_count)
                + tuple([v for v in trimmed_seq.input_acceptor_variant_count_dict.values()]))

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('trimmed')
            trnaseq_db.db.create_table('trimmed', tables.trnaseq_trimmed_table_structure, tables.trnaseq_trimmed_table_types)
        trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                 % ('trimmed', ','.join('?' * len(tables.trnaseq_trimmed_table_structure))),
                                 trimmed_table_entries)

        trimmed_seq_count = len(self.trimmed_trna_seqs)
        trnaseq_db.db.set_meta_value('num_trimmed_trna_seqs', trimmed_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()

        self.run.info("Trimmed tRNA, removing 5'/3' ends", trimmed_seq_count)


    def write_normalized_table(self):
        self.progress.new("Writing tRNA-seq database table of normalized tRNA sequences")
        self.progress.update("...")

        normalized_table_entries = []
        for normalized_seq, average_multiplicity in zip(self.normalized_trna_seqs,
                                                        self.average_multiplicities_of_normalized_seqs):
            normalized_table_entries.append(
                (normalized_seq.representative_name,
                 len(normalized_seq.trimmed_seqs),
                 normalized_seq.input_count,
                 average_multiplicity,
                 normalized_seq.trimmed_seqs_mapped_without_extra_fiveprime_count,
                 normalized_seq.input_seqs_mapped_without_extra_fiveprime_count,
                 normalized_seq.trimmed_seqs_mapped_with_extra_fiveprime_count,
                 normalized_seq.input_seqs_mapped_with_extra_fiveprime_count)
                + tuple(normalized_seq.input_acceptor_variant_count_dict.values()))

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('normalized')
            trnaseq_db.db.create_table('normalized', tables.trnaseq_normalized_table_structure, tables.trnaseq_normalized_table_types)
        trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                 % ('normalized', ','.join('?' * len(tables.trnaseq_normalized_table_structure))),
                                 normalized_table_entries)

        normalized_seq_count = len(self.normalized_trna_seqs)
        trnaseq_db.db.set_meta_value('num_normalized_trna_seqs', normalized_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()

        self.run.info("Normalized tRNA, consolidating tRNA fragments", normalized_seq_count)


    def write_modified_table(self):
        self.progress.new("Writing tRNA-seq database table of modified tRNA sequences")
        self.progress.update("...")

        modified_table_entries = []
        for modified_seq, average_multiplicity in zip(self.modified_trna_seqs,
                                                      self.average_multiplicities_of_modified_seqs):
            modified_table_entries.append(
                (modified_seq.representative_name,
                 ','.join([str(substitution_index) for substitution_index in modified_seq.substitution_indices]) + ',',
                 ','.join(sorted(str(deletion_index) for deletion_index in modified_seq.deletion_indices)) + ',',
                 modified_seq.consensus_seq_string,
                 ','.join([normalized_seq.representative_name for normalized_seq in modified_seq.all_normalized_seqs]),
                 len(modified_seq.all_normalized_seqs),
                 modified_seq.input_count,
                 average_multiplicity,
                 modified_seq.input_seqs_mapped_without_extra_fiveprime_count,
                 modified_seq.input_seqs_mapped_with_extra_fiveprime_count)
                + tuple(modified_seq.input_acceptor_variant_count_dict.values()))

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('modified')
            trnaseq_db.db.create_table('modified', tables.trnaseq_modified_table_structure, tables.trnaseq_modified_table_types)
        trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                 % ('modified', ','.join('?' * len(tables.trnaseq_modified_table_structure))),
                                 modified_table_entries)

        modified_seq_count = len(self.modified_trna_seqs)
        trnaseq_db.db.set_meta_value('num_modified_trna_seqs', modified_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()

        self.run.info("Modified tRNA", modified_seq_count)


    def write_uniqued_nontrna_supplement(self):
        self.progress.new("Writing a file of sequences not identified as tRNA.")
        self.run.info("Output non-tRNA file", self.uniqued_nontrna_path)

        with open(self.uniqued_nontrna_path, 'w') as nontrna_file:
            nontrna_file.write("\t".join(self.UNIQUED_NONTRNA_HEADER) + "\n")
            for unique_seq in self.unique_nontrna_seqs:
                nontrna_file.write(unique_seq.representative_name + "\t"
                                   + str(unique_seq.input_count) + "\t"
                                   + unique_seq.seq_string + "\n")

        self.progress.end()


    def write_trimmed_supplement(self):
        self.progress.new("Writing a file showing how trimmed tRNA sequences were formed from unique sequences")
        self.progress.update("...")

        with open(self.trimmed_ends_path, 'w') as trimmed_file:
            trimmed_file.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for trimmed_seq in sorted(self.trimmed_trna_seqs,
                                      key=lambda trimmed_seq: -trimmed_seq.input_count):
                representative_name = trimmed_seq.representative_name
                for unique_seq in sorted(trimmed_seq.unique_seqs,
                                         key=lambda unique_seq: (-unique_seq.extra_fiveprime_length,
                                                                 -unique_seq.acceptor_length)):
                    trimmed_file.write(representative_name + "\t"
                                       + unique_seq.representative_name + "\t"
                                       + unique_seq.seq_string[: unique_seq.extra_fiveprime_length] + "\t"
                                       + unique_seq.seq_string[len(unique_seq.seq_string) - unique_seq.acceptor_length: ] + "\t"
                                       + str(unique_seq.input_count) + "\n")

        self.progress.end()

        self.run.info("Output trimmed tRNA file", self.trimmed_ends_path)


    def write_normalized_supplement(self):
        self.progress.new("Writing a file showing how normalized tRNA sequences were formed from trimmed sequences")
        self.progress.update("...")

        with open(self.normalized_fragments_path, 'w') as normalized_file:
            normalized_file.write("\t".join(self.NORMALIZED_FRAGMENTS_HEADER) + "\n")
            for normalized_seq in sorted(self.normalized_trna_seqs,
                                         key=lambda normalized_seq: -normalized_seq.input_count):
                representative_name = normalized_seq.representative_name
                for trimmed_seq, start_position, end_position in sorted(
                    zip(normalized_seq.trimmed_seqs, normalized_seq.start_positions, normalized_seq.end_positions),
                    key=lambda t: (t[1], -t[2])):
                    normalized_file.write(representative_name + "\t"
                                          + trimmed_seq.representative_name + "\t"
                                          + str(start_position) + "\t"
                                          + str(end_position) + "\n")

        self.progress.end()

        self.run.info("Output normalized tRNA file", self.normalized_fragments_path)


    def process(self):
        """The main method of TRNASeqDataset, called from `anvi-trnaseq`

        Checkpoint loading and saving occurs in this method.
        """
        self.sanity_check()

        # The first checkpoints are after tRNA profiling.
        if not self.load_checkpoint:
            self.create_trnaseq_db()

            # Profile each input sequence for tRNA features.
            if self.feature_param_path:
                self.progress.new("Setting tRNA feature parameters from user file")
                self.progress.update("...")

                trnaidentifier.TRNAFeature.set_params_from_file(self.feature_param_path)

                self.progress.end()

            self.profile_trna(self.get_unique_input_seqs())

            if self.write_checkpoints:
                self.progress.new("Writing intermediate files for the \"profile\" checkpoint")
                self.progress.update("...")

                if os.path.exists(self.profile_unique_trna_seqs_path):
                    overwrote_profile_unique_trna_seqs_path = True
                else:
                    overwrote_profile_unique_trna_seqs_path = False
                with open(self.profile_unique_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.unique_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.profile_unique_nontrna_seqs_path):
                    overwrote_profile_unique_nontrna_seqs_path = True
                else:
                    overwrote_profile_unique_nontrna_seqs_path = False
                with open(self.profile_unique_nontrna_seqs_path, 'wb') as f:
                    pkl.dump(self.unique_nontrna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                self.progress.end()

                self.run.info("%s\"profile\" checkpoint intermediate file of unique tRNA"
                              % ("Overwritten " if overwrote_profile_unique_trna_seqs_path else ""),
                              self.profile_unique_trna_seqs_path)
                self.run.info("%s\"profile\" checkpoint intermediate file of unique non-tRNA"
                              % ("Overwritten " if overwrote_profile_unique_nontrna_seqs_path else ""),
                              self.profile_unique_nontrna_seqs_path)

        if self.load_checkpoint == 'profile':
            self.progress.new("Loading intermediate files at the checkpoint, \"profile\"")
            self.progress.update("...")

            with open(self.profile_unique_trna_seqs_path, 'rb') as f:
                self.unique_trna_seqs = pkl.load(f)
            with open(self.profile_unique_nontrna_seqs_path, 'rb') as f:
                self.unique_nontrna_seqs = pkl.load(f)

            self.progress.end()

        if self.load_checkpoint == 'threeprime_normalization':
            self.progress.new("Loading intermediate files at the checkpoint, \"threeprime_normalization\"")
            self.progress.update("...")

            with open(self.threeprime_normalization_unique_trna_seqs_path, 'rb') as f:
                self.unique_trna_seqs = pkl.load(f)
            with open(self.threeprime_normalization_unique_nontrna_seqs_path, 'rb') as f:
                self.unique_nontrna_seqs = pkl.load(f)
            with open(self.threeprime_normalization_trimmed_trna_seqs_path, 'rb') as f:
                self.trimmed_trna_seqs = pkl.load(f)
            with open(self.threeprime_normalization_normalized_trna_seqs_path, 'rb') as f:
                self.normalized_trna_seqs = pkl.load(f)

            self.progress.end()
        elif self.load_checkpoint == 'fragment_mapping':
            pass
        else:
            # Trim 5' and 3' ends of profiled tRNA.
            self.trim_ends(self.unique_trna_seqs)

            # Consolidate 3' fragments of longer tRNA sequences.
            self.dereplicate_threeprime()

            if self.write_checkpoints:
                self.progress.new("Writing intermediate files for the \"threeprime_normalization\" checkpoint")
                self.progress.update("...")

                if os.path.exists(self.threeprime_normalization_unique_trna_seqs_path):
                    overwrote_threeprime_normalization_unique_trna_seqs_path = True
                else:
                    overwrote_threeprime_normalization_unique_trna_seqs_path = False
                with open(self.threeprime_normalization_unique_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.unique_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.threeprime_normalization_unique_nontrna_seqs_path):
                    overwrote_threeprime_normalization_unique_nontrna_seqs_path = True
                else:
                    overwrote_threeprime_normalization_unique_nontrna_seqs_path = False
                with open(self.threeprime_normalization_unique_nontrna_seqs_path, 'wb') as f:
                    pkl.dump(self.unique_nontrna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.threeprime_normalization_trimmed_trna_seqs_path):
                    overwrote_threeprime_normalization_trimmed_trna_seqs_path = True
                else:
                    overwrote_threeprime_normalization_trimmed_trna_seqs_path = False
                with open(self.threeprime_normalization_trimmed_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.trimmed_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.threeprime_normalization_normalized_trna_seqs_path):
                    overwrote_threeprime_normalization_normalized_trna_seqs_path = True
                else:
                    overwrote_threeprime_normalization_normalized_trna_seqs_path = False
                with open(self.threeprime_normalization_normalized_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.normalized_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                self.progress.end()

                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of unique tRNA"
                              % ("Overwritten " if overwrote_threeprime_normalization_unique_trna_seqs_path else ""),
                              self.threeprime_normalization_unique_trna_seqs_path)
                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of unique non-tRNA"
                              % ("Overwritten " if overwrote_threeprime_normalization_unique_nontrna_seqs_path else ""),
                              self.threeprime_normalization_unique_nontrna_seqs_path)
                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of trimmed tRNA"
                              % ("Overwritten " if overwrote_threeprime_normalization_trimmed_trna_seqs_path else ""),
                              self.threeprime_normalization_trimmed_trna_seqs_path)
                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of normalized tRNA"
                              % ("Overwritten " if overwrote_threeprime_normalization_normalized_trna_seqs_path else ""),
                              self.threeprime_normalization_normalized_trna_seqs_path)

        if self.load_checkpoint == 'fragment_mapping':
            self.progress.new("Loading intermediate files at the checkpoint, \"fragment_mapping\"")
            self.progress.update("...")

            with open(self.fragment_mapping_unique_trna_seqs_path, 'rb') as f:
                self.unique_trna_seqs = pkl.load(f)
            with open(self.fragment_mapping_unique_nontrna_seqs_path, 'rb') as f:
                self.unique_nontrna_seqs = pkl.load(f)
            with open(self.fragment_mapping_trimmed_trna_seqs_path, 'rb') as f:
                self.trimmed_trna_seqs = pkl.load(f)
            with open(self.fragment_mapping_normalized_trna_seqs_path, 'rb') as f:
                self.normalized_trna_seqs = pkl.load(f)

            self.progress.end()
        else:
            # Map fragments derived from the interior and 5' end of tRNA.
            self.map_fragments()

            if self.write_checkpoints:
                self.progress.new("Writing intermediate files for the \"fragment_mapping\" checkpoint")
                self.progress.update("...")

                if os.path.exists(self.fragment_mapping_unique_trna_seqs_path):
                    overwrote_fragment_mapping_unique_trna_seqs_path = True
                else:
                    overwrote_fragment_mapping_unique_trna_seqs_path = False
                with open(self.fragment_mapping_unique_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.unique_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.fragment_mapping_unique_nontrna_seqs_path):
                    overwrote_fragment_mapping_unique_nontrna_seqs_path = True
                else:
                    overwrote_fragment_mapping_unique_nontrna_seqs_path = False
                with open(self.fragment_mapping_unique_nontrna_seqs_path, 'wb') as f:
                    pkl.dump(self.unique_nontrna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.fragment_mapping_trimmed_trna_seqs_path):
                    overwrote_fragment_mapping_trimmed_trna_seqs_path = True
                else:
                    overwrote_fragment_mapping_trimmed_trna_seqs_path = False
                with open(self.fragment_mapping_trimmed_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.trimmed_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.fragment_mapping_normalized_trna_seqs_path):
                    overwrote_fragment_mapping_normalized_trna_seqs_path = True
                else:
                    overwrote_fragment_mapping_normalized_trna_seqs_path = False
                with open(self.fragment_mapping_normalized_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.normalized_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                self.progress.end()

                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of unique tRNA"
                              % ("Overwritten " if overwrote_fragment_mapping_unique_trna_seqs_path else ""),
                              self.threeprime_normalization_unique_trna_seqs_path)
                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of unique non-tRNA"
                              % ("Overwritten " if overwrote_fragment_mapping_unique_nontrna_seqs_path else ""),
                              self.fragment_mapping_unique_nontrna_seqs_path)
                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of trimmed tRNA"
                              % ("Overwritten " if overwrote_fragment_mapping_trimmed_trna_seqs_path else ""),
                              self.fragment_mapping_trimmed_trna_seqs_path)
                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of normalized tRNA"
                              % ("Overwritten " if overwrote_fragment_mapping_normalized_trna_seqs_path else ""),
                              self.fragment_mapping_normalized_trna_seqs_path)

        # Find modified nucleotides, grouping sequences into modified sequences.
        self.find_modifications()

        # Calculate some statistics.
        self.calc_normalization_stats()

        # Write more tables to the database.
        self.write_trimmed_table()
        self.write_normalized_table()
        self.write_modified_table()

        # Write supplementary text files.
        self.write_uniqued_nontrna_supplement()
        self.write_trimmed_supplement()
        self.write_normalized_supplement()


def profile_worker(input_queue, output_queue):
    while True:
        seq_string, seq_name = input_queue.get()
        output_queue.put(trnaidentifier.Profile(seq_string, name=seq_name))
