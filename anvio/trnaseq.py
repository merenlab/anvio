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

from itertools import combinations, product
from collections import OrderedDict, deque, defaultdict

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


unambiguous_nts = ('A', 'C', 'G', 'T')
# The next value is used for counting nucleotides in alignments: there is one bin (value 0) for end gaps in the alignment.
num_nt_bins = len(unambiguous_nts) + 1
nt_int_dict = {nt: i for i, nt in enumerate(unambiguous_nts, start=1)}
int_nt_dict = {i: nt for i, nt in enumerate(unambiguous_nts, start=1)}


class UniqueSeq:
    __slots__ = ('seq_string',
                 'represent_name',
                 'has_complete_feature_set',
                 'read_count',
                 'id_method',
                 'acceptor_length',
                 'extra_fiveprime_length')

    def __init__(self,
                 seq_string,
                 represent_name,
                 read_count,
                 id_method=None,
                 acceptor_length=None,
                 has_complete_feature_set=None,
                 extra_fiveprime_length=None):
        """A dereplicated tRNA-seq read, with information from tRNA feature profiling"""

        self.seq_string = seq_string
        self.represent_name = represent_name
        self.read_count = read_count
        # If dealing with tRNA, identification method 0 = profiled, 1 = mapped
        self.id_method = id_method
        self.acceptor_length = acceptor_length
        self.has_complete_feature_set = has_complete_feature_set
        self.extra_fiveprime_length = extra_fiveprime_length


class TrimmedSeq:
    __slots__ = ('seq_string',
                 'uniq_seqs',
                 'has_complete_feature_set',
                 'read_count',
                 'uniq_with_extra_fiveprime_count',
                 'read_with_extra_fiveprime_count',
                 'represent_name',
                 'read_acceptor_variant_count_dict',
                 'id_method',
                 'norm_seq_count')

    def __init__(self, seq_string, uniq_seqs, skip_init=False):
        """A tRNA sequence with bases trimmed 5' of the acceptor stem and 3' of the discriminator"""

        self.seq_string = seq_string
        self.uniq_seqs = uniq_seqs # list of UniqueSeq objects
        self.has_complete_feature_set = uniq_seqs[0].has_complete_feature_set
        self.norm_seq_count = 0

        if skip_init:
            self.read_count = None
            self.uniq_with_extra_fiveprime_count = None
            self.read_with_extra_fiveprime_count = None
            self.represent_name = None
            self.read_acceptor_variant_count_dict = None
            self.id_method = None
        else:
            self.init()


    def init(self):
        """Set attributes representative of a final set of input `UniqueSeq` objects"""

        self.read_count = sum([uniq_seq.read_count for uniq_seq in self.uniq_seqs])

        self.uniq_with_extra_fiveprime_count = sum([1 if uniq_seq.extra_fiveprime_length else 0
                                                    for uniq_seq in self.uniq_seqs])
        self.read_with_extra_fiveprime_count = sum([uniq_seq.read_count if uniq_seq.extra_fiveprime_length else 0
                                                    for uniq_seq in self.uniq_seqs])

        # The representative name is chosen as follows:
        # 1. Most abundant full-length tRNA (no extra 5' bases), ignoring acceptor sequence
        # 2. Most abundant longer-than-full-length tRNA
        # 3. Most abundant fragmentary tRNA
        # Sort such that the first sequence is the most abundant longest and the last is the least abundant shortest.
        uniq_seqs = sorted(self.uniq_seqs,
                           key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, -uniq_seq.read_count))

        if uniq_seqs[0].extra_fiveprime_length > 0:
            # If there is also a unique sequence that was ultimately trimmed down
            # to the same sequence as the sequence with extra 5' bases, it must be a full-length sequence.
            if uniq_seqs[-1].extra_fiveprime_length == 0:
                # Sort such that the last sequence is the most abundant shortest.
                represent_name = sorted(uniq_seqs,
                                        key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length,
                                                              uniq_seq.read_count))[-1].represent_name
            else:
                represent_name = uniq_seqs[0].represent_name
        else:
            # ALL unique sequences are EITHER full-length OR a fragment.
            represent_name = uniq_seqs[0].represent_name

        self.represent_name = represent_name

        read_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0)
                                                        for threeprime_variant in THREEPRIME_VARIANTS])
        for uniq_seq in self.uniq_seqs:
            if uniq_seq.acceptor_length: # unique_seq need not have an acceptor
                acceptor_seq_string = uniq_seq.seq_string[-uniq_seq.acceptor_length: ]
                read_acceptor_variant_count_dict[acceptor_seq_string] += uniq_seq.read_count
        self.read_acceptor_variant_count_dict = read_acceptor_variant_count_dict

        id_methods = set(uniq_seq.id_method for uniq_seq in self.uniq_seqs)
        if len(id_methods) == 1:
            self.id_method = id_methods.pop()
        else:
            raise ConfigError("A TrimmedSeq should not be made from UniqueSeq objects "
                              "with different identification methods. "
                              "Trimmed tRNA sequences will EITHER be formed from "
                              "\"profiled\" tRNA sequences or \"mapped\" tRNA sequences, "
                              "because they are of different lengths and are fragments from different parts of the tRNA.")


class NormalizedSeq:
    __slots__ = ('trimmed_seqs',
                 'represent_name',
                 'seq_string',
                 'has_complete_feature_set',
                 'start_positions',
                 'end_positions',
                 'read_count',
                 'read_with_extra_fiveprime_count',
                 'read_acceptor_variant_count_dict',
                 'trimmed_seqs_mapped_without_extra_fiveprime_count',
                 'reads_mapped_without_extra_fiveprime_count',
                 'trimmed_seqs_mapped_with_extra_fiveprime_count',
                 'reads_mapped_with_extra_fiveprime_count',
                 'specific_covs',
                 'nonspecific_covs',
                 'mod_seqs')

    def __init__(self, trimmed_seqs, start_positions=None, end_positions=None, skip_init=False):
        """A longer tRNA sequence consolidated from shorter tRNA fragments"""

        self.trimmed_seqs = trimmed_seqs # list of TrimmedSeq objects
        for trimmed_seq in trimmed_seqs:
            trimmed_seq.norm_seq_count += 1
        self.represent_name = trimmed_seqs[0].represent_name
        self.seq_string = trimmed_seqs[0].seq_string
        self.has_complete_feature_set = trimmed_seqs[0].has_complete_feature_set
        if start_positions and end_positions:
            self.start_positions = start_positions
            self.end_positions = end_positions
        elif (not start_positions) and (not end_positions):
            # Trimmed seqs were dereplicated from the 3' end of the normalized sequence.
            norm_seq_length = len(self.seq_string)
            self.start_positions = [norm_seq_length - len(trimmed_seq.seq_string)
                                    for trimmed_seq in self.trimmed_seqs]
            self.end_positions = [norm_seq_length] * len(trimmed_seqs)
        else:
            self.start_positions = None
            self.end_positions = None

        # It is useful to know which modified sequences, if any, encompass this normalized sequence.
        # A normalized sequence without modification-induced deletions can only be assigned to one modified sequence,
        # but a normalized sequence with deletions can be assigned to more than one modified sequence.
        self.mod_seqs = []

        if skip_init:
            self.read_count = None
            self.read_with_extra_fiveprime_count = None
            self.read_acceptor_variant_count_dict = None
            self.trimmed_seqs_mapped_without_extra_fiveprime_count = None
            self.reads_mapped_without_extra_fiveprime_count = None
            self.trimmed_seqs_mapped_with_extra_fiveprime_count = None
            self.reads_mapped_with_extra_fiveprime_count = None
            self.specific_covs = None
            self.nonspecific_covs = None
        else:
            self.init()


    def init(self):
        """Set the attributes representative of a finalized list of `TrimmedSeq` objects"""

        self.read_count = sum([trimmed_seq.read_count for trimmed_seq in self.trimmed_seqs])

        read_acceptor_variant_count_dict = OrderedDict([(threeprime_variant, 0)
                                                        for threeprime_variant in THREEPRIME_VARIANTS])
        for trimmed_seq in self.trimmed_seqs:
            for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                if read_count > 0:
                    read_acceptor_variant_count_dict[acceptor_seq_string] += read_count
        self.read_acceptor_variant_count_dict = read_acceptor_variant_count_dict

        read_with_extra_fiveprime_count = 0
        trimmed_seqs_mapped_without_extra_fiveprime_count = 0
        reads_mapped_without_extra_fiveprime_count = 0
        trimmed_seqs_mapped_with_extra_fiveprime_count = 0
        reads_mapped_with_extra_fiveprime_count = 0
        specific_covs = np.zeros(len(self.seq_string), dtype=int)
        nonspecific_covs = np.zeros(len(self.seq_string), dtype=int)
        for trimmed_seq, start_pos, end_pos in zip(self.trimmed_seqs,
                                                   self.start_positions,
                                                   self.end_positions):
            if trimmed_seq.id_method == 1: # 1 => mapped
                # TrimmedSeqs are comprised of EITHER profiled (0) OR mapped (1) UniqueSeqs.
                if trimmed_seq.uniq_with_extra_fiveprime_count == 0:
                    trimmed_seqs_mapped_without_extra_fiveprime_count += 1
                    reads_mapped_without_extra_fiveprime_count += trimmed_seq.read_count
                else:
                    read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count
                    trimmed_seqs_mapped_with_extra_fiveprime_count += 1
                    reads_mapped_with_extra_fiveprime_count += trimmed_seq.read_count
            else:
                read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count

            if trimmed_seq.norm_seq_count == 1:
                specific_covs[start_pos: end_pos] += trimmed_seq.read_count
            else:
                nonspecific_covs[start_pos: end_pos] += trimmed_seq.read_count
        self.trimmed_seqs_mapped_without_extra_fiveprime_count = trimmed_seqs_mapped_without_extra_fiveprime_count
        self.reads_mapped_without_extra_fiveprime_count = reads_mapped_without_extra_fiveprime_count
        self.trimmed_seqs_mapped_with_extra_fiveprime_count = trimmed_seqs_mapped_with_extra_fiveprime_count
        self.reads_mapped_with_extra_fiveprime_count = reads_mapped_with_extra_fiveprime_count
        self.specific_covs = specific_covs
        self.nonspecific_covs = nonspecific_covs


class ModifiedSeq:
    __slots__ = ('norm_seqs_without_dels',
                 'sub_positions',
                 'represent_name',
                 'norm_seqs_with_dels',
                 'del_configs',
                 'specific_covs',
                 'nonspecific_covs',
                 'specific_sub_covs',
                 'nonspecific_sub_covs',
                 'specific_del_covs',
                 'nonspecific_del_covs',
                 'specific_read_count',
                 'nonspecific_read_count',
                 'count_of_specific_reads_with_extra_fiveprime',
                 'count_of_nonspecific_reads_with_extra_fiveprime',
                 'specific_mapped_read_count',
                 'nonspecific_mapped_read_count',
                 'consensus_seq_string')

    def __init__(self, norm_seqs_without_dels, sub_positions, init=False):
        """A tRNA sequence with sites of predicted modification-induced substitutions and deletions

        Parameters
        ==========
        norm_seqs_without_dels : list
            NormalizedSeq objects representing sequences with distinct modification-induced substitutions
            The first sequence in the list should be longest or tied for longest.
            Sequences with modification-induced deletions must be added later.

        sub_positions : list
            Positions of modification-induced substitutions in the modified sequence
            Positions in the modified sequence are equivalent to positions in the longest normalized sequence with substitutions.

        init : bool, False
            Triggers the analysis of added normalized sequences, finding nucleotide coverage and other information
            This should be set to True when normalized sequences with deletions are not going to be added.
        """

        self.norm_seqs_without_dels = norm_seqs_without_dels
        self.sub_positions = sub_positions
        # A normalized sequence without modification-induced deletions can only be assigned to one modified sequence,
        # but a normalized sequence with deletions can be assigned to more than one modified sequence.
        for norm_seq in norm_seqs_without_dels:
            norm_seq.mod_seqs.append(self)
        self.represent_name = norm_seqs_without_dels[0].represent_name

        if init:
            self.init()
        else:
            self.norm_seqs_with_dels = []
            self.del_configs = []
            self.specific_covs = None
            self.nonspecific_covs = None
            self.specific_sub_covs = None
            self.nonspecific_sub_covs = None
            self.specific_del_covs = None
            self.nonspecific_del_covs = None
            self.specific_read_count = None
            self.nonspecific_read_count = None
            self.count_of_specific_reads_with_extra_fiveprime = None
            self.count_of_nonspecific_reads_with_extra_fiveprime = None
            self.specific_mapped_read_count = None
            self.nonspecific_mapped_read_count = None
            self.consensus_seq_string = None


    def get_seqs_with_dels(self, possible_del_starts=(-2, -1, 0), possible_del_stops=(0, 1), max_del_sites=2):
        """Generate in silico modified sequences with deletions at and/or around substitution sites.

        Parameters
        ==========
        possible_del_starts : tuple, (-2, -1, 0)
            Where deletions can start relative to substitution sites.
            By default, deletions can start at the substitution site or at the two preceding 5' nucleotides.

        possible_del_stops : tuple, (0, 1)
            Where deletions can stop relative to substitution sites.
            By default, deletions can stop at the substitution site or at the preceding 5' nucleotide.
            In conjunction with the default deletion starts, deletions can be 1-3 nucleotides long.

        max_del_sites : int, 2
            The maximum number of substitution sites at which deletions can be introduced.
            For example, if this parameter is set to 2, and there are 3 substitution positions in the sequence,
            then sequences will be produced containing deletions only at the first position;
            other sequences will be produced containing deletions at the first and second positions;
            other sequences will be produced containing deletions only at the second position, etc.

        Returns
        =======
        del_set : set
            A set of tuples with 2 elements.
            The first element is a sequence string containing deletions.
            The second element is a tuple of the indices of these deletions in the input sequence.
        """

        # Make template sequences with different nucleotides at substitution sites.
        # Only consider the observed substitution configurations.
        seq_strings_without_dels = set()
        longest_norm_seq_string = self.norm_seqs_without_dels[0].seq_string
        mod_seq_length = len(longest_norm_seq_string)
        sub_positions = self.sub_positions
        for norm_seq in self.norm_seqs_without_dels:
            norm_seq_string = norm_seq.seq_string
            norm_seq_start_in_mod_seq = mod_seq_length - len(norm_seq_string)
            altered_seq_string = longest_norm_seq_string
            for sub_pos in sub_positions:
                if sub_pos < norm_seq_start_in_mod_seq:
                    # This normalized sequence is shorter, lacking the substitution position.
                    continue
                nt = norm_seq.seq_string[sub_pos - norm_seq_start_in_mod_seq]
                altered_seq_string = (altered_seq_string[: sub_pos]
                                      + nt
                                      + altered_seq_string[sub_pos + 1: ])
            seq_strings_without_dels.add(altered_seq_string)

        # Determine the different deletion sizes/positions relative to a substitution.
        del_ranges = []
        for del_start in possible_del_starts:
            for del_stop in possible_del_stops:
                if del_start < del_stop:
                    del_ranges.append(range(del_start, del_stop))

        # Introduce deletions into each template sequence, potentially producing a number of new sequences.
        del_dict = {}
        for seq_string in seq_strings_without_dels:
            del_dict_for_seq = self.introduce_dels(seq_string, sub_positions, del_ranges, max_del_sites)
            for seq_string_with_del, del_positions in del_dict_for_seq.items():
                try:
                    prior_del_pos_sum = sum(del_dict[seq_string_with_del])
                except KeyError:
                    del_dict[seq_string_with_del] = del_positions
                    continue
                if sum(del_positions) < prior_del_pos_sum:
                    del_dict[seq_string_with_del] = del_positions
        del_set = set([(seq_string_with_del, del_positions)
                       for seq_string_with_del, del_positions in del_dict.items()])

        return del_set


    @staticmethod
    def introduce_dels(seq_string, sub_positions, del_ranges, max_del_sites):
        """Generate in silico sequences with deletions at and/or around substitution sites in the input sequence.

        Parameters
        ==========
        seq_string : str
            The sequence in which deletions will be introduced

        sub_positions : list-like
            Where substitutions are located in the input sequence

        del_ranges : list-like
            A list of ranges representing the locations of deletions in relation to substitution sites
            For example, [range(-1, 0), range(-1, 1), range(0, 1)]
            allows three types of deletions to be introduced at a substitution position:
            a one-nucleotide deletion of the adjacent 5' nucleotide
            a two-nucleotide deletion of the adjacent 5' nucleotide and the nucleotide at the substitution position,
            and a one-nucleotide deletion of the nucleotide at the substitution position.

        max_del_sites : int
            The maximum number of substitution sites at which deletions can be introduced.
            For example, if this parameter is set to 2, and there are 3 substitution positions in the sequence,
            then sequences with deletions will be produced containing deletions only at the first position,
            other sequences with deletions will be produced containing deletions at the first and second positions,
            other sequences with deletions will be produced containing deletions only at the second position, etc.

        Returns
        =======
        del_dict : dict
            Each dict key is a sequence string containing deletions.
            Each dict value is a tuple of the indices of these deletions in the input sequence.
        """

        # Find all the ways deletions can be introduced into the sequence given the parameterization.
        del_pos_configs = set()
        for num_del_sites in range(1, max_del_sites + 1):
            for del_locus_config in combinations(sub_positions, num_del_sites):
                for del_range_config in product(*[del_ranges for _ in range(num_del_sites)]):
                    del_positions = set()
                    for i, del_range in enumerate(del_range_config):
                        sub_pos = del_locus_config[i]
                        for del_pos_relative_to_sub in del_range:
                            del_pos = sub_pos + del_pos_relative_to_sub
                            if del_pos >= 0:
                                del_positions.add(del_pos)
                    if del_positions:
                        del_positions = sorted(del_positions)

                        # Remove any nominal deletions at the 5' end,
                        # as these could rightly be interpreted as unseen nucleotides preceding a fragment.
                        fiveprime_del_pos = -1
                        unsupported_del_positions = []
                        for d in del_positions:
                            if d == fiveprime_del_pos + 1:
                                unsupported_del_positions.append(d)
                                fiveprime_del_pos += 1
                            else:
                                break
                        for d in unsupported_del_positions[::-1]:
                            del_positions.pop(d)

                        del_pos_configs.add(tuple(del_positions))

        # It is possible to generate the same sequence with deletions given different deletion sites.
        # For example, ACCG can become ACG by deleting either C.
        # We resolve this complication by choosing the most 5' deletion indices.
        del_dict = {}
        for del_positions in del_pos_configs:
            seq_string_with_dels = seq_string
            for del_pos in del_positions[::-1]:
                seq_string_with_dels = seq_string_with_dels[: del_pos] + seq_string_with_dels[del_pos + 1: ]
            try:
                prior_del_pos_sum = sum(del_dict[seq_string_with_dels])
            except KeyError:
                del_dict[seq_string_with_dels] = del_positions
                continue
            if sum(del_positions) < prior_del_pos_sum:
                del_dict[seq_string_with_dels] = del_positions

        return del_dict


    def init(self):
        """Analyze sequences comprising the modified sequence to find coverages and other attributes."""

        norm_seqs_without_dels = self.norm_seqs_without_dels
        norm_seqs_with_dels = self.norm_seqs_with_dels
        all_norm_seqs = norm_seqs_without_dels + norm_seqs_with_dels
        del_configs = self.del_configs
        mod_seq_len = len(norm_seqs_without_dels[0].seq_string)
        norm_seq_specific_covs = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        norm_seq_nonspecific_covs = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        num_subs = len(self.sub_positions)
        self.specific_sub_covs = specific_sub_covs = np.zeros((num_subs, len(unambiguous_nts)), dtype=int)
        self.nonspecific_sub_covs = nonspecific_sub_covs = np.zeros((num_subs, len(unambiguous_nts)), dtype=int)
        del_positions = sorted(set([i for del_config in del_configs for i in del_config]))
        self.specific_del_covs = specific_del_covs = np.zeros(len(del_positions), dtype=int)
        self.nonspecific_del_covs = nonspecific_del_covs = np.zeros(len(del_positions), dtype=int)
        specific_read_count = 0
        nonspecific_read_count = 0
        count_of_specific_reads_with_extra_fiveprime = 0
        count_of_nonspecific_reads_with_extra_fiveprime = 0
        specific_mapped_read_count = 0
        nonspecific_mapped_read_count = 0

        # Make an array of aligned nucleotide positions in all normalized sequences.
        norm_seq_array = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        for n, norm_seq in enumerate(norm_seqs_without_dels):
            norm_seq_array[n, mod_seq_len - len(norm_seq.seq_string): ] += [nt_int_dict[nt]
                                                                            for nt in norm_seq.seq_string]

        nt_positions_covered_by_norm_seqs_with_dels = []
        n = len(norm_seqs_without_dels)
        for norm_seq, del_config in zip(norm_seqs_with_dels, del_configs):
            aligned_seq = [nt_int_dict[nt] for nt in norm_seq.seq_string]
            # Insert a 0 (no nucleotide) at each deletion position in the alignment.
            for del_pos in del_config:
                aligned_seq.insert(del_pos, 0)
            norm_seq_start_in_mod_seq = mod_seq_len - len(aligned_seq)
            norm_seq_array[n, norm_seq_start_in_mod_seq: ] += aligned_seq

            covered_nt_positions = []
            for i, nt_int in enumerate(aligned_seq):
                if nt_int != 0:
                    covered_nt_positions.append(norm_seq_start_in_mod_seq + i)
            nt_positions_covered_by_norm_seqs_with_dels.append(covered_nt_positions)
            n += 1

        processed_trimmed_seq_names = []

        # First handle sequences without deletions.
        for n, norm_seq in enumerate(norm_seqs_without_dels):
            norm_seq_start_in_mod_seq = mod_seq_len - len(norm_seq.seq_string)

            for trimmed_seq, trimmed_seq_start_in_norm_seq, trimmed_seq_stop_in_norm_seq in zip(norm_seq.trimmed_seqs,
                                                                                                norm_seq.start_positions,
                                                                                                norm_seq.end_positions):
                if trimmed_seq.represent_name in processed_trimmed_seq_names:
                    continue

                # Determine whether the reads constituting the trimmed sequence are specific to the modified sequence.
                if trimmed_seq.norm_seq_count == 1:
                    is_trimmed_seq_specific_to_mod_seq = True
                else:
                    # The trimmed sequence is specific to the modified sequence
                    # if it is unique to a set of normalized sequences specific to the modified sequence.
                    # Coverage information for such trimmed sequences will be recorded
                    # in the row of the array for the first normalized sequence in which it was found.
                    num_specific_norm_seqs_containing_trimmed_seq = 1
                    trimmed_seq_name = trimmed_seq.represent_name
                    for p, other_norm_seq in enumerate(all_norm_seqs):
                        if n == p:
                            continue
                        for other_trimmed_seq in other_norm_seq.trimmed_seqs:
                            if trimmed_seq_name == other_trimmed_seq.represent_name:
                                num_specific_norm_seqs_containing_trimmed_seq += 1
                                break

                    if num_specific_norm_seqs_containing_trimmed_seq == trimmed_seq.norm_seq_count:
                        is_trimmed_seq_specific_to_mod_seq = True
                    elif num_specific_norm_seqs_containing_trimmed_seq < trimmed_seq.norm_seq_count:
                        is_trimmed_seq_specific_to_mod_seq = False
                    else:
                        raise ConfigError("The number of normalized sequences containing the trimmed sequence was somehow miscalculated.")

                trimmed_seq_start_in_mod_seq = norm_seq_start_in_mod_seq + trimmed_seq_start_in_norm_seq
                trimmed_seq_stop_in_mod_seq = trimmed_seq_start_in_mod_seq + len(trimmed_seq.seq_string)
                if is_trimmed_seq_specific_to_mod_seq:
                    specific_read_count += trimmed_seq.read_count
                    count_of_specific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    # Trimmed sequences are comprised of either profiled or mapped unique sequences.
                    if trimmed_seq.id_method == 1: # 1 => mapped
                        specific_mapped_read_count += trimmed_seq.read_count
                    norm_seq_specific_covs[n, trimmed_seq_start_in_mod_seq: trimmed_seq_stop_in_mod_seq] += trimmed_seq.read_count
                else:
                    nonspecific_read_count += trimmed_seq.read_count
                    count_of_nonspecific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    if trimmed_seq.id_method == 1:
                        nonspecific_mapped_read_count += trimmed_seq.read_count
                    norm_seq_nonspecific_covs[n, trimmed_seq_start_in_mod_seq: trimmed_seq_stop_in_mod_seq] += trimmed_seq.read_count

                processed_trimmed_seq_names.append(trimmed_seq.represent_name)

        # Handle normalized sequences with deletions.
        n = len(norm_seqs_without_dels)
        for norm_seq, del_config, nt_positions_covered_by_norm_seq in zip(norm_seqs_with_dels,
                                                                          del_configs,
                                                                          nt_positions_covered_by_norm_seqs_with_dels):
            norm_seq_start_in_mod_seq = mod_seq_len - len(norm_seq.seq_string) - len(del_config)
            # Normalized sequences with deletions can be found in multiple modified sequences,
            # unlike normalized sequences without deletions.
            num_mod_seqs_containing_norm_seq = len(norm_seq.mod_seqs)

            for trimmed_seq, trimmed_seq_start_in_norm_seq, trimmed_seq_stop_in_norm_seq in zip(norm_seq.trimmed_seqs,
                                                                                                norm_seq.start_positions,
                                                                                                norm_seq.end_positions):
                if trimmed_seq.represent_name in processed_trimmed_seq_names:
                    continue

                # Determine whether the reads constituting the trimmed sequence are specific to the modified sequence.
                if num_mod_seqs_containing_norm_seq > 1:
                    is_trimmed_seq_specific_to_mod_seq = False
                else:
                    # The trimmed sequence is specific to the modified sequence
                    # if it is unique to a set of normalized sequences specific to the modified sequence.
                    num_specific_norm_seqs_containing_trimmed_seq = 1
                    trimmed_seq_name = trimmed_seq.represent_name
                    is_trimmed_seq_specific_to_mod_seq = True # initial value
                    for p, other_norm_seq in enumerate(norm_seqs_with_dels, start=len(norm_seqs_without_dels)):
                        if n == p:
                            continue
                        for other_trimmed_seq in other_norm_seq.trimmed_seqs:
                            if trimmed_seq_name == other_trimmed_seq.represent_name:
                                if len(other_norm_seq.mod_seqs) > 1:
                                    is_trimmed_seq_specific_to_mod_seq = False
                                else:
                                    num_specific_norm_seqs_containing_trimmed_seq += 1
                                break
                        if not is_trimmed_seq_specific_to_mod_seq:
                            # The trimmed sequence was found in a normalized sequence that is in multiple modified sequences.
                            break
                    else:
                        if num_specific_norm_seqs_containing_trimmed_seq == trimmed_seq.norm_seq_count:
                            is_trimmed_seq_specific_to_mod_seq = True
                        elif num_specific_norm_seqs_containing_trimmed_seq < trimmed_seq.norm_seq_count:
                            # The trimmed sequence was found in normalized sequences that are not in the modified sequence.
                            is_trimmed_seq_specific_to_mod_seq = False
                        else:
                            raise ConfigError("The number of normalized sequences containing the trimmed sequence was somehow miscalculated.")

                nt_positions_covered_by_trimmed_seq = nt_positions_covered_by_norm_seq[
                    trimmed_seq_start_in_norm_seq: trimmed_seq_start_in_norm_seq + len(trimmed_seq.seq_string)
                ]

                trimmed_seq_read_count = trimmed_seq.read_count
                if is_trimmed_seq_specific_to_mod_seq:
                    specific_read_count += trimmed_seq_read_count
                    count_of_specific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    # Trimmed sequences are comprised of either profiled or mapped unique sequences.
                    if trimmed_seq.id_method == 1: # 1 => mapped
                        specific_mapped_read_count += trimmed_seq_read_count
                    norm_seq_specific_covs[n, nt_positions_covered_by_trimmed_seq] += trimmed_seq_read_count
                    for del_pos in del_config:
                        specific_del_covs[del_positions.index(del_pos)] += trimmed_seq_read_count
                else:
                    nonspecific_read_count += trimmed_seq_read_count
                    count_of_nonspecific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    if trimmed_seq.id_method == 1:
                        nonspecific_mapped_read_count += trimmed_seq_read_count
                    norm_seq_nonspecific_covs[n, nt_positions_covered_by_trimmed_seq] += trimmed_seq_read_count
                    for del_pos in del_config:
                        nonspecific_del_covs[del_positions.index(del_pos)] += trimmed_seq_read_count

            processed_trimmed_seq_names.append(trimmed_seq.represent_name)
            n += 1

        self.specific_covs = norm_seq_specific_covs.sum(0)
        self.nonspecific_covs = norm_seq_nonspecific_covs.sum(0)
        self.specific_read_count = specific_read_count
        self.nonspecific_read_count = nonspecific_read_count
        self.count_of_specific_reads_with_extra_fiveprime = count_of_specific_reads_with_extra_fiveprime
        self.count_of_nonspecific_reads_with_extra_fiveprime = count_of_nonspecific_reads_with_extra_fiveprime
        self.specific_mapped_read_count = specific_mapped_read_count
        self.nonspecific_mapped_read_count = nonspecific_mapped_read_count

        # For each substitution position, record the coverage of A, C, G, and T.
        for s, sub_pos in enumerate(self.sub_positions):
            aligned_nts = norm_seq_array[:, sub_pos]
            nt_counts = np.bincount(aligned_nts, minlength=num_nt_bins)[1: ]
            for nt_int, nt_count in enumerate(nt_counts, start=1):
                if nt_count > 0:
                    norm_seq_rows_with_nt = (aligned_nts == nt_int).nonzero()[0]
                    specific_sub_covs[s, nt_int - 1] = norm_seq_specific_covs[norm_seq_rows_with_nt, sub_pos].sum()
                    nonspecific_sub_covs[s, nt_int - 1] = norm_seq_nonspecific_covs[norm_seq_rows_with_nt, sub_pos].sum()

        # Set a consensus sequence from the nucleotides with the highest specific coverage at each position.
        consensus_seq_string = norm_seqs_without_dels[0].seq_string
        for sub_pos, covs in zip(self.sub_positions, specific_sub_covs):
            max_pos = covs.argmax()
            nt_int = max_pos + 1
            consensus_seq_string = (consensus_seq_string[: sub_pos]
                                    + int_nt_dict[nt_int]
                                    + consensus_seq_string[sub_pos + 1: ])
        self.consensus_seq_string = consensus_seq_string


class TRNASeqDataset:
    # Column headers for supplementary tables written to text files
    UNIQ_NONTRNA_HEADER = ["represent_name", "read_count", "sequence"]
    TRIMMED_ENDS_HEADER = ["represent_name", "unique_name", "fiveprime_sequence", "threeprime_sequence", "read_count"]
    NORM_FRAG_HEADER = ["represent_name", "trimmed_name", "start", "end"]


    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        """Class for processing a tRNA-seq dataset"""

        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # Argument group 1A: MANDATORY
        self.input_fasta_path = A('fasta_file')
        self.project_name = A('project_name')
        self.out_dir = os.path.abspath(A('output_dir')) if A('output_dir') else None

        # Argument group 1B: MUNDANE
        self.overwrite_out_dest = A('overwrite_output_destinations')
        self.descrip_path = os.path.abspath(A('description')) if A('description') else None

        # Argument group 1C: PERFORMANCE
        self.num_threads = A('num_threads')
        self.skip_fasta_check = A('skip_fasta_check')
        self.write_checkpoints = A('write_checkpoints')
        self.load_checkpoint = A('load_checkpoint')
        self.write_buffer_size = A('write_buffer_size')
        self.alignment_target_chunk_size = A('alignment_target_chunk_size')
        self.frag_mapping_query_chunk_length = A('fragment_mapping_query_chunk_length')

        # Argument group 1D: ADVANCED
        self.feature_param_path = os.path.abspath(A('feature_param_file')) if A('feature_param_file') else None
        self.min_trna_frag_size = A('min_trna_fragment_size')
        self.agglom_max_mismatch_freq = A('agglomeration_max_mismatch_freq')
        self.min_modification_count = A('min_modification_count')
        self.min_modification_fraction = A('min_modification_fraction')
        self.max_del_size = A('max_deletion_size')

        # Argument group 1E: MINUTIAE
        self.alignment_progress_interval = A('alignment_progress_interval')
        self.agglom_progress_interval = A('agglomeration_progress_interval')

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using --fasta-file or -f.")
        if not self.project_name:
            raise ConfigError("Please set a project name using --project-name or -n.")
        if not self.out_dir:
            raise ConfigError("Please provide an output directory using --output-dir or -o.")

        self.trnaseq_db_path = os.path.join(self.out_dir, self.project_name + "-TRNASEQ.db")

        # Supplementary text file paths
        self.uniq_nontrna_path = os.path.join(self.out_dir, self.project_name + "-UNIQUED_NONTRNA.txt")
        self.uniq_trna_path = os.path.join(self.out_dir, self.project_name + "-UNIQUED_TRNA.txt")
        self.trimmed_ends_path = os.path.join(self.out_dir, self.project_name + "-TRIMMED_ENDS.txt")
        self.norm_frag_path = os.path.join(self.out_dir, self.project_name + "-NORMALIZED_FRAGMENTS.txt")

        # Intermediate pickle file paths
        self.profile_uniq_trna_seqs_path = os.path.join(self.out_dir, "UNIQUE_TRNA_SEQS-PROFILE_CHECKPOINT.pkl")
        self.profile_uniq_nontrna_seqs_path = os.path.join(self.out_dir, "UNIQUE_NONTRNA_SEQS-PROFILE_CHECKPOINT.pkl")
        self.threeprime_norm_uniq_trna_seqs_path = os.path.join(self.out_dir, "UNIQUE_TRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.threeprime_norm_uniq_nontrna_seqs_path = os.path.join(self.out_dir, "UNIQUE_NONTRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.threeprime_norm_trimmed_trna_seqs_path = os.path.join(self.out_dir, "TRIMMED_TRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.threeprime_norm_norm_trna_seqs_path = os.path.join(self.out_dir, "NORMALIZED_TRNA_SEQS-THREEPRIME_NORMALIZATION_CHECKPOINT.pkl")
        self.frag_map_uniq_trna_seqs_path = os.path.join(self.out_dir, "UNIQUE_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
        self.frag_map_uniq_nontrna_seqs_path = os.path.join(self.out_dir, "UNIQUE_NONTRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
        self.frag_map_trimmed_trna_seqs_path = os.path.join(self.out_dir, "TRIMMED_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
        self.frag_map_norm_trna_seqs_path = os.path.join(self.out_dir, "NORMALIZED_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")

        self.uniq_nontrna_seqs = []
        self.uniq_trna_seqs = []
        self.trimmed_trna_seqs = []
        self.norm_trna_seqs = []
        self.mod_trna_seqs = []

        self.counts_of_norm_seqs_containing_trimmed_seqs = [] # same length as self.trimmed_trna_seqs
        # "Multiplicity" of a trimmed tRNA sequence
        # = number of normalized sequences that the sequence is in * number of input sequences represented by the trimmed sequence
        self.multiplicities_of_trimmed_seqs_among_norm_seqs = [] # same length as self.trimmed_trna_seqs
        self.mean_multiplicities_of_norm_seqs = [] # same length as self.norm_trna_seqs
        self.mean_multiplicities_of_mod_seqs = []


    def sanity_check(self):
        """Check user inputs before proceeding."""

        if os.path.exists(self.out_dir):
            if len(os.listdir(self.out_dir)) == 0:
                # There is nothing in the output directory.
                # `anvi-run-workflow` creates the output directory even before `anvi-trnaseq` is called.
                pass
            elif self.overwrite_out_dest:
                if self.load_checkpoint:
                    raise ConfigError("You cannot use `--load-checkpoint` in conjunction with `--overwrite-output-destinations`. "
                                      "Starting at a checkpoint requires loading intermediate files "
                                      "written to the output directory in a previous `anvi-trnaseq` run, "
                                      "but this directory would be removed with `--overwrite-output-destinations`.")
                shutil.rmtree(self.out_dir)
            else:
                if not self.load_checkpoint:
                    raise ConfigError("The directory that was specified by --output-dir or -o, %s, already exists. "
                                      "Use the flag --overwrite-output-destinations to overwrite this directory." % self.out_dir)

        # Check that needed intermediate pickle files exist when loading from a checkpoint.
        missing_intermed_files = False
        if self.load_checkpoint == 'profile':
            if (not os.path.exists(self.profile_uniq_trna_seqs_path)
                or not os.path.exists(self.profile_uniq_nontrna_seqs_path)):
                missing_intermed_files = True
        elif self.load_checkpoint == 'threeprime_normalization':
            if (not os.path.exists(self.threeprime_norm_uniq_trna_seqs_path)
                or not os.path.exists(self.threeprime_norm_uniq_nontrna_seqs_path)
                or not os.path.exists(self.threeprime_norm_trimmed_trna_seqs_path)
                or not os.path.exists(self.threeprime_norm_norm_trna_seqs_path)):
                missing_intermed_files = True
        elif self.load_checkpoint == 'fragment_mapping':
            if (not os.path.exists(self.frag_map_uniq_trna_seqs_path)
                or not os.path.exists(self.frag_map_uniq_nontrna_seqs_path)
                or not os.path.exists(self.frag_map_trimmed_trna_seqs_path)
                or not os.path.exists(self.frag_map_norm_trna_seqs_path)):
                missing_intermed_files = True
        else:
            if not os.path.exists(self.out_dir):
                os.mkdir(self.out_dir)
        if missing_intermed_files:
            raise ConfigError("Intermediate files needed for running `anvi-trnaseq` with `--load-checkpoint %s` are missing. "
                              "You should probably run `anvi-trnaseq` from the beginning without `--load-checkpoint`. "
                              "To generate necessary intermediate files for future use of `--load-checkpoint`, use the flag `--write-checkpoints`."
                              % self.load_checkpoint)

        filesnpaths.is_output_dir_writable(self.out_dir)

        if self.descrip_path:
            filesnpaths.is_file_plain_text(self.descrip_path)
            self.descrip = self.descrip_path.read()
        else:
            self.descrip = None
        self.run.info("Description", self.descrip_path if self.descrip_path else "No description given")

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
                       'description': self.descrip if self.descrip else '_No description is provided_'}
        TRNASeqDatabase(self.trnaseq_db_path, quiet=False).create(meta_values)


    def unique_reads(self):
        """Dereplicate input reads."""

        self.progress.new("Finding replicate reads")
        self.progress.update("Loading reads")

        fasta = fastalib.SequenceSource(self.input_fasta_path)
        names = []
        seqs = []
        read_count = 0
        while next(fasta):
            names.append(fasta.id)
            seqs.append(fasta.seq)
            read_count += 1
        fasta.close()
        self.read_count = read_count

        clusters = Dereplicator(names, seqs, progress=self.progress).full_length_dereplicate()

        uniq_reads = [UniqueSeq(cluster.representative_seq_string,
                                cluster.member_names[0],
                                len(cluster.member_names))
                      for cluster in clusters]

        self.progress.end()

        return uniq_reads


    def profile_trna(self, uniq_reads):
        """Profile tRNA features in reads.

        Appends UniqueSeq objects representing profiled tRNA sequences to `self.uniq_trna_seqs`
        Appends leftover UniqueSeq objects representing unprofiled tRNA sequences to `self.uniq_nontrna_seqs`

        Parameters
        ==========
        uniq_reads : list
            List of UniqueSeq objects
        """

        self.progress.new("Profiling tRNA features in reads")
        self.progress.update("...")

        processed_read_count = 0
        processed_seq_count = 0
        trna_read_count = 0
        uniq_trna_count = 0
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
                                     for i in range(len(uniq_reads) // self.write_buffer_size)]
                                    + [len(uniq_reads), None])
        write_point = next(write_point_iterator)
        fetched_profile_count = 0
        uniq_reads_to_write_dict = {}
        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        for uniq_read in uniq_reads:
            input_queue.put((uniq_read.seq_string, uniq_read.represent_name))
            uniq_reads_to_write_dict[uniq_read.represent_name] = uniq_read
            processed_read_count += uniq_read.read_count
            processed_seq_count += 1

            if processed_seq_count == write_point:
                # Write a chunk of sequence results.
                write_point = next(write_point_iterator)

                # List of entries for each tRNA-seq database table
                trnaseq_sequences_table_entries = []
                trnaseq_feature_table_entries = []
                trnaseq_unconserved_table_entries = []
                trnaseq_unpaired_table_entries = []

                while fetched_profile_count < processed_seq_count:
                    trna_profile = output_queue.get()
                    fetched_profile_count += 1
                    output_name = trna_profile.name
                    output_seq = trna_profile.input_seq

                    uniq_seq = uniq_reads_to_write_dict.pop(output_name)

                    if not trna_profile.is_predicted_trna:
                        self.uniq_nontrna_seqs.append(uniq_seq)
                        continue

                    uniq_seq.id_method = 0 # 0 => profiled
                    uniq_seq.acceptor_length = len(trna_profile.acceptor_variant_string)
                    uniq_seq.has_complete_feature_set = trna_profile.has_complete_feature_set
                    uniq_seq.extra_fiveprime_length = trna_profile.num_extra_fiveprime
                    self.uniq_trna_seqs.append(uniq_seq)

                    uniq_trna_count += 1
                    output_seq_length = len(output_seq)

                    num_replicates = uniq_seq.read_count
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

                    trnaseq_feature_table_entries.append(
                        (output_name,
                         trna_profile.has_complete_feature_set,
                         trna_profile.anticodon_seq,
                         trna_profile.anticodon_aa,
                         output_seq_length,
                         # Zero-based start position of identified tRNA features within the read.
                         output_seq_length - len(trna_profile.profiled_seq),
                         # Stop position of features (real stop position, not Pythonic stop index for slicing).
                         output_seq_length - trna_profile.num_extra_threeprime - 1,
                         trna_profile.num_conserved,
                         trna_profile.num_unconserved,
                         trna_profile.num_paired,
                         trna_profile.num_unpaired,
                         trna_profile.num_in_extrapolated_fiveprime_feature,
                         trna_profile.num_extra_fiveprime,
                         trna_profile.num_extra_threeprime)
                        # When tRNA features were not found at the 5' end of the read,
                        # their start and stop positions also were not found.
                        + tuple(['?' * 2 for _ in range((len(TRNA_FEATURE_NAMES) - len(trna_profile.features)))]) * 2
                        + tuple(itertools.chain(*zip(
                            [str(feature.start_pos) if hasattr(feature, 'start_pos')
                             else ','.join(map(str, feature.start_positions))
                             for feature in trna_profile.features],
                            # Convert Pythonic stop position for slicing to real stop position of feature.
                            [str(feature.stop_pos - 1) if hasattr(feature, 'stop_pos')
                             else ','.join(map(str, [stop_pos - 1 for stop_pos in feature.stop_positions]))
                             for feature in trna_profile.features])))
                        # The alpha and beta sections of the D loop are not full-fledged features, but "subfeatures,"
                        # so add them as columns to the table after the features.
                        + (alpha_start,
                           alpha_stop,
                           beta_start,
                           beta_stop)
                    )

                    for unconserved_tuple in unconserved_info:
                        trnaseq_unconserved_table_entries.append((output_name, ) + unconserved_tuple)

                    for unpaired_tuple in unpaired_info:
                        trnaseq_unpaired_table_entries.append((output_name, ) + unpaired_tuple)

                if len(trnaseq_sequences_table_entries) > 0:
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('sequences',
                                                ','.join('?' * len(tables.trnaseq_sequences_table_structure))),
                                             trnaseq_sequences_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('feature',
                                                ','.join('?' * len(tables.trnaseq_feature_table_structure))),
                                             trnaseq_feature_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('feature_unconserved_nucleotides',
                                                ','.join('?' * len(tables.trnaseq_unconserved_table_structure))),
                                             trnaseq_unconserved_table_entries)
                    trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                             % ('feature_unpaired_nucleotides',
                                                ','.join('?' * len(tables.trnaseq_unpaired_table_structure))),
                                             trnaseq_unpaired_table_entries)

                self.progress.update("%d of %d unique sequences have been profiled"
                                     % (fetched_profile_count, len(uniq_reads)))

        for p in processes:
            p.terminate()
            p.join()

        # Profiled seqs were added to the output queue as they were processed, so sort by name.
        self.uniq_trna_seqs.sort(key=lambda uniq_seq: uniq_seq.represent_name)
        self.uniq_nontrna_seqs.sort(key=lambda uniq_seq: uniq_seq.represent_name)

        trnaseq_db.db.set_meta_value('reads_processed', processed_read_count)
        trnaseq_db.db.set_meta_value('unique_reads_processed', processed_seq_count)
        trnaseq_db.db.set_meta_value('trna_reads', trna_read_count)
        trnaseq_db.db.set_meta_value('unique_trna_seqs', uniq_trna_count)
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
        self.run.info("Unique profiled tRNA sequences", uniq_trna_count)
        self.run.info("Profiled reads with anticodon", trna_containing_anticodon_read_count)
        self.run.info("Profiled reads spanning acceptor stem", full_length_trna_read_count)
        self.run.info("Profiled reads with 1-3 extra 5' bases", trna_with_one_to_three_extra_fiveprime_bases_read_count)
        self.run.info("Profiled reads with >3 extra 5' bases", trna_with_more_than_three_extra_fiveprime_bases_read_count)
        self.run.info("Profiled reads with extrapolated 5' feature", trna_with_extrapolated_fiveprime_feature_read_count)
        self.run.info("Profiled reads ending in 3'-CCA", trna_with_threeprime_cca_read_count)
        self.run.info("Profiled reads ending in 3'-CC", trna_with_threeprime_cc_read_count)
        self.run.info("Profiled reads ending in 3'-C", trna_with_threeprime_c_read_count)
        self.run.info("Profiled reads ending in 3'-CCAN/CCANN", trna_with_threeprime_ccan_ccann_read_count)


    def trim_ends(self, uniq_trna_seqs):
        """Trim any nucleotides 5' of the acceptor stem and 3' of the discriminator.

        Appends TrimmedSeq objects formed from input UniqueSeq objects to `self.trimmed_trna_seqs`

        Parameters
        ==========
        uniq_trna_seqs : list
            List of UniqueSeq objects
        """

        self.progress.new("Trimming the 3' and 5' ends of sequences")
        self.progress.update("...")

        represent_names = [uniq_seq.represent_name for uniq_seq in uniq_trna_seqs]
        trimmed_seq_strings = [
            uniq_seq.seq_string[uniq_seq.extra_fiveprime_length: len(uniq_seq.seq_string) - uniq_seq.acceptor_length]
            for uniq_seq in uniq_trna_seqs
        ]

        clusters = Dereplicator(represent_names,
                                trimmed_seq_strings,
                                extras=uniq_trna_seqs,
                                progress=self.progress).full_length_dereplicate()

        trimmed_seqs = [TrimmedSeq(cluster.representative_seq_string, cluster.member_extras)
                        for cluster in clusters]

        self.trimmed_trna_seqs.extend(trimmed_seqs)

        self.trimmed_trna_seqs.sort(key=lambda trimmed_seq: trimmed_seq.represent_name)

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

        represent_names = [trimmed_seq.represent_name for trimmed_seq in self.trimmed_trna_seqs]
        # Reverse sequence orientation to dereplicate from the 3' end.
        reversed_seq_strings = [trimmed_seq.seq_string[::-1] for trimmed_seq in self.trimmed_trna_seqs]
        clusters = Dereplicator(represent_names,
                                reversed_seq_strings,
                                extras=self.trimmed_trna_seqs,
                                progress=self.progress).prefix_dereplicate()

        # Skip initialization of NormalizedSeq objects,
        # as additional TrimmedSeq constituents are added after mapping unprofiled fragments to profiled tRNA.
        self.norm_trna_seqs = [NormalizedSeq(cluster.member_extras, skip_init=True) for cluster in clusters]

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
        max_query_length = max(map(len, [seq.seq_string for seq in self.uniq_nontrna_seqs]))
        # By default, avoid sequences shorter than 25 nucleotides, the minimum length of a profiled 3' fragment of tRNA.
        interval_start_length = self.min_trna_frag_size
        interval_stop_length = interval_start_length + self.frag_mapping_query_chunk_length
        if interval_stop_length + self.frag_mapping_query_chunk_length > max_query_length:
            # Group the last chunk in with the second to last chunk if it is less than full size.
            interval_stop_length = max_query_length + 1
        query_length_intervals.append((interval_start_length, interval_stop_length))
        query_names = []
        query_seqs = []
        query_name_chunks = [query_names]
        query_seq_chunks = [query_seqs]
        for nontrna_index, seq in sorted([t for t in zip(range(len(self.uniq_nontrna_seqs)), self.uniq_nontrna_seqs)
                                          if len(t[1].seq_string) >= self.min_trna_frag_size],
                                         key=lambda t: len(t[1].seq_string)):
            if len(seq.seq_string) < interval_stop_length:
                query_names.append((seq.represent_name, nontrna_index))
                query_seqs.append(seq.seq_string)
            else:
                interval_start_length = interval_stop_length
                interval_stop_length = interval_stop_length + self.frag_mapping_query_chunk_length
                if interval_stop_length + self.frag_mapping_query_chunk_length > max_query_length:
                    interval_stop_length = max_query_length + 1
                query_length_intervals.append((interval_start_length, interval_stop_length))
                query_names = [(seq.represent_name, nontrna_index)]
                query_seqs = [seq.seq_string]
                query_name_chunks.append(query_names)
                query_seq_chunks.append(query_seqs)

        self.progress.update("Retrieving targets from profiled tRNA")
        # Leftover non-tRNA sequences are mapped to normalized tRNA sequences
        # with extra 5' bases added when present in underlying unique tRNA sequences.
        # Multiple targets for each normalized sequence are therefore produced for different 5' sequence extensions.
        target_names = []
        target_seqs = []
        for norm_seq_index, norm_seq in enumerate(self.norm_trna_seqs):
            norm_name = norm_seq.represent_name
            norm_seq_string = norm_seq.seq_string
            # The longest trimmed sequence (the first in the list) is by design
            # the only one of the profiled trimmed sequences forming the normalized sequence that may have extra 5' bases.
            longest_trimmed_seq = norm_seq.trimmed_seqs[0]
            if longest_trimmed_seq.uniq_with_extra_fiveprime_count > 0:
                fiveprime_seq_string_set = set()
                for uniq_seq in longest_trimmed_seq.uniq_seqs:
                    if uniq_seq.extra_fiveprime_length > 0:
                        fiveprime_seq_string_set.add(uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length])

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
                    target_names.append((norm_seq_index, len(fiveprime_seq_string), fiveprime_index))
                    target_seqs.append(fiveprime_seq_string + norm_seq_string)
            else:
                target_names.append((norm_seq_index, 0, 0)) # no extra 5' bases
                target_seqs.append(norm_seq_string)

        self.progress.end()


        interval_index = 0
        nontrna_indices = []
        for query_names, query_seqs in zip(query_name_chunks, query_seq_chunks):
            self.progress.new("Mapping %d unprofiled sequences of length %d-%d to profiled tRNA"
                              % (len(query_names),
                                 query_length_intervals[interval_index][0],
                                 query_length_intervals[interval_index][1] - 1))

            aligned_query_dict, aligned_target_dict = Aligner( # aligned_target_dict is not used for anything
                query_names,
                query_seqs,
                target_names,
                target_seqs,
                num_threads=self.num_threads,
                progress=self.progress
            ).align(max_mismatch_freq=0,
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
                    ref_alignment_start = alignment.target_start
                    ref_alignment_end = alignment.target_start + alignment.alignment_length

                    norm_seq_index, ref_fiveprime_length, _ = alignment.aligned_target.name # extra 5' index doesn't matter now

                    norm_end_pos = ref_alignment_end - ref_fiveprime_length
                    if norm_end_pos < 0:
                        # Ignore queries that align entirely to extra 5' bases.
                        continue

                    norm_seq = self.norm_trna_seqs[norm_seq_index]

                    if not trimmed_seq:
                        nontrna_indices.append(nontrna_index)

                        uniq_mapped_seq = self.uniq_nontrna_seqs[nontrna_index]
                        uniq_mapped_seq.id_method = 1 # 1 => mapped
                        uniq_mapped_seq.acceptor_length = 0
                        uniq_mapped_seq.has_complete_feature_set = False

                        # Assume that 5' extensions are the same for the query regardless of the reference.
                        # This could be false in the unlikely cases of
                        # 1. tRNA profiling erroneously identifying the end of the acceptor stem
                        # or 2. the query mapping to different places at the end of the acceptor stem in different tRNAs.
                        if ref_fiveprime_length - ref_alignment_start > 0:
                            uniq_mapped_seq.extra_fiveprime_length = ref_fiveprime_length - ref_alignment_start
                            norm_start_pos = 0
                        else:
                            uniq_mapped_seq.extra_fiveprime_length = 0
                            norm_start_pos = ref_alignment_start - ref_fiveprime_length

                        trimmed_seq = TrimmedSeq(uniq_mapped_seq.seq_string[uniq_mapped_seq.extra_fiveprime_length: ],
                                                 [uniq_mapped_seq])
                        self.trimmed_trna_seqs.append(trimmed_seq)

                        norm_seq.trimmed_seqs.append(trimmed_seq)
                        trimmed_seq.norm_seq_count += 1
                        norm_seq.start_positions.append(norm_start_pos)
                        norm_seq.end_positions.append(norm_end_pos)
                    else:
                        for prev_trimmed_seq in norm_seq.trimmed_seqs[::-1]:
                            # Ensure that the trimmed sequence maps to the normalized sequence only once.
                            # Multiple targets can be created from the same normalized sequence for different 5' extensions.
                            if prev_trimmed_seq.id_method == 0:
                                norm_seq.trimmed_seqs.append(trimmed_seq)
                                trimmed_seq.norm_seq_count += 1
                                if ref_fiveprime_length - ref_alignment_start > 0:
                                    norm_start_pos = 0
                                else:
                                    norm_start_pos = ref_alignment_start - ref_fiveprime_length
                                norm_seq.start_positions.append(norm_start_pos)
                                norm_seq.end_positions.append(norm_end_pos)
                                break
                            if trimmed_seq.represent_name == prev_trimmed_seq.represent_name:
                                break
            interval_index += 1

            del aligned_query_dict
            gc.collect()

            self.progress.end()

        for norm_seq in self.norm_trna_seqs:
            norm_seq.init()

        interior_mapped_count = 0
        fiveprime_mapped_count = 0
        for nontrna_index in sorted(nontrna_indices, reverse=True):
            uniq_mapped_seq = self.uniq_nontrna_seqs.pop(nontrna_index)

            if uniq_mapped_seq.extra_fiveprime_length > 0:
                fiveprime_mapped_count += uniq_mapped_seq.read_count
            else:
                interior_mapped_count += uniq_mapped_seq.read_count

        self.run.info("Mapped reads without extra 5' tRNA bases", interior_mapped_count)
        self.run.info("Mapped reads with extra 5' tRNA bases", fiveprime_mapped_count)


    def find_modifications(self):
        self.progress.new("Finding modifications")

        # Cluster normalized tRNA sequences.
        # Clusters agglomerate sequences that differ from at least one other sequence in the cluster
        # by no more than 2 substitutions per 71 aligned positions (by default) in a gapless end-to-end alignment.
        agglomerator = Agglomerator([seq.represent_name for seq in self.norm_trna_seqs],
                                    [seq.seq_string for seq in self.norm_trna_seqs],
                                    num_threads=self.num_threads,
                                    progress=self.progress)
        # Provide a priority function for seeding clusters
        # that favors fully profiled tRNA over "longer" tRNA without a full set of profiled features.
        # Such incompletely profiled longer tRNA includes tRNA-tRNA chimeras,
        # and some of these have a long 5' section that is a long 3' fragment of tRNA,
        # which can cause other shorter normalized sequences to agglomerate by aligning to the 5' section of the chimera.
        full_length_trna_dict = {seq.represent_name: seq.has_complete_feature_set
                                 for seq in self.norm_trna_seqs}
        agglomerator.agglomerate(max_mismatch_freq=self.agglom_max_mismatch_freq,
                                 priority_function=lambda aligned_ref: (-full_length_trna_dict[aligned_ref.name],
                                                                        -len(aligned_ref.seq_string),
                                                                        -len(aligned_ref.alignments),
                                                                        aligned_ref.name),
                                 alignment_target_chunk_size=self.alignment_target_chunk_size,
                                 alignment_progress_interval=self.alignment_progress_interval,
                                 agglom_progress_interval=self.agglom_progress_interval)

        agglom_aligned_ref_dict = agglomerator.agglom_aligned_ref_dict

        self.progress.update("Separating modification-induced substitutions from \"inter-strain\" variants")

        norm_seq_dict = {seq.represent_name: seq for seq in self.norm_trna_seqs}
        names_of_norm_seqs_assigned_to_mod_seqs = []
        for ref_name, aligned_ref in agglom_aligned_ref_dict.items():
            # A modification requires at least 3 different nucleotides to be detected,
            # and each normalized sequence differs by at least 1 nucleotide (substitution or gap),
            # so for a cluster to form a modified sequence, it must contain at least 3 normalized sequences.
            if len(aligned_ref.alignments) < 2:
                continue

            aligned_ref_length = len(aligned_ref.seq_string)

            valid_aligned_queries = []
            for alignment in aligned_ref.alignments:
                # Normalized tRNA sequences should only align at the 3' end.
                # Alignments to the interior of the sequence can occur when the reference is a tRNA-tRNA chimera.
                if aligned_ref_length != alignment.target_start + alignment.alignment_length:
                    continue

                query_name = alignment.aligned_query.name
                # The normalized sequence query may have agglomerated with another reference as well.
                # If the query formed a modified sequence,
                # it would form the same modified sequence when starting with this agglomeration.
                if query_name in names_of_norm_seqs_assigned_to_mod_seqs:
                    continue

                valid_aligned_queries.append(norm_seq_dict[query_name])

            # Confirm that at > 1 query passed the filters,
            # so at least 3 normalized sequences are still in the cluster.
            if len(valid_aligned_queries) < 2:
                continue

            seq_array = np.zeros((len(aligned_ref.alignments) + 1, aligned_ref_length), dtype=int)
            # Rather than using the ASCII representation of each character,
            # which saves some time in converting the sequence string to a numpy array,
            # constrain the integer representation to the smallest possible range of integers
            # to speed up the bincount method used to determine the number of unique nucleotides at an alignment position.
            seq_array[0, :] += [nt_int_dict[nt] for nt in aligned_ref.seq_string]
            for i, aligned_query in enumerate(valid_aligned_queries, start=1):
                seq_array[i, aligned_ref_length - len(aligned_query.seq_string): ] += [nt_int_dict[nt]
                                                                                       for nt in aligned_query.seq_string]

            norm_seqs = np.array([norm_seq_dict[ref_name]] + valid_aligned_queries)

            # Find positions in the alignment with nucleotide variability.
            alignment_pos_uniq_nt_counts = (
                np.bincount((seq_array + np.arange(aligned_ref_length, dtype=int) * num_nt_bins).ravel(),
                            minlength=aligned_ref_length * num_nt_bins).reshape(-1, num_nt_bins)[:, 1:] != 0
            ).sum(axis=1)
            three_four_nt_alignment_positions = (alignment_pos_uniq_nt_counts > 2).nonzero()[0]

            # Modification sites must have at least 3 nucleotides.
            if not three_four_nt_alignment_positions.size:
                continue

            two_nt_alignment_positions = (alignment_pos_uniq_nt_counts == 2).nonzero()[0]
            clusters = deque(((seq_array, norm_seqs, three_four_nt_alignment_positions), ))
            for alignment_pos in two_nt_alignment_positions:
                next_clusters = deque() # Make a new object with each iteration rather than clearing it.

                while clusters:
                    seq_array, norm_seqs, three_four_nt_alignment_positions = clusters.pop()

                    # A modification requires at least 3 different nucleotides to be detected,
                    # and each normalized sequence differs by at least 1 nucleotide (substitution or gap),
                    # so for a cluster to form a modified sequence, it must contain at least 3 normalized sequences.
                    if norm_seqs.size < 3:
                        continue

                    aligned_nts = seq_array[:, alignment_pos]
                    nt_counts = np.bincount(aligned_nts, minlength=num_nt_bins)[1: ]

                    if (nt_counts != 0).sum() < 2:
                        # There are now < 2 nucleotides at the alignment position in the (derived) cluster under consideration.
                        # 2 different nucleotides are needed to distinguish single nucleotide variants.
                        next_clusters.appendleft((seq_array, norm_seqs, three_four_nt_alignment_positions))
                        continue

                    # Add a new cluster for each nucleotide variant to the stack of clusters to process
                    # if the new cluster contains at least 3 sequences.
                    represented_nts = nt_counts.nonzero()[0] + 1
                    for nt in represented_nts:
                        split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]
                        if split_cluster_seq_indices.size > 2:
                            next_clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                      norm_seqs[split_cluster_seq_indices],
                                                      three_four_nt_alignment_positions))
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
                seq_array, norm_seqs, three_four_nt_alignment_positions = clusters.pop()
                candidates_to_remove = []

                for i, alignment_pos in enumerate(three_four_nt_alignment_positions):
                    aligned_nts = seq_array[:, alignment_pos]
                    nt_counts = np.bincount(aligned_nts, minlength=num_nt_bins)[1: ]
                    # At least 3 different nucleotides are needed at a position to predict a modification.
                    represented_nts = nt_counts.nonzero()[0] + 1
                    if represented_nts.size < 2:
                        candidates_to_remove.append(i)
                    elif represented_nts.size == 2:
                        candidates_to_remove.append(i)
                        split_three_four_nt_alignment_positions = np.delete(three_four_nt_alignment_positions,
                                                                            candidates_to_remove)
                        for nt in represented_nts:
                            split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]
                            # At least 3 normalized sequences are needed to form a modified sequence.
                            if split_cluster_seq_indices.size > 2:
                                clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                     norm_seqs[split_cluster_seq_indices],
                                                     split_three_four_nt_alignment_positions))
                        # Reevaluate previous alignment positions in the split clusters.
                        break
                else:
                    # At least 1 position was discounted as no longer having 3-4 different nucleotides,
                    # but these positions had fewer than 2 nucleotides,
                    # and so did not cause the cluster to be split into new clusters.
                    # Therefore, do not cycle through the remaining positions again to find those with fewer than 3 nucleotides.
                    if candidates_to_remove:
                        next_clusters.appendleft(
                            (norm_seqs, np.delete(three_four_nt_alignment_positions, candidates_to_remove))
                        )
                    else:
                        next_clusters.appendleft((norm_seqs, three_four_nt_alignment_positions))

            if not next_clusters:
                continue
            clusters = next_clusters

            while clusters:
                norm_seqs, mod_positions = clusters.pop()
                norm_seqs = sorted(norm_seqs, key=lambda seq: (-len(seq.seq_string), seq.represent_name)) # Turn the `norm_seqs` array into a list.
                represent_norm_seq_start_in_array = aligned_ref_length - len(norm_seqs[0].seq_string)
                mod_positions -= represent_norm_seq_start_in_array
                mod_seq = ModifiedSeq(norm_seqs, mod_positions.tolist())
                for norm_seq in norm_seqs:
                    names_of_norm_seqs_assigned_to_mod_seqs.append(norm_seq.represent_name)
                self.mod_trna_seqs.append(mod_seq)

        self.progress.update("Finding sequences with modification-induced deletions")

        # Search for normalized sequences with modification-induced deletions,
        # which are added to uninitialized modified sequences.
        # Exclude normalized sequences that have already been added to modified sequences from the search.
        norm_seq_targets = [norm_seq for norm_seq in self.norm_trna_seqs
                            if norm_seq.represent_name not in names_of_norm_seqs_assigned_to_mod_seqs]
        # Speed up dict lookup by searching modified sequences against normalized sequences of the same length.
        norm_seq_target_dict = defaultdict(dict)
        for norm_seq in norm_seq_targets:
            norm_seq_len = len(norm_seq.seq_string)
            norm_seq_target_dict[norm_seq_len][norm_seq.seq_string] = norm_seq

        query_dict = defaultdict(list)
        for mod_seq in self.mod_trna_seqs:
            for seq_string_with_del, del_config in mod_seq.get_seqs_with_dels():
                query_dict[seq_string_with_del].append((del_config, mod_seq))

        for seq_string_with_del, mod_seq_items in query_dict.items():
            d = norm_seq_target_dict[len(seq_string_with_del)]
            try:
                norm_seq = d[seq_string_with_del]
            except KeyError:
                continue
            for del_config, mod_seq in mod_seq_items:
                norm_seq.mod_seqs.append(mod_seq)
                mod_seq.norm_seqs_with_dels.append(norm_seq)
                mod_seq.del_configs.append(del_config)

        for mod_seq in self.mod_trna_seqs:
            mod_seq.init()

        self.progress.end()


    def calc_normalization_stats(self):
        self.progress.new("Calculating normalized tRNA stats")
        self.progress.update("...")

        # Count the normalized sequences containing each trimmed sequence.
        norm_count_dict = OrderedDict([(trimmed_seq.represent_name, 0)
                                       for trimmed_seq in self.trimmed_trna_seqs])
        for norm_seq in self.norm_trna_seqs:
            for trimmed_seq in norm_seq.trimmed_seqs:
                norm_count_dict[trimmed_seq.represent_name] += 1
        self.counts_of_norm_seqs_containing_trimmed_seqs = [norm_count for norm_count
                                                            in norm_count_dict.values()]

        # Find the "multiplicity" of trimmed sequences among normalized sequences.
        multiplicity_dict = OrderedDict()
        for trimmed_seq, norm_count_item in zip(self.trimmed_trna_seqs, norm_count_dict.items()):
            trimmed_represent_name, norm_count = norm_count_item
            multiplicity_dict[trimmed_represent_name] = trimmed_seq.read_count * norm_count
        self.multiplicities_of_trimmed_seqs_among_norm_seqs = [multiplicity for multiplicity
                                                               in multiplicity_dict.values()]

        # Find the "mean multiplicity" of normalized and modified sequences.
        for norm_seq in self.norm_trna_seqs:
            multiplicity_sum = 0
            for trimmed_seq in norm_seq.trimmed_seqs:
                multiplicity_sum += multiplicity_dict[trimmed_seq.represent_name]
            self.mean_multiplicities_of_norm_seqs.append(round(multiplicity_sum / norm_seq.read_count, 1))

        for mod_seq in self.mod_trna_seqs:
            multiplicity_sum = 0
            for norm_seq in mod_seq.norm_seqs_without_dels + mod_seq.norm_seqs_with_dels:
                for trimmed_seq in norm_seq.trimmed_seqs:
                    multiplicity_sum += multiplicity_dict[trimmed_seq.represent_name]
            self.mean_multiplicities_of_mod_seqs.append(
                round(multiplicity_sum / (mod_seq.specific_read_count + mod_seq.nonspecific_read_count), 1)
            )

        self.progress.end()


    def write_trimmed_table(self):
        self.progress.new("Writing tRNA-seq database table of trimmed tRNA sequences")
        self.progress.update("...")

        trimmed_table_entries = []
        for trimmed_seq, norm_seq_count in zip(self.trimmed_trna_seqs,
                                               self.counts_of_norm_seqs_containing_trimmed_seqs):
            trimmed_table_entries.append(
                (trimmed_seq.represent_name,
                 len(trimmed_seq.uniq_seqs),
                 trimmed_seq.read_count,
                 trimmed_seq.seq_string,
                 norm_seq_count,
                 trimmed_seq.uniq_with_extra_fiveprime_count,
                 trimmed_seq.read_with_extra_fiveprime_count)
                + tuple([v for v in trimmed_seq.read_acceptor_variant_count_dict.values()])
            )

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

        norm_table_entries = []
        for norm_seq, mean_multiplicity in zip(self.norm_trna_seqs,
                                               self.mean_multiplicities_of_norm_seqs):
            norm_table_entries.append(
                (norm_seq.represent_name,
                 len(norm_seq.trimmed_seqs),
                 norm_seq.read_count,
                 norm_seq.specific_covs,
                 norm_seq.nonspecific_covs,
                 mean_multiplicity,
                 norm_seq.trimmed_seqs_mapped_without_extra_fiveprime_count,
                 norm_seq.reads_mapped_without_extra_fiveprime_count,
                 norm_seq.trimmed_seqs_mapped_with_extra_fiveprime_count,
                 norm_seq.reads_mapped_with_extra_fiveprime_count)
                + tuple(norm_seq.read_acceptor_variant_count_dict.values())
            )

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('normalized')
            trnaseq_db.db.create_table('normalized',
                                       tables.trnaseq_normalized_table_structure,
                                       tables.trnaseq_normalized_table_types)
        trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                 % ('normalized', ','.join('?' * len(tables.trnaseq_normalized_table_structure))),
                                 norm_table_entries)

        norm_seq_count = len(self.norm_trna_seqs)
        trnaseq_db.db.set_meta_value('num_normalized_trna_seqs', norm_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()

        self.run.info("Normalized tRNA, consolidating tRNA fragments", norm_seq_count)


    def write_modified_table(self):
        self.progress.new("Writing tRNA-seq database table of modified tRNA sequences")
        self.progress.update("...")

        mod_table_entries = []
        for mod_seq, mean_multiplicity in zip(self.mod_trna_seqs,
                                              self.mean_multiplicities_of_mod_seqs):
            mod_table_entries.append(
                (mod_seq.represent_name,
                 ','.join([str(sub_pos) for sub_pos in mod_seq.sub_positions]) + ',')
                + tuple([','.join(map(str, mod_seq.specific_sub_covs[:, i - 1])) + ',' for i in int_nt_dict])
                + tuple([','.join(map(str, mod_seq.nonspecific_sub_covs[:, i - 1])) + ',' for i in int_nt_dict])
                + (';'.join(','.join(map(str, del_config)) for del_config in mod_seq.del_configs) + ',', )
                + (','.join(str(del_cov) for del_cov in mod_seq.specific_del_covs) + ',',
                   ','.join(str(del_cov) for del_cov in mod_seq.nonspecific_del_covs) + ',')
                + (mod_seq.consensus_seq_string,
                   len(mod_seq.norm_seqs_without_dels),
                   ','.join([norm_seq.represent_name for norm_seq in mod_seq.norm_seqs_without_dels]),
                   len(mod_seq.norm_seqs_with_dels),
                   ','.join([norm_seq.represent_name for norm_seq in mod_seq.norm_seqs_with_dels]),
                   mod_seq.specific_read_count,
                   mod_seq.nonspecific_read_count,
                   mean_multiplicity,
                   mod_seq.count_of_specific_reads_with_extra_fiveprime,
                   mod_seq.count_of_nonspecific_reads_with_extra_fiveprime,
                   mod_seq.specific_mapped_read_count,
                   mod_seq.nonspecific_mapped_read_count)
            )

        trnaseq_db = TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('modified')
            trnaseq_db.db.create_table('modified',
                                       tables.trnaseq_modified_table_structure,
                                       tables.trnaseq_modified_table_types)
        trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                 % ('modified', ','.join('?' * len(tables.trnaseq_modified_table_structure))),
                                 mod_table_entries)

        mod_seq_count = len(self.mod_trna_seqs)
        trnaseq_db.db.set_meta_value('num_mod_trna_seqs', mod_seq_count)
        trnaseq_db.disconnect()

        self.progress.end()

        self.run.info("Modified tRNA", mod_seq_count)


    def write_uniq_nontrna_supplement(self):
        self.progress.new("Writing a file of sequences not identified as tRNA.")
        self.run.info("Output non-tRNA file", self.uniq_nontrna_path)

        with open(self.uniq_nontrna_path, 'w') as nontrna_file:
            nontrna_file.write("\t".join(self.UNIQ_NONTRNA_HEADER) + "\n")
            for uniq_seq in self.uniq_nontrna_seqs:
                nontrna_file.write(uniq_seq.represent_name + "\t"
                                   + str(uniq_seq.read_count) + "\t"
                                   + uniq_seq.seq_string + "\n")

        self.progress.end()


    def write_trimmed_supplement(self):
        self.progress.new("Writing a file showing how trimmed tRNA sequences were formed from unique sequences")
        self.progress.update("...")

        with open(self.trimmed_ends_path, 'w') as trimmed_file:
            trimmed_file.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for trimmed_seq in sorted(self.trimmed_trna_seqs,
                                      key=lambda trimmed_seq: -trimmed_seq.read_count):
                represent_name = trimmed_seq.represent_name
                for uniq_seq in sorted(trimmed_seq.uniq_seqs,
                                       key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length,
                                                             -uniq_seq.acceptor_length)):
                    trimmed_file.write(represent_name + "\t"
                                       + uniq_seq.represent_name + "\t"
                                       + uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length] + "\t"
                                       + uniq_seq.seq_string[len(uniq_seq.seq_string) - uniq_seq.acceptor_length: ] + "\t"
                                       + str(uniq_seq.read_count) + "\n")

        self.progress.end()

        self.run.info("Output trimmed tRNA file", self.trimmed_ends_path)


    def write_normalized_supplement(self):
        self.progress.new("Writing a file showing how normalized tRNA sequences were formed from trimmed sequences")
        self.progress.update("...")

        with open(self.norm_frag_path, 'w') as norm_file:
            norm_file.write("\t".join(self.NORM_FRAG_HEADER) + "\n")
            for norm_seq in sorted(self.norm_trna_seqs, key=lambda norm_seq: -norm_seq.read_count):
                represent_name = norm_seq.represent_name
                for trimmed_seq, start_pos, end_pos in sorted(
                    zip(norm_seq.trimmed_seqs, norm_seq.start_positions, norm_seq.end_positions),
                    key=lambda t: (t[1], -t[2])):
                    norm_file.write(represent_name + "\t"
                                    + trimmed_seq.represent_name + "\t"
                                    + str(start_pos) + "\t"
                                    + str(end_pos) + "\n")

        self.progress.end()

        self.run.info("Output normalized tRNA file", self.norm_frag_path)


    def process(self):
        """The main method of TRNASeqDataset, called from `anvi-trnaseq`

        Checkpoint loading and saving occurs in this method.
        """
        self.sanity_check()

        # The first checkpoints are after tRNA profiling.
        if not self.load_checkpoint:
            self.create_trnaseq_db()

            # Profile each read for tRNA features.
            if self.feature_param_path:
                self.progress.new("Setting tRNA feature parameters from user file")
                self.progress.update("...")

                trnaidentifier.TRNAFeature.set_params_from_file(self.feature_param_path)

                self.progress.end()
            self.profile_trna(self.unique_reads())

            if self.write_checkpoints:
                self.progress.new("Writing intermediate files for the \"profile\" checkpoint")
                self.progress.update("...")

                if os.path.exists(self.profile_uniq_trna_seqs_path):
                    overwrote_profile_uniq_trna_seqs_path = True
                else:
                    overwrote_profile_uniq_trna_seqs_path = False
                with open(self.profile_uniq_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.uniq_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.profile_uniq_nontrna_seqs_path):
                    overwrote_profile_uniq_nontrna_seqs_path = True
                else:
                    overwrote_profile_uniq_nontrna_seqs_path = False
                with open(self.profile_uniq_nontrna_seqs_path, 'wb') as f:
                    pkl.dump(self.uniq_nontrna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                self.progress.end()

                self.run.info("%s\"profile\" checkpoint intermediate file of unique tRNA"
                              % ("Overwritten " if overwrote_profile_uniq_trna_seqs_path else ""),
                              self.profile_uniq_trna_seqs_path)
                self.run.info("%s\"profile\" checkpoint intermediate file of unique non-tRNA"
                              % ("Overwritten " if overwrote_profile_uniq_nontrna_seqs_path else ""),
                              self.profile_uniq_nontrna_seqs_path)

        if self.load_checkpoint == 'profile':
            self.progress.new("Loading intermediate files at the checkpoint, \"profile\"")
            self.progress.update("...")

            with open(self.profile_uniq_trna_seqs_path, 'rb') as f:
                self.uniq_trna_seqs = pkl.load(f)
            with open(self.profile_uniq_nontrna_seqs_path, 'rb') as f:
                self.uniq_nontrna_seqs = pkl.load(f)

            self.progress.end()

        if self.load_checkpoint == 'threeprime_normalization':
            self.progress.new("Loading intermediate files at the checkpoint, \"threeprime_normalization\"")
            self.progress.update("...")

            with open(self.threeprime_norm_uniq_trna_seqs_path, 'rb') as f:
                self.uniq_trna_seqs = pkl.load(f)
            with open(self.threeprime_norm_uniq_nontrna_seqs_path, 'rb') as f:
                self.uniq_nontrna_seqs = pkl.load(f)
            with open(self.threeprime_norm_trimmed_trna_seqs_path, 'rb') as f:
                self.trimmed_trna_seqs = pkl.load(f)
            with open(self.threeprime_norm_norm_trna_seqs_path, 'rb') as f:
                self.norm_trna_seqs = pkl.load(f)

            self.progress.end()
        elif self.load_checkpoint == 'fragment_mapping':
            pass
        else:
            # Trim 5' and 3' ends of profiled tRNA.
            self.trim_ends(self.uniq_trna_seqs)

            # Consolidate 3' fragments of longer tRNA sequences.
            self.dereplicate_threeprime()

            if self.write_checkpoints:
                self.progress.new("Writing intermediate files for the \"threeprime_normalization\" checkpoint")
                self.progress.update("...")

                if os.path.exists(self.threeprime_norm_uniq_trna_seqs_path):
                    overwrote_threeprime_norm_uniq_trna_seqs_path = True
                else:
                    overwrote_threeprime_norm_uniq_trna_seqs_path = False
                with open(self.threeprime_norm_uniq_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.uniq_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.threeprime_norm_uniq_nontrna_seqs_path):
                    overwrote_threeprime_norm_uniq_nontrna_seqs_path = True
                else:
                    overwrote_threeprime_norm_uniq_nontrna_seqs_path = False
                with open(self.threeprime_norm_uniq_nontrna_seqs_path, 'wb') as f:
                    pkl.dump(self.uniq_nontrna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.threeprime_norm_trimmed_trna_seqs_path):
                    overwrote_threeprime_norm_trimmed_trna_seqs_path = True
                else:
                    overwrote_threeprime_norm_trimmed_trna_seqs_path = False
                with open(self.threeprime_norm_trimmed_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.trimmed_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.threeprime_norm_norm_trna_seqs_path):
                    overwrote_threeprime_norm_norm_trna_seqs_path = True
                else:
                    overwrote_threeprime_norm_norm_trna_seqs_path = False
                with open(self.threeprime_norm_norm_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.norm_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                self.progress.end()

                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of unique tRNA"
                              % ("Overwritten " if overwrote_threeprime_norm_uniq_trna_seqs_path else ""),
                              self.threeprime_norm_uniq_trna_seqs_path)
                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of unique non-tRNA"
                              % ("Overwritten " if overwrote_threeprime_norm_uniq_nontrna_seqs_path else ""),
                              self.threeprime_norm_uniq_nontrna_seqs_path)
                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of trimmed tRNA"
                              % ("Overwritten " if overwrote_threeprime_norm_trimmed_trna_seqs_path else ""),
                              self.threeprime_norm_trimmed_trna_seqs_path)
                self.run.info("%s\"threeprime_normalization\" checkpoint intermediate file of normalized tRNA"
                              % ("Overwritten " if overwrote_threeprime_norm_norm_trna_seqs_path else ""),
                              self.threeprime_norm_norm_trna_seqs_path)

        if self.load_checkpoint == 'fragment_mapping':
            self.progress.new("Loading intermediate files at the checkpoint, \"fragment_mapping\"")
            self.progress.update("...")

            with open(self.frag_map_uniq_trna_seqs_path, 'rb') as f:
                self.uniq_trna_seqs = pkl.load(f)
            with open(self.frag_map_uniq_nontrna_seqs_path, 'rb') as f:
                self.uniq_nontrna_seqs = pkl.load(f)
            with open(self.frag_map_trimmed_trna_seqs_path, 'rb') as f:
                self.trimmed_trna_seqs = pkl.load(f)
            with open(self.frag_map_norm_trna_seqs_path, 'rb') as f:
                self.norm_trna_seqs = pkl.load(f)

            self.progress.end()
        else:
            # Map fragments derived from the interior and 5' end of tRNA.
            self.map_fragments()

            if self.write_checkpoints:
                self.progress.new("Writing intermediate files for the \"fragment_mapping\" checkpoint")
                self.progress.update("...")

                if os.path.exists(self.frag_map_uniq_trna_seqs_path):
                    overwrote_frag_map_uniq_trna_seqs_path = True
                else:
                    overwrote_frag_map_uniq_trna_seqs_path = False
                with open(self.frag_map_uniq_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.uniq_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.frag_map_uniq_nontrna_seqs_path):
                    overwrote_frag_map_uniq_nontrna_seqs_path = True
                else:
                    overwrote_frag_map_uniq_nontrna_seqs_path = False
                with open(self.frag_map_uniq_nontrna_seqs_path, 'wb') as f:
                    pkl.dump(self.uniq_nontrna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.frag_map_trimmed_trna_seqs_path):
                    overwrote_frag_map_trimmed_trna_seqs_path = True
                else:
                    overwrote_frag_map_trimmed_trna_seqs_path = False
                with open(self.frag_map_trimmed_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.trimmed_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                if os.path.exists(self.frag_map_norm_trna_seqs_path):
                    overwrote_frag_map_norm_trna_seqs_path = True
                else:
                    overwrote_frag_map_norm_trna_seqs_path = False
                with open(self.frag_map_norm_trna_seqs_path, 'wb') as f:
                    pkl.dump(self.norm_trna_seqs, f, protocol=pkl.HIGHEST_PROTOCOL)

                self.progress.end()

                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of unique tRNA"
                              % ("Overwritten " if overwrote_frag_map_uniq_trna_seqs_path else ""),
                              self.threeprime_norm_uniq_trna_seqs_path)
                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of unique non-tRNA"
                              % ("Overwritten " if overwrote_frag_map_uniq_nontrna_seqs_path else ""),
                              self.frag_map_uniq_nontrna_seqs_path)
                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of trimmed tRNA"
                              % ("Overwritten " if overwrote_frag_map_trimmed_trna_seqs_path else ""),
                              self.frag_map_trimmed_trna_seqs_path)
                self.run.info("%s\"fragment_mapping\" checkpoint intermediate file of normalized tRNA"
                              % ("Overwritten " if overwrote_frag_map_norm_trna_seqs_path else ""),
                              self.frag_map_norm_trna_seqs_path)

        # Find modified nucleotides, grouping sequences into modified sequences.
        self.find_modifications()

        # Calculate some statistics.
        self.calc_normalization_stats()

        # Write more tables to the database.
        self.write_trimmed_table()
        self.write_normalized_table()
        self.write_modified_table()

        # Write supplementary text files.
        self.write_uniq_nontrna_supplement()
        self.write_trimmed_supplement()
        self.write_normalized_supplement()


def profile_worker(input_queue, output_queue):
    while True:
        seq_string, seq_name = input_queue.get()
        output_queue.put(trnaidentifier.Profile(seq_string, name=seq_name))
