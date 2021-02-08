# -*- coding: utf-8
# pylint: disable=line-too-long
"""Library for tRNA-seq dataset operations bin/anvi-trnaseq and bin/anvi-convert-trnaseq-database
are the default clients using this module. anvi-trnaseq instantiates a TRNASeqDataset.
anvi-convert-trnaseq-database instantiates a DatabaseConverter. The scripts call the objects'
process() methods to start the analytic workflows.

Each sequence library in an experiment is processed separately as a TRNASeqDataset, storing an
information-rich anvi'o tRNA-seq database. DatabaseConverter finds reference seed sequences from a
set of tRNA-seq databases, storing seeds in an anvi'o contigs database and coverage patterns for
each dataset in anvi'o profile and auxiliary databases. Contigs and profile databases interface with
a range of other tools in the anvi'o platform.
"""

import gc
import os
import sys
import time
import shutil
import random
import hashlib
import argparse
import itertools
import numpy as np
import pandas as pd
import pickle as pkl
import multiprocessing as mp

from functools import partial
from itertools import combinations, product
from collections import OrderedDict, deque, defaultdict

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.tables as tables
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.tables.miscdata as miscdata
import anvio.trnaidentifier as trnaidentifier
import anvio.auxiliarydataops as auxiliarydataops

from anvio.errors import ConfigError
from anvio.agglomeration import Agglomerator
from anvio.tables.views import TablesForViews
from anvio.sequence import Aligner, Dereplicator
from anvio.tables.miscdata import TableForLayerOrders


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


pp = terminal.pretty_print

THREEPRIME_VARIANTS = constants.THREEPRIME_VARIANTS

UNAMBIG_NTS = ('A', 'C', 'G', 'T')
NT_INT_DICT = {nt: i for i, nt in enumerate(UNAMBIG_NTS, start=1)}
INT_NT_DICT = {i: nt for i, nt in enumerate(UNAMBIG_NTS, start=1)}
# NUM_NT_BINS is used in counting the number of distinct nucleotides (1, 2, 3, or 4) at positions in
# internally ungapped alignments: there is one bin (value 0) for end gaps in the alignment.
NUM_NT_BINS = len(UNAMBIG_NTS) + 1

ANTICODON_TO_AA = constants.anticodon_to_AA


class UniqueSeq(object):
    """A dereplicated tRNA-seq read that can contain additional information from tRNA feature
    profiling.
    """

    __slots__ = (
        'seq_string',
        'represent_name',
        'read_count',
        'id_method',
        'feature_start_indices',
        'feature_stop_indices',
        'has_complete_feature_set',
        'has_his_g',
        'alpha_start_index',
        'alpha_stop_index',
        'beta_start_index',
        'beta_stop_index',
        'acceptor_length',
        'anticodon_string',
        'anticodon_aa',
        'contains_anticodon',
        'num_conserved',
        'num_unconserved',
        'num_paired',
        'num_unpaired',
        'unconserved_info',
        'unpaired_info',
        'profiled_seq_length',
        'num_extrapolated_fiveprime_nts',
        'extra_fiveprime_length',
        'extra_threeprime_length',
        'trunc_profile_index',
        'trunc_profile_recovered_by_derep',
        'trunc_profile_recovered_by_del_analysis'
    )

    def __init__(self, seq_string, represent_name, read_count):
        self.seq_string = seq_string
        self.represent_name = represent_name
        self.read_count = read_count
        self.id_method = None # If dealing with tRNA, identification method 0 = profiled, 1 = mapped
        self.feature_start_indices = None
        self.feature_stop_indices = None
        self.has_complete_feature_set = None
        self.has_his_g = None
        self.alpha_start_index = None
        self.alpha_stop_index = None
        self.beta_start_index = None
        self.beta_stop_index = None
        self.acceptor_length = None
        self.anticodon_string = None
        self.anticodon_aa = None
        self.contains_anticodon = None
        self.num_conserved = None
        self.num_unconserved = None
        self.num_paired = None
        self.num_unpaired = None
        self.unconserved_info = None
        self.unpaired_info = None
        self.profiled_seq_length = None
        self.num_extrapolated_fiveprime_nts = None
        self.extra_fiveprime_length = None
        self.extra_threeprime_length = None
        self.trunc_profile_index = None # Index at which feature profile is truncated (None maintained if not truncated)
        self.trunc_profile_recovered_by_derep = None
        self.trunc_profile_recovered_by_del_analysis = None


class TrimmedSeq(object):
    """A tRNA sequence with bases trimmed 5' of the acceptor stem (or 5'-G in the case of tRNA-His)
    and 3' of the discriminator.

    The purpose of trimming is to collapse non-biological variability prevalent at the ends of
    reads.

    EXAMPLE 1:
    E. coli tRNA-Ala-GGC-1-1
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA

    This collapses to the following trimmed sequence, removing the acceptor (the acceptor happens to
    be genomic rather than post-transcriptionally added in this example, but it doesn't matter):
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA

    Examples of possible profiled reads that collapse to this sequence:
    AGGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCAC

    EXAMPLE 2:
    3' fragment of the same tRNA, ending in 3'-CC rather than canonical 3'-CCA
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA

    This collapses to the following trimmed sequence, removing 3'-CC:
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA

    Examples of possible profiled reads that collapse to this sequence:
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCC
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTC

    Unique MAPPED tRNA fragments, including tRNA reads that only lack 3'-CCA/CC/C, each generate a
    TrimmedSeq object (in TRNASeqDataset.map_fragments and
    TRNASeqDataset.threeprime_dereplicate_acceptorless_sequences). Mapped fragments with extra 5'
    nucleotides beyond the acceptor stem are not trimmed. The 5' extension may represent all but a
    small number of nucleotides in the sequence, so it is best not to dereplicate mapped sequences
    identical in the non-5' section by grouping them as the same trimmed sequence.
    """

    __slots__ = (
        'seq_string',
        'uniq_seqs',
        'feature_start_indices',
        'feature_stop_indices',
        'has_complete_feature_set',
        'has_his_g',
        'read_count',
        'uniq_with_extra_fiveprime_count',
        'read_with_extra_fiveprime_count',
        'represent_name',
        'long_fiveprime_extension_dict',
        'read_acceptor_variant_count_dict',
        'id_method',
        'contains_anticodon',
        'has_extrapolated_fiveprime_nts',
        'has_trunc_profile',
        'trunc_profile_recovered_by_derep',
        'trunc_profile_recovered_by_del_analysis',
        'norm_seq_count'
    )

    # The user can specify what defines a long (biological vs. non-templated) 5' extension.
    min_length_of_long_fiveprime_extension = 4

    def __init__(self, seq_string, uniq_seqs, skip_init=False):
        self.seq_string = seq_string
        self.uniq_seqs = uniq_seqs # list of UniqueSeq objects
        represent_uniq_seq = uniq_seqs[0]
        # Assume that the feature profile indices of the representative unique sequence are the same
        # as the other unique sequences. The acceptor is the last feature in the profile and not
        # part of the trimmed sequence.
        self.feature_start_indices = represent_uniq_seq.feature_start_indices[:-1] if represent_uniq_seq.feature_start_indices else None
        self.feature_stop_indices = represent_uniq_seq.feature_stop_indices[:-1] if represent_uniq_seq.feature_stop_indices else None
        # Assume that if the representative unique sequence has a complete feature set, then so do
        # the other unique sequences.
        self.has_complete_feature_set = represent_uniq_seq.has_complete_feature_set
        self.has_his_g = represent_uniq_seq.has_his_g
        self.trunc_profile_recovered_by_derep = None
        self.trunc_profile_recovered_by_del_analysis = None
        self.norm_seq_count = 0

        if skip_init:
            self.read_count = None
            self.uniq_with_extra_fiveprime_count = None
            self.read_with_extra_fiveprime_count = None
            self.represent_name = None
            self.long_fiveprime_extension_dict = None
            self.read_acceptor_variant_count_dict = None
            self.id_method = None
            self.contains_anticodon = None
            self.has_extrapolated_fiveprime_nts = None
            self.has_trunc_profile = None
        else:
            self.init()


    def init(self):
        """Set attributes for a "finalized" set of input `UniqueSeq` objects."""
        self.read_count = sum([uniq_seq.read_count for uniq_seq in self.uniq_seqs])

        self.uniq_with_extra_fiveprime_count = sum(
            [1 if uniq_seq.extra_fiveprime_length else 0 for uniq_seq in self.uniq_seqs])
        self.read_with_extra_fiveprime_count = sum(
            [uniq_seq.read_count if uniq_seq.extra_fiveprime_length else 0 for uniq_seq in self.uniq_seqs])

        self.represent_name = self.get_represent_name()

        self.long_fiveprime_extension_dict, self.read_acceptor_variant_count_dict = self.get_end_counts()

        id_method_set = set(uniq_seq.id_method for uniq_seq in self.uniq_seqs)
        if len(id_method_set) == 1:
            self.id_method = id_method_set.pop()
        else:
            raise ConfigError("A TrimmedSeq should not be made from UniqueSeq objects with different identification methods. "
                              "Trimmed tRNA sequences will EITHER be formed from \"profiled\" tRNA sequences or \"mapped\" "
                              "tRNA sequences, because they are of different lengths and are fragments from different parts "
                              "of the tRNA.")

        self.contains_anticodon = self.uniq_seqs[0].contains_anticodon
        self.has_extrapolated_fiveprime_nts = True if self.uniq_seqs[0].num_extrapolated_fiveprime_nts else False

        uniq_seq_trunc_profile_statuses = list(set([0 if uniq_seq.trunc_profile_index is None else 1 for uniq_seq in self.uniq_seqs]))
        if len(uniq_seq_trunc_profile_statuses) == 1:
            self.has_trunc_profile = True if uniq_seq_trunc_profile_statuses[0] else False
        else:
            raise ConfigError("A TrimmedSeq should not be made from some UniqueSeq objects "
                              "with a truncated feature profile and some without.")


    def get_represent_name(self):
        """The representative name of the trimmed sequence is chosen as follows:
        1. Most abundant full-length tRNA (no extra 5' bases), ignoring acceptor sequence, OR
        2. Most abundant longer-than-full-length tRNA, OR
        3. Most abundant fragmentary tRNA
        """
        # Sort unique sequences such that the first sequence is the most abundant+longest and the
        # last is the least abundant+shortest.
        uniq_seqs = sorted(
            self.uniq_seqs, key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, -uniq_seq.read_count))

        if uniq_seqs[0].extra_fiveprime_length > 0:
            if uniq_seqs[-1].extra_fiveprime_length == 0:
                # If the first unique sequence has extra 5' nucleotides and the last has none, then
                # the last sequence and others without extra 5' nucleotides must be a full-length
                # tRNA (ignoring the acceptor). Therefore, select the most abundant of these
                # full-length tRNAs as the representative sequence.
                represent_name = sorted(
                    uniq_seqs, key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, uniq_seq.read_count)
                )[-1].represent_name
            else:
                represent_name = uniq_seqs[0].represent_name
        else:
            # In this case, ALL unique sequences are EITHER full-length OR a fragment.
            represent_name = uniq_seqs[0].represent_name

        return represent_name


    def get_end_counts(self):
        """Get the number of reads containing each unique 5' extension and sequence variant 3' of
        the discriminator that is collapsed into the trimmed sequence.

        There is a known number of possible 3' variants. Each variant is counted though it may not
        be represented in the trimmed sequence.
        """
        long_fiveprime_extension_dict = {}
        read_acceptor_variant_count_dict = OrderedDict(
            [(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        for uniq_seq in self.uniq_seqs:
            if uniq_seq.extra_fiveprime_length > self.min_length_of_long_fiveprime_extension:
                long_fiveprime_extension_dict[
                    uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length]
                ] = uniq_seq.read_count

            if uniq_seq.acceptor_length: # UniqueSeq need not have an acceptor
                acceptor_seq_string = uniq_seq.seq_string[-uniq_seq.acceptor_length: ]
                read_acceptor_variant_count_dict[acceptor_seq_string] += uniq_seq.read_count
        self.long_fiveprime_extension_dict = long_fiveprime_extension_dict
        self.read_acceptor_variant_count_dict = read_acceptor_variant_count_dict

        return long_fiveprime_extension_dict, read_acceptor_variant_count_dict


class NormalizedSeq(object):
    """A longer tRNA sequence consolidated from shorter tRNA fragments.

    Normalized sequences are derived from trimmed tRNA sequences. Trimmed profiled tRNA are first
    prefix-dereplicated from the 3' end by the method, `TRNASeqDataset.dereplicate_threeprime`. The
    longest sequence in a cluster of dereplicated sequences becomes the representative normalized
    sequence. `TRNASeqDataset.map_fragments` subsequently maps unprofiled reads to the set of
    normalized sequences. Mapped tRNA fragments are added as trimmed sequences, and the `init`
    method is called to finalize attributes of each normalized sequence.

    EXAMPLE:
    Consider the full-length and fragmentary E. coli tRNA-Ala-GGC-1-1 used in TrimmedSeq examples.
    GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA
                               TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA

    The sequences collapse into a single normalized sequence when dereplicated from the 3' end. The
    normalized sequence is represented by the longer sequence, which should be the first in the list
    of trimmed sequences added to the normalized sequence.
    """

    __slots__ = (
        'trimmed_seqs',
        'represent_name',
        'seq_string',
        'feature_start_indices',
        'feature_stop_indices',
        'has_complete_feature_set',
        'has_his_g',
        'has_trunc_profile',
        'start_positions',
        'stop_positions',
        'specific_read_count',
        'nonspecific_read_count',
        'count_of_specific_reads_with_extra_fiveprime',
        'count_of_nonspecific_reads_with_extra_fiveprime',
        'specific_mapped_read_count',
        'nonspecific_mapped_read_count',
        'specific_long_fiveprime_extension_dict',
        'nonspecific_long_fiveprime_extension_dict',
        'specific_read_acceptor_variant_count_dict',
        'nonspecific_read_acceptor_variant_count_dict',
        'specific_covs',
        'nonspecific_covs',
        'mean_specific_cov',
        'mean_nonspecific_cov',
        'mod_seqs',
        'profile_changed_by_del_analysis',
        'trunc_profile_recovered_by_del_analysis'
    )

    def __init__(self,
                 trimmed_seqs,
                 start_positions=None,
                 stop_positions=None,
                 skip_init=False):
        # In practice, a list of prefix-dereplicated, trimmed, profiled tRNA sequences is added when
        # a normalized sequence is instantiated, and mapped unprofiled sequences are appended later.
        self.trimmed_seqs = trimmed_seqs
        for trimmed_seq in trimmed_seqs:
            trimmed_seq.norm_seq_count += 1
        represent_trimmed_seq = trimmed_seqs[0]
        self.represent_name = represent_trimmed_seq.represent_name
        self.seq_string = represent_trimmed_seq.seq_string
        self.feature_start_indices = represent_trimmed_seq.feature_start_indices
        self.feature_stop_indices = represent_trimmed_seq.feature_stop_indices
        self.has_complete_feature_set = represent_trimmed_seq.has_complete_feature_set
        self.has_his_g = represent_trimmed_seq.has_his_g
        self.has_trunc_profile = represent_trimmed_seq.has_trunc_profile

        if start_positions and stop_positions:
            self.start_positions = start_positions
            self.stop_positions = stop_positions
        elif (not start_positions) and (not stop_positions):
            # Trimmed seqs were dereplicated from the 3' end of the normalized sequence. This is how
            # initialization occurs in the workflow, as mapped tRNA sequences that are not aligned
            # at the 3' end of the normalized sequence are added later.
            norm_seq_length = len(self.seq_string)
            self.start_positions = [norm_seq_length - len(trimmed_seq.seq_string)
                                    for trimmed_seq in self.trimmed_seqs]
            self.stop_positions = [norm_seq_length] * len(trimmed_seqs)
        else:
            self.start_positions = None
            self.stop_positions = None

        # It is useful to know which modified sequences, if any, contain this normalized sequence. A
        # normalized sequence with deletions, unlike a normalized sequence without deletions, can be
        # assigned to more than one modified sequence.
        self.mod_seqs = []
        self.profile_changed_by_del_analysis = None
        self.trunc_profile_recovered_by_del_analysis = None

        if skip_init:
            self.specific_read_count = None
            self.nonspecific_read_count = None
            self.count_of_specific_reads_with_extra_fiveprime = None
            self.count_of_nonspecific_reads_with_extra_fiveprime = None
            self.specific_mapped_read_count = None
            self.nonspecific_mapped_read_count = None
            self.specific_long_fiveprime_extension_dict = None
            self.nonspecific_long_fiveprime_extension_dict = None
            self.specific_read_acceptor_variant_count_dict = None
            self.nonspecific_read_acceptor_variant_count_dict = None
            self.specific_covs = None
            self.nonspecific_covs = None
            self.mean_specific_cov = None
            self.mean_nonspecific_cov = None
        else:
            self.init()


    def init(self):
        """Set attributes for a "finalized" set of input `TrimmedSeq` objects."""
        # Specific reads are those that are only contained in this normalized sequence.
        specific_read_count = 0
        nonspecific_read_count = 0
        count_of_specific_reads_with_extra_fiveprime = 0
        count_of_nonspecific_reads_with_extra_fiveprime = 0
        specific_mapped_read_count = 0
        nonspecific_mapped_read_count = 0
        specific_long_fiveprime_extension_dict = defaultdict(int)
        nonspecific_long_fiveprime_extension_dict = defaultdict(int)
        specific_read_acceptor_variant_count_dict = OrderedDict(
            [(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        nonspecific_read_acceptor_variant_count_dict = OrderedDict(
            [(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        specific_covs = np.zeros(len(self.seq_string), dtype=int)
        nonspecific_covs = np.zeros(len(self.seq_string), dtype=int)

        for trimmed_seq, start_pos, stop_pos in zip(
            self.trimmed_seqs, self.start_positions, self.stop_positions):
            if trimmed_seq.norm_seq_count == 1:
                specific_read_count += trimmed_seq.read_count

                if trimmed_seq.uniq_with_extra_fiveprime_count > 0:
                    count_of_specific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count

                if trimmed_seq.id_method == 1: # 1 => mapped
                    # TrimmedSeqs are comprised of EITHER profiled (0) OR mapped (1) UniqueSeqs.
                    specific_mapped_read_count += trimmed_seq.read_count

                for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                    specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

                for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                    if read_count > 0:
                        specific_read_acceptor_variant_count_dict[acceptor_seq_string] += read_count

                specific_covs[start_pos: stop_pos] += trimmed_seq.read_count
            else:
                nonspecific_read_count += trimmed_seq.read_count

                if trimmed_seq.uniq_with_extra_fiveprime_count > 0:
                    count_of_nonspecific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count

                if trimmed_seq.id_method == 1:
                    nonspecific_mapped_read_count += trimmed_seq.read_count

                    # Only nonspecific MAPPED sequences can have 5' sequence extensions, as profiled
                    # sequences with 5' extensions would span the length of the normalized sequence
                    # and would thus be specific to it.
                    for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                        nonspecific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

                for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                    if read_count > 0:
                        nonspecific_read_acceptor_variant_count_dict[acceptor_seq_string] += read_count

                nonspecific_covs[start_pos: stop_pos] += trimmed_seq.read_count

        self.specific_read_count = specific_read_count
        self.nonspecific_read_count = nonspecific_read_count
        self.count_of_specific_reads_with_extra_fiveprime = count_of_specific_reads_with_extra_fiveprime
        self.count_of_nonspecific_reads_with_extra_fiveprime = count_of_nonspecific_reads_with_extra_fiveprime
        self.specific_mapped_read_count = specific_mapped_read_count
        self.nonspecific_mapped_read_count = nonspecific_mapped_read_count
        self.specific_long_fiveprime_extension_dict = specific_long_fiveprime_extension_dict
        self.nonspecific_long_fiveprime_extension_dict = nonspecific_long_fiveprime_extension_dict
        self.specific_read_acceptor_variant_count_dict = specific_read_acceptor_variant_count_dict
        self.nonspecific_read_acceptor_variant_count_dict = nonspecific_read_acceptor_variant_count_dict
        self.specific_covs = specific_covs
        self.nonspecific_covs = nonspecific_covs
        self.mean_specific_cov = specific_covs.mean()
        self.mean_nonspecific_cov = nonspecific_covs.mean()


class ModifiedSeq(object):
    """A tRNA sequence with sites of predicted modification-induced substitutions and deletions,
    formed from normalized sequences with distinct patterns of modification-induced mutations.

    The workflow aggregates similar normalized sequences and processes the resulting clusters to
    produce clusters of normalized sequences distinguished by potential modification-induced
    substitutions (3-4 different nucleotides at one or more aligned positions). A modified sequence
    is initialized by a list of clustered normalized sequences, with the first sequence in the list
    being the longest, or tied for longest, and substitution positions indexed relative to this
    longest sequence. The workflow later finds sequences with modification-induced deletions (which
    occur in the vicinity of substitutions) and adds them to the modified sequences. After adding
    sequences with deletions, the workflow calls `init` to calculate coverages and other
    information.

    EXAMPLE:
    Consider E. coli tRNA-Ala-GGC-1-1, with detected modifications at positions 17 and 46. As seen
    in the normalized sequence example, the first sequence is the normalized sequence with unmutated
    nucleotides. The next set of sequences are other possible normalized sequences with
    modification-induced substitutions. The last set of sequences are possible normalized sequences
    with modification-induced deletions.

                     |                              |
    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA

    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG A TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG C TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC C GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA
              CTCAGC G GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA
                                    CATGGCATGCAAGAG T TCAGCGGTTCGATCCCGCTTAGCTCCA

    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG - TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGA- - TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGA- G TCAGCGGTTCGATCCCGCTTAGCTCCA
    """

    __slots__ = (
        'norm_seqs_without_dels',
        'sub_positions',
        'represent_name',
        'norm_seqs_with_dels',
        'del_configs',
        'specific_read_count',
        'nonspecific_read_count',
        'count_of_specific_reads_with_extra_fiveprime',
        'count_of_nonspecific_reads_with_extra_fiveprime',
        'specific_mapped_read_count',
        'nonspecific_mapped_read_count',
        'specific_long_fiveprime_extension_dict',
        'nonspecific_long_fiveprime_extension_dict',
        'specific_read_acceptor_variant_count_dict',
        'nonspecific_read_acceptor_variant_count_dict',
        'specific_covs',
        'nonspecific_covs',
        'mean_specific_cov',
        'mean_nonspecific_cov',
        'specific_sub_covs',
        'nonspecific_sub_covs',
        'specific_del_covs',
        'nonspecific_del_covs',
        'consensus_seq_string'
    )

    def __init__(self, norm_seqs_without_dels, sub_positions, skip_init=True):
        self.norm_seqs_without_dels = norm_seqs_without_dels
        self.sub_positions = sub_positions
        # A normalized sequence without modification-induced deletions can only be assigned to one
        # modified sequence, but a normalized sequence with deletions can be assigned to more than
        # one.
        for norm_seq in norm_seqs_without_dels:
            norm_seq.mod_seqs.append(self)
        self.represent_name = norm_seqs_without_dels[0].represent_name

        if skip_init:
            self.norm_seqs_with_dels = []
            # The deletion configurations of normalized sequences with deletions are recorded in the
            # modified sequence rather than the normalized sequence, because the positions of the
            # deletions correspond to nucleotides in the modified sequence.
            self.del_configs = []
            self.specific_read_count = None
            self.nonspecific_read_count = None
            self.count_of_specific_reads_with_extra_fiveprime = None
            self.count_of_nonspecific_reads_with_extra_fiveprime = None
            self.specific_mapped_read_count = None
            self.nonspecific_mapped_read_count = None
            self.specific_long_fiveprime_extension_dict = None
            self.nonspecific_long_fiveprime_extension_dict = None
            self.specific_read_acceptor_variant_count_dict = None
            self.nonspecific_read_acceptor_variant_count_dict = None
            self.specific_covs = None
            self.nonspecific_covs = None
            self.mean_specific_cov = None
            self.mean_nonspecific_cov = None
            self.specific_sub_covs = None
            self.nonspecific_sub_covs = None
            self.specific_del_covs = None
            self.nonspecific_del_covs = None
            self.consensus_seq_string = None
        else:
            self.init()


    def init(self):
        """Set attributes for a "finalized" set of input `NormalizedSeq` objects, some of which may
        contain deletions."""
        norm_seqs_without_dels = self.norm_seqs_without_dels
        norm_seqs_with_dels = self.norm_seqs_with_dels
        all_norm_seqs = norm_seqs_without_dels + norm_seqs_with_dels
        del_configs = self.del_configs
        specific_read_count = 0
        nonspecific_read_count = 0
        count_of_specific_reads_with_extra_fiveprime = 0
        count_of_nonspecific_reads_with_extra_fiveprime = 0
        specific_mapped_read_count = 0
        nonspecific_mapped_read_count = 0
        specific_long_fiveprime_extension_dict = defaultdict(int)
        nonspecific_long_fiveprime_extension_dict = defaultdict(int)
        specific_read_acceptor_variant_count_dict = OrderedDict(
            [(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        nonspecific_read_acceptor_variant_count_dict = OrderedDict(
            [(threeprime_variant, 0) for threeprime_variant in THREEPRIME_VARIANTS])
        mod_seq_len = len(norm_seqs_without_dels[0].seq_string)
        norm_seq_specific_covs = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        norm_seq_nonspecific_covs = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        num_subs = len(self.sub_positions)
        self.specific_sub_covs = specific_sub_covs = np.zeros((num_subs, len(UNAMBIG_NTS)), dtype=int)
        self.nonspecific_sub_covs = nonspecific_sub_covs = np.zeros((num_subs, len(UNAMBIG_NTS)), dtype=int)
        del_positions = sorted(set([i for del_config in del_configs for i in del_config]))
        self.specific_del_covs = specific_del_covs = np.zeros(len(del_positions), dtype=int)
        self.nonspecific_del_covs = nonspecific_del_covs = np.zeros(len(del_positions), dtype=int)

        # Make an array of aligned nucleotide positions in all normalized sequences.
        norm_seq_array = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        for n, norm_seq in enumerate(norm_seqs_without_dels):
            norm_seq_array[n, mod_seq_len - len(norm_seq.seq_string): ] += [
                NT_INT_DICT[nt] for nt in norm_seq.seq_string]

        nt_positions_covered_by_norm_seqs_with_dels = []
        n = len(norm_seqs_without_dels)
        for norm_seq, del_config in zip(norm_seqs_with_dels, del_configs):
            aligned_seq = [NT_INT_DICT[nt] for nt in norm_seq.seq_string]
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

            for trimmed_seq, trimmed_seq_start_in_norm_seq, trimmed_seq_stop_in_norm_seq in zip(
                norm_seq.trimmed_seqs, norm_seq.start_positions, norm_seq.stop_positions):
                if trimmed_seq.represent_name in processed_trimmed_seq_names:
                    continue

                # Determine whether the reads constituting the trimmed sequence are specific to the
                # modified sequence.
                if trimmed_seq.norm_seq_count == 1:
                    is_trimmed_seq_specific_to_mod_seq = True
                else:
                    # The trimmed sequence is specific to the modified sequence if it is unique to a
                    # set of normalized sequences specific to the modified sequence. Coverage
                    # information for such trimmed sequences will be recorded in the row of the
                    # array for the first normalized sequence in which it was found.
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
                        raise ConfigError("The number of normalized sequences containing the trimmed sequence was "
                                          "somehow miscalculated.")

                trimmed_seq_start_in_mod_seq = norm_seq_start_in_mod_seq + trimmed_seq_start_in_norm_seq
                trimmed_seq_stop_in_mod_seq = trimmed_seq_start_in_mod_seq + len(trimmed_seq.seq_string)
                if is_trimmed_seq_specific_to_mod_seq:
                    specific_read_count += trimmed_seq.read_count
                    count_of_specific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    # Trimmed sequences are comprised of either profiled or mapped unique sequences.
                    if trimmed_seq.id_method == 1: # 1 => mapped
                        specific_mapped_read_count += trimmed_seq.read_count

                    for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                        specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

                    for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                        if read_count > 0:
                            specific_read_acceptor_variant_count_dict[acceptor_seq_string] += read_count

                    norm_seq_specific_covs[n, trimmed_seq_start_in_mod_seq: trimmed_seq_stop_in_mod_seq] += trimmed_seq.read_count
                else:
                    nonspecific_read_count += trimmed_seq.read_count
                    count_of_nonspecific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    if trimmed_seq.id_method == 1:
                        nonspecific_mapped_read_count += trimmed_seq.read_count

                        # Only mapped sequences can have 5' sequence extensions, as profiled
                        # sequences with 5' extensions would span the length of the normalized
                        # sequence and would thus be specific to it.
                        for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                            nonspecific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

                    for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                        if read_count > 0:
                            nonspecific_read_acceptor_variant_count_dict[acceptor_seq_string] += read_count

                    norm_seq_nonspecific_covs[n, trimmed_seq_start_in_mod_seq: trimmed_seq_stop_in_mod_seq] += trimmed_seq.read_count

                processed_trimmed_seq_names.append(trimmed_seq.represent_name)

        # Handle normalized sequences with deletions.
        n = len(norm_seqs_without_dels)
        for norm_seq, del_config, nt_positions_covered_by_norm_seq in zip(
            norm_seqs_with_dels, del_configs, nt_positions_covered_by_norm_seqs_with_dels):
            norm_seq_start_in_mod_seq = mod_seq_len - len(norm_seq.seq_string) - len(del_config)
            # Normalized sequences with deletions can be found in multiple modified sequences,
            # unlike normalized sequences without deletions.
            num_mod_seqs_containing_norm_seq = len(norm_seq.mod_seqs)

            for trimmed_seq, trimmed_seq_start_in_norm_seq, trimmed_seq_stop_in_norm_seq in zip(
                norm_seq.trimmed_seqs, norm_seq.start_positions, norm_seq.stop_positions):
                if trimmed_seq.represent_name in processed_trimmed_seq_names:
                    continue

                # Determine whether the reads constituting the trimmed sequence are specific to the
                # modified sequence.
                if num_mod_seqs_containing_norm_seq > 1:
                    is_trimmed_seq_specific_to_mod_seq = False
                else:
                    # The trimmed sequence is specific to the modified sequence if it is unique to a
                    # set of normalized sequences specific to the modified sequence.
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
                            # The trimmed sequence was found in a normalized sequence that is in
                            # multiple modified sequences.
                            break
                    else:
                        if num_specific_norm_seqs_containing_trimmed_seq == trimmed_seq.norm_seq_count:
                            is_trimmed_seq_specific_to_mod_seq = True
                        elif num_specific_norm_seqs_containing_trimmed_seq < trimmed_seq.norm_seq_count:
                            # The trimmed sequence was found in normalized sequences that are not in
                            # the modified sequence.
                            is_trimmed_seq_specific_to_mod_seq = False
                        else:
                            raise ConfigError("The number of normalized sequences containing the trimmed sequence "
                                              "was somehow miscalculated.")

                nt_positions_covered_by_trimmed_seq = nt_positions_covered_by_norm_seq[
                    trimmed_seq_start_in_norm_seq: trimmed_seq_start_in_norm_seq + len(trimmed_seq.seq_string)]

                trimmed_seq_read_count = trimmed_seq.read_count
                if is_trimmed_seq_specific_to_mod_seq:
                    specific_read_count += trimmed_seq_read_count
                    count_of_specific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    # Trimmed sequences are comprised of either profiled or mapped unique sequences.
                    if trimmed_seq.id_method == 1: # 1 => mapped
                        specific_mapped_read_count += trimmed_seq_read_count

                    # Long 5' extensions are not counted for sequences with deletions. Recording the
                    # 5' extensions would occur here if the complex task were undertaken.

                    for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                        if read_count > 0:
                            specific_read_acceptor_variant_count_dict[acceptor_seq_string] += read_count

                    norm_seq_specific_covs[n, nt_positions_covered_by_trimmed_seq] += trimmed_seq_read_count
                    for del_pos in del_config:
                        specific_del_covs[del_positions.index(del_pos)] += trimmed_seq_read_count
                else:
                    nonspecific_read_count += trimmed_seq_read_count
                    count_of_nonspecific_reads_with_extra_fiveprime += trimmed_seq.read_with_extra_fiveprime_count
                    if trimmed_seq.id_method == 1:
                        nonspecific_mapped_read_count += trimmed_seq_read_count

                        # Long 5' extensions are not counted for sequences with deletions. Recording
                        # the 5' extensions would occur here if the complex task were undertaken.

                    for acceptor_seq_string, read_count in trimmed_seq.read_acceptor_variant_count_dict.items():
                        if read_count > 0:
                            nonspecific_read_acceptor_variant_count_dict[acceptor_seq_string] += read_count

                    norm_seq_nonspecific_covs[n, nt_positions_covered_by_trimmed_seq] += trimmed_seq_read_count
                    for del_pos in del_config:
                        nonspecific_del_covs[del_positions.index(del_pos)] += trimmed_seq_read_count

            processed_trimmed_seq_names.append(trimmed_seq.represent_name)
            n += 1

        self.specific_read_count = specific_read_count
        self.nonspecific_read_count = nonspecific_read_count
        self.count_of_specific_reads_with_extra_fiveprime = count_of_specific_reads_with_extra_fiveprime
        self.count_of_nonspecific_reads_with_extra_fiveprime = count_of_nonspecific_reads_with_extra_fiveprime
        self.specific_mapped_read_count = specific_mapped_read_count
        self.nonspecific_mapped_read_count = nonspecific_mapped_read_count
        self.specific_long_fiveprime_extension_dict = specific_long_fiveprime_extension_dict
        self.nonspecific_long_fiveprime_extension_dict = nonspecific_long_fiveprime_extension_dict
        self.specific_read_acceptor_variant_count_dict = specific_read_acceptor_variant_count_dict
        self.nonspecific_read_acceptor_variant_count_dict = nonspecific_read_acceptor_variant_count_dict
        self.specific_covs = norm_seq_specific_covs.sum(0)
        self.nonspecific_covs = norm_seq_nonspecific_covs.sum(0)
        self.mean_specific_cov = self.specific_covs.mean()
        self.mean_nonspecific_cov = self.nonspecific_covs.mean()

        # For each substitution position, record the coverage of A, C, G, and T.
        for s, sub_pos in enumerate(self.sub_positions):
            aligned_nts = norm_seq_array[:, sub_pos]
            nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]
            for nt_int, nt_count in enumerate(nt_counts, start=1):
                if nt_count > 0:
                    norm_seq_rows_with_nt = (aligned_nts == nt_int).nonzero()[0]
                    specific_sub_covs[s, nt_int - 1] = norm_seq_specific_covs[norm_seq_rows_with_nt, sub_pos].sum()
                    nonspecific_sub_covs[s, nt_int - 1] = norm_seq_nonspecific_covs[norm_seq_rows_with_nt, sub_pos].sum()

        # Set a consensus sequence from the nucleotides with the highest specific coverage at each
        # position.
        consensus_seq_string = norm_seqs_without_dels[0].seq_string
        for sub_pos, covs in zip(self.sub_positions, specific_sub_covs):
            max_pos = covs.argmax()
            nt_int = max_pos + 1
            consensus_seq_string = consensus_seq_string[: sub_pos] + INT_NT_DICT[nt_int] + consensus_seq_string[sub_pos + 1: ]
        self.consensus_seq_string = consensus_seq_string


class TRNASeqDataset(object):
    """Processes reads from a tRNA-seq library. `bin/anvi-trnaseq` is the client."""

    TRNA_FEATURE_NAMES = constants.TRNA_FEATURE_NAMES

    # Column headers for supplementary tables written to text files
    UNIQ_NONTRNA_HEADER = [
        "representative_name",
        "read_count",
        "truncated_profile_index",
        "sequence"
    ]
    TRIMMED_ENDS_HEADER = [
        "representative_name",
        "unique_name",
        "fiveprime_sequence",
        "threeprime_sequence",
        "read_count"
    ]

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # Argument group 1A: MANDATORY
        self.input_fasta_path = A('trnaseq_fasta')
        self.sample_id = A('sample_name')
        self.out_dir = os.path.abspath(A('output_dir')) if A('output_dir') else None

        # Argument group 1B: EXTRAS
        self.treatment = A('treatment')
        self.overwrite_out_dest = A('overwrite_output_destinations')
        self.descrip_path = os.path.abspath(A('description')) if A('description') else None

        # Argument group 1C: ADVANCED
        self.write_checkpoints = A('write_checkpoints')
        self.load_checkpoint = A('load_checkpoint')
        self.feature_param_path = os.path.abspath(A('feature_param_file')) if A('feature_param_file') else None
        self.min_length_of_long_fiveprime_extension = A('min_length_long_fiveprime')
        TrimmedSeq.min_length_of_long_fiveprime_extension = self.min_length_of_long_fiveprime_extension
        self.min_trna_frag_size = A('min_trna_fragment_size')
        self.agglom_max_mismatch_freq = A('agglomeration_max_mismatch_freq')
        self.skip_INDEL_profiling = A('skip_INDEL_profiling')
        self.fiveprimemost_del_start = A('fiveprimemost_deletion_start')
        self.threeprimemost_del_start = A('threeprimemost_deletion_start')
        self.fiveprimemost_del_stop = A('fiveprimemost_deletion_stop')
        self.threeprimemost_del_stop = A('threeprimemost_deletion_stop')
        self.max_distinct_dels = A('max_distinct_deletions')

        # Argument group 1D: PERFORMANCE
        self.num_threads = A('num_threads')
        self.skip_fasta_check = A('skip_fasta_check')
        self.alignment_target_chunk_size = A('alignment_target_chunk_size')
        self.frag_mapping_query_chunk_length = A('fragment_mapping_query_chunk_length')

        # Argument group 1E: PROGRESS
        self.profiling_progress_interval = A('profiling_progress_interval')
        self.alignment_progress_interval = A('alignment_progress_interval')
        self.agglom_progress_interval = A('agglomeration_progress_interval')

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using --fasta-file or -f.")
        if not self.sample_id:
            raise ConfigError("Please provide a sample name using --sample-name or -S.")
        if not self.out_dir:
            raise ConfigError("Please provide an output directory using --output-dir or -o.")

        self.descrip = None

        get_out_dir_path = partial(os.path.join, self.out_dir)

        self.trnaseq_db_path = get_out_dir_path(self.sample_id + "-TRNASEQ.db")

        self.analysis_summary_path = get_out_dir_path(self.sample_id + "-ANALYSIS_SUMMARY.txt")

        # Supplementary text file paths
        self.uniq_nontrna_path = get_out_dir_path(self.sample_id + "-UNIQUED_NONTRNA.txt")
        self.trimmed_ends_path = get_out_dir_path(self.sample_id + "-TRIMMED_ENDS.txt")

        # Intermediate pickle file paths
        self.intermed_file_path_dict = {
            'profile': {
                'uniq_trna_seqs': get_out_dir_path("UNIQUE_TRNA_SEQS-PROFILE_CHECKPOINT.pkl"),
                'uniq_trunc_seqs': get_out_dir_path("UNIQUE_TRUNCATED_SEQS-PROFILE_CHECKPOINT.pkl"),
                'uniq_nontrna_seqs': get_out_dir_path("UNIQUE_NONTRNA_SEQS-PROFILE_CHECKPOINT.pkl"),
                'trimmed_trna_seqs': get_out_dir_path("TRIMMED_TRNA_SEQS-PROFILE_CHECKPOINT.pkl"),
                'trimmed_trunc_seqs': get_out_dir_path("TRIMMED_TRUNC_SEQS-PROFILE_CHECKPOINT.pkl"),
                'norm_trna_seqs': get_out_dir_path("NORMALIZED_TRNA_SEQS-PROFILE_CHECKPOINT.pkl")
            },
            'fragment_mapping': {
                'uniq_trna_seqs': get_out_dir_path("UNIQUE_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl"),
                'uniq_trunc_seqs': get_out_dir_path("UNIQUE_TRUNCATED_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl"),
                'uniq_nontrna_seqs': get_out_dir_path("UNIQUE_NONTRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl"),
                'trimmed_trna_seqs': get_out_dir_path("TRIMMED_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl"),
                'trimmed_trunc_seqs': get_out_dir_path("TRIMMED_TRUNC_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl"),
                'norm_trna_seqs': get_out_dir_path("NORMALIZED_TRNA_SEQS-FRAGMENT_MAPPING_CHECKPOINT.pkl")
            }
        }
        self.intermed_file_label_dict = {
            'uniq_trna_seqs': 'unique tRNA',
            'uniq_trunc_seqs': 'unique sequences with a truncated feature profile',
            'uniq_nontrna_seqs': 'unique non-tRNA',
            'trimmed_trna_seqs': 'trimmed tRNA',
            'trimmed_trunc_seqs': 'trimmed sequences with a truncated feature profile',
            'norm_trna_seqs': 'normalized tRNA'
        }

        self.uniq_nontrna_seqs = []
        self.uniq_trunc_seqs = []
        self.uniq_trna_seqs = []
        self.trimmed_trna_seqs = []
        self.trimmed_trunc_seqs = []
        self.norm_trna_seqs = []
        self.norm_trunc_seqs = []
        self.mod_trna_seqs = []

        self.possible_del_starts = None
        self.possible_del_stops = None
        # Ranges representing the locations of deletions in relation to substitution sites. For
        # example, (range(-1, 0), range(-1, 1), range(0, 1)) allows three types of deletions to be
        # introduced at a substitution position: a 1 nucleotide deletion of the adjacent 5'
        # nucleotide a 2 nucleotide deletion of the adjacent 5' nucleotide and the nucleotide at the
        # substitution position, and a 1 nucleotide deletion of the nucleotide at the substitution
        # position.
        self.del_ranges = None


    def process(self):
        """The entry method of TRNASeqDataset, called from `anvi-trnaseq`.

        Checkpoint loading and saving occurs in this method.
        """
        total_time_start = time.time()
        self.sanity_check()

        # The first checkpoint occurs after tRNA profiling, tRNA trimming, and the recovery of tRNA
        # sequences with truncated feature profiles.
        if not self.load_checkpoint:
            self.create_trnaseq_db()

            # Profile each read for tRNA features.
            if self.feature_param_path:
                trnaidentifier.TRNAFeature.set_params_from_file(self.feature_param_path)
            self.report_profiling_params()
            self.report_fragment_mapping_params()
            self.report_modification_analysis_params()
            self.profile_trna(self.unique_reads())

            # Trim 5' and 3' ends of profiled tRNA.
            self.trim_trna_ends()
            # Trim 3' ends of sequences with truncated tRNA profile.
            self.trim_truncated_profile_ends()

            # Consolidate 3' fragments of longer profiled tRNA sequences.
            self.threeprime_dereplicate_trna()

            # Write tRNA feature profile information to the database.
            self.write_feature_table()
            self.write_unconserved_table()
            self.write_unpaired_table()

            # Recover tRNA sequences with truncated feature profiles by comparing to normalized
            # sequences.
            self.threeprime_dereplicate_truncated_sequences()

            # Write intermediate files at the "profile" checkpoint
            self.write_checkpoint_files('profile')

        if self.load_checkpoint == 'profile':
            self.load_checkpoint_files('profile')
            self.report_fragment_mapping_params()
            self.report_modification_analysis_params()

        if self.load_checkpoint == 'fragment_mapping':
            self.load_checkpoint_files('fragment_mapping')
            self.report_modification_analysis_params()
        else:
            # Reach this point either by starting from the beginning of the workflow
            # or loading the 'profile' checkpoint.

            # Recover 3' tRNA sequences lacking an acceptor sequence.
            self.threeprime_dereplicate_acceptorless_sequences()

            # Map fragments derived from the interior and 5' end of tRNA.
            self.map_fragments()

            if self.write_checkpoints:
                self.write_checkpoint_files('fragment_mapping')

        # Find modified nucleotides, grouping sequences into modified sequences.
        self.find_substitutions()
        if not self.skip_INDEL_profiling:
            self.find_deletions()
        for mod_seq in self.mod_trna_seqs:
            mod_seq.init()

        self.report_stats()

        # Write more tables to the database.
        self.write_sequences_table()
        self.write_trimmed_table()
        self.write_normalized_table()
        self.write_modified_table()

        # Write supplementary text files.
        self.write_uniq_nontrna_supplement()
        self.write_trimmed_supplement()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Total time elapsed (min)",
                                          time.time() - total_time_start,
                                          is_time_value=True))
            # Write an empty line to separate this run from any subsequent run starting from a
            # checkpoint writing to the same summary file.
            f.write("\n")


    def sanity_check(self):
        """Check user inputs before proceeding."""
        if os.path.exists(self.out_dir):
            if len(os.listdir(self.out_dir)) == 0:
                # There is nothing in the output directory. This may be caused by
                # `anvi-run-workflow`, which creates the output directory even before `anvi-trnaseq`
                # is called.
                pass
            elif self.overwrite_out_dest:
                if self.load_checkpoint:
                    raise ConfigError("You cannot use `--load-checkpoint` in conjunction with `--overwrite-output-destinations`. "
                                      "Starting at a checkpoint requires loading intermediate files written to the output "
                                      "directory in a previous `anvi-trnaseq` run, but this directory would be removed with "
                                      "`--overwrite-output-destinations`.")
                shutil.rmtree(self.out_dir)
            else:
                if not self.load_checkpoint:
                    raise ConfigError("The directory that was specified by --output-dir or -o, %s, already exists. "
                                      "Use the flag --overwrite-output-destinations to overwrite this directory." % self.out_dir)

        if self.load_checkpoint:
            # Check that needed intermediate pickle files exist when loading from a checkpoint.
            missing_intermed_files = []
            for intermed_file_path in self.intermed_file_path_dict[self.load_checkpoint].values():
                if not os.path.exists(intermed_file_path):
                    missing_intermed_files.append(intermed_file_path)
            if missing_intermed_files:
                raise ConfigError("Intermediate files needed for running `anvi-trnaseq` with `--load-checkpoint %s` are missing: %s. "
                                  "You should probably run `anvi-trnaseq` from the beginning without `--load-checkpoint`. "
                                  "To generate necessary intermediate files for future use of `--load-checkpoint`, use the flag `--write-checkpoints`."
                                  % (self.load_checkpoint, ', '.join(self.missing_intermed_files)))
        else:
            if not os.path.exists(self.out_dir):
                os.mkdir(self.out_dir)

        filesnpaths.is_output_dir_writable(self.out_dir)

        if self.descrip_path:
            filesnpaths.is_file_plain_text(self.descrip_path)
            self.descrip = open(self.descrip_path).read()

        if not 1 <= self.num_threads <= mp.cpu_count():
            raise ConfigError("The number of threads to use must be a positive integer less than or equal to %d. "
                              "Try again!" % mp.cpu_count())

        if (self.fiveprimemost_del_start > self.threeprimemost_del_start
            or self.fiveprimemost_del_start > self.fiveprimemost_del_stop
            or self.threeprimemost_del_start > self.threeprimemost_del_stop
            or self.fiveprimemost_del_stop > self.threeprimemost_del_stop):
            raise ConfigError("The following relations of deletion start and stop boundaries must be maintained: "
                              "1. 5'-most deletion start <= 3'-most deletion start, "
                              "2. 5'-most deletion start <= 5'-most deletion stop, "
                              "3. 3'-most deletion start <= 3'-most deletion stop, "
                              "4. 5'-most deletion stop <= 3'-most deletion stop. "
                              f"You provided a 5'-most deletion start of {self.fiveprimemost_del_start}, "
                              f"a 3'-most deletion start of {self.threeprimemost_del_start}, "
                              f"a 5'-most deletion stop of {self.fiveprimemost_del_stop}, and "
                              f"a 3'-most deletion stop of {self.threeprimemost_del_stop}.")
        self.possible_del_starts = tuple(range(self.fiveprimemost_del_start, self.threeprimemost_del_start + 1))
        self.possible_del_stops = tuple(range(self.fiveprimemost_del_stop, self.threeprimemost_del_stop + 1))
        # Determine the spectrum of different deletion sizes/positions relative to a substitution.
        del_ranges = []
        for del_start in self.possible_del_starts:
            for del_stop in self.possible_del_stops:
                if del_start < del_stop:
                    del_ranges.append(range(del_start, del_stop))
        self.del_ranges = tuple(del_ranges)

        self.run.info("Input FASTA file", self.input_fasta_path)

        if not self.skip_fasta_check and not self.load_checkpoint:
            self.progress.new("Checking input FASTA defline format")
            self.progress.update("...")

            utils.check_fasta_id_formatting(self.input_fasta_path)

            self.progress.end()

            self.run.info_single("FASTA deflines were found to be anvi'o-compliant", mc='green')


    def create_trnaseq_db(self):
        """Create an empty tRNA-seq database."""
        meta_values = {'sample_id': self.sample_id,
                       'treatment': self.treatment,
                       'description': self.descrip if self.descrip else '_No description is provided_',
                       'INDELs_profiled': not self.skip_INDEL_profiling}
        dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=False).create(meta_values)


    def get_summary_line(self, label, value, is_time_value=False, padding=68):
        """Return a string formatted to be written to the summary statistics file."""
        # Report elapsed time in seconds in minutes.
        if is_time_value:
            value = "%.2f" % round(value / 60, 2)
        return '%s%s\t%s\n' % (label, ' ' + '.' * (padding - len(label)), value)


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

        uniq_reads = [UniqueSeq(cluster.member_seqs[0], cluster.member_names[0], len(cluster.member_names))
                      for cluster in clusters]

        self.progress.end()
        return uniq_reads


    def report_profiling_params(self):
        """Add profiling parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        parameterizer = trnaidentifier.TRNAFeatureParameterizer()
        for param_tuple in parameterizer.list_accessible_param_tuples():
            trnaseq_db.db.set_meta_value(param_tuple[0], param_tuple[1])
        trnaseq_db.db.set_meta_value('min_length_of_long_fiveprime_extension', self.min_length_of_long_fiveprime_extension)
        trnaseq_db.disconnect()

        with open(self.analysis_summary_path, 'a') as f:
            for param_tuple in parameterizer.list_accessible_param_tuples(pretty=True):
                f.write(self.get_summary_line(param_tuple[0], param_tuple[1]))
            f.write(self.get_summary_line("Min length of \"long\" 5' extension", self.min_length_of_long_fiveprime_extension))


    def report_fragment_mapping_params(self):
        """Add fragment mapping parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db.set_meta_value('min_mapped_trna_fragment_size', self.min_trna_frag_size)
        trnaseq_db.disconnect()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Min length of mapped tRNA fragment", self.min_trna_frag_size))


    def report_modification_analysis_params(self):
        """Add modification analysis parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db.set_meta_value('agglomeration_max_mismatch_freq', self.agglom_max_mismatch_freq)
        trnaseq_db.db.set_meta_value('fiveprimemost_deletion_start', self.fiveprimemost_del_start)
        trnaseq_db.db.set_meta_value('threeprimemost_deletion_start', self.threeprimemost_del_start)
        trnaseq_db.db.set_meta_value('fiveprimemost_deletion_stop', self.fiveprimemost_del_stop)
        trnaseq_db.db.set_meta_value('threeprimemost_deletion_stop', self.threeprimemost_del_stop)
        trnaseq_db.db.set_meta_value('max_distinct_deletions', self.max_distinct_dels)
        trnaseq_db.disconnect()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Agglomeration max mismatch frequency", self.agglom_max_mismatch_freq))
            f.write(self.get_summary_line("INDELs profiled", not self.skip_INDEL_profiling))
            f.write(self.get_summary_line("Fiveprime-most deletion start", self.fiveprimemost_del_start))
            f.write(self.get_summary_line("Threeprime-most deletion start", self.threeprimemost_del_start))
            f.write(self.get_summary_line("Fiveprime-most deletion start", self.threeprimemost_del_start))
            f.write(self.get_summary_line("Threeprime-most deletion stop", self.threeprimemost_del_start))
            f.write(self.get_summary_line("Max distinct deletions", self.max_distinct_dels))


    def profile_trna(self, uniq_reads):
        """Profile tRNA features in reads.

        Appends UniqueSeq objects representing profiled tRNA sequences to `self.uniq_trna_seqs`.
        Appends UniqueSeq objects representing sequences with a truncated tRNA profile to
        `self.uniq_trunc_seqs`. Appends leftover UniqueSeq objects representing unprofiled tRNA
        sequences to `self.uniq_nontrna_seqs`.

        Parameters
        ==========
        uniq_reads : list
            List of UniqueSeq objects
        """
        start_time = time.time()
        self.progress.new("Profiling tRNA features in reads")
        self.progress.update("...")

        queued_read_count = 0
        queued_seq_count = 0

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        profiler = trnaidentifier.Profiler()
        processes = [mp.Process(target=profile_worker, args=(input_queue, output_queue, profiler))
                     for _ in range(self.num_threads)]
        for p in processes:
            p.start()

        fetched_profile_count = 0
        uniq_seqs_dict = {}
        profiling_progress_interval = self.profiling_progress_interval
        for uniq_read in uniq_reads:
            input_queue.put((uniq_read.seq_string, uniq_read.represent_name))
            uniq_seqs_dict[uniq_read.represent_name] = uniq_read
            queued_read_count += uniq_read.read_count
            queued_seq_count += 1

            while fetched_profile_count < queued_seq_count:
                trna_profile = output_queue.get()
                fetched_profile_count += 1
                processed_seq_name = trna_profile.name
                processed_seq_string = trna_profile.input_seq

                uniq_seq = uniq_seqs_dict.pop(processed_seq_name)

                if not trna_profile.is_predicted_trna:
                    if trna_profile.trunc_profile_index: # This cannot be 0, but will be None when the profile is not truncated
                        uniq_seq.id_method = 0 # Given choice of profiled (0) and mapped (1), choose profiled
                        uniq_seq.feature_start_indices = [feature.start_pos if hasattr(feature, 'start_pos') else feature.start_positions
                                                          for feature in trna_profile.features]
                        uniq_seq.feature_stop_indices = [feature.stop_pos if hasattr(feature, 'stop_pos') else feature.stop_positions
                                                         for feature in trna_profile.features]
                        uniq_seq.has_complete_feature_set = False # trna_profile.has_complete_feature_set should always be False here
                        uniq_seq.has_his_g = False
                        uniq_seq.alpha_start_index = trna_profile.alpha_start
                        uniq_seq.alpha_stop_index = trna_profile.alpha_stop
                        uniq_seq.beta_start_index = trna_profile.beta_start
                        uniq_seq.beta_stop_index = trna_profile.beta_stop
                        uniq_seq.anticodon_string = anticodon = trna_profile.anticodon_seq
                        uniq_seq.anticodon_aa = trna_profile.anticodon_aa if trna_profile.anticodon_aa else None
                        uniq_seq.contains_anticodon = True if anticodon else False
                        uniq_seq.acceptor_length = len(trna_profile.acceptor_variant_string)
                        uniq_seq.extra_fiveprime_length = 0 # trna_profile.extra_fiveprime_length should always be None here
                        uniq_seq.extra_threeprime_length = trna_profile.num_extra_threeprime
                        uniq_seq.profiled_seq_length = len(trna_profile.profiled_seq)
                        uniq_seq.trunc_profile_index = trna_profile.trunc_profile_index
                        self.uniq_trunc_seqs.append(uniq_seq)
                    else:
                        self.uniq_nontrna_seqs.append(uniq_seq)

                    if fetched_profile_count % profiling_progress_interval == 0:
                        self.progress.update("%s of %s unique sequences have been profiled"
                                             % (pp(fetched_profile_count), pp(len(uniq_reads))))
                    continue

                uniq_seq.id_method = 0
                uniq_seq.feature_start_indices = [feature.start_pos if hasattr(feature, 'start_pos') else feature.start_positions
                                                  for feature in trna_profile.features]
                uniq_seq.feature_stop_indices = [feature.stop_pos if hasattr(feature, 'stop_pos') else feature.stop_positions
                                                 for feature in trna_profile.features]
                uniq_seq.has_complete_feature_set = trna_profile.has_complete_feature_set
                uniq_seq.has_his_g = True if trna_profile.features[0].name == 'tRNA-His position 0' else False
                uniq_seq.alpha_start_index = trna_profile.alpha_start
                uniq_seq.alpha_stop_index = trna_profile.alpha_stop
                uniq_seq.beta_start_index = trna_profile.beta_start
                uniq_seq.beta_stop_index = trna_profile.beta_stop
                uniq_seq.anticodon_string = anticodon = trna_profile.anticodon_seq
                uniq_seq.anticodon_aa = trna_profile.anticodon_aa if trna_profile.anticodon_aa else None
                uniq_seq.contains_anticodon = True if anticodon else False
                uniq_seq.acceptor_length = len(trna_profile.acceptor_variant_string)
                uniq_seq.num_extrapolated_fiveprime_nts = trna_profile.num_in_extrapolated_fiveprime_feature
                uniq_seq.extra_fiveprime_length = trna_profile.num_extra_fiveprime
                uniq_seq.extra_threeprime_length = trna_profile.num_extra_threeprime
                uniq_seq.profiled_seq_length = len(trna_profile.profiled_seq)
                uniq_seq.num_conserved = trna_profile.num_conserved
                uniq_seq.num_unconserved = trna_profile.num_unconserved
                uniq_seq.num_paired = trna_profile.num_paired
                uniq_seq.num_unpaired = trna_profile.num_unpaired
                # Recover nucleotides that did not fit expectation, either by not being the
                # expected nucleotide or type of nucleotide or by not base pairing in a stem.
                uniq_seq.unpaired_info = trna_profile.unpaired_info
                uniq_seq.unconserved_info = trna_profile.unconserved_info
                self.uniq_trna_seqs.append(uniq_seq)

                if fetched_profile_count % profiling_progress_interval == 0:
                    self.progress.update("%s of %s unique sequences have been profiled"
                                         % (pp(fetched_profile_count), pp(len(uniq_reads))))

        for p in processes:
            p.terminate()
            p.join()

        # Profiled seqs were added to the output queue as they were processed, so sort by name.
        self.uniq_trna_seqs.sort(key=lambda uniq_seq: uniq_seq.represent_name)
        self.uniq_trunc_seqs.sort(key=lambda uniq_seq: uniq_seq.represent_name)
        self.uniq_nontrna_seqs.sort(key=lambda uniq_seq: uniq_seq.represent_name)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed profiling tRNA (min)",
                                          time.time() - start_time,
                                          is_time_value=True))
            f.write(self.get_summary_line("Reads processed", queued_read_count))
            f.write(self.get_summary_line("Unique sequences processed", queued_seq_count))

        self.progress.end()

        self.run.info("Reads processed", queued_read_count)
        self.run.info("Unique sequences processed", queued_seq_count)


    def trim_trna_ends(self):
        """Trim any nucleotides 5' of the acceptor stem and 3' of the discriminator from profiled
        tRNA. Appends TrimmedSeq objects formed from input UniqueSeq objects to
        `self.trimmed_trna_seqs`.
        """
        start_time = time.time()
        self.progress.new("Trimming the 3' and 5' ends of tRNA")
        self.progress.update("...")

        represent_names = [uniq_seq.represent_name for uniq_seq in self.uniq_trna_seqs]
        trimmed_seq_strings = [
            uniq_seq.seq_string[uniq_seq.extra_fiveprime_length: len(uniq_seq.seq_string) - uniq_seq.acceptor_length]
            for uniq_seq in self.uniq_trna_seqs
        ]

        clusters = Dereplicator(represent_names,
                                trimmed_seq_strings,
                                extras=self.uniq_trna_seqs,
                                progress=self.progress).full_length_dereplicate()

        trimmed_seqs = [TrimmedSeq(cluster.member_seqs[0], cluster.member_extras) for cluster in clusters]
        self.trimmed_trna_seqs.extend(trimmed_seqs)
        self.trimmed_trna_seqs.sort(key=lambda trimmed_seq: trimmed_seq.represent_name)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed trimming profiled sequences (min)",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def trim_truncated_profile_ends(self):
        """Trim any nucleotides 3' of the discriminator from sequences with a truncated tRNA
        profile. Appends TrimmedSeq objects formed from input UniqueSeq objects to
        `self.trimmed_trunc_seqs`."""
        start_time = time.time()
        self.progress.new("Trimming the 3' ends of sequences with truncated tRNA profiles")
        self.progress.update("...")

        represent_names = [uniq_seq.represent_name for uniq_seq in self.uniq_trunc_seqs]
        trimmed_seq_strings = [uniq_seq.seq_string[: len(uniq_seq.seq_string) - uniq_seq.acceptor_length]
                               for uniq_seq in self.uniq_trunc_seqs]

        clusters = Dereplicator(represent_names,
                                trimmed_seq_strings,
                                extras=self.uniq_trunc_seqs,
                                progress=self.progress).full_length_dereplicate()

        trimmed_seqs = [TrimmedSeq(cluster.member_seqs[0], cluster.member_extras) for cluster in clusters]
        self.trimmed_trunc_seqs.extend(trimmed_seqs)
        self.trimmed_trunc_seqs.sort(key=lambda trimmed_seq: trimmed_seq.represent_name)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed trimming sequences with truncated feature profile (min)",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def threeprime_dereplicate_trna(self):
        """Dereplicate trimmed profiled tRNA sequences from the 3' end of longer trimmed sequences.

        EXAMPLE:
        Normalized tRNA (trimmed tRNA 1): TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Trimmed tRNA 2                  :                       AATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Trimmed tRNA 3                  :                                     GCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        """
        start_time = time.time()
        self.progress.new("Dereplicating trimmed tRNA sequences from the 3' end")
        self.progress.update("...")

        represent_names = [trimmed_seq.represent_name for trimmed_seq in self.trimmed_trna_seqs]
        # Reverse sequence orientation to dereplicate from the 3' end.
        reversed_seq_strings = [trimmed_seq.seq_string[::-1] for trimmed_seq in self.trimmed_trna_seqs]
        clusters = Dereplicator(represent_names,
                                reversed_seq_strings,
                                extras=self.trimmed_trna_seqs,
                                progress=self.progress).prefix_dereplicate()

        # Profiling may have found multiple sequences that would here be 3'-dereplicated as having
        # complete, but different, feature profiles. This can be caused by the "accommodation" of
        # extra 5' nucleotides (due to a reverse transcriptase artifact or pre-tRNA leader sequence)
        # in an erroneous profile with an elongated variable loop. Retain the shortest "completely
        # profiled" sequence in the cluster, discarding any longer sequences from the list of
        # trimmed tRNA sequences, and their constituent profiled sequences from the list of unique
        # tRNA sequences.

        # Similarly, the longest sequence in the cluster may have an erroneous "incomplete profile."
        # If there is a shorter sequence in the cluster with a complete profile, then any longer
        # sequences with an incomplete profile can be discarded.

        # We do not check for feature-by-feature agreement among clustered profiles here, as 3' tRNA
        # fragments often have some incorrect feature positions due to the paucity of sequence
        # information available in profiling. It is conceivable that the seed sequence of the
        # cluster is wrongly completely profiled, but all of the shorter sequences in the cluster
        # are correctly incompletely profiled. It is also possible that the seed sequence is
        # incompletely profiled, but has the wrong profile, and the shorter sequences are correctly
        # incompletely profiled.

        # This step, which likely removes a tiny number of trimmed sequences, helps reduce the
        # number of extra, wrong nucleotides at the 5' end of seed sequences.
        self.progress.update("Inspecting normalized sequence clusters")
        norm_trna_seqs = self.norm_trna_seqs
        trimmed_trna_seq_represent_names = [trimmed_seq.represent_name for trimmed_seq in self.trimmed_trna_seqs]
        uniq_trna_seq_represent_names = [uniq_seq.represent_name for uniq_seq in self.uniq_trna_seqs]
        trimmed_seq_indices_to_remove = []
        uniq_seq_indices_to_remove = []
        for cluster in clusters:
            # Skip initialization of NormalizedSeq objects, as additional TrimmedSeq members are
            # later added to the objects after dereplicating sequences with truncated tRNA profiles
            # and mapping unprofiled tRNA fragments.
            if len(cluster.member_extras) == 1:
                norm_trna_seqs.append(NormalizedSeq(cluster.member_extras, skip_init=True))
                continue

            # Check that there are no shorter sequences in the cluster with a "complete profile".
            complete_profile_indices = []
            for trimmed_seq_index, trimmed_seq in enumerate(cluster.member_extras):
                if trimmed_seq.has_complete_feature_set:
                    complete_profile_indices.append(trimmed_seq_index)

            if not complete_profile_indices or complete_profile_indices == [0]:
                norm_trna_seqs.append(NormalizedSeq(cluster.member_extras, skip_init=True))
                continue

            # Reaching this point means that there are multiple sequences with "complete
            # profiles" in the cluster.

            # If the two shortest sequences with complete feature profiles differ by the
            # post-transcriptionally added 5'-G of tRNA-His, then they should both be allowed.
            if cluster.member_extras[complete_profile_indices[-2]].has_his_g:
                if cluster.member_extras[complete_profile_indices[-1]].seq_string == cluster.member_extras[complete_profile_indices[-2]].seq_string[1:]:
                    norm_trna_seqs.append(NormalizedSeq(cluster.member_extras[complete_profile_indices[-2]:], skip_init=True))
                    continue

            norm_trna_seqs.append(NormalizedSeq(cluster.member_extras[complete_profile_indices[-1]:], skip_init=True))
            for trimmed_seq in cluster.member_extras[:complete_profile_indices[-1]]:
                trimmed_seq_indices_to_remove.append(trimmed_trna_seq_represent_names.index(trimmed_seq.represent_name))
                for uniq_seq in trimmed_seq.uniq_seqs:
                    uniq_seq_indices_to_remove.append(uniq_trna_seq_represent_names.index(uniq_seq.represent_name))

        # Trimmed and unique tRNA sequences should already be sorted by representative name.
        trimmed_seq_indices_to_remove.sort(reverse=True)
        for trimmed_seq_index in trimmed_seq_indices_to_remove:
            self.trimmed_trna_seqs.pop(trimmed_seq_index)
        uniq_seq_indices_to_remove.sort(reverse=True)
        for uniq_seq_index in uniq_seq_indices_to_remove:
            self.uniq_trna_seqs.pop(uniq_seq_index)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed 3'-dereplicating trimmed profiled sequences",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def write_feature_table(self):
        self.progress.new("Writing tRNA-seq database table of profiled tRNA features")
        self.progress.update("...")

        feature_table_entries = []
        for uniq_seq in self.uniq_trna_seqs:
            feature_table_entries.append(
                (uniq_seq.represent_name,
                 uniq_seq.has_complete_feature_set,
                 uniq_seq.anticodon_string,
                 uniq_seq.anticodon_aa,
                 len(uniq_seq.seq_string),
                 # Zero-based start position of identified tRNA features within the read.
                 len(uniq_seq.seq_string) - uniq_seq.profiled_seq_length,
                 uniq_seq.num_conserved,
                 uniq_seq.num_unconserved,
                 uniq_seq.num_paired,
                 uniq_seq.num_unpaired,
                 uniq_seq.num_extrapolated_fiveprime_nts,
                 uniq_seq.extra_fiveprime_length,
                 uniq_seq.extra_threeprime_length)
                # When tRNA features are not found at the 5' end of the read,
                # the start and stop positions of these features also are not found.
                + tuple([None for _ in range((len(self.TRNA_FEATURE_NAMES) - len(uniq_seq.feature_start_indices)))]) * 2
                + tuple(itertools.chain(*zip(
                    [str(start) if isinstance(start, int) else ','.join(map(str, start))
                     for start in uniq_seq.feature_start_indices],
                    # Convert Pythonic stop position to real stop position of feature.
                    [str(stop - 1) if isinstance(stop, int) else ','.join(str(i - 1) for i in stop)
                     for stop in uniq_seq.feature_stop_indices])))
                # The alpha and beta sections of the D loop are "subfeatures," not "features,"
                # so add them as columns to the table after the features.
                + (uniq_seq.alpha_start_index,
                   uniq_seq.alpha_stop_index - 1 if uniq_seq.alpha_stop_index else None,
                   uniq_seq.beta_start_index,
                   uniq_seq.beta_stop_index - 1 if uniq_seq.beta_stop_index else None)
            )

        if feature_table_entries:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('feature', ','.join('?' * len(tables.trnaseq_feature_table_structure))),
                                     feature_table_entries)
            trnaseq_db.disconnect()

        self.progress.end()


    def write_unconserved_table(self):
        self.progress.new("Writing tRNA-seq database table of unconserved nucleotides in profiled tRNA")
        self.progress.update("...")

        unconserved_table_entries = []
        for uniq_seq in self.uniq_trna_seqs:
            for unconserved_tuple in uniq_seq.unconserved_info:
                unconserved_table_entries.append((uniq_seq.represent_name, ) + unconserved_tuple)

        if unconserved_table_entries:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('unconserved', ','.join('?' * len(tables.trnaseq_unconserved_table_structure))),
                                     unconserved_table_entries)
            trnaseq_db.disconnect()

        self.progress.end()


    def write_unpaired_table(self):
        self.progress.new("Writing tRNA-seq database table of unpaired nucleotides in profiled tRNA")
        self.progress.update("...")

        unpaired_table_entries = []
        for uniq_seq in self.uniq_trna_seqs:
            for unpaired_tuple in uniq_seq.unpaired_info:
                unpaired_table_entries.append((uniq_seq.represent_name, ) + unpaired_tuple)

        if unpaired_table_entries:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('unpaired', ','.join('?' * len(tables.trnaseq_unpaired_table_structure))),
                                     unpaired_table_entries)
            trnaseq_db.disconnect()

        self.progress.end()


    def threeprime_dereplicate_truncated_sequences(self):
        """Try to recover sequences with truncated tRNA profiles that are 3' subsequences of
        profiled tRNA and thus legitimate 3' tRNA fragments. These trimmed sequences are folded into
        the normalized tRNA sequences in self.norm_trna_seqs. Unrecovered trimmed sequences are
        themselves 3'-dereplicated, forming another pool of normalized sequences,
        self.norm_trunc_seqs."""
        start_time = time.time()
        self.progress.new("Dereplicating trimmed sequences with a truncated feature profile")
        self.progress.update("...")

        trimmed_trunc_seq_represent_names = [trimmed_trunc_seq.represent_name for trimmed_trunc_seq in self.trimmed_trunc_seqs]
        norm_trna_seq_represent_names = [norm_trna_seq.represent_name for norm_trna_seq in self.norm_trna_seqs]
        # Reverse sequence orientation to dereplicate from the 3' end.
        reversed_trimmed_trunc_seq_strings = [trimmed_trunc_seq.seq_string[::-1] for trimmed_trunc_seq in self.trimmed_trunc_seqs]
        reversed_norm_trna_seq_strings = [norm_trna_seq.seq_string[::-1] for norm_trna_seq in self.norm_trna_seqs]
        clusters = Dereplicator(trimmed_trunc_seq_represent_names + norm_trna_seq_represent_names,
                                reversed_trimmed_trunc_seq_strings + reversed_norm_trna_seq_strings,
                                extras=self.trimmed_trunc_seqs + self.norm_trna_seqs,
                                progress=self.progress).prefix_dereplicate()

        # Associate each truncated sequence with any normalized tRNA sequences that contain it as a
        # 3'-subsequence. Do not bother with the complex task of reconstructing a feature profile of
        # the truncated sequence (entries for each unique truncated sequence would appear in the
        # Features table of the tRNA-seq database). An attempt to do this revealed that it is a very
        # deep rabbit hole.

        # Clusters cannot contain more than one normalized tRNA sequence, as they have already been
        # 3'-dereplicated. Normalized sequences can seed clusters and also be members of clusters
        # seed by truncated sequences. In the latter case of clusters containing a normalized tRNA
        # sequence as a member but not as the seed, truncated sequences in the cluster that are
        # shorter than the normalized sequence (3' subsequences) are incorporated as members of the
        # normalized sequence.

        # Consider three types of cluster: 1. clusters consisting of a single normalized sequence
        # (ignore), 2. clusters containing a normalized sequence as seed or member with shorter
        # truncated sequence members, and 3. clusters seeded by truncated sequences. If a truncated
        # sequence is found in group 2 (part of one or more longer normalized sequences) then ignore
        # it in group 3 (do not include it in normalized truncated sequences formed from group 3
        # clusters).
        norm_trna_seq_match_dict = {}
        trimmed_trunc_seq_match_dict = {}
        for cluster in clusters:
            if len(cluster.member_names) == 1:
                if isinstance(cluster.member_extras[0], NormalizedSeq):
                    continue

            norm_trna_seq = None
            trimmed_trunc_seq_seed = cluster.member_extras[0] if isinstance(cluster.member_extras[0], TrimmedSeq) else None
            for seq in cluster.member_extras:
                # Members of each cluster are pre-sorted in descending order of sequence length.
                # There cannot be a normalized tRNA sequence and a truncated sequence of the same
                # length (they would be the same sequence).
                if isinstance(seq, NormalizedSeq):
                    norm_trna_seq = seq
                    continue

                if norm_trna_seq:
                    if isinstance(seq, TrimmedSeq):
                        try:
                            norm_trna_seq_match_dict[seq.represent_name][1].append(norm_trna_seq)
                        except KeyError:
                            norm_trna_seq_match_dict[seq.represent_name] = (seq, [norm_trna_seq])
                        continue
                    else:
                        raise ConfigError("Each cluster should only contain zero or one normalized sequence, "
                                          "but it appears that this cluster contains more than one.")

                try:
                    trimmed_trunc_seq_match_dict[trimmed_trunc_seq_seed.represent_name][1].append(seq)
                except KeyError:
                    trimmed_trunc_seq_match_dict[trimmed_trunc_seq_seed.represent_name] = (trimmed_trunc_seq_seed, [seq])

        # Add truncated sequences to matching normalized tRNA sequences.

        # To determine the count of truncated sequences with an anticodon, the location of the
        # anticodon in a matching normalized sequence must first be found.
        relative_anticodon_loop_index = self.TRNA_FEATURE_NAMES.index('anticodon_loop') - len(self.TRNA_FEATURE_NAMES) + 1
        norm_trna_seq_anticodon_dict = {}

        trimmed_trunc_seq_represent_names = [trimmed_trunc_seq.represent_name for trimmed_trunc_seq in self.trimmed_trunc_seqs]
        uniq_trunc_seq_represent_names = [uniq_trunc_seq.represent_name for uniq_trunc_seq in self.uniq_trunc_seqs]
        trimmed_trunc_seq_indices_to_remove = []
        uniq_trunc_seq_indices_to_remove = []
        for trimmed_trunc_seq_represent_name, entry in norm_trna_seq_match_dict.items():
            trimmed_trunc_seq, norm_trna_seqs = entry
            trimmed_trunc_seq_length = len(trimmed_trunc_seq.seq_string)

            trimmed_trunc_seq_indices_to_remove.append(trimmed_trunc_seq_represent_names.index(trimmed_trunc_seq_represent_name))
            for uniq_trunc_seq in trimmed_trunc_seq.uniq_seqs:
                uniq_trunc_seq_indices_to_remove.append(uniq_trunc_seq_represent_names.index(uniq_trunc_seq.represent_name))

            for norm_seq in norm_trna_seqs:
                norm_seq.trimmed_seqs.append(trimmed_trunc_seq)
                norm_seq_length = len(norm_seq.seq_string)
                norm_seq.start_positions.append(norm_seq_length - trimmed_trunc_seq_length)
                norm_seq.stop_positions.append(norm_seq_length)
                trimmed_trunc_seq.norm_seq_count += 1

                # Determine from the first normalized sequence in which the truncated sequence is
                # found whether the truncated sequence contains the anticodon.
                if not trimmed_trunc_seq.contains_anticodon:
                    try:
                        anticodon_start_relative_to_acceptor = norm_trna_seq_anticodon_dict[norm_seq.represent_name]
                    except KeyError:
                        anticodon_loop_start = norm_seq.trimmed_seqs[0].feature_start_indices[relative_anticodon_loop_index]
                        if anticodon_loop_start > -1:
                            anticodon_start = anticodon_loop_start + 2
                            anticodon_start_relative_to_acceptor = anticodon_start - norm_seq.trimmed_seqs[0].uniq_seqs[0].feature_start_indices[-1]
                        else:
                            anticodon_start_relative_to_acceptor = 1 # A positive number is used to mean the anticodon was not profiled
                    if trimmed_trunc_seq_length + anticodon_start_relative_to_acceptor >= 0:
                        trimmed_trunc_seq.contains_anticodon = True
                        for uniq_seq in trimmed_trunc_seq.uniq_seqs:
                            uniq_seq.contains_anticodon = True

        # Consolidate trimmed truncated sequences that don't match normalized tRNA sequences into
        # normalized truncated sequences.
        for trimmed_trunc_seq_target_represent_name, entry in trimmed_trunc_seq_match_dict.items():
            trimmed_trunc_seq_target, trimmed_trunc_seq_queries = entry
            if trimmed_trunc_seq_target_represent_name in norm_trna_seq_match_dict:
                # All shorter trimmed truncated sequences in the cluster will also have matched a
                # normalized tRNA sequence.
                continue

            for i, trimmed_trunc_seq_query in enumerate(trimmed_trunc_seq_queries):
                if trimmed_trunc_seq_query.represent_name in norm_trna_seq_match_dict:
                    trimmed_trunc_seq_queries = trimmed_trunc_seq_queries[: i]
                    break

            self.norm_trunc_seqs.append(NormalizedSeq([trimmed_trunc_seq_target] + trimmed_trunc_seq_queries, skip_init=True))

        # Transfer recovered truncated sequences into the tRNA lists.
        trimmed_trunc_seq_indices_to_remove.sort(reverse=True)
        for trimmed_trunc_seq_index in trimmed_trunc_seq_indices_to_remove:
            trimmed_trunc_seq = self.trimmed_trunc_seqs.pop(trimmed_trunc_seq_index)
            trimmed_trunc_seq.trunc_profile_recovered_by_derep = True
            self.trimmed_trna_seqs.append(trimmed_trunc_seq)
        uniq_trunc_seq_indices_to_remove.sort(reverse=True)
        for uniq_trunc_seq_index in uniq_trunc_seq_indices_to_remove:
            uniq_trunc_seq = self.uniq_trunc_seqs.pop(uniq_trunc_seq_index)
            uniq_trunc_seq.trunc_profile_recovered_by_derep = True
            self.uniq_trna_seqs.append(uniq_trunc_seq)

        # All trimmed sequences have now been added to normalized truncated sequences, so they can
        # be initialized.
        for norm_trunc_seq in self.norm_trunc_seqs:
            norm_trunc_seq.init()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed recovering tRNA with truncated feature profile (min)",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def write_checkpoint_files(self, checkpoint_name):
        if not self.write_checkpoints:
            return

        self.progress.new(f"Writing intermediate files for the \"{checkpoint_name}\" checkpoint")
        self.progress.update("...")

        overwrote_dict = {}
        for intermed_file_key, intermed_file_path in self.intermed_file_path_dict[checkpoint_name].items():
            if os.path.exists(intermed_file_path):
                overwrote_dict[self.intermed_file_label_dict[intermed_file_key]] = intermed_file_path
            else:
                overwrote_dict[self.intermed_file_label_dict[intermed_file_key]] = None
            with open(intermed_file_path, 'wb') as f:
                # The key, e.g., "uniq_trna_seqs", corresponds to the attribute to be saved to file.
                pkl.dump(getattr(self, intermed_file_key), f, protocol=pkl.HIGHEST_PROTOCOL)

        self.progress.end()

        for intermed_file_label, intermed_file_path in overwrote_dict.items():
            # Example: "Overwrote profile checkpoint intermediate file of unique tRNA"
            self.run.info("%s\"checkpoint_name\" checkpoint intermediate file of %s"
                          % ("Overwrote " if intermed_file_path else "", intermed_file_label),
                          intermed_file_path)


    def load_checkpoint_files(self, checkpoint_name):
        self.progress.new(f"Loading intermediate files at the checkpoint, \"{checkpoint_name}\"")
        self.progress.update("...")

        with open(self.analysis_summary_path, 'a') as f:
            f.write(f"\nAnalysis restarted from the checkpoint, \"{checkpoint_name}\"\n")

        for intermed_file_key, intermed_file_path in self.intermed_file_path_dict[checkpoint_name].items():
            with open(intermed_file_path, 'rb') as f:
                setattr(self, intermed_file_key, pkl.load(f))

        self.progress.end()


    def threeprime_dereplicate_acceptorless_sequences(self):
        """Find tRNA sequences missing a 3'-acceptor sequence variant that are 3' subsequences of
        normalized tRNA sequences. Sequences required an acceptor sequence variant to be profiled as
        tRNA. tRNA sequences salvaged here are therefore counted as `mapped` rather than `profiled`
        in the `id_method` attribute of their `UniqueSeq` and `TrimmedSeq` objects (these sequences
        form independent `TrimmedSeq` objects rather than being added to existing ones with an
        identical sequence)."""
        start_time = time.time()
        self.progress.new("Dereplicating acceptorless tRNA sequences")
        self.progress.update("...")

        represent_names = []
        reversed_seq_strings = []
        extras = []
        for uniq_nontrna_index, uniq_nontrna_seq in enumerate(self.uniq_nontrna_seqs):
            if len(uniq_nontrna_seq.seq_string) >= self.min_trna_frag_size:
                represent_names.append(uniq_nontrna_seq.represent_name)
                reversed_seq_strings.append(uniq_nontrna_seq.seq_string[::-1])
                extras.append((uniq_nontrna_index, uniq_nontrna_seq))
        for norm_seq in self.norm_trna_seqs:
            represent_names.append(norm_seq.represent_name)
            reversed_seq_strings.append(norm_seq.seq_string[::-1])
            extras.append((-1, norm_seq))
        clusters = Dereplicator(represent_names, reversed_seq_strings, extras=extras, progress=self.progress).prefix_dereplicate()

        uniq_nontrna_seq_indices_to_remove = []
        uniq_seq_norm_seqs_dict = defaultdict(list)
        for cluster in clusters:
            if len(cluster.member_seqs) == 1:
                continue

            # Check that there is a normalized tRNA sequence in the cluster -- there cannot be more
            # than one.
            cluster_norm_seq_index = None
            norm_seq = None
            norm_seq_length = None
            norm_seq_has_complete_feature_set = None
            for member_index, member_extra in enumerate(cluster.member_extras):
                if member_extra[0] == -1:
                    cluster_norm_seq_index = member_index
                    norm_seq = member_extra[1]
                    norm_seq_length = len(norm_seq.seq_string)
                    if member_index > 0:
                        norm_seq_has_complete_feature_set = norm_seq.has_complete_feature_set
                    break
            else:
                continue

            for uniq_nontrna_index, uniq_seq in cluster.member_extras[:cluster_norm_seq_index]:
                # To be longer than a normalized sequence, the sequence must only be found in this
                # one normalized-sequence-containing cluster, as to be found in multiple of these
                # clusters would mean that the normalized sequences in those clusters must have
                # 3'-subsequence relationships, which is impossible by the very existence of the
                # normalized sequences. If the normalized sequence has a complete feature profile,
                # than the overhanging 5' bases in the salvaged sequence can be trimmed as "extra"
                # 5' bases. Otherwise, it is unclear if the overhanging 5' bases are somehow part of
                # an artifact, so take the conservative option of ignoring the sequence.
                if not norm_seq_has_complete_feature_set:
                    break

                uniq_nontrna_seq_indices_to_remove.append(uniq_nontrna_index)
                uniq_seq.id_method = 1 # mapped
                uniq_seq.has_complete_feature_set = False
                uniq_seq.acceptor_length = 0
                uniq_seq.extra_fiveprime_length = len(uniq_seq.seq_string) - norm_seq_length

                trimmed_seq = TrimmedSeq(uniq_seq.seq_string[uniq_seq.extra_fiveprime_length: ], [uniq_seq])
                self.trimmed_trna_seqs.append(trimmed_seq)

                norm_seq.trimmed_seqs.append(trimmed_seq)
                trimmed_seq.norm_seq_count += 1
                norm_seq.start_positions.append(0)
                norm_seq.stop_positions.append(norm_seq_length)

            for uniq_nontrna_index, _ in cluster.member_extras[cluster_norm_seq_index + 1:]:
                uniq_nontrna_seq_indices_to_remove.append(uniq_nontrna_index)
                uniq_seq_norm_seqs_dict[uniq_nontrna_index].append(norm_seq)

        for seq_index in sorted(set(uniq_nontrna_seq_indices_to_remove), reverse=True):
            uniq_seq = self.uniq_nontrna_seqs.pop(seq_index)
            if uniq_seq.id_method == 1:
                # This sequence was already handle above, save removal from the list of non-tRNA.
                continue

            self.uniq_trna_seqs.append(uniq_seq)
            uniq_seq.id_method = 1
            uniq_seq.has_complete_feature_set = False
            uniq_seq.acceptor_length = 0
            uniq_seq.extra_fiveprime_length = 0

            trimmed_seq = TrimmedSeq(uniq_seq.seq_string, [uniq_seq])
            self.trimmed_trna_seqs.append(trimmed_seq)

            uniq_seq_length = len(uniq_seq.seq_string)
            for norm_seq in uniq_seq_norm_seqs_dict[seq_index]:
                norm_seq.trimmed_seqs.append(trimmed_seq)
                trimmed_seq.norm_seq_count += 1
                norm_seq_length = len(norm_seq.seq_string)
                norm_seq.start_positions.append(norm_seq_length - uniq_seq_length)
                norm_seq.stop_positions.append(norm_seq_length)

        self.progress.end()


    def map_fragments(self):
        """Map unprofiled tRNA fragments to longer profiled tRNA sequences.

        If the specified minimum fragment length is shorter than the minimum length of profiled tRNA
        fragments -- which does not happen with the default settings, as both lengths are 25, but
        can happen when the user adjusts `--min-trna-fragment-size` downward -- then mapped
        fragments may occur at the 3' end of a tRNA, but will not include 3' acceptor variants (CCA,
        CC, C, etc.), as these 3' extensions were trimmed off the mapping targets.

        EXAMPLE:
        Normalized tRNA:                 (GT)TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Mapped tRNA 1 (extra 5' bases) :   T TCCGTGATAGTTTAATGGTCAGAATGG
        Mapped tRNA 2 (interior)  :                 TAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGG
        """
        start_time = time.time()
        self.progress.new("Mapping unprofiled reads to profiled tRNA")

        self.progress.update("Retrieving queries from unprofiled reads")
        query_length_intervals = []
        max_query_length = max(map(len, [seq.seq_string for seq in self.uniq_nontrna_seqs]))
        # By default, avoid sequences shorter than 25 nucleotides, the default minimum length of a
        # profiled 3' fragment of tRNA.
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
        # Leftover non-tRNA sequences are mapped to normalized tRNA sequences with extra 5' bases
        # added when present in underlying unique tRNA sequences. Multiple targets for each
        # normalized sequence are therefore produced for different 5' sequence extensions.
        target_names = []
        target_seqs = []
        for norm_seq_index, norm_seq in enumerate(self.norm_trna_seqs):
            norm_seq_string = norm_seq.seq_string
            # The longest trimmed sequence (the first in the list) is by design the only one of the
            # profiled trimmed sequences forming the normalized sequence that may have extra 5'
            # bases.
            longest_trimmed_seq = norm_seq.trimmed_seqs[0]
            if longest_trimmed_seq.uniq_with_extra_fiveprime_count > 0:
                fiveprime_seq_string_set = set()
                for uniq_seq in longest_trimmed_seq.uniq_seqs:
                    if uniq_seq.extra_fiveprime_length > 0:
                        fiveprime_seq_string_set.add(uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length])

                # Avoid creating superfluous target sequences that are subsequences of other target
                # sequences due to a 5' extension of a normalized sequence being a subsequence of a
                # longer 5' extension of the same normalized sequence.
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
                    # Use an index to distinguish otherwise equivalent targets with different 5'
                    # extensions of the same length.
                    target_names.append((norm_seq_index, len(fiveprime_seq_string), fiveprime_index))
                    target_seqs.append(fiveprime_seq_string + norm_seq_string)
            else:
                target_names.append((norm_seq_index, 0, 0)) # No extra 5' bases
                target_seqs.append(norm_seq_string)

        self.progress.end()


        interval_index = 0
        nontrna_indices = []
        for query_names, query_seqs in zip(query_name_chunks, query_seq_chunks):
            self.progress.new("Mapping %s unprofiled reads of length %d-%d to profiled tRNA"
                              % (pp(len(query_names)),
                                 query_length_intervals[interval_index][0],
                                 query_length_intervals[interval_index][1] - 1))

            aligned_query_dict, aligned_target_dict = Aligner( # aligned_target_dict is not used for anything and can be big
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
                    ref_alignment_stop = alignment.target_start + alignment.alignment_length

                    norm_seq_index, ref_fiveprime_length, _ = alignment.aligned_target.name # extra 5' index doesn't matter now

                    norm_stop_pos = ref_alignment_stop - ref_fiveprime_length
                    if norm_stop_pos <= 0:
                        # Ignore queries that align entirely to extra 5' bases. Sequences mapping
                        # exclusively to the 5' extension that are long enough to fulfill the
                        # minimum length requirement are often mapping to an artifactual chimeric
                        # sequence.
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

                        # There can be multiple trimmed sequence objects representing the same
                        # sequence that remains after the 5' extension has been trimmed. The 5'
                        # extension may represent all but a small number of nucleotides in the
                        # sequence, so it is best not to dereplicate mapped trimmed sequences
                        # identical in the non-5' section, as is done for trimmed sequences formed
                        # from profiled tRNA that are grouped into the same normalized sequence.
                        trimmed_seq = TrimmedSeq(uniq_mapped_seq.seq_string[uniq_mapped_seq.extra_fiveprime_length: ],
                                                 [uniq_mapped_seq])
                        self.trimmed_trna_seqs.append(trimmed_seq)

                        norm_seq.trimmed_seqs.append(trimmed_seq)
                        trimmed_seq.norm_seq_count += 1
                        norm_seq.start_positions.append(norm_start_pos)
                        norm_seq.stop_positions.append(norm_stop_pos)
                    else:
                        for prev_trimmed_seq in norm_seq.trimmed_seqs[::-1]:
                            # Ensure that the trimmed sequence maps to the normalized sequence only
                            # once. Multiple targets can be created from the same normalized
                            # sequence for different 5' extensions.
                            if prev_trimmed_seq.id_method == 0:
                                norm_seq.trimmed_seqs.append(trimmed_seq)
                                trimmed_seq.norm_seq_count += 1
                                if ref_fiveprime_length - ref_alignment_start > 0:
                                    norm_start_pos = 0
                                else:
                                    norm_start_pos = ref_alignment_start - ref_fiveprime_length
                                norm_seq.start_positions.append(norm_start_pos)
                                norm_seq.stop_positions.append(norm_stop_pos)
                                break
                            if trimmed_seq.represent_name == prev_trimmed_seq.represent_name:
                                break
            interval_index += 1

            del aligned_query_dict
            gc.collect()

            self.progress.end()

        for norm_seq in self.norm_trna_seqs:
            norm_seq.init()

        for nontrna_index in sorted(nontrna_indices, reverse=True):
            self.uniq_trna_seqs.append(self.uniq_nontrna_seqs.pop(nontrna_index))

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed mapping tRNA fragments (min)",
                                          time.time() - start_time,
                                          is_time_value=True))


    def find_substitutions(self):
        """Find potential modification-induced substitutions."""
        start_time = time.time()
        self.progress.new("Finding modification-induced substitutions")

        # Cluster normalized tRNA sequences. Clusters agglomerate sequences that differ from at
        # least one other sequence in the cluster by no more than 2 substitutions per 71 aligned
        # positions (by default) in a gapless end-to-end alignment.
        agglomerator = Agglomerator([seq.represent_name for seq in self.norm_trna_seqs],
                                    [seq.seq_string for seq in self.norm_trna_seqs],
                                    num_threads=self.num_threads,
                                    progress=self.progress)
        # Provide a priority function for seeding clusters that favors fully profiled tRNA over
        # "longer" tRNA without a full set of profiled features. Such incompletely profiled longer
        # tRNA includes tRNA-tRNA chimeras -- some of these have a long 5' section that is a long 3'
        # fragment of tRNA, which can cause other shorter normalized sequences to agglomerate by
        # aligning to the 5' section of the chimera.
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
            # A modification requires at least 3 different nucleotides to be detected, and each
            # normalized sequence differs by at least 1 nucleotide (substitution or gap), so for a
            # cluster to form a modified sequence, it must contain at least 3 normalized sequences.
            if len(aligned_ref.alignments) < 2:
                continue

            aligned_ref_length = len(aligned_ref.seq_string)

            valid_aligned_queries = []
            for alignment in aligned_ref.alignments:
                # Normalized tRNA sequences should only align at the 3' end. Alignments to the
                # interior of the sequence can occur when the reference is a tRNA-tRNA chimera.
                if aligned_ref_length != alignment.target_start + alignment.alignment_length:
                    continue

                query_name = alignment.aligned_query.name
                # The normalized sequence query may have agglomerated with another reference as
                # well. If the query formed a modified sequence, it would form the same modified
                # sequence when starting with this agglomeration.
                if query_name in names_of_norm_seqs_assigned_to_mod_seqs:
                    continue

                valid_aligned_queries.append(norm_seq_dict[query_name])

            # Confirm that 2 or more queries passed the filters, so at least 3 normalized sequences
            # are still in the cluster.
            if len(valid_aligned_queries) < 2:
                continue

            seq_array = np.zeros((len(valid_aligned_queries) + 1, aligned_ref_length), dtype=int)
            # Rather than using the ASCII representation of each character, which saves some time in
            # converting the sequence string to a numpy array, constrain the integer representation
            # to the smallest possible range of integers to speed up the bincount method used to
            # determine the number of unique nucleotides at an alignment position.
            seq_array[0, :] += [NT_INT_DICT[nt] for nt in aligned_ref.seq_string]
            for i, aligned_query in enumerate(valid_aligned_queries, start=1):
                seq_array[i, aligned_ref_length - len(aligned_query.seq_string): ] += [NT_INT_DICT[nt]
                                                                                       for nt in aligned_query.seq_string]

            norm_seqs = np.array([norm_seq_dict[ref_name]] + valid_aligned_queries)

            # Find positions in the alignment with nucleotide variability.
            alignment_pos_uniq_nt_counts = (
                np.bincount(
                    (seq_array + np.arange(aligned_ref_length, dtype=int) * NUM_NT_BINS).ravel(),
                    minlength=aligned_ref_length * NUM_NT_BINS
                ).reshape(-1, NUM_NT_BINS)[:, 1:] != 0
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

                    # A modification requires at least 3 different nucleotides to be detected, and
                    # each normalized sequence differs by at least 1 nucleotide (substitution or
                    # gap), so for a cluster to form a modified sequence, it must contain at least 3
                    # normalized sequences.
                    if norm_seqs.size < 3:
                        continue

                    aligned_nts = seq_array[:, alignment_pos]
                    nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]

                    if (nt_counts != 0).sum() < 2:
                        # There are now < 2 nucleotides at the alignment position in the (derived)
                        # cluster under consideration. 2 different nucleotides are needed to
                        # distinguish single nucleotide variants.
                        next_clusters.appendleft((seq_array, norm_seqs, three_four_nt_alignment_positions))
                        continue

                    # Add a new cluster for each nucleotide variant to the stack of clusters to
                    # process if the new cluster contains at least 3 sequences.
                    represented_nts = nt_counts.nonzero()[0] + 1
                    for nt in represented_nts:
                        split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]
                        if split_cluster_seq_indices.size > 2:
                            next_clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                      norm_seqs[split_cluster_seq_indices],
                                                      three_four_nt_alignment_positions.copy()))
                if next_clusters:
                    clusters = next_clusters
                else:
                    break
            if not clusters:
                continue

            # Check alignment positions previously found to have 3-4 nucleotides. Further split
            # (derived) clusters when positions now have 2 nucleotides.
            next_clusters = deque()
            while clusters:
                seq_array, norm_seqs, three_four_nt_alignment_positions = clusters.pop()
                candidates_to_remove = []

                for i, alignment_pos in enumerate(three_four_nt_alignment_positions):
                    aligned_nts = seq_array[:, alignment_pos]
                    nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]
                    # At least 3 different nucleotides are needed at a position to predict a
                    # modification.
                    represented_nts = nt_counts.nonzero()[0] + 1
                    if represented_nts.size < 2:
                        candidates_to_remove.append(i)
                    elif represented_nts.size == 2:
                        candidates_to_remove.append(i)
                        for nt in represented_nts:
                            split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]
                            # At least 3 normalized sequences are needed to form a modified
                            # sequence.
                            if split_cluster_seq_indices.size > 2:
                                clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                     norm_seqs[split_cluster_seq_indices],
                                                     np.delete(three_four_nt_alignment_positions, candidates_to_remove)))
                        # Reevaluate previous alignment positions in the split clusters.
                        break
                else:
                    # At least 1 position was discounted as no longer having 3-4 different
                    # nucleotides, but these positions had fewer than 2 nucleotides, and so did not
                    # cause the cluster to be split into new clusters. Therefore, do not cycle
                    # through the remaining positions again to find those with fewer than 3
                    # nucleotides.
                    if candidates_to_remove:
                        next_clusters.appendleft((norm_seqs, np.delete(three_four_nt_alignment_positions, candidates_to_remove)))
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

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed finding modification-induced substitutions (min)",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def find_deletions(self):
        """Find potential modification-induced deletions. First, in silico "test" deletions are
        introduced at and around substitution sites in sequences with potential modification-induced
        substitutions. These query sequences are searched against two pools of sequence targets: 1.
        normalized sequences not assigned to modified sequences and 2. normalized "non-tRNA"
        sequences with truncated tRNA profiles. Why are these two specific pools used? 1. Why are
        normalized tRNA sequences considered at all, when they have been successfully profiled, and
        therefore, presumably, do not have deletions interrupting the profile? Deletions can be
        erroneously accommodated by flexibility in feature lengths. For example, a deletion
        associated with a modification in the D loop can cause the variable-length alpha or beta
        sections of the D loop to be assigned one fewer nucleotide than is correct to optimize the
        profile. 2. Why not all "non-tRNAs," why just those with a truncated feature profile?
        Deletions can cause the truncation of the profile. "Non-tRNAs" without even a truncated
        profile (sequences that were not profiled past the minimum length threshold of the T arm)
        also have fewer opportunities for modification-induced mutations.

        Normalized sequences rather than trimmed or unique sequences are searched for the sake of
        speed. Ideally, normalized sequences found to have deletions would be further analyzed,
        finding which constituent nonspecific trimmed sequences have the deletions in the "middle"
        (deletions at the 5' end of the trimmed sequence could be confused with nontemplated
        nucleotides). These nonspecific sequences are also in other normalized sequences, so the
        deletions would be identified in those as well. Since nonspecific sequences tend to be
        shorter, it may not be possible to confidentally assign deletions in the middle of these
        sequences."""
        start_time = time.time()
        self.progress.new("Finding sequences with modification-induced deletions")
        self.progress.update("...")

        # Since query and target sequences must be the same length, grouping targets by length in
        # the dictionary speeds lookup.
        norm_trna_seq_target_dict = defaultdict(dict)
        for norm_trna_seq in [norm_trna_seq for norm_trna_seq in self.norm_trna_seqs if not norm_trna_seq.mod_seqs]:
            norm_trna_seq_len = len(norm_trna_seq.seq_string)
            norm_trna_seq_target_dict[norm_trna_seq_len][norm_trna_seq.seq_string] = norm_trna_seq

        norm_trunc_seq_target_dict = defaultdict(dict)
        for norm_trunc_seq_index, norm_trunc_seq in enumerate(self.norm_trunc_seqs):
            norm_trunc_seq_len = len(norm_trunc_seq.seq_string)
            norm_trunc_seq_target_dict[norm_trunc_seq_len][norm_trunc_seq.seq_string] = (norm_trunc_seq, norm_trunc_seq_index)

        query_dict = defaultdict(list)
        for mod_seq in self.mod_trna_seqs:
            for seq_string_with_del, del_config in self.get_seqs_with_dels(mod_seq):
                query_dict[seq_string_with_del].append((del_config, mod_seq))

        norm_trunc_seq_indices_to_remove = []
        trimmed_trunc_seq_represent_names = [trimmed_trunc_seq.represent_name for trimmed_trunc_seq in self.trimmed_trunc_seqs]
        trimmed_trunc_seq_indices_to_remove = []
        uniq_trunc_seq_represent_names = [uniq_trunc_seq.represent_name for uniq_trunc_seq in self.uniq_trunc_seqs]
        uniq_trunc_seq_indices_to_remove = []
        for query_seq_string, mod_seq_items in query_dict.items():
            if len(mod_seq_items) > 1:
                # Ignore sequences that can be produced by introducing deletions in different
                # modified sequences.
                continue

            # Search the normalized tRNA sequence targets.
            inner_dict = norm_trna_seq_target_dict[len(query_seq_string)]
            norm_trna_seq = None
            try:
                norm_trna_seq = inner_dict[query_seq_string]
            except KeyError:
                pass
            if norm_trna_seq:
                for del_config, mod_seq in mod_seq_items:
                    norm_trna_seq.mod_seqs.append(mod_seq)
                    # The feature profiles of many of the constituent trimmed (and unique) sequences
                    # are also changed, but since these may be shorter sequences not containing the
                    # deletion or nonspecific sequences found in other normalized sequences, the
                    # sequences are not marked as having changed profiles.
                    norm_trna_seq.profile_changed_by_del_analysis = True
                    # Normalized sequences with potential modification-induced deletions are added
                    # to uninitialized modified sequences.
                    mod_seq.norm_seqs_with_dels.append(norm_trna_seq)
                    mod_seq.del_configs.append(del_config)
                continue # The target pools are non-overlapping, so avoid searching the next one

            # Search the normalized truncated sequence targets.
            inner_dict = norm_trunc_seq_target_dict[len(query_seq_string)]
            norm_trunc_seq = None
            try:
                norm_trunc_seq, norm_trunc_seq_index = inner_dict[query_seq_string]
            except KeyError:
                pass
            if norm_trunc_seq:
                norm_trunc_seq_indices_to_remove.append(norm_trunc_seq_index)
                for trimmed_trunc_seq in norm_trunc_seq.trimmed_seqs:
                    trimmed_trunc_seq_indices_to_remove.append(trimmed_trunc_seq_represent_names.index(trimmed_trunc_seq.represent_name))
                    for uniq_trunc_seq in trimmed_trunc_seq.uniq_seqs:
                        uniq_trunc_seq_indices_to_remove.append(uniq_trunc_seq_represent_names.index(uniq_trunc_seq.represent_name))
                for del_config, mod_seq in mod_seq_items:
                    norm_trunc_seq.trunc_profile_recovered_by_del_analysis = True
                    norm_trunc_seq.mod_seqs.append(mod_seq)
                    mod_seq.norm_seqs_with_dels.append(norm_trunc_seq)
                    mod_seq.del_configs.append(del_config)

        for norm_trunc_seq_index in sorted(norm_trunc_seq_indices_to_remove, reverse=True):
            self.norm_trna_seqs.append(self.norm_trunc_seqs.pop(norm_trunc_seq_index))
        # The recovery of trimmed and unique sequences here is a little confusing, for the reason
        # mentioned in the docstring. The normalized sequence is identified as having deletions (and
        # thus being tRNA), but constituent trimmed sequences may be nonspecific, and therefore
        # found in other normalized sequences that are not recovered by this method.
        for trimmed_trunc_seq_index in sorted(set(trimmed_trunc_seq_indices_to_remove), reverse=True):
            trimmed_seq = self.trimmed_trunc_seqs.pop(trimmed_trunc_seq_index)
            trimmed_seq.trunc_profile_recovered_by_del_analysis = True
            self.trimmed_trna_seqs.append(trimmed_seq)
        for uniq_trunc_seq_index in sorted(set(uniq_trunc_seq_indices_to_remove), reverse=True):
            uniq_seq = self.uniq_trunc_seqs.pop(uniq_trunc_seq_index)
            uniq_seq.trunc_profile_recovered_by_del_analysis = True
            self.uniq_trna_seqs.append(uniq_seq)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed finding modification-induced deletions (min)",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def get_seqs_with_dels(self, mod_seq):
        """Generate in silico modified sequences with deletions at and/or around substitution sites.

        This method introduces in silico substitutions at substitution sites to create template
        sequences. Substitutions correspond to the nucleotides observed in the normalized sequences
        comprising the modified sequence, i.e., if 3 nucleotides are observed at a substitution
        site, then 3 are introduced in silico, rather than the maximum 4. For each template,
        self.introduce_dels is called to introduce deletions of different lengths at the possible
        configurations of substitution sites. This method then removes redundant sequences that may
        be produced from the different templates.

        Returns
        =======
        del_set : set
            A set of tuples with 2 elements. The first element is a sequence string containing
            deletions. The second element is a tuple of the indices of these deletions in the input
            sequence.
        """
        # Make template sequences with different nucleotides at substitution sites. Only consider
        # the observed substitution configurations.
        seq_strings_without_dels = set()
        longest_norm_seq_string = mod_seq.norm_seqs_without_dels[0].seq_string
        mod_seq_length = len(longest_norm_seq_string)
        sub_positions = mod_seq.sub_positions
        for norm_seq in mod_seq.norm_seqs_without_dels:
            norm_seq_string = norm_seq.seq_string
            norm_seq_start_in_mod_seq = mod_seq_length - len(norm_seq_string)
            altered_seq_string = longest_norm_seq_string
            for sub_pos in sub_positions:
                if sub_pos < norm_seq_start_in_mod_seq:
                    # This normalized sequence is shorter, lacking the substitution position.
                    continue
                nt = norm_seq.seq_string[sub_pos - norm_seq_start_in_mod_seq]
                altered_seq_string = altered_seq_string[: sub_pos] + nt + altered_seq_string[sub_pos + 1: ]
            seq_strings_without_dels.add(altered_seq_string)

        # Introduce deletions into each template sequence, potentially producing a number of new
        # sequences.
        del_dict = {}
        for seq_string in seq_strings_without_dels:
            del_dict_for_seq = self.introduce_dels(seq_string, sub_positions)
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


    def introduce_dels(self, seq_string, sub_positions):
        """Generate in silico sequences with deletions at and/or around substitution sites in the
        input sequence.

        This method is called by 'get_seqs_with_dels'.

        Parameters
        ==========
        seq_string : str
            The sequence in which deletions will be introduced

        sub_positions : list-like
            Where substitutions are located in the input sequence

        Returns
        =======
        del_dict : dict
            Each dict key is a sequence string containing deletions. Each dict value is a tuple of
            the indices of these deletions in the input sequence.
        """
        # Find all the ways deletions can be introduced into the sequence given the
        # parameterization.
        del_pos_configs = set()
        for num_del_sites in range(1, self.max_distinct_dels + 1):
            for del_locus_config in combinations(sub_positions, num_del_sites):
                for del_range_config in product(*[self.del_ranges for _ in range(num_del_sites)]):
                    del_positions = set()
                    for i, del_range in enumerate(del_range_config):
                        sub_pos = del_locus_config[i]
                        for del_pos_relative_to_sub in del_range:
                            del_pos = sub_pos + del_pos_relative_to_sub
                            if del_pos >= 0:
                                del_positions.add(del_pos)
                    if del_positions:
                        del_positions = sorted(del_positions)

                        # Remove any nominal deletions at the 5' end, as these could rightly be
                        # interpreted as unseen nucleotides preceding a fragment.
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

        # It is possible to generate the same sequence with deletions given different deletion
        # sites. For example, ACCG can become ACG by deleting either C. We resolve this complication
        # by choosing the most 5' deletion indices.
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


    def report_stats(self):
        """Add run statistics to the database, write them to the summary file, and write less
        technical statistics to the terminal."""

        profiled_trna_reads = 0
        trna_reads_with_threeprime_cca = 0
        trna_reads_with_threeprime_cc = 0
        trna_reads_with_threeprime_c = 0
        trna_reads_with_threeprime_ccan_ccann = 0
        profiled_trna_reads_containing_anticodon = 0
        full_length_trna_reads = 0
        trna_reads_with_extrapolated_fiveprime_feature = 0
        min_length_of_long_fiveprime_extension = self.min_length_of_long_fiveprime_extension
        profiled_trna_reads_with_short_fiveprime_extension = 0
        profiled_trna_reads_with_long_fiveprime_extension = 0
        trna_reads_with_trunc_profile_recovered_by_derep = 0
        trna_reads_with_trunc_profile_recovered_by_del_analysis = 0
        uniq_profiled_trna_seqs = 0
        uniq_trna_seqs_with_threeprime_cca = 0
        uniq_trna_seqs_with_threeprime_cc = 0
        uniq_trna_seqs_with_threeprime_c = 0
        uniq_trna_seqs_with_threeprime_ccan_ccann = 0
        uniq_trna_seqs_containing_anticodon = 0
        full_length_uniq_seqs = 0
        uniq_seqs_with_extrapolated_fiveprime_feature = 0
        uniq_seqs_with_short_fiveprime_extension = 0
        uniq_seqs_with_long_fiveprime_extension = 0
        uniq_seqs_with_trunc_profile_recovered_by_derep = 0
        uniq_seqs_with_trunc_profile_recovered_by_del_analysis = 0
        interior_mapped_reads = 0
        fiveprime_mapped_reads = 0
        interior_mapped_uniq_seqs = 0
        fiveprime_mapped_uniq_seqs = 0
        for uniq_trna_seq in self.uniq_trna_seqs:
            read_count = uniq_trna_seq.read_count
            if uniq_trna_seq.id_method == 0: # profiled
                profiled_trna_reads += read_count
                uniq_profiled_trna_seqs += 1

                if uniq_trna_seq.acceptor_length == 3:
                    trna_reads_with_threeprime_cca += read_count
                    uniq_trna_seqs_with_threeprime_cca += 1
                elif uniq_trna_seq.acceptor_length == 2:
                    trna_reads_with_threeprime_cc += read_count
                    uniq_trna_seqs_with_threeprime_cc += 1
                elif uniq_trna_seq.acceptor_length == 1:
                    trna_reads_with_threeprime_c += read_count
                    uniq_trna_seqs_with_threeprime_c += 1
                elif uniq_trna_seq.acceptor_length == 4 or uniq_trna_seq.acceptor_length == 5:
                    trna_reads_with_threeprime_ccan_ccann += read_count
                    uniq_trna_seqs_with_threeprime_ccan_ccann += 1
                else:
                    raise ConfigError("Acceptor variants are expected to be 1-5 nucleotides in length.")

                if uniq_trna_seq.contains_anticodon:
                    profiled_trna_reads_containing_anticodon += read_count
                    uniq_trna_seqs_containing_anticodon += 1

                if uniq_trna_seq.has_complete_feature_set:
                    full_length_trna_reads += read_count
                    full_length_uniq_seqs += 1

                if uniq_trna_seq.num_extrapolated_fiveprime_nts:
                    trna_reads_with_extrapolated_fiveprime_feature += read_count
                    uniq_seqs_with_extrapolated_fiveprime_feature += 1

                if uniq_trna_seq.extra_fiveprime_length:
                    if uniq_trna_seq.extra_fiveprime_length >= min_length_of_long_fiveprime_extension:
                        profiled_trna_reads_with_long_fiveprime_extension += read_count
                        uniq_seqs_with_long_fiveprime_extension += 1
                    else:
                        profiled_trna_reads_with_short_fiveprime_extension += read_count
                        uniq_seqs_with_short_fiveprime_extension += 1

                if uniq_trna_seq.trunc_profile_index:
                    if uniq_trna_seq.trunc_profile_recovered_by_derep:
                        trna_reads_with_trunc_profile_recovered_by_derep += read_count
                        uniq_seqs_with_trunc_profile_recovered_by_derep += 1
                    elif uniq_trna_seq.trunc_profile_recovered_by_del_analysis:
                        trna_reads_with_trunc_profile_recovered_by_del_analysis += read_count
                        uniq_seqs_with_trunc_profile_recovered_by_del_analysis += 1
                    else:
                        ConfigError("It should not be possible for a UniqueSeq object listed as tRNA in self.uniq_trna_seqs "
                                    "to have a value of trunc_profile_index but not show "
                                    "that it was recovered by dereplication or deletion analysis.")
            elif uniq_trna_seq.id_method == 1: # mapped
                if uniq_trna_seq.extra_fiveprime_length:
                    fiveprime_mapped_reads += read_count
                    fiveprime_mapped_uniq_seqs += 1
                else:
                    interior_mapped_reads += read_count
                    interior_mapped_uniq_seqs += 1
            else:
                raise ConfigError("The only recognized tRNA identification methods are 0 (profiled) and 1 (mapped).")
        total_trna_reads = profiled_trna_reads + interior_mapped_reads + fiveprime_mapped_reads
        total_uniq_trna_seqs = uniq_profiled_trna_seqs + interior_mapped_uniq_seqs + fiveprime_mapped_uniq_seqs
        mean_profiled_reads_per_uniq_seq = profiled_trna_reads / uniq_profiled_trna_seqs
        mean_mapped_reads_per_uniq_seq = (fiveprime_mapped_reads + interior_mapped_reads) / (fiveprime_mapped_uniq_seqs + interior_mapped_uniq_seqs)

        # The TOTAL count of trimmed sequence objects combines apples and oranges (profiled and
        # mapped sequences) and so should not be reported. Separate trimmed sequences are formed for
        # each unique mapped sequence, because the 5' extension of a mapped sequence can be most of
        # its length, so consolidation of unique mapped sequences with different 5' extensions would
        # result in very short trimmed sequences of only a handful of nucleotides.
        trimmed_profiled_seqs_containing_anticodon = 0
        full_length_trimmed_seqs = 0
        trimmed_seqs_with_extrapolated_fiveprime_feature = 0
        trimmed_seqs_with_trunc_profile_recovered_by_derep = 0
        trimmed_seqs_with_trunc_profile_recovered_by_del_analysis = 0
        for trimmed_trna_seq in self.trimmed_trna_seqs:
            represent_uniq_seq = trimmed_trna_seq.uniq_seqs[0]
            if trimmed_trna_seq.contains_anticodon:
                trimmed_profiled_seqs_containing_anticodon += 1

            if trimmed_trna_seq.has_complete_feature_set:
                full_length_trimmed_seqs += 1

            if represent_uniq_seq.num_extrapolated_fiveprime_nts:
                trimmed_seqs_with_extrapolated_fiveprime_feature += 1

            if trimmed_trna_seq.has_trunc_profile:
                if trimmed_trna_seq.trunc_profile_recovered_by_derep:
                    trimmed_seqs_with_trunc_profile_recovered_by_derep += 1
                elif trimmed_trna_seq.trunc_profile_recovered_by_del_analysis:
                    trimmed_seqs_with_trunc_profile_recovered_by_del_analysis += 1
                else:
                    raise ConfigError("It should not be possible for a TrimmedSeq object listed as tRNA in self.trimmed_trna_seqs "
                                      "to have a truncated profile but not show "
                                      "that it was recovered by dereplication or deletion analysis.")

        norm_trna_seqs = len(self.norm_trna_seqs)
        norm_trna_seqs_without_mods = 0
        norm_trna_seqs_containing_anticodon = 0
        full_length_norm_seqs = 0
        norm_seqs_with_profile_changed_by_del_analysis = 0
        norm_seqs_with_trunc_profile_recovered_by_del_analysis = 0
        for norm_trna_seq in self.norm_trna_seqs:
            if not norm_trna_seq.mod_seqs:
                norm_trna_seqs_without_mods += 1

            if norm_trna_seq.trimmed_seqs[0].contains_anticodon:
                norm_trna_seqs_containing_anticodon += 1

            if norm_trna_seq.has_complete_feature_set:
                full_length_norm_seqs += 1

            if norm_trna_seq.profile_changed_by_del_analysis:
                norm_seqs_with_profile_changed_by_del_analysis += 1

            if norm_trna_seq.trunc_profile_recovered_by_del_analysis:
                norm_seqs_with_trunc_profile_recovered_by_del_analysis += 1
        norm_seqs_with_dels = norm_seqs_with_profile_changed_by_del_analysis + norm_seqs_with_trunc_profile_recovered_by_del_analysis

        mod_trna_seqs = len(self.mod_trna_seqs)
        mod_seqs_with_dels = 0
        for mod_trna_seq in self.mod_trna_seqs:
            if mod_trna_seq.del_configs:
                mod_seqs_with_dels += 1

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db.set_meta_value('total_trna_reads', total_trna_reads)
        trnaseq_db.db.set_meta_value('profiled_trna_reads', profiled_trna_reads)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_cca', trna_reads_with_threeprime_cca)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_cc', trna_reads_with_threeprime_cc)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_c', trna_reads_with_threeprime_c)
        trnaseq_db.db.set_meta_value('trna_reads_with_threeprime_ccan_ccann', trna_reads_with_threeprime_ccan_ccann)
        trnaseq_db.db.set_meta_value('profiled_trna_reads_containing_anticodon', profiled_trna_reads_containing_anticodon)
        trnaseq_db.db.set_meta_value('full_length_trna_reads', full_length_trna_reads)
        trnaseq_db.db.set_meta_value('trna_reads_with_extrapolated_fiveprime_feature', trna_reads_with_extrapolated_fiveprime_feature)
        trnaseq_db.db.set_meta_value('profiled_trna_reads_with_short_fiveprime_extension', profiled_trna_reads_with_short_fiveprime_extension)
        trnaseq_db.db.set_meta_value('profiled_trna_reads_with_long_fiveprime_extension', profiled_trna_reads_with_long_fiveprime_extension)
        trnaseq_db.db.set_meta_value('trna_reads_with_truncated_profile_recovered_by_dereplication', trna_reads_with_trunc_profile_recovered_by_derep)
        trnaseq_db.db.set_meta_value('trna_reads_with_truncated_profile_recovered_by_deletion_analysis', trna_reads_with_trunc_profile_recovered_by_del_analysis)
        trnaseq_db.db.set_meta_value('interior_mapped_reads', interior_mapped_reads)
        trnaseq_db.db.set_meta_value('fiveprime_mapped_reads', fiveprime_mapped_reads)
        trnaseq_db.db.set_meta_value('total_unique_trna_seqs', total_uniq_trna_seqs)
        trnaseq_db.db.set_meta_value('unique_profiled_trna_seqs', uniq_profiled_trna_seqs)
        trnaseq_db.db.set_meta_value('unique_trna_seqs_with_threeprime_cca', uniq_trna_seqs_with_threeprime_cca)
        trnaseq_db.db.set_meta_value('unique_trna_seqs_with_threeprime_cc', uniq_trna_seqs_with_threeprime_cc)
        trnaseq_db.db.set_meta_value('unique_trna_seqs_with_threeprime_c', uniq_trna_seqs_with_threeprime_c)
        trnaseq_db.db.set_meta_value('unique_trna_seqs_with_threeprime_ccan_ccann', uniq_trna_seqs_with_threeprime_ccan_ccann)
        trnaseq_db.db.set_meta_value('unique_trna_seqs_containing_anticodon', uniq_trna_seqs_containing_anticodon)
        trnaseq_db.db.set_meta_value('full_length_unique_seqs', full_length_uniq_seqs)
        trnaseq_db.db.set_meta_value('unique_seqs_with_extrapolated_fiveprime_feature', uniq_seqs_with_extrapolated_fiveprime_feature)
        trnaseq_db.db.set_meta_value('unique_seqs_with_short_fiveprime_extension', uniq_seqs_with_short_fiveprime_extension)
        trnaseq_db.db.set_meta_value('unique_seqs_with_long_fiveprime_extension', uniq_seqs_with_long_fiveprime_extension)
        trnaseq_db.db.set_meta_value('unique_seqs_with_truncated_profile_recovered_by_dereplication', uniq_seqs_with_trunc_profile_recovered_by_derep)
        trnaseq_db.db.set_meta_value('unique_seqs_with_truncated_profile_recovered_by_deletion_analysis', uniq_seqs_with_trunc_profile_recovered_by_del_analysis)
        trnaseq_db.db.set_meta_value('interior_mapped_uniq_seqs', interior_mapped_uniq_seqs)
        trnaseq_db.db.set_meta_value('fiveprime_mapped_uniq_seqs', fiveprime_mapped_uniq_seqs)
        trnaseq_db.db.set_meta_value('trimmed_profiled_seqs_containing_anticodon', trimmed_profiled_seqs_containing_anticodon)
        trnaseq_db.db.set_meta_value('full_length_trimmed_seqs', full_length_trimmed_seqs)
        trnaseq_db.db.set_meta_value('trimmed_seqs_with_extrapolated_fiveprime_feature', trimmed_seqs_with_extrapolated_fiveprime_feature)
        trnaseq_db.db.set_meta_value('trimmed_seqs_with_truncated_profile_recovered_by_dereplication', trimmed_seqs_with_trunc_profile_recovered_by_derep)
        trnaseq_db.db.set_meta_value('trimmed_seqs_with_truncated_profile_recovered_by_deletion_analysis', trimmed_seqs_with_trunc_profile_recovered_by_del_analysis)
        trnaseq_db.db.set_meta_value('normalized_trna_seqs', norm_trna_seqs)
        trnaseq_db.db.set_meta_value('normalized_trna_seqs_without_potential_modifications', norm_trna_seqs_without_mods)
        trnaseq_db.db.set_meta_value('normalized_trna_seqs_containing_anticodon', norm_trna_seqs_containing_anticodon)
        trnaseq_db.db.set_meta_value('full_length_normalized_seqs', full_length_norm_seqs)
        trnaseq_db.db.set_meta_value('normalized_seqs_with_profile_changed_by_deletion_analysis', norm_seqs_with_profile_changed_by_del_analysis)
        trnaseq_db.db.set_meta_value('normalized_seqs_with_truncated_profile_recovered_by_deletion_analysis', norm_seqs_with_trunc_profile_recovered_by_del_analysis)
        trnaseq_db.db.set_meta_value('normalized_seqs_with_deletions', norm_seqs_with_dels)
        trnaseq_db.db.set_meta_value('potentially_modified_seqs', mod_trna_seqs)
        trnaseq_db.db.set_meta_value('potentially_modified_seqs_with_deletions', mod_seqs_with_dels)
        trnaseq_db.disconnect()

        # Use debug flag to write additional statistics of technical but not general interest to the
        # summary text file.
        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Total tRNA reads (profiled and mapped)", total_trna_reads))
            f.write(self.get_summary_line("Reads profiled as tRNA", profiled_trna_reads))
            f.write(self.get_summary_line("Profiled reads ending in 3'-CCA", trna_reads_with_threeprime_cca))
            f.write(self.get_summary_line("Profiled reads ending in 3'-CC", trna_reads_with_threeprime_cc))
            f.write(self.get_summary_line("Profiled reads ending in 3'-C", trna_reads_with_threeprime_c))
            f.write(self.get_summary_line("Profiled reads ending in 3'-CCAN/CCANN", trna_reads_with_threeprime_ccan_ccann))
            f.write(self.get_summary_line("Profiled reads containing anticodon", profiled_trna_reads_containing_anticodon))
            f.write(self.get_summary_line("Profiled reads spanning acceptor stem", full_length_trna_reads))
            if anvio.DEBUG:
                f.write(self.get_summary_line("Profiled reads with extrapolated 5' feature", trna_reads_with_extrapolated_fiveprime_feature))
            f.write(self.get_summary_line("Profiled reads with 1-%d extra 5' bases" % (self.min_length_of_long_fiveprime_extension - 1),
                                          profiled_trna_reads_with_short_fiveprime_extension))
            f.write(self.get_summary_line("Profiled reads with >%d extra 5' bases" % (self.min_length_of_long_fiveprime_extension - 1),
                                          profiled_trna_reads_with_long_fiveprime_extension))
            if anvio.DEBUG:
                f.write(self.get_summary_line("Profiled reads with truncated profile recovered by 3'-dereplication",
                                              trna_reads_with_trunc_profile_recovered_by_derep))
                f.write(self.get_summary_line("Profiled reads with truncated profile recovered by deletion analysis",
                                              trna_reads_with_trunc_profile_recovered_by_del_analysis))
            f.write(self.get_summary_line("Reads mapped to tRNA interior", interior_mapped_reads))
            f.write(self.get_summary_line("Reads mapped to tRNA with extra 5' bases", fiveprime_mapped_reads))
            f.write(self.get_summary_line("Total unique tRNA sequences", total_uniq_trna_seqs))
            f.write(self.get_summary_line("Unique profiled tRNA sequences", uniq_profiled_trna_seqs))
            f.write(self.get_summary_line("Unique tRNA sequences ending in 3'-CCA", uniq_trna_seqs_with_threeprime_cca))
            f.write(self.get_summary_line("Unique tRNA sequences ending in 3'-CC", uniq_trna_seqs_with_threeprime_cc))
            f.write(self.get_summary_line("Unique tRNA sequences ending in 3'-C", uniq_trna_seqs_with_threeprime_c))
            f.write(self.get_summary_line("Unique tRNA sequences ending in 3'-CCAN/CCANNs", uniq_trna_seqs_with_threeprime_ccan_ccann))
            f.write(self.get_summary_line("Unique tRNA sequences containing anticodon", uniq_trna_seqs_containing_anticodon))
            f.write(self.get_summary_line("Unique tRNA sequences spanning acceptor stem", full_length_uniq_seqs))
            if anvio.DEBUG:
                f.write(self.get_summary_line("Unique tRNA sequences with extrapolated 5' feature", uniq_seqs_with_extrapolated_fiveprime_feature))
            f.write(self.get_summary_line("Unique tRNA sequences with 1-%d extra 5' bases" % (self.min_length_of_long_fiveprime_extension - 1),
                                          uniq_seqs_with_short_fiveprime_extension))
            f.write(self.get_summary_line("Unique tRNA sequences with >%d extra 5' bases" % (self.min_length_of_long_fiveprime_extension - 1),
                                          uniq_seqs_with_long_fiveprime_extension))
            if anvio.DEBUG:
                f.write(self.get_summary_line("Unique tRNA sequences with truncated profile recovered by 3'-dereplication",
                                              uniq_seqs_with_trunc_profile_recovered_by_derep))
                f.write(self.get_summary_line("Unique tRNA sequences with truncated profile recovered by deletion analysis",
                                              uniq_seqs_with_trunc_profile_recovered_by_del_analysis))
            f.write(self.get_summary_line("Mean profiled reads per unique sequence", mean_profiled_reads_per_uniq_seq))
            f.write(self.get_summary_line("Unique sequences mapped to tRNA interior", interior_mapped_uniq_seqs))
            f.write(self.get_summary_line("Unique sequences mapped to tRNA with extra 5' bases", fiveprime_mapped_uniq_seqs))
            f.write(self.get_summary_line("Mean mapped reads per unique sequence", mean_mapped_reads_per_uniq_seq))
            f.write(self.get_summary_line("Trimmed tRNA sequences containing anticodon", trimmed_profiled_seqs_containing_anticodon))
            f.write(self.get_summary_line("Trimmed tRNA sequences spanning acceptor stem", full_length_trimmed_seqs))
            if anvio.DEBUG:
                f.write(self.get_summary_line("Trimmed tRNA sequences with extrapolated 5' feature", trimmed_seqs_with_extrapolated_fiveprime_feature))
                f.write(self.get_summary_line("Trimmed tRNA sequences with truncated profile recovered by 3'-dereplication",
                                              trimmed_seqs_with_trunc_profile_recovered_by_derep))
                f.write(self.get_summary_line("Trimmed tRNA sequences with truncated profile recovered by deletion analysis",
                                              trimmed_seqs_with_trunc_profile_recovered_by_del_analysis))
            f.write(self.get_summary_line("Normalized tRNA sequences", norm_trna_seqs))
            f.write(self.get_summary_line("Normalized tRNA sequences without potential modifications", norm_trna_seqs_without_mods))
            f.write(self.get_summary_line("Normalized tRNA sequences containing anticodon", norm_trna_seqs_containing_anticodon))
            f.write(self.get_summary_line("Normalized tRNA sequences spanning acceptor stem", full_length_norm_seqs))
            if anvio.DEBUG:
                f.write(self.get_summary_line("Normalized tRNA sequences with profile changed by deletion analysis",
                                              norm_seqs_with_profile_changed_by_del_analysis))
                f.write(self.get_summary_line("Normalized tRNA sequences with truncated profile recovered by deletion analysis",
                                              norm_seqs_with_trunc_profile_recovered_by_del_analysis))
            f.write(self.get_summary_line("Normalized sequences with deletions", norm_seqs_with_dels))
            f.write(self.get_summary_line("Potentially modified sequences", mod_trna_seqs))
            f.write(self.get_summary_line("Potentially modified sequences with deletions", mod_seqs_with_dels))


    def write_sequences_table(self):
        self.progress.new("Writing tRNA-seq database table of unique tRNA sequences")
        self.progress.update("...")

        sequences_table_entries = []
        for uniq_seq in self.uniq_trna_seqs:
            sequences_table_entries.append(
                (uniq_seq.represent_name,
                 uniq_seq.read_count,
                 uniq_seq.id_method,
                 0 if uniq_seq.trunc_profile_index is None else 1,
                 uniq_seq.seq_string)
            )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('sequences')
            trnaseq_db.db.create_table('sequences',
                                       tables.trnaseq_sequences_table_structure,
                                       tables.trnaseq_sequences_table_types)
        if sequences_table_entries:
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                    % ('sequences', ','.join('?' * len(tables.trnaseq_sequences_table_structure))),
                                    sequences_table_entries)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_trimmed_table(self):
        self.progress.new("Writing tRNA-seq database table of trimmed tRNA sequences")
        self.progress.update("...")

        trimmed_table_entries = []
        for trimmed_seq in self.trimmed_trna_seqs:
            trimmed_table_entries.append(
                (trimmed_seq.represent_name,
                 len(trimmed_seq.uniq_seqs),
                 trimmed_seq.read_count,
                 trimmed_seq.id_method,
                 trimmed_seq.has_trunc_profile,
                 trimmed_seq.trunc_profile_recovered_by_derep,
                 trimmed_seq.seq_string,
                 trimmed_seq.norm_seq_count,
                 trimmed_seq.uniq_with_extra_fiveprime_count,
                 trimmed_seq.read_with_extra_fiveprime_count)
                + tuple([v for v in trimmed_seq.read_acceptor_variant_count_dict.values()])
            )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('trimmed')
            trnaseq_db.db.create_table('trimmed',
                                       tables.trnaseq_trimmed_table_structure,
                                       tables.trnaseq_trimmed_table_types)
        if trimmed_table_entries:
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                    % ('trimmed', ','.join('?' * len(tables.trnaseq_trimmed_table_structure))),
                                    trimmed_table_entries)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_normalized_table(self):
        self.progress.new("Writing tRNA-seq database table of fragment-dereplicated tRNA sequences")
        self.progress.update("...")

        norm_table_entries = []
        for norm_seq in self.norm_trna_seqs:
            specific_long_fiveprime_extensions = ''
            specific_long_fiveprime_extension_read_counts = ''
            for fiveprime_extension_string, read_count in sorted(norm_seq.specific_long_fiveprime_extension_dict.items(),
                                                                 key=lambda item: -len(item[0])):
                specific_long_fiveprime_extensions += fiveprime_extension_string + ','
                specific_long_fiveprime_extension_read_counts += str(read_count) + ','

            nonspecific_long_fiveprime_extensions = ''
            nonspecific_long_fiveprime_extension_read_counts = ''
            for fiveprime_extension_string, read_count in sorted(norm_seq.nonspecific_long_fiveprime_extension_dict.items(),
                                                                 key=lambda item: -len(item[0])):
                nonspecific_long_fiveprime_extensions += fiveprime_extension_string + ','
                nonspecific_long_fiveprime_extension_read_counts += str(read_count) + ','

            norm_table_entries.append(
                (norm_seq.represent_name,
                 len(norm_seq.trimmed_seqs),
                 norm_seq.mean_specific_cov,
                 norm_seq.mean_nonspecific_cov,
                 ','.join(map(str, norm_seq.specific_covs)) + ',',
                 ','.join(map(str, norm_seq.nonspecific_covs)) + ',',
                 len(norm_seq.mod_seqs),
                 norm_seq.specific_read_count,
                 norm_seq.nonspecific_read_count,
                 norm_seq.profile_changed_by_del_analysis,
                 norm_seq.trunc_profile_recovered_by_del_analysis,
                 norm_seq.count_of_specific_reads_with_extra_fiveprime,
                 norm_seq.count_of_nonspecific_reads_with_extra_fiveprime,
                 norm_seq.specific_mapped_read_count,
                 norm_seq.nonspecific_mapped_read_count,
                 specific_long_fiveprime_extensions,
                 specific_long_fiveprime_extension_read_counts,
                 nonspecific_long_fiveprime_extensions,
                 nonspecific_long_fiveprime_extension_read_counts)
                + tuple(norm_seq.specific_read_acceptor_variant_count_dict.values())
                + tuple(norm_seq.nonspecific_read_acceptor_variant_count_dict.values())
            )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('normalized')
            trnaseq_db.db.create_table('normalized',
                                       tables.trnaseq_normalized_table_structure,
                                       tables.trnaseq_normalized_table_types)
        if norm_table_entries:
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                    % ('normalized', ','.join('?' * len(tables.trnaseq_normalized_table_structure))),
                                    norm_table_entries)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_modified_table(self):
        self.progress.new("Writing tRNA-seq database table of modified tRNA sequences")
        self.progress.update("...")

        mod_table_entries = []
        for mod_seq in self.mod_trna_seqs:
            specific_long_fiveprime_extensions = ''
            specific_long_fiveprime_extension_read_counts = ''
            for fiveprime_extension_string, read_count in sorted(mod_seq.specific_long_fiveprime_extension_dict.items(),
                                                                 key=lambda item: -len(item[0])):
                specific_long_fiveprime_extensions += fiveprime_extension_string + ','
                specific_long_fiveprime_extension_read_counts += str(read_count) + ','

            nonspecific_long_fiveprime_extensions = ''
            nonspecific_long_fiveprime_extension_read_counts = ''
            for fiveprime_extension_string, read_count in sorted(mod_seq.nonspecific_long_fiveprime_extension_dict.items(),
                                                                 key=lambda item: -len(item[0])):
                nonspecific_long_fiveprime_extensions += fiveprime_extension_string + ','
                nonspecific_long_fiveprime_extension_read_counts += str(read_count) + ','

            mod_table_entries.append(
                (mod_seq.represent_name,
                 mod_seq.mean_specific_cov,
                 mod_seq.mean_nonspecific_cov,
                 ','.join(map(str, mod_seq.specific_covs)) + ',',
                 ','.join(map(str, mod_seq.nonspecific_covs)) + ',',
                 ','.join([str(sub_pos) for sub_pos in mod_seq.sub_positions]) + ',')
                + tuple([','.join(map(str, mod_seq.specific_sub_covs[:, i - 1])) + ',' for i in INT_NT_DICT])
                + tuple([','.join(map(str, mod_seq.nonspecific_sub_covs[:, i - 1])) + ',' for i in INT_NT_DICT])
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
                   mod_seq.count_of_specific_reads_with_extra_fiveprime,
                   mod_seq.count_of_nonspecific_reads_with_extra_fiveprime,
                   mod_seq.specific_mapped_read_count,
                   mod_seq.nonspecific_mapped_read_count,
                   specific_long_fiveprime_extensions,
                   specific_long_fiveprime_extension_read_counts,
                   nonspecific_long_fiveprime_extensions,
                   nonspecific_long_fiveprime_extension_read_counts)
                + tuple(mod_seq.specific_read_acceptor_variant_count_dict.values())
                + tuple(mod_seq.nonspecific_read_acceptor_variant_count_dict.values())
            )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('modified')
            trnaseq_db.db.create_table('modified',
                                       tables.trnaseq_modified_table_structure,
                                       tables.trnaseq_modified_table_types)
        if mod_table_entries:
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                    % ('modified', ','.join('?' * len(tables.trnaseq_modified_table_structure))),
                                    mod_table_entries)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_uniq_nontrna_supplement(self):
        self.progress.new("Writing a file of unique sequences not identified as tRNA")
        self.progress.update("...")

        with open(self.uniq_nontrna_path, 'w') as nontrna_file:
            nontrna_file.write("\t".join(self.UNIQ_NONTRNA_HEADER) + "\n")
            for uniq_seq in self.uniq_nontrna_seqs:
                nontrna_file.write(uniq_seq.represent_name + "\t"
                                   + str(uniq_seq.read_count) + "\t"
                                   + "\t"
                                   + uniq_seq.seq_string + "\n")
            for uniq_seq in self.uniq_trunc_seqs:
                nontrna_file.write(uniq_seq.represent_name + "\t"
                                   + str(uniq_seq.read_count) + "\t"
                                   + str(uniq_seq.trunc_profile_index) + "\t"
                                   + uniq_seq.seq_string + "\n")

        self.progress.end()

        self.run.info("Unique non-tRNA supplement", self.uniq_nontrna_path)


    def write_trimmed_supplement(self):
        """Write a supplementary file showing the spectrum of 5'/3' extensions of core trimmed tRNA
        sequences.

        Mapped, as opposed to profiled, trimmed sequences are not considered, as each unique mapped
        sequence with potentially varying 5' ends generates its own trimmed sequence object. The 5'
        extension of a mapped sequence may represent all but a small number of nucleotides in the
        sequence, so sequences identical in the non-5' section are not dereplicated into a single
        trimmed sequence object.
        """
        self.progress.new("Writing a file showing the 5'/3' ends of each trimmed profiled tRNA sequence")
        self.progress.update("...")

        with open(self.trimmed_ends_path, 'w') as trimmed_file:
            trimmed_file.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for trimmed_seq in sorted(self.trimmed_trna_seqs, key=lambda trimmed_seq: -trimmed_seq.read_count):
                if trimmed_seq.id_method == 1:
                    continue

                represent_name = trimmed_seq.represent_name
                for uniq_seq in sorted(trimmed_seq.uniq_seqs,
                                       key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, -uniq_seq.acceptor_length)):
                    trimmed_file.write(represent_name + "\t"
                                       + uniq_seq.represent_name + "\t"
                                       + uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length] + "\t"
                                       + uniq_seq.seq_string[len(uniq_seq.seq_string) - uniq_seq.acceptor_length: ] + "\t"
                                       + str(uniq_seq.read_count) + "\n")

        self.progress.end()

        self.run.info("Trimmed tRNA supplement", self.trimmed_ends_path)


def profile_worker(input_queue, output_queue, profiler):
    """This client for `trnaidentifier.Profiler.profile` is located outside the `TRNASeqDataset`
    class to allow multiprocessing."""
    while True:
        seq_string, seq_name = input_queue.get()
        output_queue.put(profiler.profile(seq_string, name=seq_name))


class NormalizedSeqSummary(object):
    """Relevant data from normalized sequences stored in anvi'o tRNA-seq databases are reloaded into
    these objects."""

    __slots__ = (
        'name',
        'sample_id',
        'seq_string',
        'threshold_feature_start',
        'anticodon_seq_string',
        'feature_threshold_start',
        'mean_specific_cov',
        'specific_covs',
        'nonspecific_covs',
        'specific_nt_covs_dict',
        'nonspecific_nt_covs_dict',
        'mod_seq_summary'
    )

    def __init__(self):
        self.name = None
        self.sample_id = None
        self.seq_string = None
        self.threshold_feature_start = None
        self.anticodon_seq_string = None
        self.feature_threshold_start = None
        self.mean_specific_cov = None
        self.specific_covs = None
        self.nonspecific_covs = None
        self.specific_nt_covs_dict = None
        self.nonspecific_nt_covs_dict = None
        self.mod_seq_summary = None


class ModifiedSeqSummary(object):
    """Relevant data from modified sequences stored in anvi'o tRNA-seq databases are reloaded into
    these objects."""

    __slots__ = (
        'name',
        'sample_id',
        'consensus_seq_string',
        'sub_positions',
        'specific_nt_covs_dict',
        'nonspecific_nt_covs_dict',
        'specific_covs',
        'nonspecific_covs',
        'specific_del_covs',
        'nonspecific_del_covs',
        'norm_seq_summaries'
    )

    def __init__(self):
        name = None
        sample_id = None
        consensus_seq_string = None
        sub_positions = None
        specific_nt_covs_dict = None
        nonspecific_nt_covs_dict = None
        specific_covs = None
        nonspecific_covs = None
        specific_del_covs = None
        nonspecific_del_covs = None
        norm_seq_summaries = None


class SeedSeq(object):

    __slots__ = (
        'name',
        'seq_string',
        'meets_feature_threshold',
        'unmod_norm_seq_summaries',
        'mod_seq_summaries',
        'anticodon_seq_string',
        'total_specific_covs',
        'total_nonspecific_covs',
        'total_mean_specific_cov',
        'total_mean_nonspecific_cov',
        'sample_specific_covs_dict',
        'sample_nonspecific_covs_dict',
        'sample_summed_covs_dict',
        'sample_specific_nt_covs_dict',
        'sample_nonspecific_nt_covs_dict',
        'sample_summed_nt_covs_dict',
        'sample_mean_specific_cov_dict',
        'sample_mean_nonspecific_cov_dict',
        'sample_mean_summed_cov_dict',
        'sample_std_specific_cov_dict',
        'sample_std_nonspecific_cov_dict',
        'sample_std_summed_cov_dict',
        'sample_specific_abundances_dict',
        'sample_nonspecific_abundances_dict',
        'sample_summed_abundances_dict',
        'sample_specific_relative_abundances_dict',
        'sample_nonspecific_relative_abundances_dict',
        'sample_summed_relative_abundances_dict',
        'sample_specific_detection_dict',
        'sample_nonspecific_detection_dict',
        'sample_summed_detection_dict',
        'sample_mean_Q2Q3_specific_cov_dict',
        'sample_mean_Q2Q3_nonspecific_cov_dict',
        'sample_mean_Q2Q3_summed_cov_dict',
        'sample_normalized_mean_Q2Q3_specific_cov_dict',
        'sample_normalized_mean_Q2Q3_nonspecific_cov_dict',
        'sample_normalized_mean_Q2Q3_summed_cov_dict',
        'sample_specific_max_normalized_ratio_dict',
        'sample_nonspecific_max_normalized_ratio_dict',
        'sample_summed_max_normalized_ratio_dict',
        'gc_fraction',
        'sample_sub_positions_dict',
        'total_mod_positions',
        'sample_mod_positions_dict',
        'sample_variability_dict',
        'sample_specific_dels_dict',
        'sample_nonspecific_dels_dict'
    )

    def __init__(self):
        self.name = None
        self.seq_string = None
        self.meets_feature_threshold = None
        self.unmod_norm_seq_summaries = None
        self.mod_seq_summaries = None
        self.anticodon_seq_string = None
        self.total_specific_covs = None
        self.total_nonspecific_covs = None
        self.total_mean_specific_cov = None
        self.total_mean_nonspecific_cov = None
        self.sample_specific_covs_dict = None
        self.sample_nonspecific_covs_dict = None
        self.sample_summed_covs_dict = None
        self.sample_specific_nt_covs_dict = None
        self.sample_nonspecific_nt_covs_dict = None
        self.sample_summed_nt_covs_dict = None
        self.sample_mean_specific_cov_dict = None
        self.sample_mean_nonspecific_cov_dict = None
        self.sample_mean_summed_cov_dict = None
        self.sample_std_specific_cov_dict = None
        self.sample_std_nonspecific_cov_dict = None
        self.sample_std_summed_cov_dict = None
        self.sample_specific_abundances_dict = None
        self.sample_nonspecific_abundances_dict = None
        self.sample_summed_abundances_dict = None
        self.sample_specific_relative_abundances_dict = None
        self.sample_nonspecific_relative_abundances_dict = None
        self.sample_summed_relative_abundances_dict = None
        self.sample_specific_detection_dict = None
        self.sample_nonspecific_detection_dict = None
        self.sample_summed_detection_dict = None
        self.sample_mean_Q2Q3_specific_cov_dict = None
        self.sample_mean_Q2Q3_nonspecific_cov_dict = None
        self.sample_mean_Q2Q3_summed_cov_dict = None
        self.sample_normalized_mean_Q2Q3_specific_cov_dict = None
        self.sample_normalized_mean_Q2Q3_nonspecific_cov_dict = None
        self.sample_normalized_mean_Q2Q3_summed_cov_dict = None
        self.sample_specific_max_normalized_ratio_dict = None
        self.sample_nonspecific_max_normalized_ratio_dict = None
        self.sample_summed_max_normalized_ratio_dict = None
        self.gc_fraction = None
        self.sample_sub_positions_dict = None
        self.total_mod_positions = None
        self.sample_mod_positions_dict = None
        self.sample_variability_dict = None
        self.sample_specific_dels_dict = None
        self.sample_nonspecific_dels_dict = None


class DatabaseConverter(object):
    """Converts tRNA-seq database(s) into contigs, auxiliary, and profile databases.

    "Contigs" in this context are tRNA seed sequences representing tRNA identified in the collection
    of samples.
    """

    # The columns needed from tables of a tRNA-seq database.
    FEATURE_TABLE_COLS_OF_INTEREST = [
        'name',
        'anticodon_sequence'
    ]
    TRIMMED_TABLE_COLS_OF_INTEREST = [
        'name',
        'sequence'
    ]
    NORM_TABLE_COLS_OF_INTEREST = [
        'name',
        'mean_specific_coverage',
        'specific_coverages',
        'nonspecific_coverages'
    ]
    MOD_TABLE_COLS_OF_INTEREST = [
        'name',
        'specific_coverages',
        'nonspecific_coverages',
        'names_of_normalized_seqs_without_dels',
        'names_of_normalized_seqs_with_dels',
        'substitution_positions',
        'substitution_A_specific_coverage',
        'substitution_C_specific_coverage',
        'substitution_G_specific_coverage',
        'substitution_T_specific_coverage',
        'substitution_A_nonspecific_coverage',
        'substitution_C_nonspecific_coverage',
        'substitution_G_nonspecific_coverage',
        'substitution_T_nonspecific_coverage',
        'deletion_positions',
        'deletion_specific_coverage',
        'deletion_nonspecific_coverage',
        'consensus_sequence'
    ]

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        # Argument group A: MANDATORY
        self.trnaseq_db_paths = A('input')
        self.out_dir = A('output_dir')
        self.project_name = A('project_name')

        # Argument group B: EXTRAS
        self.num_threads = A('num_threads')
        self.seed_seq_limit = A('max_reported_trna_seeds')
        self.overwrite_out_dest = A('overwrite_output_destinations')
        self.descrip_path = os.path.abspath(A('description')) if A('description') else None

        # Argument group C: ADVANCED
        self.feature_threshold = A('feature_threshold')
        self.preferred_treatment = A('preferred_treatment')
        self.nonspecific_output = A('nonspecific_output')
        self.min_variation = A('min_variation')
        self.min_third_fourth_nt = A('min_third_fourth_nt')
        self.min_del_fraction = A('min_del_fraction')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default

        if not self.project_name:
            raise ConfigError("Please specify a name for the collection of input tRNA-seq databases "
                              "using --project-name or -n.")
        if not self.out_dir:
            raise ConfigError("Please provide an output directory using --output-dir or -o.")

        self.contigs_db_path = None
        self.contigs_db_hash = None

        self.specific_out_dir = None
        self.specific_profile_db_path = None
        self.specific_auxiliary_db_path = None

        self.nonspecific_out_dir = None
        self.nonspecific_profile_db_path = None
        self.nonspecific_auxiliary_db_path = None
        self.summed_out_dir = None
        self.summed_profile_db_path = None
        self.summed_auxiliary_db_path = None
        self.combined_out_dir = None
        self.combined_profile_db_path = None
        self.combined_auxiliary_db_path = None

        self.descrip = None

        self.preferred_trnaseq_db_sample_ids = None
        self.preferred_trnaseq_db_nums = None

        self.trnaseq_dbs_info_dict = OrderedDict()
        self.num_trnaseq_dbs = None
        self.trnaseq_db_sample_ids = None
        self.unmod_norm_seq_summaries_dict = OrderedDict()
        self.mod_seq_summaries_dict = OrderedDict()

        self.seed_seqs = None
        self.total_seed_length = None

        self.sample_total_specific_cov_dict = None
        self.sample_total_nonspecific_cov_dict = None
        self.sample_total_summed_cov_dict = None
        self.sample_overall_mean_specific_cov_dict = None
        self.sample_mean_nonspecific_cov_dict = None
        self.sample_mean_summed_cov_dict = None
        self.sample_normalization_multiplier_dict = None

        self.overall_mean_specific_cov = None
        self.overall_mean_nonspecific_cov = None

        self.variable_nts_table_entries = None
        self.specific_indels_table_entries = None
        self.nonspecific_indels_table_entries = None
        self.summed_indels_table_entries = None


    def process(self):
        """Orchestrate the steps needed to create contigs, profile and auxiliary databases."""
        self.sanity_check()
        filesnpaths.gen_output_directory(self.out_dir, delete_if_exists=self.overwrite_out_dest)

        self.load_trnaseq_dbs()

        self.form_seeds()

        filesnpaths.gen_output_directory(self.specific_out_dir)
        self.gen_contigs_db()
        self.gen_auxiliary_db('specific')

        self.set_sample_total_covs()
        self.set_sample_overall_mean_covs()
        self.set_sample_mean_covs()
        self.set_sample_std_covs()
        self.set_sample_abundances()
        self.set_sample_normalization_multipliers()
        self.set_sample_normalized_mean_Q2Q3_coverages()
        self.set_sample_detections()
        self.set_sample_relative_abundances()
        self.set_sample_max_normalized_ratios()
        self.set_variable_nts_table_entries()
        self.set_indels_table_entries()
        self.gen_profile_db('specific')

        if self.nonspecific_out_dir:
            filesnpaths.gen_output_directory(self.nonspecific_out_dir, delete_if_exists=self.overwrite_out_dest)
            self.gen_auxiliary_db('nonspecific')
            self.gen_profile_db('nonspecific')
        if self.combined_out_dir:
            filesnpaths.gen_output_directory(self.combined_out_dir, delete_if_exists=self.overwrite_out_dest)
            self.gen_auxiliary_db('combined')
            self.gen_profile_db('combined')
        if self.summed_out_dir:
            filesnpaths.gen_output_directory(self.summed_out_dir, delete_if_exists=self.overwrite_out_dest)
            self.gen_auxiliary_db('summed')
            self.gen_profile_db('summed')


    def sanity_check(self):
        """Check user inputs before proceeding."""
        for trnaseq_db_path in self.trnaseq_db_paths:
            utils.is_trnaseq_db(trnaseq_db_path)
        self.populate_trnaseq_dbs_info_dict()
        self.trnaseq_db_sample_ids = [inner_dict['sample_id'] for inner_dict in self.trnaseq_dbs_info_dict.values()]
        if len(self.trnaseq_dbs_info_dict) != len(set(self.trnaseq_db_sample_ids)):
            raise ConfigError("Sample IDs in each input tRNA-seq database must be unique. This is "
                              "not the case with your input. Here are the sample names so you can "
                              "see which ones occur more than once: '%s'" % (", ".join(self.trnaseq_db_sample_ids)))
        self.num_trnaseq_dbs = len(self.trnaseq_db_sample_ids)

        self.check_trnaseq_db_versions()

        self.out_dir = filesnpaths.check_output_directory(self.out_dir,
                                                          ok_if_exists=self.overwrite_out_dest)
        self.out_dir = os.path.abspath(self.out_dir)

        self.contigs_db_path = os.path.join(self.out_dir, 'CONTIGS.db')
        self.contigs_db_hash = 'hash' + str('%08x' % random.randrange(16**8))

        self.specific_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'SPECIFIC_COVERAGE'),
                                                                   ok_if_exists=self.overwrite_out_dest)
        self.specific_profile_db_path = os.path.join(self.specific_out_dir, 'PROFILE.db')
        self.specific_auxiliary_db_path = os.path.join(self.specific_out_dir, 'AUXILIARY-DATA.db')

        if not 1 <= self.num_threads <= mp.cpu_count():
            raise ConfigError("The number of threads to use must be a positive integer less than or equal to %d. "
                              "Try again!" % mp.cpu_count())

        self.set_treatment_preference()

        self.set_nonspecific_db_info()

        utils.check_sample_id(self.project_name)

        if self.descrip_path:
            filesnpaths.is_file_plain_text(self.descrip_path)
            self.descrip_path = os.path.abspath(self.descrip_path)
            self.descrip = open(self.descrip_path).read()

        if self.seed_seq_limit == -1:
            self.seed_seq_limit = sys.maxsize
        elif self.seed_seq_limit < 1:
            raise ConfigError(f"{self.seed_seq_limit} is an invalid value for `--max-reported-seed-seqs`. "
                              "To remove the limit on tRNA seeds reported to the contigs database, "
                              "provide a value of -1. Otherwise provide an integer greater than 0.")

        self.run.info("Input tRNA-seq databases", ", ".join(self.trnaseq_db_paths))
        if self.preferred_treatment:
            self.run.info("Databases preferred for seed formation",
                          ", ".join([trnaseq_db_path for trnaseq_db_num, trnaseq_db_path
                                     in enumerate(self.trnaseq_db_paths)
                                     if trnaseq_db_num in self.preferred_trnaseq_db_nums]))
        self.run.info("Output directory", self.out_dir)


    def populate_trnaseq_dbs_info_dict(self):
        """Get the meta-data from the input tRNA-seq databases."""
        for trnaseq_db_path in self.trnaseq_db_paths:
            trnaseq_db = dbops.TRNASeqDatabase(trnaseq_db_path)
            self.trnaseq_dbs_info_dict[trnaseq_db_path] = trnaseq_db.meta


    def check_trnaseq_db_versions(self):
        if len(set([inner_dict['version'] for inner_dict in self.trnaseq_dbs_info_dict.values()])) > 1:
            trnaseq_db_version_report = "\n".join([trnaseq_db_path + " : " + inner_dict['version']
                                                   for trnaseq_db_path, inner_dict in self.trnaseq_dbs_info_dict.items()])
            if anvio.FORCE:
                self.run.warning("Not all input tRNA-seq databases have the same version number, "
                                 "but since you have used the `--force` flag, `anvi-convert-trnaseq-database` "
                                 "will proceed though this is dangerous and may lead to errors. "
                                 f"Here is the version number of each database:\n{trnaseq_db_version_report}")
            else:
                raise ConfigError("Not all input tRNA-seq databases have the same version number. "
                                  f"Here is the version number of each database:\n{trnaseq_db_version_report}")


    def set_treatment_preference(self):
        if not self.preferred_treatment:
            return

        input_treatments = [inner_dict['treatment'] for inner_dict in self.trnaseq_dbs_info_dict.values()]
        self.preferred_trnaseq_db_sample_ids = []
        self.preferred_trnaseq_db_nums = []
        if self.preferred_treatment not in input_treatments:
            raise ConfigError("You provided a preferred treatment type, %s, "
                              "but it was not found in any of the input databases, "
                              "which were found to have the following treatments: %s."
                              % (self.preferred_treatment, ', '.join(input_treatments)))
        for trnaseq_db_num, treatment in enumerate(input_treatments):
            if self.preferred_treatment == treatment:
                self.preferred_trnaseq_db_sample_ids.append(self.trnaseq_db_sample_ids[trnaseq_db_num])
                self.preferred_trnaseq_db_nums.append(trnaseq_db_num)


    def set_nonspecific_db_info(self):
        self.nonspecific_db_types = self.nonspecific_output.split(',')
        for nonspecific_db_type in self.nonspecific_db_types:
            if nonspecific_db_type not in ['nonspecific_db', 'combined_db', 'summed_db']:
                raise ConfigError("The nonspecific profile database types provided by `--nonspecific-output` are not recognized. "
                                  "The database types must be comma separated without spaces, "
                                  f"e.g., 'nonspecific_db,combined_db,summed_db'. Your argument was: {self.nonspecific_output}'")

        if 'nonspecific_db' in self.nonspecific_db_types:
            self.nonspecific_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'NONSPECIFIC_COVERAGE'),
                                                                          ok_if_exists=self.overwrite_out_dest)
            self.nonspecific_profile_db_path = os.path.join(self.nonspecific_out_dir, 'PROFILE.db')
            self.nonspecific_auxiliary_db_path = os.path.join(self.nonspecific_out_dir, 'AUXILIARY-DATA.db')

        if 'combined_db' in self.nonspecific_db_types:
            self.combined_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'COMBINED_COVERAGE'),
                                                                       ok_if_exists=self.overwrite_out_dest)
            self.combined_profile_db_path = os.path.join(self.combined_out_dir, 'PROFILE.db')
            self.combined_auxiliary_db_path = os.path.join(self.combined_out_dir, 'AUXILIARY-DATA.db')

        if 'summed_db' in self.nonspecific_db_types:
            self.summed_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'SUMMED_COVERAGE'),
                                                                     ok_if_exists=self.overwrite_out_dest)
            self.summed_profile_db_path = os.path.join(self.summed_out_dir, 'PROFILE.db')
            self.summed_auxiliary_db_path = os.path.join(self.summed_out_dir, 'AUXILIARY-DATA.db')


    def load_trnaseq_dbs(self):
        loaded_db_count = 0
        num_trnaseq_db_paths = len(self.trnaseq_db_paths)
        self.progress.new("Loading sequence information from tRNA-seq databases")
        self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} databases loaded")

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        processes = [mp.Process(target=trnaseq_db_loader, args=(input_queue, output_queue, self))
                     for _ in range(self.num_threads)]
        for p in processes:
            p.start()

        for trnaseq_db_path in self.trnaseq_db_paths:
            input_queue.put(trnaseq_db_path)

        while loaded_db_count < len(self.trnaseq_db_paths):
            trnaseq_db_path, unmod_norm_seq_summaries, mod_seq_summaries = output_queue.get()
            self.unmod_norm_seq_summaries_dict[trnaseq_db_path] = unmod_norm_seq_summaries
            self.mod_seq_summaries_dict[trnaseq_db_path] = mod_seq_summaries
            loaded_db_count += 1
            self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} databases loaded")

        for p in processes:
            p.terminate()
            p.join()

        self.progress.end()


    def load_trnaseq_db_seq_info(self, trnaseq_db_path):
        """Load necessary tRNA sequence data from the input tRNA-seq database.

        Unmodified normalized sequences and "modified" sequences, comprising clustered normalized
        sequences, are stored in distinct data structures.
        """
        trnaseq_db_num = list(self.trnaseq_dbs_info_dict.keys()).index(trnaseq_db_path)
        sample_id = self.trnaseq_db_sample_ids[trnaseq_db_num]

        trnaseq_db = dbops.TRNASeqDatabase(trnaseq_db_path)

        norm_seq_summary_dict = {} # Used to link normalized to modified sequence summaries

        # Store normalized sequence strings and feature information.
        seq_string_and_feature_df = pd.DataFrame(
            trnaseq_db.db.get_some_columns_from_table('feature', ', '.join(self.FEATURE_TABLE_COLS_OF_INTEREST + [self.feature_threshold + '_start'])),
            columns=self.FEATURE_TABLE_COLS_OF_INTEREST + [self.feature_threshold + '_start']
        ).set_index('name')
        if 'stem' in self.feature_threshold:
            # The starts of both strands of the stem are recorded, so pick the start of the 5' strand.
            seq_string_and_feature_df[self.feature_threshold + '_start'] = [int(entry.split(',')[0]) if isinstance(entry, str) else entry for entry
                                                                            in seq_string_and_feature_df[self.feature_threshold + '_start'].fillna(-1).tolist()]
        else:
            seq_string_and_feature_df[self.feature_threshold + '_start'] = seq_string_and_feature_df[self.feature_threshold + '_start'].fillna(-1)
        seq_string_and_feature_df['anticodon_sequence'] = seq_string_and_feature_df['anticodon_sequence'].fillna('')
        seq_string_and_feature_df = pd.merge(
            pd.DataFrame(
                trnaseq_db.db.get_some_columns_from_table('trimmed', ', '.join(self.TRIMMED_TABLE_COLS_OF_INTEREST)),
                columns=self.TRIMMED_TABLE_COLS_OF_INTEREST
            ).set_index('name'),
            seq_string_and_feature_df,
            left_index=True,
            right_index=True)

        for norm_seq_info in trnaseq_db.db.get_some_columns_from_table('normalized', ', '.join(self.NORM_TABLE_COLS_OF_INTEREST)):
            (name,
             mean_specific_cov,
             specific_covs_string,
             nonspecific_covs_string) = norm_seq_info

            norm_seq_summary = NormalizedSeqSummary()
            norm_seq_summary.name = name
            norm_seq_summary.sample_id = sample_id
            norm_seq_summary.mean_specific_cov = mean_specific_cov
            # There is always a trailing comma in the coverage strings.
            norm_seq_summary.specific_covs = np.fromiter(map(int, specific_covs_string.split(',')[: -1]), int)
            norm_seq_summary.nonspecific_covs = np.fromiter(map(int, nonspecific_covs_string.split(',')[: -1]), int)
            try:
                (norm_seq_summary.seq_string,
                 norm_seq_summary.anticodon_seq_string,
                 norm_seq_summary.feature_threshold_start) = seq_string_and_feature_df.loc[norm_seq_summary.name, ['sequence', 'anticodon_sequence', self.feature_threshold + '_start']]
            except KeyError:
                # Normalized sequences with truncated profiles that were discovered and recovered
                # through anvi-trnaseq deletion analysis do not have entries in the Feature table,
                # and so are ignored.
                continue

            norm_seq_summary_dict[norm_seq_summary.name] = norm_seq_summary

        mod_seq_summaries = []
        for mod_seq_info in trnaseq_db.db.get_some_columns_from_table('modified', ', '.join(self.MOD_TABLE_COLS_OF_INTEREST)):
            (name,
             specific_covs,
             nonspecific_covs,
             names_of_norm_seqs_without_dels,
             names_of_norm_seqs_with_dels,
             sub_positions,
             sub_A_specific_covs,
             sub_C_specific_covs,
             sub_G_specific_covs,
             sub_T_specific_covs,
             sub_A_nonspecific_covs,
             sub_C_nonspecific_covs,
             sub_G_nonspecific_covs,
             sub_T_nonspecific_covs,
             del_positions,
             del_specific_covs,
             del_nonspecific_covs,
             consensus_seq_string) = mod_seq_info

            mod_seq_summary = ModifiedSeqSummary()
            mod_seq_summary.name = name
            mod_seq_summary.sample_id = sample_id
            # There is always a trailing comma in the substitution and deletion coverage and
            # position strings.
            mod_seq_summary.sub_positions = np.fromiter(map(int, sub_positions.split(',')[: -1]), int)
            mod_seq_summary.consensus_seq_string = consensus_seq_string

            # Make nucleotide variability arrays covering every position in the sequence. Start with
            # arrays of overall specific/nonspecific coverage with nonzero values for the
            # nucleotides found in the consensus sequence, and then correct the variable positions.
            seq_length = len(consensus_seq_string)
            specific_nt_covs_dict = {nt: np.zeros(seq_length, int) for nt in UNAMBIG_NTS}
            nonspecific_nt_covs_dict = {nt: np.zeros(seq_length, int) for nt in UNAMBIG_NTS}
            pos = 0
            for nt, specific_cov, nonspecific_cov in zip(consensus_seq_string,
                                                         specific_covs.split(',')[: -1],
                                                         nonspecific_covs.split(',')[: -1]):
                specific_nt_covs_dict[nt][pos] = specific_cov
                nonspecific_nt_covs_dict[nt][pos] = nonspecific_cov
                pos += 1
            for (sub_pos,
                 specific_A_cov,
                 specific_C_cov,
                 specific_G_cov,
                 specific_T_cov,
                 nonspecific_A_cov,
                 nonspecific_C_cov,
                 nonspecific_G_cov,
                 nonspecific_T_cov) in zip(map(int, sub_positions.split(',')[: -1]),
                                           map(int, sub_A_specific_covs.split(',')[: -1]),
                                           map(int, sub_C_specific_covs.split(',')[: -1]),
                                           map(int, sub_G_specific_covs.split(',')[: -1]),
                                           map(int, sub_T_specific_covs.split(',')[: -1]),
                                           map(int, sub_A_nonspecific_covs.split(',')[: -1]),
                                           map(int, sub_C_nonspecific_covs.split(',')[: -1]),
                                           map(int, sub_G_nonspecific_covs.split(',')[: -1]),
                                           map(int, sub_T_nonspecific_covs.split(',')[: -1])):
                     specific_nt_covs_dict['A'][sub_pos] = specific_A_cov
                     specific_nt_covs_dict['C'][sub_pos] = specific_C_cov
                     specific_nt_covs_dict['G'][sub_pos] = specific_G_cov
                     specific_nt_covs_dict['T'][sub_pos] = specific_T_cov
                     nonspecific_nt_covs_dict['A'][sub_pos] = nonspecific_A_cov
                     nonspecific_nt_covs_dict['C'][sub_pos] = nonspecific_C_cov
                     nonspecific_nt_covs_dict['G'][sub_pos] = nonspecific_G_cov
                     nonspecific_nt_covs_dict['T'][sub_pos] = nonspecific_T_cov
            mod_seq_summary.specific_nt_covs_dict = specific_nt_covs_dict
            mod_seq_summary.nonspecific_nt_covs_dict = nonspecific_nt_covs_dict

            specific_del_covs = np.zeros(seq_length, int)
            nonspecific_del_covs = np.zeros(seq_length, int)
            if del_positions != ',':
                for (del_pos,
                     specific_del_cov,
                     nonspecific_del_cov) in zip(map(int, del_positions.replace(';', ',').split(',')[: -1]),
                                                 map(int, del_specific_covs.split(',')[: -1]),
                                                 map(int, del_nonspecific_covs.split(',')[: -1])):
                    specific_del_covs[del_pos] = specific_del_cov
                    nonspecific_del_covs[del_pos] = nonspecific_del_cov
            mod_seq_summary.specific_del_covs = specific_del_covs
            mod_seq_summary.nonspecific_del_covs = nonspecific_del_covs

            norm_seq_names = names_of_norm_seqs_without_dels.split(',')
            if names_of_norm_seqs_with_dels:
                norm_seq_names.extend(names_of_norm_seqs_with_dels.split(','))
            mod_seq_summary.norm_seq_summaries = []
            for norm_seq_name in norm_seq_names:
                try:
                    norm_seq_summary = norm_seq_summary_dict[norm_seq_name]
                except KeyError:
                    # Normalized sequences with truncated profiles that were discovered and recovered
                    # through anvi-trnaseq deletion analysis are ignored.
                    continue

                # Cross-reference the modified sequence summary and constituent modified normalized
                # sequence summary objects.
                norm_seq_summary.mod_seq_summary = mod_seq_summary
                mod_seq_summary.norm_seq_summaries.append(norm_seq_summary)

                # Ensure that all of the constituent modified normalized sequences have coverage
                # arrays flush with the modified sequence.
                if len(norm_seq_summary.seq_string) < len(mod_seq_summary.consensus_seq_string):
                    fiveprime_extension = np.zeros(len(mod_seq_summary.consensus_seq_string) - len(norm_seq_summary.seq_string), int)
                    self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
            mod_seq_summaries.append(mod_seq_summary)

        unmod_norm_seq_summaries = [norm_seq for norm_seq in norm_seq_summary_dict.values()
                                    if not norm_seq.mod_seq_summary]

        return unmod_norm_seq_summaries, mod_seq_summaries


    def extend_norm_seq_fiveprime_end(self, norm_seq_summary, fiveprime_extension):
        """Seed sequences can be longer than the normalized sequences from the individual samples,
        requiring addition of empty positions in the normalized sequence coverage arrays at the 5'
        end, as normalized (and modified) sequences are aligned from the 3' end."""
        norm_seq_summary.specific_covs = np.concatenate([fiveprime_extension, norm_seq_summary.specific_covs])
        norm_seq_summary.nonspecific_covs = np.concatenate([fiveprime_extension, norm_seq_summary.nonspecific_covs])


    def form_seeds(self):
        """Form tRNA seed sequences through comparison of the input samples.

        Seed sequences are formed through comparison of the samples' normalized sequences (both
        unmodified normalized sequences and normalized sequences underlying modified sequences).
        Seed sequences need not be found in every sample.

        Modification-induced mutations complicate seed formation. (anvi-trnaseq is capable of
        finding substitutions -- the main type of mutation -- and deletions, though insertions can
        also occur but go undetected.) Modified sequences are derived from clusters of normalized
        sequences; unmodified normalized sequences and the normalized sequences underlying modified
        sequences are here compared between samples to find seeds. If a normalized sequence is
        shared identically (not as a subsequence) between samples, then the normalized seed
        sequences and any modified sequences that the normalized sequences are part of are combined
        into a single seed sequence.

        It is a heuristic to exactly match normalized sequences, rather than to check whether one
        normalized sequence is a 3' subsequence of the other, or to instead compare underlying
        trimmed sequences making up normalized sequences. This heuristic should not distort sample
        merging for the more abundant tRNA species, in particular, as these are most likely to be
        represented by reads spanning the full length of the tRNA, producing the same normalized
        sequences."""
        self.progress.new("Forming seed sequences from input samples")

        norm_seq_string_seed_seq_dict = {}
        for trnaseq_db_num, trnaseq_db_path in enumerate(self.trnaseq_db_paths):
            sample_id = self.trnaseq_db_sample_ids[trnaseq_db_num]
            self.progress.update(f"Adding {sample_id}")

            unmod_norm_seq_summaries = self.unmod_norm_seq_summaries_dict[trnaseq_db_path]
            mod_seq_summaries = self.mod_seq_summaries_dict[trnaseq_db_path]

            for norm_seq_summary in unmod_norm_seq_summaries:
                # Process normalized sequences without any detected potential modifications.
                norm_seq_string = norm_seq_summary.seq_string
                try:
                    # The normalized sequence has already been found in another dataset.
                    seed_seq = norm_seq_string_seed_seq_dict[norm_seq_string]
                except KeyError:
                    # Create a new seed sequence based on the normalized sequence.
                    seed_seq = SeedSeq()
                    seed_seq.name = norm_seq_summary.name + '_' + sample_id
                    seed_seq.seq_string = norm_seq_string
                    seed_seq.meets_feature_threshold = True if norm_seq_summary.feature_threshold_start >= 0 else False
                    seed_seq.unmod_norm_seq_summaries = [norm_seq_summary]
                    seed_seq.mod_seq_summaries = []
                    norm_seq_string_seed_seq_dict[norm_seq_string] = seed_seq
                    continue

                if len(norm_seq_string) < len(seed_seq.seq_string):
                    # The normalized sequence string is shorter than the seed sequence string,
                    # implying that the seed string is a longer modified sequence in another
                    # dataset. Extend the normalized sequence coverage arrays the needed amount at
                    # the 5' end. (Note: It is impossible here for the normalized sequence string to
                    # be longer than the seed sequence string.)
                    fiveprime_extension = np.zeros(len(seed_seq.seq_string) - len(norm_seq_string), int)
                    self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
                seed_seq.unmod_norm_seq_summaries.append(norm_seq_summary)

            for mod_seq_summary in mod_seq_summaries:
                # Find seed sequences from other datasets containing any of the normalized sequences
                # forming the modified sequence under consideration. If more than one seed sequence
                # is identified, they are merged.
                seed_seq_dict = {}
                for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                    try:
                        # The normalized sequence is represented in another dataset.
                        seed_seq = norm_seq_string_seed_seq_dict[norm_seq_summary.seq_string]
                    except KeyError:
                        continue
                    seed_seq_dict[seed_seq.name] = seed_seq

                if not seed_seq_dict:
                    # Create a new seed sequence based on the modified sequence.
                    seed_seq = SeedSeq()
                    seed_seq.name = mod_seq_summary.name + '_' + sample_id
                    seed_seq.seq_string = mod_seq_summary.consensus_seq_string
                    seed_seq.unmod_norm_seq_summaries = []
                    seed_seq.mod_seq_summaries = []
                    seed_seq.mod_seq_summaries.append(mod_seq_summary)
                    for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                        if norm_seq_summary.feature_threshold_start >= 0:
                            seed_seq.meets_feature_threshold = True
                            break
                    else:
                        seed_seq.meets_feature_threshold = False
                    for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                        norm_seq_string_seed_seq_dict[norm_seq_summary.seq_string] = seed_seq
                    continue

                if len(seed_seq_dict) == 1:
                    # The modified sequence shares one or more normalized sequences with one seed sequence.
                    seed_seq_name, seed_seq = seed_seq_dict.popitem()
                    if len(mod_seq_summary.consensus_seq_string) < len(seed_seq.seq_string):
                        # The modified sequence is shorter than the seed sequence, so its coverage
                        # arrays must be extended with zeros at the 5' end.
                        fiveprime_extension = np.zeros(len(seed_seq.seq_string) - len(mod_seq_summary.consensus_seq_string), int)
                        self.extend_mod_seq_fiveprime_end(mod_seq_summary, fiveprime_extension)
                        for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                            self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
                    elif len(mod_seq_summary.consensus_seq_string) > len(seed_seq.seq_string):
                        # The modified sequence is longer than the seed sequence, so the coverage
                        # arrays of the sequences forming the seed sequence must be extended with
                        # zeros at the 5' end.
                        fiveprime_extension = np.zeros(len(mod_seq_summary.consensus_seq_string) - len(seed_seq.seq_string), int)
                        for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                            self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
                        for other_mod_seq_summary in seed_seq.mod_seq_summaries:
                            self.extend_mod_seq_fiveprime_end(other_mod_seq_summary, fiveprime_extension)
                        seed_seq.name = mod_seq_summary.name + '_' + sample_id
                        seed_seq.seq_string = mod_seq_summary.consensus_seq_string
                        for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                            if norm_seq_summary.feature_threshold_start >= 0:
                                seed_seq.meets_feature_threshold = True
                                break
                        else:
                            seed_seq.meets_feature_threshold = False
                    seed_seq.mod_seq_summaries.append(mod_seq_summary)
                    for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                        norm_seq_string_seed_seq_dict[norm_seq_summary.seq_string] = seed_seq
                    continue

                # To reach this point, the modified sequence must map to more than one seed
                # sequence.
                sorted_seed_seqs = sorted([seed_seq for seed_seq in seed_seq_dict.values()],
                                          key=lambda seed_seq: -len(seed_seq.seq_string))
                max_seed_seq_length = len(sorted_seed_seqs[0].seq_string)
                mod_seq_length = len(mod_seq_summary.consensus_seq_string)
                if mod_seq_length < max_seed_seq_length:
                    # Extend coverage arrays of the modified sequence.
                    fiveprime_extension = np.zeros(max_seed_seq_length - mod_seq_length, int)
                    self.extend_mod_seq_fiveprime_end(mod_seq_summary, fiveprime_extension)

                    # Extend coverage arrays of shorter seed sequences now grouped with a longer
                    # seed sequence.
                    for seed_seq in seed_seq_dict.values():
                        if len(seed_seq.seq_string) < max_seed_seq_length:
                            fiveprime_extension = np.zeros(max_seed_seq_length - len(seed_seq.seq_string), int)
                            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                                self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
                            for other_mod_seq_summary in seed_seq.mod_seq_summaries:
                                self.extend_mod_seq_fiveprime_end(other_mod_seq_summary, fiveprime_extension)

                    new_seed_seq = SeedSeq()
                    longest_seed_seq = sorted_seed_seqs[0]
                    new_seed_seq.name = longest_seed_seq.name
                    new_seed_seq.seq_string = longest_seed_seq.seq_string
                    new_seed_seq.meets_feature_threshold = longest_seed_seq.meets_feature_threshold
                elif mod_seq_length > max_seed_seq_length:
                    # Extend coverage arrays of seed sequences.
                    for seed_seq in seed_seq_dict.values():
                        fiveprime_extension = np.zeros(mod_seq_length - len(seed_seq.seq_string), int)
                        for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                            self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
                        for other_mod_seq_summary in seed_seq.mod_seq_summaries:
                            self.extend_mod_seq_fiveprime_end(other_mod_seq_summary, fiveprime_extension)

                    new_seed_seq = SeedSeq()
                    new_seed_seq.name = mod_seq_summary.name + '_' + sample_id
                    new_seed_seq.seq_string = mod_seq_summary.consensus_seq_string
                    for seed_seq in seed_seq_dict.values():
                        if seed_seq.meets_feature_threshold:
                            new_seed_seq.meets_feature_threshold = True
                            break
                    else:
                        for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                            if norm_seq_summary.feature_threshold_start >= 0:
                                new_seed_seq.meets_feature_threshold = True
                                break
                        else:
                            new_seed_seq.meets_feature_threshold = False
                else:
                    # The modified sequence is the same length as the longest seed sequence. Extend
                    # coverage arrays of shorter seed sequences now grouped with a longer seed
                    # sequence.
                    for seed_seq in seed_seq_dict.values():
                        if len(seed_seq.seq_string) < max_seed_seq_length:
                            fiveprime_extension = np.zeros(max_seed_seq_length - len(seed_seq.seq_string), int)
                            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                                self.extend_norm_seq_fiveprime_end(norm_seq_summary, fiveprime_extension)
                            for other_mod_seq_summary in seed_seq.mod_seq_summaries:
                                self.extend_mod_seq_fiveprime_end(other_mod_seq_summary, fiveprime_extension)

                    new_seed_seq = SeedSeq()
                    new_seed_seq.name = mod_seq_summary.name + '_' + sample_id
                    new_seed_seq.seq_string = mod_seq_summary.consensus_seq_string
                    if sorted_seed_seqs[0].meets_feature_threshold:
                        new_seed_seq.meets_feature_threshold = True
                    else:
                        for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                            if norm_seq_summary.feature_threshold_start >= 0:
                                new_seed_seq.meets_feature_threshold = True
                                break
                        else:
                            new_seed_seq.meets_feature_threshold = False

                # Now that all of the coverage arrays are reconciled in length, the modified
                # sequence query and constituent sequences of the matching seeds can be added to the
                # new seed.
                new_seed_seq.unmod_norm_seq_summaries = []
                new_seed_seq.mod_seq_summaries = []
                for seed_seq in seed_seq_dict.values():
                    new_seed_seq.unmod_norm_seq_summaries += seed_seq.unmod_norm_seq_summaries
                    new_seed_seq.mod_seq_summaries += seed_seq.mod_seq_summaries
                new_seed_seq.mod_seq_summaries.append(mod_seq_summary)

                for norm_seq_summary in new_seed_seq.unmod_norm_seq_summaries:
                    norm_seq_string_seed_seq_dict[norm_seq_summary.seq_string] = new_seed_seq
                for mod_seq_summary in new_seed_seq.mod_seq_summaries:
                    for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                        norm_seq_string_seed_seq_dict[norm_seq_summary.seq_string] = new_seed_seq

        self.progress.update("Finalizing seeds")

        # The seed references in the dict need to be dereplicated.
        seed_seqs = list({seed_seq.name: seed_seq for seed_seq in norm_seq_string_seed_seq_dict.values()}.values())

        # Disregard seed sequences that do not reach the 5' feature threshold.
        seed_seqs = [seed_seq for seed_seq in seed_seqs if seed_seq.meets_feature_threshold]
        self.set_anticodon(seed_seqs)
        seed_seqs = [seed_seq for seed_seq in seed_seqs if seed_seq.anticodon_seq_string]

        # Find specific coverages of modified sequences from specific coverages of nucleotides.
        self.sum_specific_nt_covs(seed_seqs)
        # Select the top seeds by specific coverage.
        self.set_total_specific_covs(seed_seqs)
        seed_seqs = sorted(seed_seqs, key=lambda seed_seq: -seed_seq.total_mean_specific_cov)[: self.seed_seq_limit]

        self.set_nt_covs(seed_seqs)
        self.sum_nonspecific_nt_covs(seed_seqs)
        self.set_total_nonspecific_covs(seed_seqs)
        self.set_consensus_seq_string(seed_seqs)
        self.set_gc_fraction(seed_seqs)

        self.seed_seqs = seed_seqs
        self.total_seed_length = sum([len(seed_seq.seq_string) for seed_seq in seed_seqs])

        self.set_sample_covs()
        self.set_mods()
        self.set_consensus_mod_nts()

        self.set_sample_dels()

        self.progress.end()


    def extend_mod_seq_fiveprime_end(self, mod_seq_summary, fiveprime_extension):
        """Seed sequences can be longer than the normalized sequences from the individual samples,
        requiring addition of empty positions in the normalized sequence coverage arrays at the 5'
        end, as normalized (and modified) sequences are aligned from the 3' end."""
        # The positions of substitutions are recorded in the seed sequence index.
        mod_seq_summary.sub_positions += fiveprime_extension.size
        for nt in UNAMBIG_NTS:
            mod_seq_summary.specific_nt_covs_dict[nt] = np.concatenate([fiveprime_extension, mod_seq_summary.specific_nt_covs_dict[nt]])
            mod_seq_summary.nonspecific_nt_covs_dict[nt] = np.concatenate([fiveprime_extension, mod_seq_summary.nonspecific_nt_covs_dict[nt]])
        mod_seq_summary.specific_del_covs = np.concatenate([fiveprime_extension, mod_seq_summary.specific_del_covs])
        mod_seq_summary.nonspecific_del_covs = np.concatenate([fiveprime_extension, mod_seq_summary.nonspecific_del_covs])


    def set_anticodon(self, seed_seqs):
        """Assign the anticodon by comparing the mean specific coverage of the normalized sequences
        comprising the seed (a simpler, approximate substitute for extracting the anticodon coverage
        from each normalized sequence). This method is called by `DatabaseConverter.process` after
        sequences from input tRNA-seq databases have been assigned to seeds."""
        for seed_seq in seed_seqs:
            if not seed_seq.mod_seq_summaries:
                # The seed is comprised entirely of identical unmodified normalized sequences.
                seed_seq.anticodon_seq_string = seed_seq.unmod_norm_seq_summaries[0].anticodon_seq_string
                continue

            anticodon_cov_dict = defaultdict(int)
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                # Mean specific coverage does not include any 5' padding of zero coverage from the
                # formation of the seed sequence.
                anticodon_cov_dict[norm_seq_summary.anticodon_seq_string] += norm_seq_summary.mean_specific_cov
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                for norm_seq_summary in mod_seq_summary.norm_seq_summaries:
                    anticodon_cov_dict[norm_seq_summary.anticodon_seq_string] += norm_seq_summary.mean_specific_cov

            seed_seq.anticodon_seq_string = sorted(anticodon_cov_dict.items(), key=lambda t: -t[1])[0][0]


    def sum_specific_nt_covs(self, seed_seqs):
        """Coverages of the 4 nucleotides are maintained in the tRNA-seq database for modified
        sequences, since modified sequences have nucleotide variability, so this method is used to
        calculate overall specific coverages."""
        for seed_seq in seed_seqs:
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                mod_seq_summary.specific_covs = np.array([covs for covs in mod_seq_summary.specific_nt_covs_dict.values()]).sum(axis=0)


    def set_total_specific_covs(self, seed_seqs):
        """Sum specific coverages from each sequence comprising the seed sequence, and also
        calculate the mean thereof."""
        for seed_seq in seed_seqs:
            seed_seq.total_specific_covs = np.zeros(len(seed_seq.seq_string), int)
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                seed_seq.total_specific_covs += norm_seq_summary.specific_covs
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                seed_seq.total_specific_covs += mod_seq_summary.specific_covs

            seed_seq.total_mean_specific_cov = seed_seq.total_specific_covs.mean()


    def set_nt_covs(self, seed_seqs):
        """Make separate coverage arrays for A, C, G and T, which are needed in seed sequence formation."""
        for seed_seq in seed_seqs:
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                norm_seq_summary.specific_nt_covs_dict = {nt: np.zeros(norm_seq_summary.specific_covs.size, int) for nt in UNAMBIG_NTS}
                norm_seq_summary.nonspecific_nt_covs_dict = {nt: np.zeros(norm_seq_summary.nonspecific_covs.size, int) for nt in UNAMBIG_NTS}
                start_pos = norm_seq_summary.specific_covs.size - len(norm_seq_summary.seq_string)
                pos = start_pos
                for nt, specific_cov, nonspecific_cov in zip(norm_seq_summary.seq_string,
                                                             norm_seq_summary.specific_covs[start_pos: ],
                                                             norm_seq_summary.nonspecific_covs[start_pos: ]):
                    norm_seq_summary.specific_nt_covs_dict[nt][pos] = specific_cov
                    norm_seq_summary.nonspecific_nt_covs_dict[nt][pos] = nonspecific_cov
                    pos += 1


    def sum_nonspecific_nt_covs(self, seed_seqs):
        """Coverages of the 4 nucleotides are maintained in the tRNA-seq database for modified
        sequences, since modified sequences have nucleotide variability, so this method is used to
        calculate overall nonspecific coverages."""
        for seed_seq in seed_seqs:
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                mod_seq_summary.nonspecific_covs = np.array([covs for covs in mod_seq_summary.nonspecific_nt_covs_dict.values()]).sum(axis=0)


    def set_total_nonspecific_covs(self, seed_seqs):
        """Sum nonspecific coverages from each sequence forming the seed sequence, and also
        calculate the mean thereof."""
        for seed_seq in seed_seqs:
            seed_seq.total_nonspecific_covs = np.zeros(len(seed_seq.seq_string), int)
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                seed_seq.total_nonspecific_covs += norm_seq_summary.nonspecific_covs
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                nonspecific_nt_covs = np.array([covs for covs in mod_seq_summary.nonspecific_nt_covs_dict.values()])
                seed_seq.total_nonspecific_covs += nonspecific_nt_covs.sum(axis=0)

            seed_seq.total_mean_nonspecific_cov = seed_seq.total_nonspecific_covs.mean()


    def set_consensus_seq_string(self, seed_seqs):
        """The consensus sequence for the seed consists of the nucleotides with the maximum specific
        coverage summed across constituent sequences. When certain tRNA-seq treatments are preferred
        (e.g., demethylase), nucleotides with predicted modifications are called on the basis of
        sequences from preferred samples."""
        for seed_seq in seed_seqs:
            if not seed_seq.mod_seq_summaries:
                # The seed is composed entirely of identical unmodified normalized sequences.
                return

            total_nt_cov_dict = {nt: np.zeros(len(seed_seq.seq_string), int) for nt in UNAMBIG_NTS}
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                for nt in UNAMBIG_NTS:
                    total_nt_cov_dict[nt] += norm_seq_summary.specific_nt_covs_dict[nt]
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                for nt in UNAMBIG_NTS:
                    total_nt_cov_dict[nt] += mod_seq_summary.specific_nt_covs_dict[nt]

            seed_seq.seq_string = ''.join(
                [INT_NT_DICT[i + 1] for i in
                np.argmax(np.array([total_nt_cov_dict[nt] for nt in UNAMBIG_NTS]), axis=0)]
            )


    def set_gc_fraction(self, seed_seqs):
        for seed_seq in seed_seqs:
            seed_seq.gc_fraction = sum([1 for nt in seed_seq.seq_string if nt == 'C' or nt == 'G']) / len(seed_seq.seq_string)


    def set_sample_covs(self):
        """Determine sample-specific coverages of seeds. Specific, nonspecific and summed coverages
        are found for A, C, G and T, as well as overall and for Q2-Q3."""
        for seed_seq in self.seed_seqs:
            sample_specific_covs_dict = {}
            sample_nonspecific_covs_dict = {}
            sample_summed_covs_dict = {}
            sample_specific_nt_covs_dict = {}
            sample_nonspecific_nt_covs_dict = {}
            sample_summed_nt_covs_dict = {}
            seed_seq_length = len(seed_seq.seq_string)

            for sample_id in self.trnaseq_db_sample_ids:
                sample_specific_covs_dict[sample_id] = np.zeros(seed_seq_length, int)
                sample_nonspecific_covs_dict[sample_id] = np.zeros(seed_seq_length, int)

                sample_specific_nt_covs_dict[sample_id] = [np.zeros(seed_seq_length, int) for _ in UNAMBIG_NTS]
                sample_nonspecific_nt_covs_dict[sample_id] = [np.zeros(seed_seq_length, int) for _ in UNAMBIG_NTS]

            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                sample_id = norm_seq_summary.sample_id

                sample_specific_covs_dict[sample_id] += norm_seq_summary.specific_covs
                sample_nonspecific_covs_dict[sample_id] += norm_seq_summary.nonspecific_covs

                for i, nt in enumerate(UNAMBIG_NTS):
                    sample_specific_nt_covs_dict[sample_id][i] += norm_seq_summary.specific_nt_covs_dict[nt]
                    sample_nonspecific_nt_covs_dict[sample_id][i] += norm_seq_summary.nonspecific_nt_covs_dict[nt]

            for mod_seq_summary in seed_seq.mod_seq_summaries:
                sample_id = mod_seq_summary.sample_id

                sample_specific_covs_dict[sample_id] += mod_seq_summary.specific_covs
                sample_nonspecific_covs_dict[sample_id] += mod_seq_summary.nonspecific_covs

                for i, nt in enumerate(UNAMBIG_NTS):
                    sample_specific_nt_covs_dict[sample_id][i] += mod_seq_summary.specific_nt_covs_dict[nt]
                    sample_nonspecific_nt_covs_dict[sample_id][i] += mod_seq_summary.nonspecific_nt_covs_dict[nt]

            for sample_id in self.trnaseq_db_sample_ids:
                sample_summed_covs_dict[sample_id] = sample_specific_covs_dict[sample_id] + sample_nonspecific_covs_dict[sample_id]

                sample_specific_nt_covs = sample_specific_nt_covs_dict[sample_id]
                sample_nonspecific_nt_covs = sample_nonspecific_nt_covs_dict[sample_id]
                sample_summed_nt_covs_dict[sample_id] = [specific_covs + nonspecific_covs for specific_covs, nonspecific_covs
                                                         in zip(sample_specific_nt_covs, sample_nonspecific_nt_covs)]

            seed_seq.sample_specific_covs_dict = sample_specific_covs_dict
            seed_seq.sample_nonspecific_covs_dict = sample_nonspecific_covs_dict
            seed_seq.sample_summed_covs_dict = sample_summed_covs_dict
            seed_seq.sample_specific_nt_covs_dict = sample_specific_nt_covs_dict
            seed_seq.sample_nonspecific_nt_covs_dict = sample_nonspecific_nt_covs_dict
            seed_seq.sample_summed_nt_covs_dict = sample_summed_nt_covs_dict

            q = int(seed_seq_length * 0.25)
            seed_seq.sample_mean_Q2Q3_specific_cov_dict = {}
            seed_seq.sample_mean_Q2Q3_nonspecific_cov_dict = {}
            seed_seq.sample_mean_Q2Q3_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                seed_seq.sample_mean_Q2Q3_specific_cov_dict[sample_id] = np.mean(
                    sorted(sample_specific_covs_dict[sample_id])[q: -q])
                seed_seq.sample_mean_Q2Q3_nonspecific_cov_dict[sample_id] = np.mean(
                    sorted(sample_nonspecific_covs_dict[sample_id])[q: -q])
                seed_seq.sample_mean_Q2Q3_summed_cov_dict[sample_id] = np.mean(
                    sorted(sample_summed_covs_dict[sample_id])[q: -q])


    def set_mods(self):
        """Predict modified positions in the tRNA seed.

        A modification requires a certain level of third- and/or fourth-most abundant nucleotides at
        the position in one or more samples. A modification in any particular sample additionally
        requires a certain level of second- through fourth-most abundant nucleotides at the
        position.

        There is currently an idiosyncracy in how modifications are set that results in the
        retention, but potential masking, of SNVs. If the position of a potential modification does
        not meet the coverage threshold for third- and fourth-most abundant nucleotides, the seed is
        not split into separate seeds around those SNVs, as occurs in anvi-trnaseq. Instead, the
        SNVs are simply not reported. This is a downside to imposing the aforementioned coverage
        threshold."""
        # Division by zero issues a numpy warning, but we handle it immediately by converting the
        # nan result to zero, so you don't need to see the warning. Unfortunately, this is the only
        # way we have found to suppress the warning.
        np.seterr(invalid='ignore')
        min_variation = self.min_variation
        min_third_fourth_nt = self.min_third_fourth_nt
        for seed_seq in self.seed_seqs:
            sample_mod_positions_dict = {}
            sample_variability_dict = {}
            sample_specific_nt_covs_dict = seed_seq.sample_specific_nt_covs_dict
            seed_seq_length = len(seed_seq.seq_string)
            sample_variations = []
            third_fourth_variations = np.zeros(seed_seq_length)
            for sample_id, specific_nt_covs in sample_specific_nt_covs_dict.items():
                specific_nt_covs_array = np.array(specific_nt_covs)
                specific_nt_covs_array.sort(axis=0)
                first_covs = specific_nt_covs_array[-1, :]
                second_covs = specific_nt_covs_array[-2, :]
                summed_covs = specific_nt_covs_array.sum(axis=0)
                sample_variations.append(np.nan_to_num(1 - first_covs / summed_covs))
                third_fourth_variations += np.nan_to_num(1 - (first_covs + second_covs) / summed_covs) >= min_third_fourth_nt
            sample_variations = np.array(sample_variations)
            third_fourth_variations = (third_fourth_variations > 0)
            total_mod_positions = np.nonzero((sample_variations >= min_variation).any(axis=0) & third_fourth_variations)[0]
            mod_sample_variations = sample_variations[:, total_mod_positions]
            for sample_num, sample_id in enumerate(sample_specific_nt_covs_dict.keys()):
                sample_mod_positions = total_mod_positions[np.nonzero(mod_sample_variations[sample_num, :] >= min_variation)[0]]
                sample_mod_positions_dict[sample_id] = sample_mod_positions.tolist()
                sample_variability_dict[sample_id] = sample_mod_positions.size * 1000 / seed_seq_length
            seed_seq.total_mod_positions = total_mod_positions.tolist()
            seed_seq.sample_mod_positions_dict = sample_mod_positions_dict
            seed_seq.sample_variability_dict = sample_variability_dict
        np.seterr(invalid='warn')


    def set_consensus_mod_nts(self):
        """Change predicted modified nucleotides in the seed consensus sequences to the nucleotides
        supported by the "preferred" treated samples, e.g., demethylase splits, with the goal of
        increasing the accuracy of the underlying base call."""
        if not self.preferred_treatment:
            return

        for seed_seq in self.seed_seqs:
            seq_string = seed_seq.seq_string
            preferred_nt_cov_dict = {nt: np.zeros(len(seq_string), int) for nt in UNAMBIG_NTS}
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                if norm_seq_summary.sample_id in self.preferred_trnaseq_db_sample_ids:
                    for nt in UNAMBIG_NTS:
                        preferred_nt_cov_dict[nt] += norm_seq_summary.specific_nt_covs_dict[nt]
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                if mod_seq_summary.sample_id in self.preferred_trnaseq_db_sample_ids:
                    for nt in UNAMBIG_NTS:
                        preferred_nt_cov_dict[nt] += mod_seq_summary.specific_nt_covs_dict[nt]

            preferred_nt_cov_array = np.array([preferred_nt_cov_dict[nt] for nt in UNAMBIG_NTS])
            for mod_pos in seed_seq.total_mod_positions:
                mod_covs = preferred_nt_cov_array[:, mod_pos]
                if mod_covs.sum() == 0:
                    # The preferred treatments do not have specific coverage of the modified site.
                    continue
                seq_string = seq_string[: mod_pos] + INT_NT_DICT[np.argmax(mod_covs) + 1] + seq_string[mod_pos + 1: ]
            seed_seq.seq_string = seq_string


    def set_sample_dels(self):
        for seed_seq in self.seed_seqs:
            sample_specific_dels_dict = {}
            sample_nonspecific_dels_dict = {}
            seed_seq_length = len(seed_seq.seq_string)

            for sample_id in self.trnaseq_db_sample_ids:
                sample_specific_dels_dict[sample_id] = np.zeros(seed_seq_length, int)
                sample_nonspecific_dels_dict[sample_id] = np.zeros(seed_seq_length, int)

            for mod_seq_summary in seed_seq.mod_seq_summaries:
                sample_id = mod_seq_summary.sample_id

                sample_specific_dels_dict[sample_id] += mod_seq_summary.specific_del_covs
                sample_nonspecific_dels_dict[sample_id] += mod_seq_summary.nonspecific_del_covs

            seed_seq.sample_specific_dels_dict = sample_specific_dels_dict
            seed_seq.sample_nonspecific_dels_dict = sample_nonspecific_dels_dict


    def gen_contigs_db(self):
        """Generate a contigs database of tRNA seeds. The create method of dbops.ContigsDatabase is
        not used because it tries to call genes, count kmers, and do other things that are
        irrelevant to tRNA-seq reads. There are no tRNA splits, but to satisfy the structure of the
        database, call every contig a split, and maintain tables on both."""
        self.progress.new("Generating a contigs database of tRNA seeds")
        self.progress.update("...")

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        contigs_db.touch()

        # Meta-values are set like in `dbops.ContigsDatabase.create`.
        contigs_db.db.set_meta_value('db_type', 'contigs')
        contigs_db.db.set_meta_value('db_variant', 'trnaseq')
        contigs_db.db.set_meta_value('project_name', self.project_name)
        contigs_db.db.set_meta_value('description', self.descrip if self.descrip else '_No description is provided_')
        contigs_db.db.set_meta_value('contigs_db_hash', self.contigs_db_hash)
        contigs_db.db.set_meta_value('split_length', 10000) # sys.maxsize
        contigs_db.db.set_meta_value('num_contigs', len(self.seed_seqs))
        contigs_db.db.set_meta_value('num_splits', len(self.seed_seqs))
        contigs_db.db.set_meta_value('total_length', self.total_seed_length)
        contigs_db.db.set_meta_value('kmer_size', 0)
        contigs_db.db.set_meta_value('gene_level_taxonomy_source', None)
        contigs_db.db.set_meta_value('gene_function_sources', 'Transfer_RNAs')
        contigs_db.db.set_meta_value('genes_are_called', True)
        contigs_db.db.set_meta_value('external_gene_calls', True)
        contigs_db.db.set_meta_value('external_gene_amino_acid_seqs', False)
        contigs_db.db.set_meta_value('skip_predict_frame', True)
        contigs_db.db.set_meta_value('splits_consider_gene_calls', False)
        contigs_db.db.set_meta_value('scg_taxonomy_was_run', False)
        contigs_db.db.set_meta_value('scg_taxonomy_database_version', None)
        contigs_db.db.set_meta_value('trna_taxonomy_was_run', False)
        contigs_db.db.set_meta_value('trna_taxonomy_database_version', None)
        contigs_db.db.set_meta_value('creation_date', time.time())

        contigs_db.db.insert_many('contig_sequences', [(seed_seq.name, seed_seq.seq_string) for seed_seq in self.seed_seqs])
        contigs_db.db.insert_many('contigs_basic_info', self.get_contigs_basic_info_table_entries())
        contigs_db.db.insert_many('splits_basic_info', self.get_splits_basic_info_table_entries())
        contigs_db.db.insert_many('hmm_hits', self.get_hmm_hits_table_entries())
        contigs_db.db.insert_many('hmm_hits_in_splits', self.get_hmm_hits_in_splits_table_entries())
        # tRNA predictions are treated like HMM or tRNAScan-SE hits. The blank columns of the HMM
        # hits info table are 'ref', 'search_type', 'domain' and 'genes'.
        contigs_db.db.insert('hmm_hits_info', ('Transfer_RNAs', '', 'Transfer_RNAs', None, ''))
        contigs_db.db.insert_many('genes_in_contigs', self.get_genes_in_contigs_table_entries())
        contigs_db.db.insert_many('gene_amino_acid_sequences', [(i, '') for i in range(len(self.seed_seqs))])
        contigs_db.db.insert_many('genes_in_splits', self.get_genes_in_splits_table_entries())
        contigs_db.db.insert_many('gene_functions', self.get_gene_functions_table_entries())

        contigs_db.disconnect()

        self.progress.end()


    def get_contigs_basic_info_table_entries(self):
        entries = []
        for seed_seq in self.seed_seqs:
            seq_string = seed_seq.seq_string
            entries.append(
                (seed_seq.name,
                 len(seq_string),
                 seed_seq.gc_fraction,
                 1)
            )
        return entries


    def get_splits_basic_info_table_entries(self):
        entries = []
        for seed_seq in self.seed_seqs:
            entries.append(
                (seed_seq.name + '_split_00001',
                 0, # Order of split in parent contig
                 0, # Start in contig
                 len(seed_seq.seq_string), # Stop in contig
                 len(seed_seq.seq_string), # Split length
                 seed_seq.gc_fraction, # GC content of split
                 seed_seq.gc_fraction, # GC content of parent contig
                 seed_seq.name)
            )
        return entries


    def get_hmm_hits_table_entries(self):
        """tRNA seeds are analogous to tRNA gene predictions from a metagenomic contigs database."""
        entries = []
        for i, seed_seq in enumerate(self.seed_seqs):
            entries.append(
                (i, # Entry ID
                'Transfer_RNAs', # Source, à la tRNA gene prediction via tRNAScan-SE
                hashlib.sha1(seed_seq.seq_string.encode('utf-8')).hexdigest(), # "Gene unique identifier"
                i, # "Gene callers ID"
                ANTICODON_TO_AA[seed_seq.anticodon_seq_string] + '_' + seed_seq.anticodon_seq_string, # "Gene name", à la tRNA gene prediction via tRNAScan-SE
                '-', # "Gene HMM ID"
                0.0) # "HMM E-value"
            )
        return entries


    def get_hmm_hits_in_splits_table_entries(self):
        entries = []
        for i, seed_seq in enumerate(self.seed_seqs):
            entries.append(
                (i, # Entry ID
                 seed_seq.name + '_split_00001', # Split name
                 100, # Percentage of "HMM hit" in split
                 'Transfer_RNAs')
            )
        return entries


    def get_genes_in_contigs_table_entries(self):
        entries = []
        for i, seed_seq in enumerate(self.seed_seqs):
            entries.append(
                (i, # Gene callers ID
                 seed_seq.name, # Contig name
                 0, # Gene start in contig
                 len(seed_seq.seq_string), # Gene stop in contig
                 'f', # Direction of gene call on contig
                 0, # Is partial gene call: for now, say all seeds are "full tRNAs"
                 2, # Call type: 1 = coding, 2 = noncoding, 3 = unknown
                 'anvi-trnaseq', # Gene caller
                 tables.trnaseq_db_version) # Version of caller
            )
        return entries


    def get_genes_in_splits_table_entries(self):
        entries = []
        for i, seed_seq in enumerate(self.seed_seqs):
            entries.append(
                (seed_seq.name + '_split_00001',
                 i,
                 0,
                 len(seed_seq.seq_string),
                 100)
            )
        return entries


    def get_gene_functions_table_entries(self):
        entries = []
        for i, seed_seq in enumerate(self.seed_seqs):
            entries.append(
                (i,
                 'Transfer_RNAs',
                 '%s_%s_%d' % (ANTICODON_TO_AA[seed_seq.anticodon_seq_string], seed_seq.anticodon_seq_string, i),
                 'tRNA transcript',
                 0.0)
            )
        return entries


    def gen_auxiliary_db(self, db_cov_type):
        if db_cov_type == 'specific':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(
                self.specific_auxiliary_db_path, self.contigs_db_hash, create_new=True)
            for seed_seq in self.seed_seqs:
                split_name = seed_seq.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:\
                    auxiliary_db.append(split_name,
                                        sample_id,
                                        seed_seq.sample_specific_covs_dict[sample_id].tolist())
        elif db_cov_type == 'nonspecific':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(
                self.nonspecific_auxiliary_db_path, self.contigs_db_hash, create_new=True)
            for seed_seq in self.seed_seqs:
                split_name = seed_seq.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name,
                                        sample_id,
                                        seed_seq.sample_nonspecific_covs_dict[sample_id].tolist())
        elif db_cov_type == 'combined':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(
                self.combined_auxiliary_db_path, self.contigs_db_hash, create_new=True)
            for seed_seq in self.seed_seqs:
                split_name = seed_seq.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name,
                                        sample_id + '_specific',
                                        seed_seq.sample_specific_covs_dict[sample_id].tolist())
                    auxiliary_db.append(split_name,
                                        sample_id + '_nonspecific',
                                        seed_seq.sample_nonspecific_covs_dict[sample_id].tolist())
        elif db_cov_type == 'summed':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(
                self.summed_auxiliary_db_path, self.contigs_db_hash, create_new=True)
            for seed_seq in self.seed_seqs:
                split_name = seed_seq.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name,
                                        sample_id,
                                        (seed_seq.sample_specific_covs_dict[sample_id] + seed_seq.sample_nonspecific_covs_dict[sample_id]).tolist())
        else:
            raise ConfigError(f"The type of profile database provided, {db_cov_type}, is not among "
                              "those that are recognized: 'specific', 'nonspecific', 'combined', and 'summed'.")

        auxiliary_db.store()
        auxiliary_db.close()


    def set_sample_total_covs(self):
        """For each input sample, find the total specific, nonspecific and summed coverage of the
        seeds across all positions (single integers)."""
        sample_total_specific_cov_dict = {sample_id: 0 for sample_id in self.trnaseq_db_sample_ids}
        sample_total_nonspecific_cov_dict = {sample_id: 0 for sample_id in self.trnaseq_db_sample_ids}
        for seed_seq in self.seed_seqs:
            for norm_seq_summary in seed_seq.unmod_norm_seq_summaries:
                sample_id = norm_seq_summary.sample_id
                specific_cov = norm_seq_summary.specific_covs.sum()
                nonspecific_cov = norm_seq_summary.nonspecific_covs.sum()
                sample_total_specific_cov_dict[sample_id] += specific_cov
                sample_total_nonspecific_cov_dict[sample_id] += nonspecific_cov
            for mod_seq_summary in seed_seq.mod_seq_summaries:
                sample_id = mod_seq_summary.sample_id
                specific_cov = mod_seq_summary.specific_covs.sum()
                nonspecific_cov = mod_seq_summary.nonspecific_covs.sum()
                sample_total_specific_cov_dict[sample_id] += specific_cov
                sample_total_nonspecific_cov_dict[sample_id] += nonspecific_cov

        sample_total_summed_cov_dict = {}
        for sample_id in self.trnaseq_db_sample_ids:
            sample_total_summed_cov_dict[sample_id] = sample_total_specific_cov_dict[sample_id] + sample_total_nonspecific_cov_dict[sample_id]

        self.sample_total_specific_cov_dict = sample_total_specific_cov_dict
        self.sample_total_nonspecific_cov_dict = sample_total_nonspecific_cov_dict
        self.sample_total_summed_cov_dict = sample_total_summed_cov_dict


    def set_sample_overall_mean_covs(self):
        """For each input sample, find the mean specific, nonspecific and summed coverage of all
        seeds across all positions (single numbers)."""
        sample_overall_mean_specific_cov_dict = {}
        sample_overall_mean_nonspecific_cov_dict = {}
        sample_overall_mean_summed_cov_dict = {}
        for sample_id, total_specific_cov in self.sample_total_specific_cov_dict.items():
            sample_overall_mean_specific_cov_dict[sample_id] = total_specific_cov / self.total_seed_length
        for sample_id, total_nonspecific_cov in self.sample_total_nonspecific_cov_dict.items():
            sample_overall_mean_nonspecific_cov_dict[sample_id] = total_nonspecific_cov / self.total_seed_length
        for sample_id, total_summed_cov in self.sample_total_summed_cov_dict.items():
            sample_overall_mean_summed_cov_dict[sample_id] = total_summed_cov / self.total_seed_length

        self.sample_overall_mean_specific_cov_dict = sample_overall_mean_specific_cov_dict
        self.sample_overall_mean_nonspecific_cov_dict = sample_overall_mean_nonspecific_cov_dict
        self.sample_overall_mean_summed_cov_dict = sample_overall_mean_summed_cov_dict


    def set_sample_mean_covs(self):
        """Mean coverage of each seed in a sample."""
        for seed_seq in self.seed_seqs:
            sample_mean_specific_cov_dict = {}
            sample_mean_nonspecific_cov_dict = {}
            sample_mean_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_mean_specific_cov_dict[sample_id] = seed_seq.sample_specific_covs_dict[sample_id].mean()
                sample_mean_nonspecific_cov_dict[sample_id] = seed_seq.sample_nonspecific_covs_dict[sample_id].mean()
                sample_mean_summed_cov_dict[sample_id] = seed_seq.sample_summed_covs_dict[sample_id].mean()
            seed_seq.sample_mean_specific_cov_dict = sample_mean_specific_cov_dict
            seed_seq.sample_mean_nonspecific_cov_dict = sample_mean_nonspecific_cov_dict
            seed_seq.sample_mean_summed_cov_dict = sample_mean_summed_cov_dict


    def set_sample_std_covs(self):
        """Standard deviation of the coverage of each seed in a sample."""
        for seed_seq in self.seed_seqs:
            sample_std_specific_cov_dict = {}
            sample_std_nonspecific_cov_dict = {}
            sample_std_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_std_specific_cov_dict[sample_id] = seed_seq.sample_specific_covs_dict[sample_id].std()
                sample_std_nonspecific_cov_dict[sample_id] = seed_seq.sample_nonspecific_covs_dict[sample_id].std()
                sample_std_summed_cov_dict[sample_id] = seed_seq.sample_summed_covs_dict[sample_id].std()
            seed_seq.sample_std_specific_cov_dict = sample_std_specific_cov_dict
            seed_seq.sample_std_nonspecific_cov_dict = sample_std_nonspecific_cov_dict
            seed_seq.sample_std_summed_cov_dict = sample_std_summed_cov_dict


    def set_sample_abundances(self):
        """For each sample, and for specific and nonspecific coverages, abundance is defined as the
        mean coverage of the seed divided by the mean total coverage of the sample across all
        seeds."""
        for seed_seq in self.seed_seqs:
            sample_specific_abundances_dict = {}
            sample_nonspecific_abundances_dict = {}
            sample_summed_abundances_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_specific_abundances_dict[sample_id] = seed_seq.sample_mean_specific_cov_dict[sample_id] / self.sample_overall_mean_specific_cov_dict[sample_id]
                sample_nonspecific_abundances_dict[sample_id] = seed_seq.sample_mean_nonspecific_cov_dict[sample_id] / self.sample_overall_mean_nonspecific_cov_dict[sample_id]
                sample_summed_abundances_dict[sample_id] = seed_seq.sample_mean_summed_cov_dict[sample_id] / self.sample_overall_mean_summed_cov_dict[sample_id]
            seed_seq.sample_specific_abundances_dict = sample_specific_abundances_dict
            seed_seq.sample_nonspecific_abundances_dict = sample_nonspecific_abundances_dict
            seed_seq.sample_summed_abundances_dict = sample_summed_abundances_dict


    def set_sample_normalization_multipliers(self):
        """Set a normalization constant for each sample to scale their coverages, allowing the
        relative abundance of seeds in a sample to compared between samples. Normalization is based
        on the total specific coverage of each sample -- one can imagine other ways of doing this,
        including use of summed specific and nonspecific coverage, but this would require
        deconvoluting the multiple representation of nonspecific reads."""
        sample_normalization_multiplier_dict = {}
        min_total_specific_cov = min([v for v in self.sample_total_specific_cov_dict.values()])
        for sample_id, total_specific_cov in self.sample_total_specific_cov_dict.items():
            sample_normalization_multiplier_dict[sample_id] = min_total_specific_cov / total_specific_cov
        self.sample_normalization_multiplier_dict = sample_normalization_multiplier_dict


    def set_sample_normalized_mean_Q2Q3_coverages(self):
        """Scale mean coverages for comparison across samples."""
        for seed_seq in self.seed_seqs:
            sample_normalized_mean_Q2Q3_specific_cov_dict = {}
            sample_normalized_mean_Q2Q3_nonspecific_cov_dict = {}
            sample_normalized_mean_Q2Q3_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_normalized_mean_Q2Q3_specific_cov_dict[sample_id] = seed_seq.sample_mean_Q2Q3_specific_cov_dict[sample_id] * self.sample_normalization_multiplier_dict[sample_id]
                sample_normalized_mean_Q2Q3_nonspecific_cov_dict[sample_id] = seed_seq.sample_mean_Q2Q3_nonspecific_cov_dict[sample_id] * self.sample_normalization_multiplier_dict[sample_id]
                sample_normalized_mean_Q2Q3_summed_cov_dict[sample_id] = seed_seq.sample_mean_Q2Q3_summed_cov_dict[sample_id] * self.sample_normalization_multiplier_dict[sample_id]
            seed_seq.sample_normalized_mean_Q2Q3_specific_cov_dict = sample_normalized_mean_Q2Q3_specific_cov_dict
            seed_seq.sample_normalized_mean_Q2Q3_nonspecific_cov_dict = sample_normalized_mean_Q2Q3_nonspecific_cov_dict
            seed_seq.sample_normalized_mean_Q2Q3_summed_cov_dict = sample_normalized_mean_Q2Q3_summed_cov_dict


    def set_sample_detections(self):
        """Find the proportion of each seed sequence covered by reads in a sample."""
        for seed_seq in self.seed_seqs:
            seed_seq_length = len(seed_seq.seq_string)
            sample_specific_detection_dict = {}
            sample_nonspecific_detection_dict = {}
            sample_summed_detection_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_specific_detection_dict[sample_id] = seed_seq.sample_specific_covs_dict[sample_id].nonzero()[0].size / seed_seq_length
                sample_nonspecific_detection_dict[sample_id] = seed_seq.sample_nonspecific_covs_dict[sample_id].nonzero()[0].size / seed_seq_length
                sample_summed_detection_dict[sample_id] = seed_seq.sample_summed_covs_dict[sample_id].nonzero()[0].size / seed_seq_length
            seed_seq.sample_specific_detection_dict = sample_specific_detection_dict
            seed_seq.sample_nonspecific_detection_dict = sample_nonspecific_detection_dict
            seed_seq.sample_summed_detection_dict = sample_summed_detection_dict


    def set_sample_relative_abundances(self):
        """Relative abundance represents the coverage of the seed in one sample relative to the
        total coverage of the seed across samples -- relative abundances sum to one across
        samples."""
        np.seterr(invalid='ignore')
        for seed_seq in self.seed_seqs:
            sample_specific_relative_abundances_dict = {}
            sample_nonspecific_relative_abundances_dict = {}
            sample_summed_relative_abundances_dict = {}
            pansample_normalized_mean_Q2Q3_specific_cov = sum(seed_seq.sample_normalized_mean_Q2Q3_specific_cov_dict.values())
            pansample_normalized_mean_Q2Q3_nonspecific_cov = sum(seed_seq.sample_normalized_mean_Q2Q3_nonspecific_cov_dict.values())
            pansample_normalized_mean_Q2Q3_summed_cov = sum(seed_seq.sample_normalized_mean_Q2Q3_summed_cov_dict.values())
            for sample_id in self.trnaseq_db_sample_ids:
                sample_specific_relative_abundances_dict[sample_id] = seed_seq.sample_normalized_mean_Q2Q3_specific_cov_dict[sample_id] / pansample_normalized_mean_Q2Q3_specific_cov
                sample_nonspecific_relative_abundances_dict[sample_id] = seed_seq.sample_normalized_mean_Q2Q3_nonspecific_cov_dict[sample_id] / pansample_normalized_mean_Q2Q3_nonspecific_cov
                sample_summed_relative_abundances_dict[sample_id] = seed_seq.sample_normalized_mean_Q2Q3_summed_cov_dict[sample_id] / pansample_normalized_mean_Q2Q3_summed_cov
            seed_seq.sample_specific_relative_abundances_dict = sample_specific_relative_abundances_dict
            seed_seq.sample_nonspecific_relative_abundances_dict = sample_nonspecific_relative_abundances_dict
            seed_seq.sample_summed_relative_abundances_dict = sample_summed_relative_abundances_dict
        np.seterr(invalid='warn')


    def set_sample_max_normalized_ratios(self):
        """The max normalized coverage ratio represents the coverage of the seed in one sample
        relative to the maximum coverage amongst the samples -- one sample will always have a value
        equal to one."""
        for seed_seq in self.seed_seqs:
            sample_specific_max_normalized_ratio_dict = {}
            sample_nonspecific_max_normalized_ratio_dict = {}
            sample_summed_max_normalized_ratio_dict = {}
            max_normalized_mean_Q2Q3_specific_cov = max(seed_seq.sample_normalized_mean_Q2Q3_specific_cov_dict.values())
            max_normalized_mean_Q2Q3_nonspecific_cov = max(seed_seq.sample_normalized_mean_Q2Q3_nonspecific_cov_dict.values())
            max_normalized_mean_Q2Q3_summed_cov = max(seed_seq.sample_normalized_mean_Q2Q3_summed_cov_dict.values())
            for sample_id in self.trnaseq_db_sample_ids:
                sample_specific_max_normalized_ratio_dict[sample_id] = seed_seq.sample_normalized_mean_Q2Q3_specific_cov_dict[sample_id] / max_normalized_mean_Q2Q3_specific_cov if max_normalized_mean_Q2Q3_specific_cov else 0
                sample_nonspecific_max_normalized_ratio_dict[sample_id] = seed_seq.sample_normalized_mean_Q2Q3_nonspecific_cov_dict[sample_id] / max_normalized_mean_Q2Q3_nonspecific_cov if max_normalized_mean_Q2Q3_nonspecific_cov else 0
                sample_summed_max_normalized_ratio_dict[sample_id] = seed_seq.sample_normalized_mean_Q2Q3_summed_cov_dict[sample_id] / max_normalized_mean_Q2Q3_summed_cov if max_normalized_mean_Q2Q3_summed_cov else 0
            seed_seq.sample_specific_max_normalized_ratio_dict = sample_specific_max_normalized_ratio_dict
            seed_seq.sample_nonspecific_max_normalized_ratio_dict = sample_nonspecific_max_normalized_ratio_dict
            seed_seq.sample_summed_max_normalized_ratio_dict = sample_summed_max_normalized_ratio_dict


    def set_variable_nts_table_entries(self):
        """Variable nucleotides in the profile databases are nucleotides with predicted
        modifications, not single nucleotide variants. Modifications are determined from specific
        coverage, but are displayed in nonspecific and summed profile databases as well."""
        entries = []
        for sample_id in self.trnaseq_db_sample_ids:
            for i, seed_seq in enumerate(self.seed_seqs):
                specific_covs = seed_seq.sample_specific_covs_dict[sample_id]
                specific_nt_cov_arrays = seed_seq.sample_specific_nt_covs_dict[sample_id]
                for pos in seed_seq.sample_mod_positions_dict[sample_id]:
                    total_cov = specific_covs[pos]
                    specific_nt_covs = [arr[pos] for arr in specific_nt_cov_arrays]
                    max_nt_cov = max(specific_nt_covs)
                    sorted_nt_covs = sorted(zip(UNAMBIG_NTS, specific_nt_covs), key=lambda x: -x[1])
                    ref_nt = sorted_nt_covs[0][0]
                    secondary_nt = sorted_nt_covs[1][0]
                    entries.append((sample_id,
                                    seed_seq.name + '_split_00001',
                                    pos, # Position in split
                                    pos, # Position in contig
                                    i, # Corresponding gene call
                                    1, # In noncoding gene call
                                    0, # In coding gene call
                                    0, # Base position in codon (0 for noncoding gene call)
                                    -1, # Codon order in gene (-1 for noncoding gene call)
                                    total_cov,
                                    0, # Coverage outlier in split (0 or 1)
                                    0, # Coverage outlier in contig (0 or 1)
                                    1 - max_nt_cov / total_cov, # Departure from reference
                                    ref_nt + secondary_nt, # Competing nts (top 2)
                                    ref_nt,
                                    specific_nt_covs[0], # A coverage
                                    specific_nt_covs[1], # C coverage
                                    specific_nt_covs[2], # G coverage
                                    specific_nt_covs[3], # T coverage
                                    0))
        self.variable_nts_table_entries = entries


    def set_indels_table_entries(self):
        """Deletions are determined separately from specific and nonspecific coverages. Currently,
        insertions are not sought by `anvi-trnaseq`, as the addition of nontemplated nucleotides to
        cDNA in reverse transcription is less common than the deletion of nucleotides from the
        transcript."""
        specific_entries = []
        nonspecific_entries = []
        summed_entries = []
        for sample_id in self.trnaseq_db_sample_ids:
            for i, seed_seq in enumerate(self.seed_seqs):
                specific_del_covs = seed_seq.sample_specific_dels_dict[sample_id]
                specific_covs = seed_seq.sample_specific_covs_dict[sample_id]
                for pos in np.nonzero(specific_del_covs)[0]:
                    specific_del_cov = specific_del_covs[pos]
                    specific_cov = specific_covs[pos]
                    if specific_cov == 0:
                        # This only occurs when there is specific coverage of the deletion but not
                        # the reference nucleotide.
                        del_fraction = 1
                    else:
                        del_fraction = specific_del_cov / specific_cov
                    if del_fraction >= self.min_del_fraction:
                        specific_entries.append((sample_id,
                                                 seed_seq.name + '_split_00001',
                                                 pos, # Position in split
                                                 pos, # Position in contig
                                                 i, # Corresponding gene call
                                                 1, # In noncoding gene call
                                                 0, # In coding gene call
                                                 0, # Base position in codon (0 for noncoding gene call)
                                                 -1, # Codon order in gene (-1 for noncoding gene call)
                                                 0, # Coverage outlier in split (0 or 1)
                                                 0, # Coverage outlier in contig (0 or 1)
                                                 seed_seq.seq_string[pos], # Reference nt
                                                 'DEL', # Type of indel
                                                 '', # Indel sequence ('' for deletion)
                                                 1, # Indel length
                                                 specific_del_cov, # Deletion count (coverage)
                                                 specific_cov)) # Reference sequence coverage

                nonspecific_del_covs = seed_seq.sample_nonspecific_dels_dict[sample_id]
                nonspecific_covs = seed_seq.sample_nonspecific_covs_dict[sample_id]
                for pos in np.nonzero(nonspecific_del_covs)[0]:
                    nonspecific_del_cov = nonspecific_del_covs[pos]
                    nonspecific_cov = nonspecific_covs[pos]
                    if nonspecific_cov == 0:
                        del_fraction = 1
                    else:
                        del_fraction = nonspecific_del_cov / nonspecific_cov
                    if del_fraction >= self.min_del_fraction:
                        nonspecific_entries.append((sample_id,
                                                    seed_seq.name + '_split_00001',
                                                    pos,
                                                    pos,
                                                    i,
                                                    1,
                                                    0,
                                                    0,
                                                    -1,
                                                    0,
                                                    0,
                                                    seed_seq.seq_string[pos],
                                                    'DEL',
                                                    '',
                                                    1,
                                                    nonspecific_del_cov,
                                                    nonspecific_cov))

                summed_del_covs = specific_del_covs + nonspecific_del_covs
                summed_covs = specific_covs + nonspecific_covs
                for pos in np.nonzero(summed_del_covs)[0]:
                    summed_del_cov = summed_del_covs[pos]
                    summed_cov = summed_covs[pos]
                    if summed_cov == 0:
                        del_fraction = 1
                    else:
                        del_fraction = summed_del_cov / summed_cov
                    if del_fraction >= self.min_del_fraction:
                        summed_entries.append((sample_id,
                                               seed_seq.name + '_split_00001',
                                               pos,
                                               pos,
                                               i,
                                               1,
                                               0,
                                               0,
                                               -1,
                                               0,
                                               0,
                                               seed_seq.seq_string[pos],
                                               'DEL',
                                               '',
                                               1,
                                               summed_del_cov,
                                               summed_cov))
        self.specific_indels_table_entries = specific_entries
        self.nonspecific_indels_table_entries = nonspecific_entries
        self.summed_indels_table_entries = summed_entries


    def gen_profile_db(self, db_cov_type):
        self.progress.new(f"Generating {db_cov_type} profile database")
        self.progress.update("...")

        if db_cov_type == 'specific':
            profile_db_path = self.specific_profile_db_path
        elif db_cov_type == 'nonspecific':
            profile_db_path = self.nonspecific_profile_db_path
        elif db_cov_type == 'combined':
            profile_db_path = self.combined_profile_db_path
        elif db_cov_type == 'summed':
            profile_db_path = self.summed_profile_db_path
        else:
            raise ConfigError(f"The tRNA-seq coverage type, {db_cov_type}, is not recognized. "
                              "The only valid options are specific, nonspecific, summed, and combined.")

        profile_db = dbops.ProfileDatabase(profile_db_path)
        profile_db.touch()

        # Profile database meta-values are set in a parallel fashion to `merger.MultipleRuns.merge`.
        profile_db.db.set_meta_value('creation_date', time.time())
        profile_db.db.set_meta_value('db_type', 'profile')

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, quiet=True)
        profile_db.db.set_meta_value('contigs_db_hash', contigs_db.meta['contigs_db_hash'])
        profile_db.db.set_meta_value('sample_id', contigs_db.meta['project_name'])
        contigs_db.disconnect()

        if db_cov_type == 'combined':
            profile_db.db.set_meta_value('samples', ', '.join([sample_id + '_' + cov_type
                                                               for sample_id in self.trnaseq_db_sample_ids
                                                               for cov_type in ('specific', 'nonspecific')]))
        else:
            profile_db.db.set_meta_value('samples', ', '.join([sample_id for sample_id in self.trnaseq_db_sample_ids]))
        # The total number of reads mapped is not calculated, as that would require deconvoluting
        # the number of reads that mapped nonspecifically. Also, the total number of mapped reads is
        # less informative here than in metagenomics, since they vary greatly in length.
        # profile_db.db.set_meta_value('total_reads_mapped', -1)
        profile_db.db.set_meta_value('merged', True)
        profile_db.db.set_meta_value('blank', False)
        profile_db.db.set_meta_value('default_view', 'mean_coverage')
        profile_db.db.set_meta_value('min_contig_length', 1)
        profile_db.db.set_meta_value('max_contig_length', sys.maxsize)
        profile_db.db.set_meta_value('SNVs_profiled', False)
        profile_db.db.set_meta_value('SCVs_profiled', False)
        profile_db.db.set_meta_value('INDELs_profiled', False)
        profile_db.db.set_meta_value('num_contigs', len(self.seed_seqs))
        profile_db.db.set_meta_value('num_splits', len(self.seed_seqs))
        profile_db.db.set_meta_value('total_length', self.total_seed_length)
        profile_db.db.set_meta_value('min_coverage_for_variability', 1)
        profile_db.db.set_meta_value('min_indel_fraction', 0)
        profile_db.db.set_meta_value('report_variability_full', False)
        profile_db.db.set_meta_value('description', self.descrip if self.descrip else '_No description is provided_')
        # profile_db.db.set_meta_value('min_percent_identity', -1)

        # Whereas variability in metagenomics refers to SNVs, here it refers to modifications.
        # Modifications are only identified from specific coverage.
        if db_cov_type == 'specific':
            profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('variable_nucleotides', ','.join('?' * len(tables.variable_nts_table_structure))),
                                     self.variable_nts_table_entries)
            profile_db.db.commit()

        if db_cov_type == 'specific' or db_cov_type == 'nonspecific' or db_cov_type == 'summed':
            for attr, table_basename in [
                ('sample_mean_' + db_cov_type + '_cov_dict', 'mean_coverage'),
                ('sample_std_' + db_cov_type + '_cov_dict', 'std_coverage'),
                ('sample_' + db_cov_type + '_abundances_dict', 'abundance'),
                ('sample_' + db_cov_type + '_relative_abundances_dict', 'relative_abundance'),
                ('sample_' + db_cov_type + '_detection_dict', 'detection'),
                ('sample_' + db_cov_type + '_max_normalized_ratio_dict', 'max_normalized_ratio'),
                ('sample_mean_Q2Q3_' + db_cov_type + '_cov_dict', 'mean_coverage_Q2Q3')]:
                data_dict = self.get_specific_nonspecific_or_summed_data_dict(attr)
                self.create_specific_nonspecific_or_summed_contigs_and_splits_tables(profile_db_path, table_basename, data_dict)
            variability_data_dict = self.get_specific_nonspecific_or_summed_data_dict('sample_variability_dict')
            self.create_specific_nonspecific_or_summed_contigs_and_splits_tables(profile_db_path, 'variability', data_dict)

            profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('indels', ','.join('?' * len(tables.indels_table_structure))),
                                     getattr(self, db_cov_type + '_indels_table_entries'))
        elif db_cov_type == 'combined':
            for specific_attr, nonspecific_attr, table_basename in [
                ('sample_mean_specific_cov_dict', 'sample_mean_nonspecific_cov_dict', 'mean_coverage'),
                ('sample_std_specific_cov_dict', 'sample_std_nonspecific_cov_dict', 'std_coverage'),
                ('sample_specific_abundances_dict', 'sample_nonspecific_abundances_dict', 'abundance'),
                ('sample_specific_relative_abundances_dict', 'sample_nonspecific_relative_abundances_dict', 'relative_abundance'),
                ('sample_specific_detection_dict', 'sample_nonspecific_detection_dict', 'detection'),
                ('sample_specific_max_normalized_ratio_dict', 'sample_nonspecific_max_normalized_ratio_dict', 'max_normalized_ratio'),
                ('sample_mean_Q2Q3_specific_cov_dict', 'sample_mean_Q2Q3_nonspecific_cov_dict', 'mean_coverage_Q2Q3')]:
                data_dict = self.get_combined_data_dict(specific_attr, nonspecific_attr)
                self.create_combined_contigs_and_splits_tables(profile_db_path, table_basename, data_dict)
            variability_data_dict = self.get_combined_data_dict('sample_variability_dict', 'sample_variability_dict')
            self.create_combined_contigs_and_splits_tables(profile_db_path, 'variability', data_dict)

            combined_indels_table_entries = []
            for entry in self.specific_indels_table_entries:
                combined_indels_table_entries.append((entry[0] + '_specific', ) + entry[1: ])
            combined_indels_table_entries = []
            for entry in self.nonspecific_indels_table_entries:
                combined_indels_table_entries.append((entry[0] + '_nonspecific', ) + entry[1: ])
            profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('indels', ','.join('?' * len(tables.indels_table_structure))),
                                     combined_indels_table_entries)

        profile_db.db.commit()
        self.progress.end()

        # Add layers for anticodon and corresponding amino acid.
        items_additional_data_table = miscdata.MiscDataTableFactory(argparse.Namespace(profile_db=profile_db_path, target_data_table='items'))
        data_dict = {}
        for seed_seq in self.seed_seqs:
            data_dict[seed_seq.name + '_split_00001'] = {'anticodon': seed_seq.anticodon_seq_string,
                                                         'amino_acid': ANTICODON_TO_AA[seed_seq.anticodon_seq_string]}
        items_additional_data_table.add(data_dict, ['anticodon', 'amino_acid'])

        # Cluster tRNA seeds to form the central dendrogram in anvi-interactive.
        dbops.do_hierarchical_clustering_of_items(profile_db_path,
                                                  constants.clustering_configs['trnaseq'],
                                                  [seed_seq.name + '_split_00001' for seed_seq in self.seed_seqs],
                                                  {'CONTIGS.db': self.contigs_db_path, 'PROFILE.db': profile_db_path},
                                                  input_directory=os.path.dirname(profile_db_path),
                                                  default_clustering_config=constants.trnaseq_default,
                                                  distance=self.distance,
                                                  linkage=self.linkage,
                                                  run=self.run,
                                                  progress=self.progress)
        profile_db.db.set_meta_value('items_ordered', True)
        profile_db.db.disconnect()

        # Cluster samples by "view" data to find possible sample layer orderings.
        profile_db_super = dbops.ProfileSuperclass(argparse.Namespace(profile_db=profile_db_path))
        profile_db_super.load_views(omit_parent_column=True)
        layer_orders_data_dict = {}
        failed_attempts = []
        for essential_field in [
            'abundance',
            'std_coverage',
            'mean_coverage',
            'mean_coverage_Q2Q3',
            'max_normalized_ratio',
            'relative_abundance',
            'detection',
            'variability'
        ]:
            try:
                data_value = clustering.get_newick_tree_data_for_dict(profile_db_super.views[essential_field]['dict'],
                                                                      distance=self.distance,
                                                                      linkage=self.linkage,
                                                                      transpose=True)
                layer_orders_data_dict[essential_field] = {'data_value': data_value, 'data_type': 'newick'}
            except:
                failed_attempts.append(essential_field)
        if not len(layer_orders_data_dict):
            self.run.warning("This may or may not be important: anvi'o attempted to generate orders for your "
                             "samples based on the view data, however, it failed :/")
            return
        if len(failed_attempts):
            self.run.warning("While anvi'o was trying to generate clusterings of samples based on view data "
                             f"available in the {db_cov_type} profile, clustering of some of the essential data "
                             "failed. It is likely not a very big deal, but you shall be the judge of it. "
                             "Anvi'o now proceeds to store layers order information for those view items "
                             "the clustering in fact worked. Here is the list of stuff that failed: '%s'"\
                              % (', '.join(failed_attempts)))
        # Add the layer orders quietly.
        TableForLayerOrders(argparse.Namespace(profile_db=profile_db_path), r=terminal.Run(verbose=False)).add(layer_orders_data_dict)


    def get_specific_nonspecific_or_summed_data_dict(self, seed_seq_attr):
        """Get data from seed sequences to generate a table in a specific, nonspecific, or summed
        profile database."""
        data_dict = {}
        for seed_seq in self.seed_seqs:
            seed_seq_name = seed_seq.name
            seed_seq_split_name = seed_seq_name + '_split_00001'
            data_dict[seed_seq_split_name] = {'__parent__': seed_seq_name}
            for sample_id in self.trnaseq_db_sample_ids:
                data_dict[seed_seq_split_name][sample_id] = getattr(seed_seq, seed_seq_attr)[sample_id]
        return data_dict


    def get_combined_data_dict(self, specific_seed_seq_attr, nonspecific_seed_seq_attr):
        """Get data from seed sequences to generate a table in a combined profile database."""
        data_dict = {}
        for seed_seq in self.seed_seqs:
            seed_seq_name = seed_seq.name
            seed_seq_split_name = seed_seq_name + '_split_00001'
            data_dict[seed_seq_name + '_split_00001'] = {'__parent__': seed_seq_name}
            for sample_id in self.trnaseq_db_sample_ids:
                data_dict[seed_seq_split_name][sample_id + '_specific'] = getattr(seed_seq, specific_seed_seq_attr)[sample_id]
                data_dict[seed_seq_split_name][sample_id + '_nonspecific'] = getattr(seed_seq, nonspecific_seed_seq_attr)[sample_id]
        return data_dict


    def create_specific_nonspecific_or_summed_contigs_and_splits_tables(self, profile_db_path, table_basename, data_dict):
        """Create a pair of tables in a specific, nonspecific, or summed profile database. Contigs
        and splits tables contain the same information since tRNA, unlike a metagenomic contig, is
        not long enough to be split."""
        TablesForViews(profile_db_path).create_new_view(
            data_dict=data_dict,
            table_name=table_basename + '_contigs',
            table_structure=['contig'] + self.trnaseq_db_sample_ids + ['__parent__'],
            table_types=['text'] + ['numeric'] * len(self.trnaseq_db_sample_ids) + ['text'],
            view_name=None)
        TablesForViews(profile_db_path).create_new_view(
            data_dict=data_dict,
            table_name=table_basename + '_splits',
            table_structure=['contig'] + self.trnaseq_db_sample_ids + ['__parent__'],
            table_types=['text'] + ['numeric'] * len(self.trnaseq_db_sample_ids) + ['text'],
            view_name=table_basename)


    def create_combined_contigs_and_splits_tables(self, profile_db_path, table_basename, data_dict):
        """Create a pair of tables in a combined profile database. Contigs and splits tables contain
        the same information since tRNA, unlike a metagenomic contig, is not long enough to be
        split."""
        TablesForViews(profile_db_path).create_new_view(
            data_dict=data_dict,
            table_name=table_basename + '_contigs',
            table_structure=['contig'] + list(itertools.chain(*[(sample_id + '_specific', sample_id + '_nonspecific') for sample_id in self.trnaseq_db_sample_ids])) + ['__parent__'],
            table_types=['text'] + ['numeric'] * 2 * len(self.trnaseq_db_sample_ids) + ['text'],
            view_name=None)
        TablesForViews(profile_db_path).create_new_view(
            data_dict=data_dict,
            table_name=table_basename + '_splits',
            table_structure=['contig'] + list(itertools.chain(*[(sample_id + '_specific', sample_id + '_nonspecific') for sample_id in self.trnaseq_db_sample_ids])) + ['__parent__'],
            table_types=['text'] + ['numeric'] * 2 * len(self.trnaseq_db_sample_ids) + ['text'],
            view_name=table_basename)


def trnaseq_db_loader(input_queue, output_queue, db_converter):
    """This client for `DatabaseConverter.load_trnaseq_db_seq_info` is located outside the
    `DatabaseConverter` class to allow multiprocessing."""
    while True:
        trnaseq_db_path = input_queue.get()
        unmod_norm_seq_summaries, mod_seq_summaries = db_converter.load_trnaseq_db_seq_info(trnaseq_db_path)
        output_queue.put((trnaseq_db_path, unmod_norm_seq_summaries, mod_seq_summaries))
