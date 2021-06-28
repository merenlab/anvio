# -*- coding: utf-8
# pylint: disable=line-too-long
"""Library for tRNA-seq dataset operations

`bin/anvi-trnaseq` and `bin/anvi-convert-trnaseq-database` are the default clients using this
module. `anvi-trnaseq` instantiates a `TRNASeqDataset` object. `anvi-convert-trnaseq-database`
instantiates a `DatabaseConverter` object. The clients call the objects' `process` methods to start
the analytic workflows.

Each sequence library in an experiment is processed separately as a `TRNASeqDataset`, storing an
information-rich anvi'o tRNA-seq database. `DatabaseConverter` finds reference seed sequences from a
set of tRNA-seq databases, storing seeds in an anvi'o contigs database and coverage patterns for
each dataset in anvi'o profile and auxiliary databases. Contigs and profile databases interface with
a range of other tools in the anvi'o platform.


GLOSSARY of essential terms
===========================
Feature: Canonical feature or structural element of tRNA, e.g., anticodon stem

Read: Synonymous with merged paired-end tRNA-seq read oriented 5'->3'

Feature profile (profile): 5'->3' features identified de novo from the 3' end in a merged tRNA-seq
    read

Profiled sequence: Sequence with an assigned feature profile, which may or may not span the whole
    length of the sequence from the 3' end, but which at minimum includes the T arm

Full profile (tRNA profile): Profile that spans (nearly) the full length of the sequence, with a
    small number of unprofiled nucleotides allowed at the 5' end when that number is less than the
    minimum length of a missing next 5' tRNA feature

Truncated profile: Profile that does not span (nearly) the full length of the sequence (e.g., a
    sequence is a chimera of two 3' tRNA fragments and the profile covers the 3' fragment but not
    the unexpected 5' fragment)

Potential modification-induced substitution (sub): Detected as 3-4 different nucleotides at a tRNA
    position, potentially the effect of semi-random nucleotide addition at the site of a modified
    nucleotide during reverse transcription

Potential modification-induced indel (indel): Detected by alignment of tRNA sequences with and
    without potential modification-induced substitutions, indels result from reverse transcriptase
    skipping or adding extra nucleotides due to interaction with a modification (substitutions are
    generally more common than deletions, which in turn are more common than insertions)

Unique sequence (U): Set of dereplicated merged paired-end tRNA-seq reads

Nontemplated nucleotide: Reverse transcription artifact typically added to the 5' end of a read

Trimmed sequence (T): Set of unique sequences that are identical after trimming sequence extensions
    5' of the acceptor stem and 3' of the discriminator nucleotide, e.g., nontemplated 5'
    nucleotides and 3'-CCA acceptor

Normalized sequence (N): The longest of a set of trimmed sequences, with shorter sequences being
    tRNA fragment subsequences

Nonspecific sequence: In contrast to a specific sequence, a trimmed sequence (or its component
    unique sequences and reads) that occurs in multiple normalized sequences (cannot be resolved to
    a single normalized sequence) due to it being a tRNA fragment

Mapped fragment: Sequence without a feature profile that maps to a normalized sequence and may
    include extra nucleotides beyond the trimmed 5' end of the normalized sequence but not
    nucleotides in the trimmed 3' terminus of the normalized sequence

Modified sequence (M): Set of normalized sequences differing by potential modification-induced
    substitutions and optionally indels


ABBREVIATIONS
=============
M: Modified seq

N: Normalized seq
-----------------
Nf: N with full profile
Nc: N with truncated profile
Nb: N with subs but not indels
Nqf: Nf with full profile but no subs (only important as queries in finding indels)
Nq: Nqf or Nc queried against Nb targets to find indels
Ni: N with indels and optionally subs

T: Trimmed seq
--------------
Tp: T with profile
Tf: T with full profile
Tc: T with truncated profile
Tm: Mapped T
Ti: T that is part of Ni (Ti does not necessarily contain indels, though Ni does)
Tip: Ti derived from Tp
Tim: Ti derived from Tm

U: Unique seq
-------------
Un: U not found to have tRNA profile
Up: U with profile
Uc: U with truncated profile
Uf: U with full profile
Us: U with full profile transferred from another Uf
Um: Mapped U
Ui: U that is part of Ti and Ni (Ui does not necessarily contain indels, though Ni does)
Uif: Ui derived from Uf
Uim: Ui derived from Um
"""

import gc
import os
import sys
import time
import random
import shutil
import hashlib
import argparse
import itertools
import numpy as np
import pandas as pd
import pickle as pkl
import multiprocessing as mp

from hashlib import sha1
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
from anvio.drivers.vmatch import Vmatch
from anvio.agglomeration import Agglomerator
from anvio.tables.views import TablesForViews
from anvio.tables.miscdata import TableForLayerOrders
from anvio.sequence import Aligner, AlignedTarget, Cluster, Dereplicator, Kmerizer


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


pp = terminal.pretty_print

THREEPRIME_VARIANTS = constants.THREEPRIME_VARIANTS

ALL_NTS = tuple(constants.nucleotides)
UNAMBIG_NTS = ('A', 'C', 'G', 'T')
NT_INT_DICT = {nt: i for i, nt in enumerate(UNAMBIG_NTS, start=1)}
INT_NT_DICT = {i: nt for i, nt in enumerate(UNAMBIG_NTS, start=1)}
# NUM_NT_BINS is used in counting the number of distinct nucleotides (1, 2, 3, or 4) at positions in
# internally ungapped alignments: there is one bin (value 0) for end gaps in the alignment.
NUM_NT_BINS = len(UNAMBIG_NTS) + 1

ANTICODON_TO_AA = constants.anticodon_to_AA


class UniqueSequence(object):
    """A dereplicated tRNA-seq read."""

    __slots__ = ('seq_string', 'represent_name', 'read_count')

    def __init__(self, seq_string, represent_name, read_count):
        self.seq_string = seq_string
        self.represent_name = represent_name
        self.read_count = read_count


class UniqueProfileSequence(UniqueSequence):
    """A tRNA feature profile, either full or truncated, was assigned to the sequence."""

    __slots__ = (
        'feature_start_indices',
        'feature_stop_indices',
        'has_complete_feature_set',
        'has_his_g',
        'alpha_start_index',
        'alpha_stop_index',
        'beta_start_index',
        'beta_stop_index',
        'anticodon_string',
        'anticodon_aa',
        'contains_anticodon',
        'threeprime_terminus_length',
        'num_conserved',
        'num_unconserved',
        'num_paired',
        'num_unpaired',
        'unconserved_info',
        'unpaired_info',
        'profiled_seq_length',
        'trimmed_seq_represent_name'
    )

    def __init__(self, seq_string, represent_name, read_count, profile):
        super().__init__(seq_string, represent_name, read_count)

        self.feature_start_indices = [feature.start_pos if hasattr(feature, 'start_pos') else feature.start_positions
                                      for feature in profile.features]
        self.feature_stop_indices = [feature.stop_pos if hasattr(feature, 'stop_pos') else feature.stop_positions
                                     for feature in profile.features]
        self.alpha_start_index = profile.alpha_start
        self.alpha_stop_index = profile.alpha_stop
        self.beta_start_index = profile.beta_start
        self.beta_stop_index = profile.beta_stop
        self.anticodon_string = anticodon = profile.anticodon_seq
        self.anticodon_aa = profile.anticodon_aa if profile.anticodon_aa else None
        self.contains_anticodon = True if anticodon else False
        self.threeprime_terminus_length = len(profile.threeprime_terminus_seq)
        self.num_conserved = profile.num_conserved
        self.num_unconserved = profile.num_unconserved
        self.num_paired = profile.num_paired
        self.num_unpaired = profile.num_unpaired
        self.unconserved_info = profile.unconserved_info
        self.unpaired_info = profile.unpaired_info
        self.profiled_seq_length = len(profile.profiled_seq)
        self.trimmed_seq_represent_name = None


class UniqueFullProfileSequence(UniqueProfileSequence):
    """A full tRNA feature profile was assigned to the sequence."""

    __slots__ = (
        'has_complete_feature_set',
        'has_his_g',
        'num_extrapolated_fiveprime_nts',
        'extra_fiveprime_length'
    )

    def __init__(self, seq_string, represent_name, read_count, profile):
        super().__init__(seq_string, represent_name, read_count, profile)

        self.has_complete_feature_set = profile.has_complete_feature_set
        self.has_his_g = True if profile.features[0].name == 'tRNA-His position 0' else False
        self.num_extrapolated_fiveprime_nts = profile.num_in_extrapolated_fiveprime_feature
        self.extra_fiveprime_length = 0 if profile.num_extra_fiveprime is None else profile.num_extra_fiveprime


class UniqueTruncatedProfileSequence(UniqueProfileSequence):
    """A truncated tRNA feature profile was assigned to the sequence."""

    __slots__ = ('trunc_profile_index', )

    def __init__(self, seq_string, represent_name, read_count, profile):
        super().__init__(seq_string, represent_name, read_count, profile)

        self.trunc_profile_index = profile.trunc_profile_index


class UniqueTransferredProfileSequence(UniqueFullProfileSequence):
    """This object is generated as part of the determination of normalized tRNA sequences from
    trimmed tRNA sequences. This type of sequence is produced in the special circumstance that the
    profile of a shorter sequence is transferred to a longer sequence, because the longer sequence
    was originally found to have a complete profile, but a shorter 3' subsequence also had a
    complete profile; so, parsimoniously, the profile of the shorter sequence was transferred to the
    longer, and the additional 5' nucleotides of the longer reclassified as extra nucleotides beyond
    the 5' end of a mature tRNA sequence."""

    __slots__ = ('defunct_uniq_seq', )

    def __init__(self, defunct_uniq_seq, replacement_info_dict):
        UniqueSequence.__init__(self, defunct_uniq_seq.seq_string, defunct_uniq_seq.represent_name, defunct_uniq_seq.read_count)

        uniq_seq_string = defunct_uniq_seq.seq_string
        uniq_seq_length = len(uniq_seq_string)
        trimmed_seq_stop_in_uniq_seq = uniq_seq_length - uniq_seq_string[::-1].index(replacement_info_dict['trimmed_seq_string'][::-1])
        feature_start_indices = []
        for feature_start_index_from_trimmed_threeprime in replacement_info_dict['feature_start_indices_from_trimmed_threeprime']:
            if isinstance(feature_start_index_from_trimmed_threeprime, int):
                feature_start_indices.append(trimmed_seq_stop_in_uniq_seq + feature_start_index_from_trimmed_threeprime)
            else:
                feature_start_indices.append(tuple([trimmed_seq_stop_in_uniq_seq + strand_start_index_from_trimmed_threeprime for strand_start_index_from_trimmed_threeprime in feature_start_index_from_trimmed_threeprime]))
        self.feature_start_indices = feature_start_indices
        feature_stop_indices = []
        for feature_stop_index_from_trimmed_threeprime in replacement_info_dict['feature_stop_indices_from_trimmed_threeprime']:
            if isinstance(feature_stop_index_from_trimmed_threeprime, int):
                feature_stop_indices.append(trimmed_seq_stop_in_uniq_seq + feature_stop_index_from_trimmed_threeprime)
            else:
                feature_stop_indices.append(tuple([trimmed_seq_stop_in_uniq_seq + strand_stop_index_from_trimmed_threeprime for strand_stop_index_from_trimmed_threeprime in feature_stop_index_from_trimmed_threeprime]))
        self.feature_stop_indices = feature_stop_indices
        self.has_complete_feature_set = True
        self.num_extrapolated_fiveprime_nts = 0
        self.has_his_g = replacement_info_dict['has_his_g']
        self.alpha_start_index = (None if replacement_info_dict['alpha_start_index_from_trimmed_threeprime'] is None
                                  else trimmed_seq_stop_in_uniq_seq + replacement_info_dict['alpha_start_index_from_trimmed_threeprime'])
        self.alpha_stop_index = (None if replacement_info_dict['alpha_stop_index_from_trimmed_threeprime'] is None
                                 else trimmed_seq_stop_in_uniq_seq - replacement_info_dict['alpha_stop_index_from_trimmed_threeprime'])
        self.beta_start_index = (None if replacement_info_dict['beta_start_index_from_trimmed_threeprime'] is None
                                 else trimmed_seq_stop_in_uniq_seq - replacement_info_dict['beta_start_index_from_trimmed_threeprime'])
        self.beta_stop_index = (None if replacement_info_dict['beta_stop_index_from_trimmed_threeprime'] is None
                                else trimmed_seq_stop_in_uniq_seq - replacement_info_dict['beta_stop_index_from_trimmed_threeprime'])
        self.anticodon_string = replacement_info_dict['anticodon_string']
        self.anticodon_aa = replacement_info_dict['anticodon_aa']
        self.contains_anticodon = replacement_info_dict['contains_anticodon']
        self.threeprime_terminus_length = uniq_seq_length - trimmed_seq_stop_in_uniq_seq
        self.num_conserved = replacement_info_dict['num_conserved']
        self.num_unconserved = replacement_info_dict['num_unconserved']
        self.num_paired = replacement_info_dict['num_paired']
        self.num_unpaired = replacement_info_dict['num_unpaired']
        unconserved_info = []
        for unconserved_tuple in replacement_info_dict['unconserved_info_from_trimmed_threeprime']:
            unconserved_info.append((trimmed_seq_stop_in_uniq_seq + unconserved_tuple[0],
                                     unconserved_tuple[1],
                                     unconserved_tuple[2]))
        self.unconserved_info = unconserved_info
        unpaired_info = []
        for unpaired_tuple in replacement_info_dict['unpaired_info_from_trimmed_threeprime']:
            unpaired_info.append((trimmed_seq_stop_in_uniq_seq + unpaired_tuple[0],
                                  trimmed_seq_stop_in_uniq_seq + unpaired_tuple[1],
                                  unpaired_tuple[2],
                                  unpaired_tuple[3]))
        self.unpaired_info = unpaired_info
        self.profiled_seq_length = replacement_info_dict['profiled_seq_without_terminus_length'] + self.threeprime_terminus_length
        self.extra_fiveprime_length = uniq_seq_length - self.profiled_seq_length
        self.trimmed_seq_represent_name = None

        # Store the defunct profile information for posterity.
        self.defunct_uniq_seq = defunct_uniq_seq


class UniqueMappedSequence(UniqueSequence):
    """This object is generated in the identification of tRNA fragments by mapping."""

    __slots__ = ('extra_fiveprime_length', 'trimmed_seq_represent_name')

    def __init__(self, seq_string, represent_name, read_count, extra_fiveprime_length=0):
        super().__init__(seq_string, represent_name, read_count)

        self.extra_fiveprime_length = extra_fiveprime_length
        self.trimmed_seq_represent_name = None


class TrimmedSequence(object):
    """A tRNA sequence with bases trimmed 5' of the acceptor stem (or 5'-G in the case of tRNA-His)
    and 3' of the discriminator.

    The purpose of trimming is to collapse non-biological variability prevalent at the ends of
    reads.

    EXAMPLE 1:
    E. coli tRNA-Ala-GGC-1-1
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA

    This collapses to the following trimmed sequence, removing the 3' terminus (the acceptor happens
    to be genomic rather than post-transcriptionally added in this example, but it doesn't matter):
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA

    Examples of possible profiled reads that collapse to this sequence:
    AGGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCAC

    EXAMPLE 2:
    3' fragment of the same tRNA, ending in 3'-CC rather than canonical 3'-CCA
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACC

    This collapses to the following trimmed sequence, removing 3'-CC:
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA
    """

    # The user can specify what defines a long (biological vs. non-templated) 5' extension. This is
    # set by `TRNASeqDataset`.
    min_length_of_long_fiveprime_extension = 4

    __slots__ = ('seq_string', 'read_count', 'uniq_seqs', 'norm_seq_represent_names')

    def __init__(self, seq_string, uniq_seqs):
        self.seq_string = seq_string
        self.uniq_seqs = uniq_seqs
        for uniq_seq in uniq_seqs:
            if uniq_seq.trimmed_seq_represent_name is not None:
                raise ConfigError(f"The unique sequence with the representative name {uniq_seq.represent_name} "
                                  f"was already assigned to a trimmed sequence with the representative name {uniq_seq.trimmed_seq_represent_name} "
                                  "and so cannot be assigned to a new trimmed sequence.")
        self.read_count = sum([uniq_seq.read_count for uniq_seq in self.uniq_seqs])
        self.norm_seq_represent_names = []


class TrimmedFullProfileSequence(TrimmedSequence):
    """This object is formed from sequences with a full tRNA feature profile."""

    __slots__ = (
        'represent_uniq_seq',
        'represent_name',
        'feature_start_indices',
        'feature_stop_indices',
        'contains_anticodon',
        'read_threeprime_terminus_count_dict',
        'has_complete_feature_set',
        'has_his_g',
        'num_extrapolated_fiveprime_nts',
        'uniq_with_extra_fiveprime_count',
        'read_with_extra_fiveprime_count',
        'long_fiveprime_extension_dict'
    )

    def __init__(self, seq_string, uniq_seqs):
        super().__init__(seq_string, uniq_seqs)

        # The representative unique sequence is chosen as follows:
        # 1. Most abundant full-length tRNA without extra 5' bases, ignoring the 3' terminus sequence, OR
        # 2. Most abundant full-length tRNA with extra 5' bases, OR
        # 3. Most abundant 3' tRNA fragment
        # Sort unique sequences such that the first sequence is the most abundant+longest and the
        # last is the least abundant+shortest.
        uniq_seqs.sort(key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, -uniq_seq.read_count, uniq_seq.represent_name))
        if uniq_seqs[0].extra_fiveprime_length > 0:
            if uniq_seqs[-1].extra_fiveprime_length == 0:
                # If the first unique sequence has extra 5' nucleotides and the last has none, then
                # the last sequence and others without extra 5' nucleotides must be a full-length
                # tRNA (ignoring the 3' terminus). Therefore, select the most abundant of these
                # full-length tRNAs with extra 5' nucleotides as the representative sequence.
                self.represent_uniq_seq = represent_uniq_seq = sorted(uniq_seqs, key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, uniq_seq.read_count, uniq_seq.represent_name))[-1]
            else:
                self.represent_uniq_seq = represent_uniq_seq = uniq_seqs[0]
        else:
            # In this case, ALL unique sequences are EITHER full-length tRNA OR a 3' tRNA fragment.
            self.represent_uniq_seq = represent_uniq_seq = uniq_seqs[0]
        self.represent_name = represent_name = represent_uniq_seq.represent_name
        for uniq_seq in uniq_seqs:
            uniq_seq.trimmed_seq_represent_name = represent_name

        # Assume that the feature profile indices of the representative unique sequence are the same
        # as the other unique sequences. The 3' terminus is the last feature in the profile and not
        # part of the trimmed sequence.
        self.feature_start_indices = represent_uniq_seq.feature_start_indices[: -1] if represent_uniq_seq.feature_start_indices else None
        self.feature_stop_indices = represent_uniq_seq.feature_stop_indices[: -1] if represent_uniq_seq.feature_stop_indices else None
        self.contains_anticodon = represent_uniq_seq.contains_anticodon
        self.has_complete_feature_set = represent_uniq_seq.has_complete_feature_set
        self.has_his_g = represent_uniq_seq.has_his_g
        self.num_extrapolated_fiveprime_nts = represent_uniq_seq.num_extrapolated_fiveprime_nts

        read_threeprime_terminus_count_dict = defaultdict(int)
        for uniq_seq in self.uniq_seqs:
            if uniq_seq.threeprime_terminus_length:
                read_threeprime_terminus_count_dict[uniq_seq.seq_string[-uniq_seq.threeprime_terminus_length: ]] += uniq_seq.read_count
            else:
                read_threeprime_terminus_count_dict[''] += uniq_seq.read_count
        self.read_threeprime_terminus_count_dict = read_threeprime_terminus_count_dict

        uniq_with_extra_fiveprime_count = 0
        read_with_extra_fiveprime_count = 0
        # Find the number of reads containing each unique 5' extension.
        long_fiveprime_extension_dict = {}
        for uniq_seq in self.uniq_seqs:
            if uniq_seq.extra_fiveprime_length:
                uniq_with_extra_fiveprime_count += 1
                read_with_extra_fiveprime_count += uniq_seq.read_count
                if uniq_seq.extra_fiveprime_length >= self.min_length_of_long_fiveprime_extension:
                    long_fiveprime_extension_dict[uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length]] = uniq_seq.read_count
        self.uniq_with_extra_fiveprime_count = uniq_with_extra_fiveprime_count
        self.read_with_extra_fiveprime_count = read_with_extra_fiveprime_count
        self.long_fiveprime_extension_dict = long_fiveprime_extension_dict


class TrimmedTruncatedProfileSequence(TrimmedSequence):
    """This object is formed from sequences with a truncated tRNA feature profile."""

    __slots__ = (
        'represent_uniq_seq',
        'represent_name',
        'feature_start_indices',
        'feature_stop_indices',
        'contains_anticodon',
        'read_threeprime_terminus_count_dict',
        'trunc_profile_index'
    )

    def __init__(self, seq_string, uniq_seqs):
        super().__init__(seq_string, uniq_seqs)

        # Make the most abundant unique sequence the representative sequence.
        uniq_seqs.sort(key=lambda uniq_seq: (-uniq_seq.read_count, uniq_seq.represent_name))
        represent_uniq_seq = uniq_seqs[0]
        self.represent_name = represent_name = represent_uniq_seq.represent_name
        for uniq_seq in uniq_seqs:
            uniq_seq.trimmed_seq_represent_name = represent_name

        # Assume that the feature profile indices of the representative unique sequence are the same
        # as the other unique sequences. The 3' terminus is the last feature in the profile and not
        # part of the trimmed sequence.
        self.feature_start_indices = represent_uniq_seq.feature_start_indices[: -1] if represent_uniq_seq.feature_start_indices else None
        self.feature_stop_indices = represent_uniq_seq.feature_stop_indices[: -1] if represent_uniq_seq.feature_stop_indices else None
        self.contains_anticodon = represent_uniq_seq.contains_anticodon

        read_threeprime_terminus_count_dict = defaultdict(int)
        for uniq_seq in self.uniq_seqs:
            if uniq_seq.threeprime_terminus_length:
                read_threeprime_terminus_count_dict[uniq_seq.seq_string[-uniq_seq.threeprime_terminus_length: ]] += uniq_seq.read_count
            else:
                read_threeprime_terminus_count_dict[''] += uniq_seq.read_count
        self.read_threeprime_terminus_count_dict = read_threeprime_terminus_count_dict

        self.trunc_profile_index = represent_uniq_seq.trunc_profile_index


class TrimmedMappedSequence(TrimmedSequence):
    """This object is formed from a single unique sequence (`UniqueMappedSequence`) in the process
    of mapping unique unprofiled sequences to normalized sequences. It is not like the other trimmed
    sequence objects. Its purpose is to be one of the trimmed sequence objects added to a normalized
    sequence, however no 5' bases are trimmed from the unique sequence string in creating the
    trimmed sequence string (and mapped sequences do not have extra 3' bases). The reason for this
    is that the 5' extension may represent all but a small number of nucleotides in the sequence, so
    it is best not to dereplicate mapped sequences identical in the non-5' section by lumping them
    together as the same trimmed sequence."""

    __slots__ = (
        'represent_name',
        'uniq_with_extra_fiveprime_count',
        'read_with_extra_fiveprime_count',
        'long_fiveprime_extension_dict'
    )

    def __init__(self, uniq_seq):
        super().__init__(uniq_seq.seq_string, [uniq_seq])

        self.represent_name = uniq_seq.represent_name
        uniq_seq.trimmed_seq_represent_name = self.represent_name

        extra_fiveprime_length = uniq_seq.extra_fiveprime_length
        self.uniq_with_extra_fiveprime_count = 1 if extra_fiveprime_length else 0
        self.read_with_extra_fiveprime_count = uniq_seq.read_count if extra_fiveprime_length else 0
        self.long_fiveprime_extension_dict = {}
        if extra_fiveprime_length >= self.min_length_of_long_fiveprime_extension:
            self.long_fiveprime_extension_dict[uniq_seq.seq_string[: extra_fiveprime_length]] = uniq_seq.read_count


class NormalizedSequence(object):
    """A tRNA sequence that can contain shorter tRNA fragment subsequences.

    Normalized sequences are derived from trimmed sequences. Trimmed profiled tRNA are first
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
        'start_positions',
        'stop_positions',
        'specific_read_count',
        'nonspecific_read_count',
        'specific_covs',
        'nonspecific_covs',
        'mean_specific_cov',
        'mean_nonspecific_cov',
        'mod_seqs'
    )

    def __init__(self, trimmed_seqs, start_positions=None, stop_positions=None, skip_init=False):
        self.trimmed_seqs = trimmed_seqs
        represent_trimmed_seq = trimmed_seqs[0]
        self.represent_name = represent_name = represent_trimmed_seq.represent_name
        for trimmed_seq in trimmed_seqs:
            trimmed_seq.norm_seq_represent_names.append(represent_name)
        self.seq_string = represent_trimmed_seq.seq_string

        if start_positions and stop_positions:
            self.start_positions = start_positions
            self.stop_positions = stop_positions
        elif (not start_positions) and (not stop_positions):
            # This is what occurs in `anvi-trnaseq` for the instantiation of
            # `NormalizedFullProfileSequence` and `NormalizedTruncatedProfileSequence` objects. All
            # trimmed sequences provided as input are aligned from the 3' end. Mapped tRNA sequences
            # that can be aligned to other places in the normalized sequence are added later.
            norm_seq_length = len(self.seq_string)
            self.start_positions = [norm_seq_length - len(trimmed_seq.seq_string) for trimmed_seq in trimmed_seqs]
            self.stop_positions = [norm_seq_length] * len(trimmed_seqs)
        else:
            self.start_positions = None
            self.stop_positions = None

        # It is useful to know which modified sequences, if any, contain this normalized sequence. A
        # normalized sequence with deletions, unlike a normalized sequence without deletions, can
        # theoretically be assigned to more than one modified sequence.
        self.mod_seqs = []

        if skip_init:
            self.specific_read_count = None
            self.nonspecific_read_count = None
            self.specific_covs = None
            self.nonspecific_covs = None
            self.mean_specific_cov = None
            self.mean_nonspecific_cov = None
        else:
            self.init()


    def init(self):
        """Set attributes for a "finalized" set of input trimmed sequences, including those that
        were added after object instantiation, such as mapped sequences."""
        # Specific reads are those that are only contained in this normalized sequence.
        specific_read_count = 0
        nonspecific_read_count = 0
        specific_covs = np.zeros(len(self.seq_string), dtype=int)
        nonspecific_covs = np.zeros(len(self.seq_string), dtype=int)

        for trimmed_seq, start_pos, stop_pos in zip(self.trimmed_seqs, self.start_positions, self.stop_positions):
            if len(trimmed_seq.norm_seq_represent_names) == 1:
                specific_read_count += trimmed_seq.read_count
                specific_covs[start_pos: stop_pos] += trimmed_seq.read_count
            else:
                nonspecific_read_count += trimmed_seq.read_count
                nonspecific_covs[start_pos: stop_pos] += trimmed_seq.read_count

        self.specific_read_count = specific_read_count
        self.nonspecific_read_count = nonspecific_read_count
        self.specific_covs = specific_covs
        self.nonspecific_covs = nonspecific_covs
        self.mean_specific_cov = specific_covs.mean()
        self.mean_nonspecific_cov = nonspecific_covs.mean()


class NormalizedFullProfileSequence(NormalizedSequence):
    """This object is formed from `TrimmedFullProfileSequence` objects and, optionally,
    `TrimmedTruncatedProfileSequence` objects that are subsequences of the seed."""

    __slots__ = (
        'feature_start_indices',
        'feature_stop_indices',
        'has_complete_feature_set',
        'has_his_g',
        'specific_read_with_extra_fiveprime_count',
        'nonspecific_read_with_extra_fiveprime_count',
        'specific_mapped_read_count',
        'nonspecific_mapped_read_count',
        'specific_long_fiveprime_extension_dict',
        'nonspecific_long_fiveprime_extension_dict',
        'specific_read_threeprime_terminus_count_dict',
        'nonspecific_read_threeprime_terminus_count_dict'
    )

    def __init__(self, trimmed_seqs, start_positions=None, stop_positions=None):
        super().__init__(trimmed_seqs, start_positions=None, stop_positions=None, skip_init=True)

        represent_trimmed_seq = self.trimmed_seqs[0]
        self.feature_start_indices = represent_trimmed_seq.feature_start_indices
        self.feature_stop_indices = represent_trimmed_seq.feature_stop_indices
        self.has_complete_feature_set = represent_trimmed_seq.has_complete_feature_set
        self.has_his_g = represent_trimmed_seq.has_his_g

        self.specific_read_with_extra_fiveprime_count = None
        self.nonspecific_read_with_extra_fiveprime_count = None
        self.specific_mapped_read_count = None
        self.nonspecific_mapped_read_count = None
        self.specific_long_fiveprime_extension_dict = None
        self.nonspecific_long_fiveprime_extension_dict = None


    def init(self):
        """Set attributes for a "finalized" set of input trimmed sequences."""
        super().init()

        specific_mapped_read_count = 0
        nonspecific_mapped_read_count = 0
        specific_read_with_extra_fiveprime_count = 0
        nonspecific_read_with_extra_fiveprime_count = 0
        specific_long_fiveprime_extension_dict = defaultdict(int)
        nonspecific_long_fiveprime_extension_dict = defaultdict(int)
        specific_read_threeprime_terminus_count_dict = defaultdict(int)
        nonspecific_read_threeprime_terminus_count_dict = defaultdict(int)

        for trimmed_seq in self.trimmed_seqs:
            if len(trimmed_seq.norm_seq_represent_names) == 1:
                if isinstance(trimmed_seq, TrimmedMappedSequence):
                    specific_mapped_read_count += trimmed_seq.read_count
                    specific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count
                    for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                        specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count
                else:
                    for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                        specific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count

                    if not isinstance(trimmed_seq, TrimmedTruncatedProfileSequence):
                        specific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count
                        for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                            specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count
            else:
                if isinstance(trimmed_seq, TrimmedMappedSequence):
                    nonspecific_mapped_read_count += trimmed_seq.read_count
                    nonspecific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count
                else:
                    for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                        nonspecific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count

                    # Nonspecific, unlike specific, profiled sequences can have 5' sequence
                    # extensions, as profiled sequences with 5' extensions would span the length of
                    # the normalized sequence and would thus be specific to it.
                    if not isinstance(trimmed_seq, TrimmedTruncatedProfileSequence):
                        nonspecific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count
                        for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                            nonspecific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

        self.specific_mapped_read_count = specific_mapped_read_count
        self.nonspecific_mapped_read_count = nonspecific_mapped_read_count
        self.specific_read_with_extra_fiveprime_count = specific_read_with_extra_fiveprime_count
        self.nonspecific_read_with_extra_fiveprime_count = nonspecific_read_with_extra_fiveprime_count
        self.specific_long_fiveprime_extension_dict = specific_long_fiveprime_extension_dict
        self.nonspecific_long_fiveprime_extension_dict = nonspecific_long_fiveprime_extension_dict
        self.specific_read_threeprime_terminus_count_dict = specific_read_threeprime_terminus_count_dict
        self.nonspecific_read_threeprime_terminus_count_dict = nonspecific_read_threeprime_terminus_count_dict


class NormalizedTruncatedProfileSequence(NormalizedSequence):
    """This object is formed exclusively from `TrimmedTruncatedProfileSequence` objects."""

    __slots__ = (
        'feature_start_indices',
        'feature_stop_indices',
        'trunc_profile_index',
        'specific_read_threeprime_terminus_count_dict',
        'nonspecific_read_threeprime_terminus_count_dict'
    )

    def __init__(self, trimmed_seqs, start_positions=None, stop_positions=None):
        super().__init__(trimmed_seqs, start_positions=None, stop_positions=None, skip_init=False)

        represent_trimmed_seq = self.trimmed_seqs[0]
        self.feature_start_indices = represent_trimmed_seq.feature_start_indices
        self.feature_stop_indices = represent_trimmed_seq.feature_stop_indices
        self.trunc_profile_index = represent_trimmed_seq.trunc_profile_index

        specific_read_threeprime_terminus_count_dict = defaultdict(int)
        nonspecific_read_threeprime_terminus_count_dict = defaultdict(int)
        for trimmed_seq in self.trimmed_seqs:
            if len(trimmed_seq.norm_seq_represent_names) == 1:
                for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                    specific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count
            else:
                for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                    nonspecific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count
        self.specific_read_threeprime_terminus_count_dict = specific_read_threeprime_terminus_count_dict
        self.nonspecific_read_threeprime_terminus_count_dict = nonspecific_read_threeprime_terminus_count_dict


class NormalizedDeletionSequence(NormalizedSequence):
    """This object is generated in the identification of tRNAs with deletions.
    `NormalizedFullProfileSequence` and `NormalizedTruncatedProfileSequence` objects can be found to
    contain deletions, invalidating their existing feature profile attributes. This object is
    created in their stead. Note that the input defunct normalized sequence can contain
    `TrimmedMappedSequence` objects in addition to `TrimmedFullProfileSequence` and
    `TrimmedTruncatedProfileSequence` objects."""

    __slots__ = (
        'defunct_norm_seqs',
        'mod_seq_del_configs',
        'start_positions_in_mod_seqs',
        'norm_seq_del_configs',
        'mod_seqs_contain_anticodon',
        'specific_del_covs_in_mod_seqs',
        'nonspecific_del_covs_in_mod_seqs',
        'specific_read_with_extra_fiveprime_count',
        'nonspecific_read_with_extra_fiveprime_count',
        'specific_mapped_read_count',
        'nonspecific_mapped_read_count',
        'specific_long_fiveprime_extension_dict',
        'nonspecific_long_fiveprime_extension_dict',
        'specific_read_threeprime_terminus_count_dict',
        'nonspecific_read_threeprime_terminus_count_dict'
    )

    # The user can specify what defines a long (biological vs. non-templated) 5' extension. This is
    # set by `TRNASeqDataset`.
    min_length_of_long_fiveprime_extension = 4

    def __init__(self, seq_string, defunct_norm_seq, mod_seq, del_config, start_pos_in_mod_seq, norm_seq_del_config, contains_anticodon):
        self.seq_string = seq_string
        self.trimmed_seqs = defunct_norm_seq.trimmed_seqs
        self.defunct_norm_seqs = [defunct_norm_seq for _ in self.trimmed_seqs]
        # Modified sequences are found before normalized sequences with deletions, and so can be
        # recorded upon instantiation of this object. Currently, the workflow only recognizes
        # normalized sequences with deletions that derive from a single modified sequence. If
        # nonspecific assignment to modified sequences were allowed, then each modified sequence can
        # have different deletion configurations, the normalized sequence can start at different
        # places in the modified sequence, and, given the modified sequence, the normalized sequence
        # may or may not contain the anticodon.
        self.mod_seqs = [mod_seq]
        self.mod_seq_del_configs = [del_config]
        self.start_positions_in_mod_seqs = [start_pos_in_mod_seq]
        self.norm_seq_del_configs = [norm_seq_del_config]
        self.mod_seqs_contain_anticodon = [contains_anticodon]


    def init(self):
        """Set attributes for a "finalized" set of input trimmed sequences."""
        # Here is an important note regarding the trimmed sequences comprising this normalized
        # sequence. The feature profiles of sequences found to have deletions are invalidated, but
        # trimmed (and unique) sequence profiles are not changed. A key reason for this, beside
        # convenience, is that a trimmed sequence may be non-specific -- also assigned to another
        # normalized sequence -- not necessarily with a deletion in the same place. Ideally,
        # additional information on each trimmed sequence would be recorded as attributes of this
        # normalized sequence. For example, each trimmed sequence may have a different number of
        # extra 5' nucleotides in this normalized sequence. Likewise, the unique sequences and reads
        # comprising the trimmed sequence would have a different number of extra 5' nucleotides.

        # If there is a trimmed profiled sequence equivalent to the normalized sequence with
        # deletions, adopt its representative name. Otherwise, use the representative name of the
        # longest trimmed sequence.
        norm_seq_length = len(self.seq_string)
        max_trimmed_seq_length = 0
        for trimmed_seq in self.trimmed_seqs:
            if isinstance(trimmed_seq, TrimmedMappedSequence):
                continue
            trimmed_seq_length = len(trimmed_seq.seq_string)
            if trimmed_seq_length == norm_seq_length:
                represent_name = trimmed_seq.represent_name
                break
            if trimmed_seq_length > max_trimmed_seq_length:
                max_trimmed_seq_length = trimmed_seq_length
                represent_name = trimmed_seq.represent_name
        self.represent_name = represent_name

        derep_defunct_norm_seqs = []
        for norm_seq in self.defunct_norm_seqs:
            for other_norm_seq in derep_defunct_norm_seqs:
                if norm_seq.represent_name == other_norm_seq.represent_name:
                    break
            else:
                derep_defunct_norm_seqs.append(norm_seq)

        # Affiliate each trimmed sequence with the present rather than the defunct normalized
        # sequence(s). A trimmed sequence may be part of more than one defunct normalized sequence
        # converted into the present normalized sequence with deletions.
        derep_defunct_norm_seq_names = [norm_seq.represent_name for norm_seq in derep_defunct_norm_seqs]
        trimmed_seq.norm_seq_represent_names = [name for name in trimmed_seq.norm_seq_represent_names if name not in derep_defunct_norm_seq_names] + [represent_name]

        # Determine attributes of this object.
        start_positions = []
        stop_positions = []
        specific_read_count = 0
        nonspecific_read_count = 0
        specific_covs = np.zeros(norm_seq_length, dtype=int)
        nonspecific_covs = np.zeros(norm_seq_length, dtype=int)
        norm_seq_del_configs = self.norm_seq_del_configs
        specific_del_covs_in_mod_seqs = [np.zeros(len(norm_seq_del_config), dtype=int) for norm_seq_del_config in norm_seq_del_configs]
        nonspecific_del_covs_in_mod_seqs = [np.zeros(len(norm_seq_del_config), dtype=int) for norm_seq_del_config in norm_seq_del_configs]
        specific_read_with_extra_fiveprime_count = 0
        nonspecific_read_with_extra_fiveprime_count = 0
        specific_mapped_read_count = 0
        nonspecific_mapped_read_count = 0
        specific_long_fiveprime_extension_dict = defaultdict(int)
        nonspecific_long_fiveprime_extension_dict = defaultdict(int)
        specific_read_threeprime_terminus_count_dict = defaultdict(int)
        nonspecific_read_threeprime_terminus_count_dict = defaultdict(int)
        for trimmed_seq, defunct_norm_seq in zip(self.trimmed_seqs, self.defunct_norm_seqs):
            trimmed_seq_length = len(trimmed_seq.seq_string)
            trimmed_seq_read_count = trimmed_seq.read_count
            if isinstance(trimmed_seq, TrimmedMappedSequence):
                # Find where the mapped sequence was located in the defunct normalized sequence.
                trimmed_seq_name = trimmed_seq.represent_name
                for trimmed_seq_index, other_trimmed_seq in enumerate(defunct_norm_seq.trimmed_seqs):
                    if trimmed_seq_name == other_trimmed_seq.represent_name:
                        defunct_start_pos = defunct_norm_seq.start_positions[trimmed_seq_index]
                        defunct_stop_pos = defunct_norm_seq.stop_positions[trimmed_seq_index]
                        break
                # If the normalized sequence with deletions has extra 5' bases relative to the
                # defunct normalized sequence, then start and stop positions of the trimmed
                # sequences decrease.
                pos_shift = len(defunct_norm_seq.seq_string) - norm_seq_length
                if defunct_start_pos + pos_shift <= 0:
                    start_pos = 0
                else:
                    start_pos = defunct_start_pos + pos_shift
                start_positions.append(start_pos)
                stop_pos = defunct_stop_pos + pos_shift
                stop_positions.append(stop_pos)

                if len(trimmed_seq.norm_seq_represent_names) == 1:
                    specific_read_count += trimmed_seq_read_count
                    specific_covs[start_pos: stop_pos] += trimmed_seq_read_count
                    for specific_del_covs, norm_seq_del_config in zip(specific_del_covs_in_mod_seqs, norm_seq_del_configs):
                        for del_index, adjacent_fiveprime_pos in enumerate(norm_seq_del_config):
                            # A trimmed sequence covering the deletion must start at or before the
                            # 5' nucleotide adjacent to the deletion and a stop position at or after
                            # the 3' nucleotide adjacent to the deletion.
                            if start_pos <= adjacent_fiveprime_pos and stop_pos >= adjacent_fiveprime_pos + 1:
                                specific_del_covs[del_index] += trimmed_seq_read_count
                    specific_mapped_read_count += trimmed_seq_read_count
                    if defunct_start_pos + pos_shift < 0:
                        specific_read_with_extra_fiveprime_count += trimmed_seq_read_count
                        if -defunct_start_pos - pos_shift >= self.min_length_of_long_fiveprime_extension:
                            specific_long_fiveprime_extension_dict[trimmed_seq.seq_string] = trimmed_seq_read_count
                else:
                    nonspecific_read_count += trimmed_seq_read_count
                    nonspecific_covs[start_pos: stop_pos] += trimmed_seq_read_count
                    for nonspecific_del_covs, norm_seq_del_config in zip(nonspecific_del_covs_in_mod_seqs, norm_seq_del_configs):
                        for del_index, adjacent_fiveprime_pos in enumerate(norm_seq_del_config):
                            if start_pos <= adjacent_fiveprime_pos and stop_pos >= adjacent_fiveprime_pos + 1:
                                nonspecific_del_covs[del_index] += trimmed_seq_read_count
                    nonspecific_mapped_read_count += trimmed_seq_read_count
                    if defunct_start_pos + pos_shift < 0:
                        nonspecific_read_with_extra_fiveprime_count += trimmed_seq_read_count
                        if -defunct_start_pos - pos_shift >= self.min_length_of_long_fiveprime_extension:
                            nonspecific_long_fiveprime_extension_dict[trimmed_seq.seq_string] = trimmed_seq_read_count
                continue

            # Profiled trimmed sequences are still anchored to the 3' end of the normalized
            # sequence.
            if trimmed_seq_length >= norm_seq_length:
                start_pos = 0
                start_positions.append(start_pos)
            else:
                start_pos = norm_seq_length - trimmed_seq_length
                start_positions.append(start_pos)
            stop_pos = norm_seq_length
            stop_positions.append(stop_pos)

            if len(trimmed_seq.norm_seq_represent_names) == 1:
                specific_read_count += trimmed_seq_read_count
                specific_covs[start_pos: stop_pos] += trimmed_seq_read_count
                for specific_del_covs, norm_seq_del_config in zip(specific_del_covs_in_mod_seqs, norm_seq_del_configs):
                    for del_index, adjacent_fiveprime_pos in enumerate(norm_seq_del_config):
                        if start_pos <= adjacent_fiveprime_pos and stop_pos >= adjacent_fiveprime_pos + 1:
                            specific_del_covs[del_index] += trimmed_seq_read_count
                if trimmed_seq_length > norm_seq_length:
                    specific_read_with_extra_fiveprime_count += trimmed_seq_read_count
                    if trimmed_seq_length - norm_seq_length >= self.min_length_of_long_fiveprime_extension:
                        specific_long_fiveprime_extension_dict[trimmed_seq.seq_string] = trimmed_seq_read_count
                for threeprime_terminus_seq_string, threeprime_terminus_read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                    specific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += threeprime_terminus_read_count
            else:
                nonspecific_read_count += trimmed_seq_read_count
                nonspecific_covs[start_pos: stop_pos] += trimmed_seq_read_count
                for nonspecific_del_covs, norm_seq_del_config in zip(nonspecific_del_covs_in_mod_seqs, norm_seq_del_configs):
                    for del_index, adjacent_fiveprime_pos in enumerate(norm_seq_del_config):
                        if start_pos <= adjacent_fiveprime_pos and stop_pos >= adjacent_fiveprime_pos + 1:
                            nonspecific_del_covs[del_index] += trimmed_seq_read_count
                if trimmed_seq_length > norm_seq_length:
                    nonspecific_read_with_extra_fiveprime_count += trimmed_seq_read_count
                    if trimmed_seq_length - norm_seq_length >= self.min_length_of_long_fiveprime_extension:
                        nonspecific_long_fiveprime_extension_dict[trimmed_seq.seq_string] = trimmed_seq_read_count
                for threeprime_terminus_seq_string, threeprime_terminus_read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                    nonspecific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += threeprime_terminus_read_count

        # Record the new start and stop positions of the trimmed sequences.
        self.start_positions = start_positions
        self.stop_positions = stop_positions
        self.specific_read_count = specific_read_count
        self.nonspecific_read_count = nonspecific_read_count
        self.specific_covs = specific_covs
        self.nonspecific_covs = nonspecific_covs
        self.specific_del_covs_in_mod_seqs = specific_del_covs_in_mod_seqs
        self.nonspecific_del_covs_in_mod_seqs = nonspecific_del_covs_in_mod_seqs
        self.mean_specific_cov = specific_covs.mean()
        self.mean_nonspecific_cov = nonspecific_covs.mean()
        self.specific_read_with_extra_fiveprime_count = specific_read_with_extra_fiveprime_count
        self.nonspecific_read_with_extra_fiveprime_count = nonspecific_read_with_extra_fiveprime_count
        self.specific_mapped_read_count = specific_mapped_read_count
        self.nonspecific_mapped_read_count = nonspecific_mapped_read_count
        self.specific_long_fiveprime_extension_dict = specific_long_fiveprime_extension_dict
        self.nonspecific_long_fiveprime_extension_dict = nonspecific_long_fiveprime_extension_dict
        self.specific_read_threeprime_terminus_count_dict = specific_read_threeprime_terminus_count_dict
        self.nonspecific_read_threeprime_terminus_count_dict = nonspecific_read_threeprime_terminus_count_dict

        # Store a uniqued list of defunct normalized sequences.
        self.defunct_norm_seqs = derep_defunct_norm_seqs


class ModifiedSequence(object):
    """A tRNA sequence with sites of predicted modification-induced substitutions and deletions,
    formed from normalized sequences with distinct patterns of these mutations.

    The `anvi-trnaseq` workflow aggregates similar normalized sequences and processes the resulting
    clusters to produce clusters of normalized sequences distinguished by potential
    modification-induced substitutions (3-4 different nucleotides at one or more aligned positions).
    A modified sequence is initialized with a list of clustered normalized sequences, with the first
    sequence in the list being the longest, or tied for longest, and substitution positions indexed
    relative to this longest sequence. The workflow later finds sequences with modification-induced
    deletions (which occur in the vicinity of substitutions) and adds them to the modified
    sequences. After adding sequences with deletions, the workflow calls `init` to calculate
    coverages and other information.

    The workflow currently requires that a normalized sequence with deletions must be assigned to
    only one modified sequence. In other words, if the normalized sequence with deletions can be
    found through the introduction of in silico deletions in multiple modified sequences, then the
    sequence is disregarded.

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
        'specific_mapped_read_count',
        'nonspecific_mapped_read_count',
        'specific_read_with_extra_fiveprime_count',
        'nonspecific_read_with_extra_fiveprime_count',
        'specific_long_fiveprime_extension_dict',
        'nonspecific_long_fiveprime_extension_dict',
        'specific_read_threeprime_terminus_count_dict',
        'nonspecific_read_threeprime_terminus_count_dict',
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

    def __init__(self, norm_seqs_without_dels, sub_positions):
        self.norm_seqs_without_dels = norm_seqs_without_dels
        self.sub_positions = sub_positions
        # A normalized sequence without modification-induced deletions can only be assigned to one
        # modified sequence.
        for norm_seq in norm_seqs_without_dels:
            norm_seq.mod_seqs.append(self)
        self.represent_name = norm_seqs_without_dels[0].represent_name

        self.norm_seqs_with_dels = []
        # The deletion configurations of normalized sequences with deletions are recorded in the
        # modified sequence rather than the normalized sequences, because the positions of the
        # deletions correspond to nucleotides in the modified sequence.
        self.del_configs = []
        self.specific_read_count = None
        self.nonspecific_read_count = None
        self.specific_mapped_read_count = None
        self.nonspecific_mapped_read_count = None
        self.specific_read_with_extra_fiveprime_count = None
        self.nonspecific_read_with_extra_fiveprime_count = None
        self.specific_long_fiveprime_extension_dict = None
        self.nonspecific_long_fiveprime_extension_dict = None
        self.specific_read_threeprime_terminus_count_dict = None
        self.nonspecific_read_threeprime_terminus_count_dict = None
        self.specific_covs = None
        self.nonspecific_covs = None
        self.mean_specific_cov = None
        self.mean_nonspecific_cov = None
        self.specific_sub_covs = None
        self.nonspecific_sub_covs = None
        # Currently, each normalized sequence with deletions can only be derived from a single
        # modified sequence, so trimmed sequences specific to a normalized sequence with deletions
        # must also be specific to the modified sequence. This allows specific deletion coverages to
        # be determined. Nonspecific deletion coverages are not determined, because a nonspecific
        # trimmed sequence may occur in multiple normalized sequences with deletions that are
        # assigned to the same modified sequence.
        self.specific_del_covs = None
        self.nonspecific_del_covs = None
        self.consensus_seq_string = None


    def init(self):
        """Set attributes for a "finalized" set of input `NormalizedSeq` objects, some of which may
        contain deletions."""
        norm_seqs_without_dels = self.norm_seqs_without_dels
        norm_seqs_with_dels = self.norm_seqs_with_dels
        all_norm_seqs = norm_seqs_without_dels + norm_seqs_with_dels
        del_configs = self.del_configs
        specific_mapped_read_count = 0
        nonspecific_mapped_read_count = 0
        specific_read_count = 0
        nonspecific_read_count = 0
        specific_read_with_extra_fiveprime_count = 0
        nonspecific_read_with_extra_fiveprime_count = 0
        specific_long_fiveprime_extension_dict = defaultdict(int)
        nonspecific_long_fiveprime_extension_dict = defaultdict(int)
        specific_read_threeprime_terminus_count_dict = defaultdict(int)
        nonspecific_read_threeprime_terminus_count_dict = defaultdict(int)

        mod_seq_len = len(norm_seqs_without_dels[0].seq_string)
        norm_seq_specific_covs = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        norm_seq_nonspecific_covs = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        num_subs = len(self.sub_positions)
        self.specific_sub_covs = specific_sub_covs = np.zeros((num_subs, len(UNAMBIG_NTS)), dtype=int)
        self.nonspecific_sub_covs = nonspecific_sub_covs = np.zeros((num_subs, len(UNAMBIG_NTS)), dtype=int)

        del_positions = sorted(set([i for del_config in del_configs for i in del_config]))
        self.specific_del_covs = specific_del_covs = np.zeros(len(del_positions), dtype=int)
        # Nonspecific deletion coverages are not currently calculated.
        self.nonspecific_del_covs = nonspecific_del_covs = np.zeros(len(del_positions), dtype=int)

        # Trimmed sequences may be shared among the normalized sequences comprising the modified
        # sequence. Therefore, most of the processing involves the analysis of trimmed sequences.
        processed_trimmed_seq_names = []
        all_norm_seq_names = []
        for norm_seq in all_norm_seqs:
            all_norm_seq_names.append(norm_seq.represent_name)
        all_norm_seq_names = [norm_seq.represent_name for norm_seq in all_norm_seqs]

        # Process sequences without deletions.
        for i, norm_seq in enumerate(norm_seqs_without_dels):
            norm_seq_start_in_mod_seq = mod_seq_len - len(norm_seq.seq_string)

            for trimmed_seq, trimmed_seq_start_in_norm_seq, trimmed_seq_stop_in_norm_seq in zip(norm_seq.trimmed_seqs, norm_seq.start_positions, norm_seq.stop_positions):
                if trimmed_seq.represent_name in processed_trimmed_seq_names:
                    continue

                # Determine whether the reads constituting the trimmed sequence are specific to the
                # modified sequence or are found in other normalized sequences outside the modified
                # sequence.
                if len(trimmed_seq.norm_seq_represent_names) == 1:
                    is_trimmed_seq_specific_to_mod_seq = True
                else:
                    # The trimmed sequence is specific to this modified sequence if it is unique to a
                    # set of normalized sequences that are all part of this modified sequence.
                    for norm_seq_containing_trimmed_seq_name in trimmed_seq.norm_seq_represent_names:
                        if norm_seq_containing_trimmed_seq_name not in all_norm_seq_names:
                            is_trimmed_seq_specific_to_mod_seq = False
                            break
                    else:
                        is_trimmed_seq_specific_to_mod_seq = True

                trimmed_seq_start_in_mod_seq = norm_seq_start_in_mod_seq + trimmed_seq_start_in_norm_seq
                trimmed_seq_stop_in_mod_seq = trimmed_seq_start_in_mod_seq + len(trimmed_seq.seq_string)
                trimmed_seq_class_name = type(trimmed_seq).__name__
                if is_trimmed_seq_specific_to_mod_seq:
                    specific_read_count += trimmed_seq.read_count

                    if trimmed_seq_class_name == 'TrimmedMappedSequence':
                        specific_mapped_read_count += trimmed_seq.read_count
                    else:
                        for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                            specific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count

                    if trimmed_seq_class_name != 'TrimmedTruncatedProfileSequence':
                        specific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count

                        for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                            specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

                    norm_seq_specific_covs[i, trimmed_seq_start_in_mod_seq: trimmed_seq_stop_in_mod_seq] += trimmed_seq.read_count
                else:
                    nonspecific_read_count += trimmed_seq.read_count

                    if trimmed_seq_class_name == 'TrimmedMappedSequence':
                        nonspecific_mapped_read_count += trimmed_seq.read_count

                        # Only mapped nonspecific sequences can have 5' extensions, as profiled
                        # sequences with 5' extensions would span the length of the modified
                        # sequence and would thus be specific to it.
                        nonspecific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count

                        for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
                            nonspecific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count
                    else:
                        for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                            nonspecific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count

                    norm_seq_nonspecific_covs[i, trimmed_seq_start_in_mod_seq: trimmed_seq_stop_in_mod_seq] += trimmed_seq.read_count

                processed_trimmed_seq_names.append(trimmed_seq.represent_name)

        # Find deletion site coverages *assuming* that each normalized sequence with
        # deletions is only derived from this modified sequence.
        del_pos_index_dict = {del_pos: del_index for del_index, del_pos in enumerate(del_positions)}
        for norm_seq in self.norm_seqs_with_dels:
            for norm_seq_del_config, mod_seq_del_config, norm_seq_specific_del_covs in zip(norm_seq.norm_seq_del_configs, norm_seq.mod_seq_del_configs, norm_seq.specific_del_covs_in_mod_seqs):
                # The normalized sequence with deletions is aligned with the modified sequence at
                # the 3' end. If the normalized sequence is shorter than the modified sequence, then
                # it may be missing some deletion positions toward the 5' end.
                for mod_seq_del_pos, norm_seq_specific_del_cov in zip(mod_seq_del_config[-len(norm_seq_del_config): ], norm_seq_specific_del_covs):
                    del_index = del_pos_index_dict[mod_seq_del_pos]
                    specific_del_covs[del_index] += norm_seq_specific_del_cov

        # The following commented code was a sketch of how normalized sequences with deletions
        # derived from multiple mdoified sequences can be handled. Currently, normalized sequences
        # with deletions do not contribute to modified sequence coverage. Again, trimmed sequences
        # can be nonspecific to the normalized sequence. They can be found at the same time in
        # normalized sequences both with and without deletions. If a trimmed sequence contains a
        # deletion, its feature profile is invalidated, but an alternate profile is not assigned.

        # Make an array of aligned nucleotides from all normalized sequences.
        # norm_seq_array = np.zeros((len(all_norm_seqs), mod_seq_len), dtype=int)
        # for i, norm_seq in enumerate(norm_seqs_without_dels):
        #     norm_seq_array[i, mod_seq_len - len(norm_seq.seq_string): ] += [NT_INT_DICT[nt] for nt in norm_seq.seq_string]

        # nt_positions_covered_by_norm_seqs_with_dels = []
        # i = len(norm_seqs_without_dels)
        # for norm_seq, del_config in zip(norm_seqs_with_dels, del_configs):
        #     aligned_seq = [NT_INT_DICT[nt] for nt in norm_seq.seq_string]
        #     # Insert a 0 (no nucleotide) at each deletion position in the alignment.
        #     for del_pos in del_config:
        #         aligned_seq.insert(del_pos, 0)
        #     norm_seq_start_in_mod_seq = mod_seq_len - len(aligned_seq)
        #     norm_seq_array[i, norm_seq_start_in_mod_seq: ] += aligned_seq

        #     covered_nt_positions = []
        #     for j, nt_int in enumerate(aligned_seq, norm_seq_start_in_mod_seq):
        #         if nt_int != 0:
        #             covered_nt_positions.append(j)
        #     nt_positions_covered_by_norm_seqs_with_dels.append(covered_nt_positions)
        #     i += 1

        # # Process normalized sequences with deletions.
        # i = len(norm_seqs_without_dels)
        # for norm_seq, del_config, nt_positions_covered_by_norm_seq in zip(norm_seqs_with_dels, del_configs, nt_positions_covered_by_norm_seqs_with_dels):
        #     norm_seq_start_in_mod_seq = mod_seq_len - len(norm_seq.seq_string) - len(del_config)

        #     # Again, the `anvi-trnaseq` workflow currently requires that a normalized sequence with
        #     # deletions be assigned exclusively to a single modified sequence, so the following
        #     # variable must be 1.
        #     num_mod_seqs_containing_norm_seq = len(norm_seq.mod_seqs)

        #     for trimmed_seq, trimmed_seq_start_in_norm_seq, trimmed_seq_stop_in_norm_seq in zip(norm_seq.trimmed_seqs, norm_seq.start_positions, norm_seq.stop_positions):
        #         if trimmed_seq.represent_name in processed_trimmed_seq_names:
        #             continue

        #         if num_mod_seqs_containing_norm_seq > 1:
        #             is_trimmed_seq_specific_to_mod_seq = False
        #         else:
        #             for norm_seq_containing_trimmed_seq_name in trimmed_seq.norm_seq_represent_names:
        #                 try:
        #                     if len(all_norm_seqs[all_norm_seq_names.index(norm_seq_containing_trimmed_seq_name)].mod_seqs) > 1:
        #                         # The trimmed sequence is part of another normalized sequence that
        #                         # is part of multiple modified sequences.
        #                         is_trimmed_seq_specific_to_mod_seq = False
        #                         break
        #                 except ValueError:
        #                     # The trimmed sequence is part of another normalized sequence that is
        #                     # not part of the modified sequence.
        #                     is_trimmed_seq_specific_to_mod_seq = False
        #                     break
        #             else:
        #                 is_trimmed_seq_specific_to_mod_seq = True

        #         nt_positions_covered_by_trimmed_seq = nt_positions_covered_by_norm_seq[trimmed_seq_start_in_norm_seq: trimmed_seq_start_in_norm_seq + len(trimmed_seq.seq_string)]

        #         trimmed_seq_read_count = trimmed_seq.read_count
        #         trimmed_seq_class_name = type(trimmed_seq).__name__
        #         if is_trimmed_seq_specific_to_mod_seq:
        #             specific_read_count += trimmed_seq_read_count

        #             if trimmed_seq_class_name == 'TrimmedMappedSequence':
        #                 specific_mapped_read_count += trimmed_seq_read_count

        #             try:
        #                 specific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count
        #             except:
        #                 print(type(trimmed_seq))
        #                 raise Exception

        #             for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
        #                 specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

        #             for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
        #                 specific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count

        #             norm_seq_specific_covs[i, nt_positions_covered_by_trimmed_seq] += trimmed_seq_read_count
        #             for del_pos in del_config:
        #                 specific_del_covs[del_positions.index(del_pos)] += trimmed_seq_read_count
        #         else:
        #             nonspecific_read_count += trimmed_seq_read_count

        #             if isinstance(trimmed_seq, TrimmedMappedSequence):
        #                 nonspecific_mapped_read_count += trimmed_seq_read_count

        #             # Nonspecific profiled trimmed sequences from normalized sequences without
        #             # deletions cannot be both nonspecific and have 5' extensions. On the other
        #             # hand, the same full-length normalized sequence with deletions can sometimes
        #             # arise from different modified sequences, regardless of the origin of the
        #             # normalized sequence. Trimmed sequences from these nonspecific normalized
        #             # sequences may have 5' extensions.
        #             nonspecific_read_with_extra_fiveprime_count += trimmed_seq.read_with_extra_fiveprime_count

        #             for fiveprime_extension_seq_string, read_count in trimmed_seq.long_fiveprime_extension_dict.items():
        #                 specific_long_fiveprime_extension_dict[fiveprime_extension_seq_string] += read_count

        #             for threeprime_terminus_seq_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
        #                 nonspecific_read_threeprime_terminus_count_dict[threeprime_terminus_seq_string] += read_count

        #             norm_seq_nonspecific_covs[i, nt_positions_covered_by_trimmed_seq] += trimmed_seq_read_count
        #             for del_pos in del_config:
        #                 nonspecific_del_covs[del_positions.index(del_pos)] += trimmed_seq_read_count

        #     processed_trimmed_seq_names.append(trimmed_seq.represent_name)
        #     i += 1

        self.specific_read_count = specific_read_count
        self.nonspecific_read_count = nonspecific_read_count
        self.specific_mapped_read_count = specific_mapped_read_count
        self.nonspecific_mapped_read_count = nonspecific_mapped_read_count
        self.specific_read_with_extra_fiveprime_count = specific_read_with_extra_fiveprime_count
        self.nonspecific_read_with_extra_fiveprime_count = nonspecific_read_with_extra_fiveprime_count
        self.specific_long_fiveprime_extension_dict = specific_long_fiveprime_extension_dict
        self.nonspecific_long_fiveprime_extension_dict = nonspecific_long_fiveprime_extension_dict
        self.specific_read_threeprime_terminus_count_dict = specific_read_threeprime_terminus_count_dict
        self.nonspecific_read_threeprime_terminus_count_dict = nonspecific_read_threeprime_terminus_count_dict
        self.specific_covs = norm_seq_specific_covs.sum(0)
        self.nonspecific_covs = norm_seq_nonspecific_covs.sum(0)
        self.mean_specific_cov = self.specific_covs.mean()
        self.mean_nonspecific_cov = self.nonspecific_covs.mean()
        self.specific_del_covs = self.specific_del_covs

        # For each substitution position, record the coverage of A, C, G, and T. The following array
        # of aligned sequences ignores normalized sequences with deletions.
        norm_seq_array = np.zeros((len(norm_seqs_without_dels), mod_seq_len), dtype=int)
        for norm_seq_index, norm_seq in enumerate(norm_seqs_without_dels):
            norm_seq_array[norm_seq_index, mod_seq_len - len(norm_seq.seq_string): ] += [NT_INT_DICT[nt] for nt in norm_seq.seq_string]
        for sub_num, sub_pos in enumerate(self.sub_positions):
            aligned_nts = norm_seq_array[:, sub_pos]
            nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]
            for nt_int, nt_count in enumerate(nt_counts, start=1):
                if nt_count > 0:
                    norm_seq_rows_with_nt = (aligned_nts == nt_int).nonzero()[0]
                    specific_sub_covs[sub_num, nt_int - 1] = norm_seq_specific_covs[norm_seq_rows_with_nt, sub_pos].sum()
                    nonspecific_sub_covs[sub_num, nt_int - 1] = norm_seq_nonspecific_covs[norm_seq_rows_with_nt, sub_pos].sum()

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
    RELATIVE_ANTICODON_LOOP_INDEX = TRNA_FEATURE_NAMES.index('anticodon_loop') - len(TRNA_FEATURE_NAMES) + 1

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
        self.threeprime_termini_param = A('threeprime_termini')
        self.min_length_of_long_fiveprime_extension = A('min_length_long_fiveprime')
        TrimmedSequence.min_length_of_long_fiveprime_extension = self.min_length_of_long_fiveprime_extension
        NormalizedDeletionSequence.min_length_of_long_fiveprime_extension = self.min_length_of_long_fiveprime_extension
        self.min_trna_frag_size = A('min_trna_fragment_size')
        agglom_max_mismatch_freq = A('agglomeration_max_mismatch_freq')
        self.agglom_max_mismatch_freq = round(agglom_max_mismatch_freq * 100) / 100
        self.skip_INDEL_profiling = A('skip_INDEL_profiling')
        self.fiveprimemost_del_start = A('fiveprimemost_deletion_start')
        self.threeprimemost_del_start = A('threeprimemost_deletion_start')
        self.fiveprimemost_del_stop = A('fiveprimemost_deletion_stop')
        self.threeprimemost_del_stop = A('threeprimemost_deletion_stop')
        self.max_distinct_dels = A('max_distinct_deletions')
        self.min_dist_between_dels = A('min_distance_between_deletions')
        self.max_del_configs = A('max_deletion_configurations')

        # Argument group 1D: PERFORMANCE
        self.num_threads = A('num_threads')
        self.skip_fasta_check = A('skip_fasta_check')
        self.profiling_chunk_size = A('profiling_chunk_size')
        self.alignment_target_chunk_size = A('alignment_target_chunk_size')

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

        # Supplementary text file paths produced by DEBUG flag
        self.uniq_nontrna_path = get_out_dir_path(self.sample_id + "-UNIQUED_NONTRNA.txt")
        self.trimmed_ends_path = get_out_dir_path(self.sample_id + "-TRIMMED_ENDS.txt")
        self.consol_seqs_with_inconsis_profiles_path = get_out_dir_path(self.sample_id + "-CONSOLIDATED_SEQS_WITH_INCONSISTENT_PROFILES.txt")

        # Intermediate pickle file paths
        self.intermed_file_path_dict = {
            'profile': {
                'uniq_trna_seq_dict': get_out_dir_path("UNIQUE_TRNA_SEQS-PROFILE_CHECKPOINT.pkl"),
                'uniq_trunc_seq_dict': get_out_dir_path("UNIQUE_TRUNCATED_SEQS-PROFILE_CHECKPOINT.pkl"),
                'uniq_nontrna_seq_dict': get_out_dir_path("UNIQUE_NONTRNA_SEQS-PROFILE_CHECKPOINT.pkl")
            },
            'normalize': {
                'uniq_trna_seq_dict': get_out_dir_path("UNIQUE_TRNA_SEQS-NORMALIZE_CHECKPOINT.pkl"),
                'uniq_trunc_seq_dict': get_out_dir_path("UNIQUE_TRUNCATED_SEQS-NORMALIZE_CHECKPOINT.pkl"),
                'uniq_nontrna_seq_dict': get_out_dir_path("UNIQUE_NONTRNA_SEQS-NORMALIZE_CHECKPOINT.pkl"),
                'trimmed_trna_seq_dict': get_out_dir_path("TRIMMED_TRNA_SEQS-NORMALIZE_CHECKPOINT.pkl"),
                'trimmed_trunc_seq_dict': get_out_dir_path("TRIMMED_TRUNC_SEQS-NORMALIZE_CHECKPOINT.pkl"),
                'norm_trna_seq_dict': get_out_dir_path("NORMALIZED_TRNA_SEQS-NORMALIZE_CHECKPOINT.pkl"),
                'norm_trunc_seq_dict': get_out_dir_path("NORMALIZED_TRUNC_SEQS-NORMALIZE_CHECKPOINT.pkl")
            },
            'map_fragments': {
                'uniq_trna_seq_dict': get_out_dir_path("UNIQUE_TRNA_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl"),
                'uniq_trunc_seq_dict': get_out_dir_path("UNIQUE_TRUNCATED_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl"),
                'uniq_nontrna_seq_dict': get_out_dir_path("UNIQUE_NONTRNA_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl"),
                'trimmed_trna_seq_dict': get_out_dir_path("TRIMMED_TRNA_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl"),
                'trimmed_trunc_seq_dict': get_out_dir_path("TRIMMED_TRUNC_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl"),
                'norm_trna_seq_dict': get_out_dir_path("NORMALIZED_TRNA_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl"),
                'norm_trunc_seq_dict': get_out_dir_path("NORMALIZED_TRUNC_SEQS-MAP_FRAGMENTS_CHECKPOINT.pkl")
            },
            'substitutions': {
                'uniq_trna_seq_dict': get_out_dir_path("UNIQUE_TRNA_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl"),
                'uniq_trunc_seq_dict': get_out_dir_path("UNIQUE_TRUNCATED_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl"),
                'uniq_nontrna_seq_dict': get_out_dir_path("UNIQUE_NONTRNA_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl"),
                'trimmed_trna_seq_dict': get_out_dir_path("TRIMMED_TRNA_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl"),
                'trimmed_trunc_seq_dict': get_out_dir_path("TRIMMED_TRUNC_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl"),
                'norm_trna_seq_dict': get_out_dir_path("NORMALIZED_TRNA_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl"),
                'norm_trunc_seq_dict': get_out_dir_path("NORMALIZED_TRUNC_SEQS-SUBSTITUTIONS_CHECKPOINT.pkl")
            }
        }
        self.intermed_file_label_dict = {
            'uniq_trna_seq_dict': 'unique tRNA',
            'uniq_trunc_seq_dict': 'unique seqs with a truncated feature profile',
            'uniq_nontrna_seq_dict': 'unique non-tRNA',
            'trimmed_trna_seq_dict': 'trimmed tRNA',
            'trimmed_trunc_seq_dict': 'trimmed seqs with a truncated feature profile',
            'norm_trna_seq_dict': 'normalized tRNA',
            'norm_trunc_seq_dict': 'normalized seqs with a truncated feature profile'
        }

        self.uniq_nontrna_seq_dict = {}
        self.uniq_trunc_seq_dict = {}
        self.uniq_trna_seq_dict = {}
        self.trimmed_trna_seq_dict = {}
        self.trimmed_trunc_seq_dict = {}
        self.norm_trna_seq_dict = {}
        self.norm_trunc_seq_dict = {}
        self.mod_trna_seq_dict = {}
        self.norm_del_seq_dict = {}
        self.uniq_del_seq_dict = {}
        self.trimmed_del_seq_dict = {}

        self.possible_del_starts = None
        self.possible_del_stops = None
        # Ranges representing the locations of deletions in relation to substitution sites. For
        # example, (range(-1, 0), range(-1, 1), range(0, 1)) allows three types of deletions to be
        # introduced at a substitution position: a 1 nucleotide deletion of the adjacent 5'
        # nucleotide a 2 nucleotide deletion of the adjacent 5' nucleotide and the nucleotide at the
        # substitution position, and a 1 nucleotide deletion of the nucleotide at the substitution
        # position.
        self.del_ranges = None

        # The same normalized sequence with deletions may arise from multiple modified sequences.
        # For the sake of simplicity, only consider normalized sequences with deletions that can
        # come from a single modified sequence.
        self.allow_norm_seq_with_dels_from_multiple_mod_seqs = False


    def process(self):
        """The entry method of TRNASeqDataset, called by `anvi-trnaseq`."""
        total_time_start = time.time()

        self.sanity_check()

        if not self.load_checkpoint:
            # Do the steps before the "profile" checkpoint.

            self.create_trnaseq_database()

            if self.feature_param_path:
                # The user provided an optional tRNA feature parameterization file.
                trnaidentifier.TRNAFeatureParameterizer().set_params_from_file(self.feature_param_path)

            # Add the user parameterizations as meta-values in the "self" table of the tRNA-seq
            # database.
            self.report_profiling_parameters()
            self.report_fragment_mapping_parameters()
            self.report_substitution_analysis_parameters()
            self.report_deletion_analysis_parameters()

            # Profile each (unique) read for tRNA features.
            self.profile_trna()

            if self.write_checkpoints:
                self.write_checkpoint_files('profile')
        elif self.load_checkpoint == 'profile':
            self.load_checkpoint_files('profile')
            self.report_fragment_mapping_parameters()
            self.report_substitution_analysis_parameters()
            self.report_deletion_analysis_parameters()


        if (self.load_checkpoint == 'profile'
            or not self.load_checkpoint):
            # Do the steps between the "profile" and "normalize" checkpoints.

            # Trim 5' and 3' ends of profiled tRNA.
            self.trim_trna_ends()
            # Trim 3' ends of sequences with truncated tRNA profile.
            self.trim_truncated_profile_ends()

            # Consolidate 3' fragments of longer profiled tRNA sequences, forming normalized sequences.
            self.threeprime_dereplicate_trna()

            # Recover tRNA sequences with truncated feature profiles by comparing to normalized
            # sequences.
            self.threeprime_dereplicate_truncated_sequences()

            if self.write_checkpoints:
                self.write_checkpoint_files('normalize')
        elif self.load_checkpoint == 'normalize':
            self.load_checkpoint_files('normalize')
            self.report_fragment_mapping_parameters()
            self.report_substitution_analysis_parameters()
            self.report_deletion_analysis_parameters()


        if (self.load_checkpoint == 'normalize'
            or self.load_checkpoint == 'profile'
            or not self.load_checkpoint):
            # Do the steps between the "normalize" and "map_fragments" checkpoints.

            if '' not in self.threeprime_termini:
                # Recover 3' tRNA sequences lacking a 3' terminus.
                self.threeprime_dereplicate_sequences_without_terminus()

            # Map fragments derived from the interior and 5' end of tRNA.
            self.map_fragments()

            # "Finalize" normalized tRNA sequence objects now that all of the trimmed tRNA sequences
            # found through various means have been added to them.
            self.progress.new("Finalizing normalized tRNA sequences")
            self.progress.update("...")
            for norm_seq in self.norm_trna_seq_dict.values():
                norm_seq.init()
            self.progress.end()

            if self.write_checkpoints:
                self.write_checkpoint_files('map_fragments')
        elif self.load_checkpoint == 'map_fragments':
            self.load_checkpoint_files('map_fragments')
            self.report_substitution_analysis_parameters()
            self.report_deletion_analysis_parameters()


        if (self.load_checkpoint == 'map_fragments'
            or self.load_checkpoint == 'normalized'
            or self.load_checkpoint == 'profile'
            or not self.load_checkpoint):
            # Do the steps between the "map_fragments" and "substitutions" checkpoints.

            # Find modified nucleotides, grouping normalized sequences into modified sequences.
            self.find_substitutions()

            if self.write_checkpoints:
                self.write_checkpoint_files('substitutions')
        elif self.load_checkpoint == 'substitutions':
            self.load_checkpoint_files('substitutions')
            self.report_deletion_analysis_parameters()

        # Do the steps after the "substitutions" checkpoint.
        if not self.skip_INDEL_profiling:
            self.find_deletions()
        for mod_seq in self.mod_trna_seq_dict.values():
            mod_seq.init()

        self.report_statistics()

        self.write_feature_table()
        self.write_unconserved_table()
        self.write_unpaired_table()
        self.write_sequences_table()
        self.write_trimmed_table()
        self.write_normalized_table()
        self.write_modified_table()

        # Write supplementary text files.
        self.write_unique_nontrna_supplement()
        self.write_trimmed_supplement()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Total time elapsed (min)",
                                          time.time() - total_time_start,
                                          is_time_value=True))
            # Write an empty line to separate this run from any subsequent run starting from a
            # checkpoint writing to the same summary file.
            f.write("\n")


    def sanity_check(self):
        """Check `anvi-trnaseq` user inputs."""
        if os.path.exists(self.out_dir):
            self.existing_output_directory_sanity_check()

        if self.load_checkpoint:
            self.load_checkpoint_sanity_check()
        elif not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

        filesnpaths.is_output_dir_writable(self.out_dir)

        if self.descrip_path:
            filesnpaths.is_file_plain_text(self.descrip_path)
            self.descrip = open(self.descrip_path).read()

        if self.num_threads < 1:
            raise ConfigError("Surely you must be joking, Mr. Feynman! "
                              "`--num-threads` wants a value greater than 0. "
                              f"Last we checked, {self.num_threads} is not greater than 0.")

        self.threeprime_termini = self.threeprime_termini_sanity_check()
        trnaidentifier.TRNAFeatureParameterizer.set_threeprime_termini(self.threeprime_termini)

        self.deletion_start_stop_sanity_check()

        self.parameterize_deletions()

        self.run.info("Input FASTA file", self.input_fasta_path, nl_after=1)

        if not self.skip_fasta_check and not self.load_checkpoint:
            self.progress.new("Checking input FASTA defline format")
            self.progress.update("...")

            utils.check_fasta_id_formatting(self.input_fasta_path)

            self.progress.end()

            self.run.info_single("FASTA deflines were found to be anvi'o-compliant", mc='green', nl_after=1)


    def existing_output_directory_sanity_check(self):
        """Conditions must be fulfilled for the `anvi-trnaseq` output directory to already exist."""
        if len(os.listdir(self.out_dir)) == 0:
            # There is nothing in the output directory.
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
                raise ConfigError(f"The directory that was specified by --output-dir or -o, {self.out_dir}, already exists. "
                                  "Use the flag --overwrite-output-destinations to overwrite this directory.")


    def load_checkpoint_sanity_check(self):
        """Needed intermediate files must exist to load from a checkpoint."""
        missing_intermed_files = []
        for intermed_file_path in self.intermed_file_path_dict[self.load_checkpoint].values():
            if not os.path.exists(intermed_file_path):
                missing_intermed_files.append(intermed_file_path)
        if missing_intermed_files:
            raise ConfigError(f"Intermediate files needed for running `anvi-trnaseq` "
                              f"with `--load-checkpoint {self.load_checkpoint}` are missing: {', '.join(missing_intermed_files)}. "
                              "You should probably run `anvi-trnaseq` from the beginning without `--load-checkpoint`. "
                              "To generate necessary intermediate files for future use of `--load-checkpoint`, use the flag `--write-checkpoints`.")


    def threeprime_termini_sanity_check(self):
        """Check validity of provided tRNA 3' termini, returning a list of terminus strings."""
        valid_threeprime_termini = []
        invalid_threeprime_termini = []
        for threeprime_terminus in self.threeprime_termini_param.split(','):
            if threeprime_terminus == '_':
                valid_threeprime_termini.append('')
                continue

            for nt in threeprime_terminus:
                if nt not in ALL_NTS:
                    invalid_threeprime_termini.append(threeprime_terminus)
                    break
            valid_threeprime_termini.append(threeprime_terminus)

        if invalid_threeprime_termini:
            raise ConfigError(f"3' termini can consist of A, C, G, T, and N (any nucleotide) "
                              "or the discriminator nucleotide with no extension, symbolized by a single underscore, \"_\". "
                              f"The following invalid 3' sequence parameterizations were provided: {', '.join(invalid_threeprime_termini)}")

        return valid_threeprime_termini


    def deletion_start_stop_sanity_check(self):
        """Check the logic of allowed start and stop positions of in silico deletions relative to
        sites of predicted modification-induced substitutions."""
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


    def parameterize_deletions(self):
        """Parameterize in silico deletions."""
        # In silico deletions are allowed to start and stop at positions relative to sites of
        # predicted modification-induced substitutions.
        self.possible_del_starts = tuple(range(self.fiveprimemost_del_start, self.threeprimemost_del_start + 1))
        self.possible_pythonic_del_stops = tuple(range(self.fiveprimemost_del_stop + 1, self.threeprimemost_del_stop + 2))
        # In silico deletions can be different sizes, so, in combination with possible start and
        # stop positions, find all of the different deletions allowed relative to predicted
        # substitutions.
        del_ranges = []
        for del_start in self.possible_del_starts:
            for del_stop in self.possible_pythonic_del_stops:
                if del_start < del_stop:
                    del_ranges.append(range(del_start, del_stop))
        self.del_ranges = tuple(del_ranges)


    def create_trnaseq_database(self):
        """Create an empty tRNA-seq database."""
        meta_values = {'sample_id': self.sample_id,
                       'treatment': self.treatment,
                       'description': self.descrip if self.descrip else '_No description is provided_',
                       'INDELs_profiled': not self.skip_INDEL_profiling}
        dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True).create(meta_values)
        self.run.info("New tRNA-seq db", self.trnaseq_db_path, nl_after=1)


    def report_profiling_parameters(self):
        """Add profiling parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        db = trnaseq_db.db
        parameterizer = trnaidentifier.TRNAFeatureParameterizer()
        for param_tuple in parameterizer.list_accessible_param_tuples():
            db.set_meta_value(param_tuple[0], param_tuple[1])
        db.set_meta_value('min_length_of_long_fiveprime_extension', self.min_length_of_long_fiveprime_extension)
        trnaseq_db.disconnect()

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            for param_name, param_value in parameterizer.list_accessible_param_tuples(pretty=True):
                if 'Conserved nucleotides' in param_name:
                    f.write(get_summary_line(param_name.replace('Conserved nucleotides', "Conserved nts"), param_value))
                    continue
                elif 'Number allowed unconserved' in param_name:
                    f.write(get_summary_line(param_name.replace('Number allowed unconserved', "Allowed number of unconserved nts"), param_value))
                    continue
                elif 'Number allowed unpaired' in param_name:
                    f.write(get_summary_line(param_name.replace('Number allowed unpaired', "Allowed number of unpaired bps"), param_value))
                    continue
                f.write(get_summary_line(param_name, param_value))
            f.write(get_summary_line("Allowed 3' termini", ",".join(self.threeprime_termini)))
            f.write(get_summary_line("Min length of \"long\" 5' extension", self.min_length_of_long_fiveprime_extension))


    def report_fragment_mapping_parameters(self):
        """Add fragment mapping parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db.set_meta_value('min_mapped_trna_fragment_size', self.min_trna_frag_size)
        trnaseq_db.disconnect()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Min length of mapped tRNA fragment", self.min_trna_frag_size))


    def report_substitution_analysis_parameters(self):
        """Add modification-induced substitution analysis parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        db = trnaseq_db.db
        db.set_meta_value('agglomeration_max_mismatch_freq', self.agglom_max_mismatch_freq)
        trnaseq_db.disconnect()

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("Agglomeration max mismatch frequency", self.agglom_max_mismatch_freq))


    def report_deletion_analysis_parameters(self):
        """Add modification-induced deletion analysis parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        db = trnaseq_db.db
        db.set_meta_value('fiveprimemost_deletion_start', self.fiveprimemost_del_start)
        db.set_meta_value('threeprimemost_deletion_start', self.threeprimemost_del_start)
        db.set_meta_value('fiveprimemost_deletion_stop', self.fiveprimemost_del_stop)
        db.set_meta_value('threeprimemost_deletion_stop', self.threeprimemost_del_stop)
        db.set_meta_value('max_distinct_deletions', self.max_distinct_dels)
        db.set_meta_value('min_distance_between_deletions', self.min_dist_between_dels)
        db.set_meta_value('max_deletion_configurations', self.max_del_configs)
        trnaseq_db.disconnect()

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("INDELs profiled", not self.skip_INDEL_profiling))
            f.write(get_summary_line("5'-most del start", self.fiveprimemost_del_start))
            f.write(get_summary_line("3'-most del start", self.threeprimemost_del_start))
            f.write(get_summary_line("5'-most del start", self.threeprimemost_del_start))
            f.write(get_summary_line("3'-most del stop", self.threeprimemost_del_start))
            f.write(get_summary_line("Max distinct dels", self.max_distinct_dels))
            f.write(get_summary_line("Min distance between dels", self.min_dist_between_dels))
            f.write(get_summary_line("Max del configurations", self.max_del_configs))


    def get_summary_line(self, label, value, is_time_value=False, padding=68):
        """Return a string formatted to be written to the summary statistics file."""
        # Report elapsed time in seconds in minutes.
        if is_time_value:
            value = "%.2f" % round(value / 60, 2)
        return '%s%s\t%s\n' % (label, ' ' + '.' * (padding - len(label)), value)


    def profile_trna(self):
        """Profile tRNA features in reads. Add `UniqueSequence` objects representing profiled tRNA
        sequences to `self.uniq_trna_seq_dict`. Add `UniqueSequence` objects representing sequences
        with a truncated tRNA profile to `self.uniq_trunc_seq_dict`. Add leftover `UniqueSequence`
        objects representing unprofiled tRNA sequences to `self.uniq_nontrna_seq_dict`."""
        uniq_reads = self.unique_reads()

        start_time = time.time()

        pid = "Profiling tRNA features in unique reads"
        self.progress.new(pid)
        self.progress.update("...")

        # Count the number of reads and unique read sequences that have been added to the
        # multiprocessing input queue.
        total_read_count = 0
        total_uniq_count = len(uniq_reads)

        manager = mp.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        profiler = trnaidentifier.Profiler()
        processes = [mp.Process(target=profile_worker, args=(input_queue, output_queue, profiler))
                     for _ in range(self.num_threads)]
        for p in processes:
            p.start()

        # Count the number of unique sequences that have been profiled and fetched from the
        # multiprocessing output queue.
        fetched_profile_count = 0
        input_count = 0
        interval_start = 0
        profiling_chunk_size = self.profiling_chunk_size
        interval_stop = profiling_chunk_size if profiling_chunk_size < total_uniq_count else total_uniq_count
        uniq_trna_seq_dict = self.uniq_trna_seq_dict
        uniq_trunc_seq_dict = self.uniq_trunc_seq_dict
        uniq_nontrna_seq_dict = self.uniq_nontrna_seq_dict
        pp_total_uniq_count = pp(len(uniq_reads))
        while fetched_profile_count < total_uniq_count:
            self.progress.update_pid(pid)
            self.progress.update(f"{pp(input_count + 1)}-{pp(interval_stop)}/{pp_total_uniq_count}")

            while input_count < interval_stop:
                for uniq_read_info in uniq_reads[interval_start: interval_stop]:
                    input_queue.put(uniq_read_info)
                    total_read_count += uniq_read_info[2]
                    input_count += 1

            while fetched_profile_count < interval_stop:
                profile, read_count = output_queue.get()
                fetched_profile_count += 1

                name = profile.name
                if profile.is_predicted_trna:
                    uniq_trna_seq_dict[name] = UniqueFullProfileSequence(profile.input_seq, name, read_count, profile)
                else:
                    if profile.trunc_profile_index:
                        uniq_trunc_seq_dict[name] = UniqueTruncatedProfileSequence(profile.input_seq, name, read_count, profile)
                    else:
                        uniq_nontrna_seq_dict[name] = UniqueSequence(profile.input_seq, name, read_count)

            interval_start = interval_stop
            interval_stop += profiling_chunk_size if interval_stop + profiling_chunk_size < total_uniq_count else total_uniq_count - interval_stop

        for p in processes:
            p.terminate()
            p.join()

        # Profiled seqs were added to the output queue as they were processed, so sort by name.
        self.uniq_trna_seq_dict = {name: seq for name, seq in sorted(uniq_trna_seq_dict.items())}
        self.uniq_trunc_seq_dict = {name: seq for name, seq in sorted(uniq_trunc_seq_dict.items())}
        self.uniq_nontrna_seq_dict = {name: seq for name, seq in sorted(uniq_nontrna_seq_dict.items())}

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("Time elapsed profiling tRNA (min)", time.time() - start_time, is_time_value=True))
            f.write(get_summary_line("Reads processed", total_read_count))
            f.write(get_summary_line("Unique seqs processed", total_uniq_count))

        self.progress.end()

        self.run.info("Reads processed", total_read_count)
        self.run.info("Unique seqs processed", total_uniq_count, nl_after=1)


    def unique_reads(self):
        """Dereplicate input reads."""
        self.progress.new("Loading reads")
        self.progress.update("...")

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

        self.progress.end()


        self.progress.new("Dereplicating reads")
        self.progress.update("...")
        uniq_reads = []
        for cluster in Dereplicator(names, seqs).full_length_dereplicate():
            uniq_reads.append((cluster.member_seqs[0], cluster.member_names[0], len(cluster.member_names)))

        self.progress.end()
        return uniq_reads


    def trim_trna_ends(self):
        """Trim any nucleotides 5' of the acceptor stem and 3' of the discriminator from profiled
        tRNA. Add `TrimmedSequence` objects formed from input fully profiled `UniqueProfileSequence`
        objects to `self.trimmed_trna_seq_dict`."""
        start_time = time.time()
        self.progress.new("Trimming the 3' and 5' ends of tRNA")
        self.progress.update("...")

        trimmed_seqs = self.get_trimmed_seqs([uniq_seq for uniq_seq in self.uniq_trna_seq_dict.values()], TrimmedFullProfileSequence)
        trimmed_trna_seq_dict = self.trimmed_trna_seq_dict
        for trimmed_seq in sorted(trimmed_seqs, key=lambda trimmed_seq: trimmed_seq.represent_name):
            trimmed_trna_seq_dict[trimmed_seq.represent_name] = trimmed_seq

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed trimming profiled seqs (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def get_trimmed_seqs(self, uniq_seqs, trimmed_seq_class):
        """Find trimmed sequences with either a full or truncated profile
        ('TrimmedFullProfileSequence` or `TrimmedTruncatedProfileSequence`, respectively) from
        unique profiled sequences."""
        names = [uniq_seq.represent_name for uniq_seq in uniq_seqs]
        if trimmed_seq_class == TrimmedFullProfileSequence:
            trimmed_seq_strings = [uniq_seq.seq_string[uniq_seq.extra_fiveprime_length: len(uniq_seq.seq_string) - uniq_seq.threeprime_terminus_length]
                                   for uniq_seq in uniq_seqs]
        elif trimmed_seq_class == TrimmedTruncatedProfileSequence:
            trimmed_seq_strings = [uniq_seq.seq_string[: len(uniq_seq.seq_string) - uniq_seq.threeprime_terminus_length]
                                   for uniq_seq in uniq_seqs]

        clusters = Dereplicator(names, trimmed_seq_strings, extras=uniq_seqs).full_length_dereplicate()

        trimmed_seqs = [trimmed_seq_class(cluster.member_seqs[0], cluster.member_extras) for cluster in clusters]
        return trimmed_seqs


    def trim_truncated_profile_ends(self):
        """Trim any nucleotides 3' of the discriminator from sequences with a truncated tRNA
        profile. Add `TrimmedTruncatedProfileSequence` objects formed from input
        `UniqueTruncatedProfileSequence` objects to `self.trimmed_trunc_seq_dict`."""
        start_time = time.time()
        self.progress.new("Trimming the 3' ends of seqs with truncated tRNA profiles")
        self.progress.update("...")

        trimmed_seqs = self.get_trimmed_seqs([uniq_seq for uniq_seq in self.uniq_trunc_seq_dict.values()], TrimmedTruncatedProfileSequence)
        trimmed_trunc_seq_dict = self.trimmed_trunc_seq_dict
        for trimmed_seq in sorted(trimmed_seqs, key=lambda trimmed_seq: trimmed_seq.represent_name):
            trimmed_trunc_seq_dict[trimmed_seq.represent_name] = trimmed_seq

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed trimming seqs with truncated feature profile (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def threeprime_dereplicate_trna(self):
        """Dereplicate trimmed profiled tRNA sequences from the 3' end of longer trimmed sequences.

        EXAMPLE:
        Normalized tRNA (trimmed tRNA 1): TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Trimmed tRNA 2                  :                       AATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Trimmed tRNA 3                  :                                     GCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        """
        start_time = time.time()
        pid = "Dereplicating trimmed tRNA seqs from the 3' end"
        self.progress.new(pid)
        self.progress.update("...")

        # Prefix dereplicate trimmed sequences from the 3' end.
        names = []
        reversed_seq_strings = []
        trimmed_seqs = []
        for name, trimmed_seq in self.trimmed_trna_seq_dict.items():
            names.append(name)
            reversed_seq_strings.append(trimmed_seq.seq_string[::-1])
            trimmed_seqs.append(trimmed_seq)
        clusters = Dereplicator(names, reversed_seq_strings, extras=trimmed_seqs).prefix_dereplicate()

        # Profiling may have found multiple sequences that would here be 3'-dereplicated as having
        # complete, but different, feature profiles. Consider the shortest "completely profiled"
        # trimmed sequence in the cluster to be the correct profile. Reclassify discrepant 5'
        # nucleotides from the unique sequences in the longer trimmed sequences as extra 5'
        # nucleotides. Transfer the profile information from the representative unique sequence of
        # the shortest trimmed sequence to the unique sequences of the longer trimmed sequences.
        # Produce a new trimmed sequence from all of these unique sequences, updating
        # `self.trimmed_trna_seq_dict` and `self.uniq_trna_seq_dict`.

        # Similarly, the longest sequence in the cluster may have an erroneous "incomplete profile."
        # If there is a shorter sequence in the cluster with a complete profile, then any longer
        # sequences with an incomplete profile can be reevaluated in the same manner.

        # This step, which likely affects a tiny number of trimmed sequences, helps reduce the
        # number of extra, wrong nucleotides at the 5' end of downstream seed sequences. In debug
        # mode, consolidated trimmed sequences are written to a file.

        # Nota bene: A complete profile may be extrapolated at the 5' end of the 5' strand of the
        # acceptor stem. By default, only one nucleotide may be extrapolated in the stem when all
        # other nucleotides form base pairs. So it is hard to see how an extrapolated complete
        # profile is more likely to be inaccurate than the removed profiles of longer sequences in
        # the cluster.

        # Do not check for feature-by-feature agreement among clustered profiles here, as 3' tRNA
        # fragments occasionally have some incorrect feature positions due to a paucity of sequence
        # information available in fragment profiling. It is conceivable that the seed sequence of
        # the cluster is assigned a wrong complete or incomplete profile, but all of the shorter
        # sequences in the cluster are assigned correct incomplete profiles -- presumably a rare
        # inaccuracy that goes unchecked.
        self.progress.update_pid(pid)
        self.progress.update("Inspecting normalized seq clusters")

        norm_trna_seq_dict = self.norm_trna_seq_dict
        uniq_trna_seq_dict = self.uniq_trna_seq_dict

        # It is possible that trimmed sequences from multiple clusters can consolidate.
        consol_trimmed_seq_dict = {}

        # This dict is for an edge case explained below.
        his_trimmed_seq_dict = {}

        if anvio.DEBUG:
            consol_seqs_with_inconsis_profiles_file = open(self.consol_seqs_with_inconsis_profiles_path, 'w')
            consol_seqs_with_inconsis_profiles_file.write("Index\tTrimmed (0) or Unique (1)\tSequence\n")
            inconsis_profile_cluster_count = 0

        for cluster in clusters:
            # Skip initialization of `NormalizedSequence` objects, as additional `TrimmedSequence`
            # members are later added to the objects after dereplicating sequences with truncated
            # tRNA profiles and mapping unprofiled tRNA fragments.
            trimmed_seqs = cluster.member_extras
            if len(trimmed_seqs) == 1:
                norm_trna_seq_dict[trimmed_seqs[0].represent_name] = NormalizedFullProfileSequence(trimmed_seqs)
                continue

            # Check that there are no shorter sequences in the cluster with a "complete profile".
            complete_profile_indices = []
            for trimmed_seq_index, trimmed_seq in enumerate(trimmed_seqs):
                if trimmed_seq.has_complete_feature_set:
                    complete_profile_indices.append(trimmed_seq_index)

            if not complete_profile_indices:
                norm_trna_seq_dict[trimmed_seqs[0].represent_name] = NormalizedFullProfileSequence(trimmed_seqs)
                continue

            if complete_profile_indices == [0]:
                norm_trna_seq_dict[trimmed_seqs[0].represent_name] = NormalizedFullProfileSequence(trimmed_seqs)
                continue

            # Reaching this point means that there are multiple sequences with "complete
            # profiles" in the cluster.

            # If the two shortest sequences with complete feature profiles differ by the
            # post-transcriptionally added 5'-G of tRNA-His, then they should both be maintained as
            # separate trimmed sequences.
            if trimmed_seqs[complete_profile_indices[-2]].has_his_g:
                if trimmed_seqs[complete_profile_indices[-1]].seq_string == trimmed_seqs[complete_profile_indices[-2]].seq_string[1: ]:

                    if len(complete_profile_indices) == 2:
                        assert complete_profile_indices[-2] == 0
                        his_trimmed_seq_dict[trimmed_seqs[1].represent_name] = trimmed_seqs
                        continue

                    # Perhaps more than two trimmed sequences in the cluster have "complete"
                    # profiles, though this has not been checked. In this case, consolidate the
                    # sequences, retaining the profile of the shortest.

            if anvio.DEBUG:
                # Report the consolidated sequences with different complete feature profiles.
                inconsis_profile_cluster_count += 1
                for trimmed_seq in trimmed_seqs[: complete_profile_indices[-1] + 1]:
                    consol_seqs_with_inconsis_profiles_file.write(str(inconsis_profile_cluster_count) + "\t")
                    consol_seqs_with_inconsis_profiles_file.write("0\t")
                    consol_seqs_with_inconsis_profiles_file.write(trimmed_seq.seq_string + "\n")
                    for uniq_seq in trimmed_seq.uniq_seqs:
                        consol_seqs_with_inconsis_profiles_file.write(str(inconsis_profile_cluster_count) + "\t")
                        consol_seqs_with_inconsis_profiles_file.write("1\t")
                        consol_seqs_with_inconsis_profiles_file.write(uniq_seq.seq_string + "\n")

            short_trimmed_seq = trimmed_seqs[complete_profile_indices[-1]]
            if short_trimmed_seq.represent_name in consol_trimmed_seq_dict:
                # Trimmed sequences from multiple clusters consolidate with the same shorter trimmed
                # sequence with a complete profile.

                # Some of the longer trimmed sequences with rejected "complete" profiles occur in
                # multiple clusters. Clearly, the longest seed sequence, which is also being
                # consolidated, must differ between the clusters. Trimmed sequences in the clusters
                # shorter than that with the chosen complete profile must be the same in each
                # cluster.
                long_trimmed_seqs = consol_trimmed_seq_dict[short_trimmed_seq.represent_name]['long_trimmed_seqs']
                encountered_long_trimmed_seq_names = [trimmed_seq.represent_name for trimmed_seq in long_trimmed_seqs]
                for long_trimmed_seq in trimmed_seqs[: complete_profile_indices[-1]]:
                    if long_trimmed_seq.represent_name not in encountered_long_trimmed_seq_names:
                        long_trimmed_seqs.append(long_trimmed_seq)
            else:
                # This is the first time the shortest sequence with a complete profile has been
                # found in a cluster.
                consol_trimmed_seq_dict[short_trimmed_seq.represent_name] = {'short_trimmed_seq': short_trimmed_seq,
                                                                             'long_trimmed_seqs': trimmed_seqs[: complete_profile_indices[-1]],
                                                                             'norm_seq_members': trimmed_seqs[complete_profile_indices[-1] + 1: ]}
        if anvio.DEBUG:
            consol_seqs_with_inconsis_profiles_file.close()

        # Consider the following edge case. One cluster had two trimmed sequences with complete
        # profiles, T1 and T2, so the sequences were consolidated, forming normalized sequence, N1.
        # T1:      GGTGGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # T2 == N1:   GGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # Another cluster had two trimmed sequences, T3 and the same T2, only differing by a
        # supposed tRNA-His 5'-G.
        # T3:        GGGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # T2:         GGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # Rather than creating another normalized sequence, N2, with the sequence of T3,
        # consolidate the trimmed sequences from the two clusters. This avoids producing two
        # normalized sequences, one of which is the 3' subsequence of the other.
        for name, trimmed_seqs in his_trimmed_seq_dict.items():
            if name in consol_trimmed_seq_dict:
                consol_trimmed_seq_dict[name]['long_trimmed_seqs'].append(trimmed_seqs[0])
            else:
                norm_trna_seq_dict[trimmed_seqs[0].represent_name] = NormalizedFullProfileSequence(trimmed_seqs)

        for consol_trimmed_seq_subdict in consol_trimmed_seq_dict.values():
            consol_trimmed_seq = self.consolidate_trimmed_sequences(consol_trimmed_seq_subdict['short_trimmed_seq'], consol_trimmed_seq_subdict['long_trimmed_seqs'])
            norm_trna_seq_dict[consol_trimmed_seq.represent_name] = NormalizedFullProfileSequence([consol_trimmed_seq] + consol_trimmed_seq_subdict['norm_seq_members'])

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed 3'-dereplicating trimmed profiled seqs", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def consolidate_trimmed_sequences(self, short_trimmed_seq, long_trimmed_seqs):
        """Consolidate longer trimmed sequences with a shorter trimmed sequence by pooling all of
        their unique sequences, changing the profile information of these sequences, and generating
        a new trimmed sequence. Update `self.trimmed_trna_seq_dict` and `self.uniq_trna_seq_dict` in
        the process."""
        short_trimmed_seq_string = short_trimmed_seq.seq_string
        short_uniq_seq = short_trimmed_seq.uniq_seqs[0]

        # To transfer profile information from the representative unique sequence to other unique
        # sequences, determine where features are relative to the 3' end of the short trimmed
        # sequence, which is found in the other unique sequences.
        replacement_info_dict = {'trimmed_seq_string': short_trimmed_seq.seq_string}
        feature_index_adjustment = -short_uniq_seq.extra_fiveprime_length - len(short_trimmed_seq_string)
        feature_start_indices_from_trimmed_threeprime = []
        for feature_start_index in short_uniq_seq.feature_start_indices:
            if isinstance(feature_start_index, int):
                feature_start_indices_from_trimmed_threeprime.append(feature_start_index + feature_index_adjustment)
            else:
                feature_start_indices_from_trimmed_threeprime.append(tuple([strand_start_index + feature_index_adjustment for strand_start_index in feature_start_index]))
        replacement_info_dict['feature_start_indices_from_trimmed_threeprime'] = feature_start_indices_from_trimmed_threeprime
        feature_stop_indices_from_trimmed_threeprime = []
        for feature_stop_index in short_uniq_seq.feature_stop_indices:
            if isinstance(feature_stop_index, int):
                feature_stop_indices_from_trimmed_threeprime.append(feature_stop_index + feature_index_adjustment)
            else:
                feature_stop_indices_from_trimmed_threeprime.append(tuple([strand_stop_index + feature_index_adjustment for strand_stop_index in feature_stop_index]))
        replacement_info_dict['feature_stop_indices_from_trimmed_threeprime'] = feature_stop_indices_from_trimmed_threeprime
        replacement_info_dict['has_his_g'] = short_uniq_seq.has_his_g
        replacement_info_dict['alpha_start_index_from_trimmed_threeprime'] = None if short_uniq_seq.alpha_start_index is None else short_uniq_seq.alpha_start_index + feature_index_adjustment
        replacement_info_dict['alpha_stop_index_from_trimmed_threeprime'] = None if short_uniq_seq.alpha_stop_index is None else short_uniq_seq.alpha_stop_index + feature_index_adjustment
        replacement_info_dict['beta_start_index_from_trimmed_threeprime'] = None if short_uniq_seq.beta_start_index is None else short_uniq_seq.beta_start_index + feature_index_adjustment
        replacement_info_dict['beta_stop_index_from_trimmed_threeprime'] = None if short_uniq_seq.beta_stop_index is None else short_uniq_seq.beta_stop_index + feature_index_adjustment
        replacement_info_dict['anticodon_string'] = short_uniq_seq.anticodon_string
        replacement_info_dict['anticodon_aa'] = short_uniq_seq.anticodon_aa
        replacement_info_dict['contains_anticodon'] = short_uniq_seq.contains_anticodon
        replacement_info_dict['num_conserved'] = short_uniq_seq.num_conserved
        replacement_info_dict['num_unconserved'] = short_uniq_seq.num_unconserved
        replacement_info_dict['num_paired'] = short_uniq_seq.num_paired
        replacement_info_dict['num_unpaired'] = short_uniq_seq.num_unpaired
        unconserved_info_from_trimmed_threeprime = []
        for unconserved_tuple in short_uniq_seq.unconserved_info:
            unconserved_info_from_trimmed_threeprime.append((unconserved_tuple[0] + feature_index_adjustment,
                                                             unconserved_tuple[1],
                                                             unconserved_tuple[2]))
        replacement_info_dict['unconserved_info_from_trimmed_threeprime'] = unconserved_info_from_trimmed_threeprime
        unpaired_info_from_trimmed_threeprime = []
        for unpaired_tuple in short_uniq_seq.unpaired_info:
            unpaired_info_from_trimmed_threeprime.append((unpaired_tuple[0] + feature_index_adjustment,
                                                          unpaired_tuple[1] + feature_index_adjustment,
                                                          unpaired_tuple[2],
                                                          unpaired_tuple[3]))
        replacement_info_dict['unpaired_info_from_trimmed_threeprime'] = unpaired_info_from_trimmed_threeprime
        replacement_info_dict['profiled_seq_without_terminus_length'] = short_uniq_seq.profiled_seq_length - short_uniq_seq.threeprime_terminus_length

        trimmed_trna_seq_dict = self.trimmed_trna_seq_dict
        uniq_trna_seq_dict = self.uniq_trna_seq_dict
        trimmed_trna_seq_dict.pop(short_trimmed_seq.represent_name)
        uniq_seqs = short_trimmed_seq.uniq_seqs
        for uniq_seq in uniq_seqs:
            uniq_seq.trimmed_seq_represent_name = None
        for long_trimmed_seq in long_trimmed_seqs:
            trimmed_trna_seq_dict.pop(long_trimmed_seq.represent_name)

            for uniq_seq in long_trimmed_seq.uniq_seqs:
                uniq_seq.trimmed_seq_represent_name = None
                uniq_trna_seq_dict.pop(uniq_seq.represent_name)
                new_uniq_seq = UniqueTransferredProfileSequence(uniq_seq, replacement_info_dict)
                uniq_seqs.append(new_uniq_seq)
                uniq_trna_seq_dict[new_uniq_seq.represent_name] = new_uniq_seq

        consol_trimmed_seqs = self.get_trimmed_seqs(uniq_seqs, TrimmedFullProfileSequence)

        if len(consol_trimmed_seqs) > 1:
            raise ConfigError(f"Consolidation should have produced only 1 trimmed profiled tRNA sequence, not {len(consol_trimmed_seqs)}.")

        consol_trimmed_seq = consol_trimmed_seqs[0]
        trimmed_trna_seq_dict[consol_trimmed_seq.represent_name] = consol_trimmed_seq

        return consol_trimmed_seq


    def threeprime_dereplicate_truncated_sequences(self):
        """Recover sequences with truncated tRNA profiles that are found to be 3' subsequences of
        profiled tRNA and thus legitimate 3' tRNA fragments. These trimmed sequences are folded into
        the normalized tRNA sequences of `self.norm_trna_seq_dict`. Unrecovered trimmed sequences
        are themselves 3'-dereplicated, forming another pool of normalized sequences in
        `self.norm_trunc_seq_dict`."""
        start_time = time.time()
        self.progress.new("Dereplicating trimmed seqs with a truncated feature profile")
        self.progress.update("...")

        # Prefix dereplicate both trimmed sequences with truncated profiles and normalized sequences
        # with full profiles from the 3' end.
        trimmed_trunc_seq_names = []
        trimmed_trunc_reversed_seq_strings = []
        trimmed_trunc_seqs = []
        for name, trimmed_seq in self.trimmed_trunc_seq_dict.items():
            trimmed_trunc_seq_names.append(name)
            trimmed_trunc_reversed_seq_strings.append(trimmed_seq.seq_string[::-1])
            trimmed_trunc_seqs.append(trimmed_seq)
        norm_seq_names = []
        norm_reversed_seq_strings = []
        norm_seqs = []
        for name, norm_seq in self.norm_trna_seq_dict.items():
            norm_seq_names.append(name)
            norm_reversed_seq_strings.append(norm_seq.seq_string[::-1])
            norm_seqs.append(norm_seq)
        clusters = Dereplicator(trimmed_trunc_seq_names + norm_seq_names,
                                trimmed_trunc_reversed_seq_strings + norm_reversed_seq_strings,
                                extras=trimmed_trunc_seqs + norm_seqs).prefix_dereplicate()

        # Associate each truncated sequence with any normalized tRNA sequences that contain it as a
        # 3' subsequence. Since the trimmed truncated sequence can be a 3' subsequence of multiple
        # normalized sequences, do not reconstruct a feature profile for the truncated sequence.
        # (Similarly, trimmed sequences with full profiles need not have the same profile as the
        # seed trimmed sequence of a normalized sequence.) The unique sequences in the trimmed
        # sequence are therefore not included in the features table of the tRNA-seq database.

        # Clusters cannot contain more than one normalized tRNA sequence, as these have already been
        # 3'-dereplicated (by definition). Normalized sequences can seed clusters and also be
        # members of clusters seeded by a truncated sequence. In the latter case, truncated
        # sequences in the cluster that are shorter than the normalized sequence (3' subsequences)
        # are incorporated as members of the normalized sequence.

        # There are three types of cluster: 1. clusters consisting of a single normalized sequence
        # (ignore), 2. clusters containing a normalized sequence as seed or member with shorter
        # truncated sequence members, and 3. clusters seeded by truncated sequences. If a truncated
        # sequence is found in group 2 (part of one or more longer normalized sequences) then ignore
        # it in group 3 (do not include it in normalized truncated sequences formed from group 3
        # clusters). The alternatives do not make sense -- including the truncated sequence in
        # normalized sequences with truncated but not full profiles, or withholding it from both,
        # perhaps as a new category of sequence.

        # This dict relates trimmed truncated sequences to normalized sequences containing them.
        trunc_seq_norm_seq_dict = defaultdict(list)
        # This dict relates trimmed truncated sequences to other trimmed truncated sequences found
        # to be subsequences of the former.
        trunc_seq_trunc_seq_dict = {}
        for cluster in clusters:
            if len(cluster.member_names) == 1:
                if isinstance(cluster.member_extras[0], NormalizedSequence):
                    continue

            norm_trna_seq = None
            trimmed_trunc_seq_seed = cluster.member_extras[0] if isinstance(cluster.member_extras[0], TrimmedTruncatedProfileSequence) else None
            if trimmed_trunc_seq_seed:
                trunc_seed_seq_members = []
                trunc_seq_trunc_seq_dict[trimmed_trunc_seq_seed.represent_name] = (trimmed_trunc_seq_seed, trunc_seed_seq_members)
            for seq in cluster.member_extras:
                # Members of each cluster are pre-sorted in descending order of sequence length.
                # There cannot be a normalized tRNA and trimmed truncated profile sequence of the
                # same length.
                if norm_trna_seq:
                    if isinstance(seq, TrimmedTruncatedProfileSequence):
                        trunc_seq_norm_seq_dict[seq.represent_name].append(norm_trna_seq)
                        continue
                    else:
                        raise ConfigError("It appears that a cluster in the 3' dereplication "
                                          "of trimmed sequences with truncated profiles and normalized sequences with full profiles "
                                          "contains >1 normalized sequence, when it should only contain 0 or 1.")

                if isinstance(seq, NormalizedSequence):
                    norm_trna_seq = seq
                    continue

                # The cluster is seeded by a truncated sequence.
                trunc_seed_seq_members.append(seq)

        # Add truncated sequences to matching normalized tRNA sequences.
        trimmed_trunc_seq_dict = self.trimmed_trunc_seq_dict
        trimmed_trna_seq_dict = self.trimmed_trna_seq_dict
        uniq_trunc_seq_dict = self.uniq_trunc_seq_dict
        uniq_trna_seq_dict = self.uniq_trna_seq_dict
        # The count of recovered truncated sequences containing an anticodon is recorded to
        # facilitate measurement of isoacceptor abundances. A truncated profile may stop 3' of the
        # anticodon, but the anticodon may still be in the sequence. The presence of the anticodon
        # is therefore inferred from the longest trimmed tRNA sequence in the matching normalized
        # sequence.
        norm_trna_seq_anticodon_dict = {} # This saves time finding the position of the anticodon relative to the 3' terminus in normalized sequences
        for trimmed_trunc_seq_name, norm_trna_seqs in trunc_seq_norm_seq_dict.items():
            trimmed_trunc_seq = trimmed_trunc_seq_dict.pop(trimmed_trunc_seq_name)
            # The sequence with a truncated tRNA profile has been confirmed as tRNA, so transfer the
            # object between dictionaries.
            trimmed_trna_seq_dict[trimmed_trunc_seq.represent_name] = trimmed_trunc_seq
            trimmed_trunc_seq_length = len(trimmed_trunc_seq.seq_string)

            for uniq_trunc_seq in trimmed_trunc_seq.uniq_seqs:
                # The sequence with a truncated tRNA profile has been confirmed as tRNA.
                uniq_trna_seq_dict[uniq_trunc_seq.represent_name] = uniq_trunc_seq_dict.pop(uniq_trunc_seq.represent_name)

            for norm_seq in norm_trna_seqs:
                norm_seq.trimmed_seqs.append(trimmed_trunc_seq)
                norm_seq_length = len(norm_seq.seq_string)
                norm_seq.start_positions.append(norm_seq_length - trimmed_trunc_seq_length)
                norm_seq.stop_positions.append(norm_seq_length)
                trimmed_trunc_seq.norm_seq_represent_names.append(norm_seq.represent_name)

                if not trimmed_trunc_seq.contains_anticodon:
                    # Determine from the first normalized sequence in which the truncated sequence is
                    # found whether the truncated sequence contains the anticodon.
                    try:
                        # The position of the anticodon in the normalized sequence has already been found.
                        anticodon_start_relative_to_threeprime_terminus = norm_trna_seq_anticodon_dict[norm_seq.represent_name]
                    except KeyError:
                        try:
                            anticodon_loop_start = norm_seq.trimmed_seqs[0].feature_start_indices[self.RELATIVE_ANTICODON_LOOP_INDEX]
                        except IndexError:
                            # The anticodon loop was not reached in the profile.
                            anticodon_loop_start = -1
                        if anticodon_loop_start > -1:
                            # The anticodon loop was profiled.
                            anticodon_start = anticodon_loop_start + 2
                            # The position of the anticodon relative to the 3' terminus is a negative number.
                            anticodon_start_relative_to_threeprime_terminus = anticodon_start - norm_seq.trimmed_seqs[0].uniq_seqs[0].feature_start_indices[-1]
                        else:
                            # The anticodon loop was not profiled, indicated by a positive number.
                            anticodon_start_relative_to_threeprime_terminus = 1
                        norm_trna_seq_anticodon_dict[norm_seq.represent_name] = anticodon_start_relative_to_threeprime_terminus
                    if anticodon_start_relative_to_threeprime_terminus == 1:
                        continue
                    if trimmed_trunc_seq_length + anticodon_start_relative_to_threeprime_terminus >= 0:
                        trimmed_trunc_seq.contains_anticodon = True
                        for uniq_seq in trimmed_trunc_seq.uniq_seqs:
                            uniq_seq.contains_anticodon = True

        # Trimmed truncated sequences that don't match normalized tRNA sequences are grouped into
        # normalized truncated sequences.
        for trimmed_trunc_seq_seed_represent_name, entry in trunc_seq_trunc_seq_dict.items():
            trimmed_trunc_seq_seed, trimmed_trunc_seq_members = entry

            # REMOVE
            if trimmed_trunc_seq_seed_represent_name in trunc_seq_norm_seq_dict:
                raise ConfigError("How is a trimmed truncated seed a part of a normalized tRNA sequence?")

            # If a trimmed truncated sequence is a 3' subsequence of a normalized tRNA sequence,
            # then it and all shorter sequences should be excluded from the new normalized truncated
            # sequence.
            for i, trimmed_trunc_seq in enumerate(trimmed_trunc_seq_members):
                if trimmed_trunc_seq.represent_name in trunc_seq_norm_seq_dict:
                    trimmed_trunc_seq_members = trimmed_trunc_seq_members[: i]
                    break

            self.norm_trunc_seq_dict[trimmed_trunc_seq_seed_represent_name] = NormalizedTruncatedProfileSequence(trimmed_trunc_seq_members)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed recovering tRNA with truncated feature profile (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def write_checkpoint_files(self, checkpoint_name):
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
            if intermed_file_path:
                self.run.info("Overwrote \"%s\" checkpoint intermediate file of %s" % (checkpoint_name, intermed_file_label), intermed_file_path)
            else:
                self.run.info_single("Wrote \"%s\" checkpoint intermediate file of %s" % (checkpoint_name, intermed_file_label), mc='green', cut_after=200)
        print()


    def load_checkpoint_files(self, checkpoint_name):
        self.progress.new(f"Loading intermediate files at the checkpoint, \"{checkpoint_name}\"")
        self.progress.update("...")

        with open(self.analysis_summary_path, 'a') as f:
            f.write(f"\nAnalysis restarted from the checkpoint, \"{checkpoint_name}\"\n")

        for intermed_file_key, intermed_file_path in self.intermed_file_path_dict[checkpoint_name].items():
            with open(intermed_file_path, 'rb') as f:
                setattr(self, intermed_file_key, pkl.load(f))

        self.progress.end()


    def threeprime_dereplicate_sequences_without_terminus(self):
        """Find tRNA sequences missing a 3' terminus. Sequences required a 3' terminus to have been
        profiled as tRNA. Unique non-tRNA sequences are searched against normalized tRNA sequences.
        Non-tRNA is recovered as tRNA when it is a 3' subsequence of normalized tRNA or is longer
        than normalized tRNA sequences with a complete feature profile and thus is shown to have a
        5' extension. Recovered sequences each generate a `TrimmedMappedSequence` object."""
        start_time = time.time()
        self.progress.new("Dereplicating tRNA seqs ending in discriminator nt")
        self.progress.update("...")

        # 3'-dereplicate normalized tRNA sequences and unique non-tRNA sequences.
        names = []
        reversed_seq_strings = []
        extras = []
        # The normalized sequences are added first to the dereplicator so that they appear first in
        # the clusters. This allows unique non-tRNA sequences that are identical to the normalized
        # sequence (due to prior trimming of the 3' terminus from the normalized sequence) to always
        # be recovered.
        for norm_seq_represent_name, norm_seq in self.norm_trna_seq_dict.items():
            names.append(norm_seq_represent_name)
            reversed_seq_strings.append(norm_seq.seq_string[::-1])
            extras.append(norm_seq)
        for uniq_nontrna_seq_represent_name, uniq_nontrna_seq in self.uniq_nontrna_seq_dict.items():
            if len(uniq_nontrna_seq.seq_string) >= self.min_trna_frag_size:
                names.append(uniq_nontrna_seq_represent_name)
                reversed_seq_strings.append(uniq_nontrna_seq.seq_string[::-1])
                extras.append(uniq_nontrna_seq)
        clusters = Dereplicator(names, reversed_seq_strings, extras=extras).prefix_dereplicate()

        # Search clusters for normalized tRNA sequences and qualifying unique non-tRNA sequences.
        # Some unique sequences may be 3' subsequences of more than one normalized sequence. Unique
        # sequences cannot be longer than more than one normalized sequence, as normalized sequences
        # have already been dereplicated. The same normalized sequence can be found in multiple
        # clusters as subsequences of different longer unique sequences.
        uniq_seq_norm_seqs_dict = defaultdict(list)
        uniq_nontrna_seq_dict = self.uniq_nontrna_seq_dict
        uniq_trna_seq_dict = self.uniq_trna_seq_dict
        trimmed_trna_seq_dict = self.trimmed_trna_seq_dict
        for cluster in clusters:
            if len(cluster.member_seqs) == 1:
                continue

            # Check that there is a normalized tRNA sequence in the cluster -- there cannot be more
            # than one.
            cluster_norm_seq_index = None
            norm_seq = None
            norm_seq_length = None
            norm_seq_has_complete_feature_set = None
            for member_index, seq in enumerate(cluster.member_extras):
                if isinstance(seq, NormalizedFullProfileSequence):
                    cluster_norm_seq_index = member_index
                    norm_seq = seq
                    norm_seq_length = len(norm_seq.seq_string)
                    if member_index > 0:
                        norm_seq_has_complete_feature_set = norm_seq.has_complete_feature_set
                    break
            else:
                continue

            # To reach this point, a normalized tRNA sequence must have been found in the cluster.
            # Now process any longer unique sequences.
            for uniq_nontrna_seq in cluster.member_extras[: cluster_norm_seq_index]:
                # If the normalized sequence has a complete feature profile (is a full-length tRNA),
                # then the overhanging 5' bases in the recovered unique sequence can be trimmed as
                # "extra" 5' bases. Otherwise, it is possible the overhanging bases are part of an
                # artifact joined to the 5' end of a tRNA fragment, so conservatively ignore the
                # unique sequence.
                if not norm_seq_has_complete_feature_set:
                    break

                uniq_seq_represent_name = uniq_nontrna_seq.represent_name
                try:
                    uniq_nontrna_seq_dict.pop(uniq_seq_represent_name)
                except KeyError:
                    # The unique sequence has already been recovered in another cluster, necessarily
                    # as a supersequence of the same normalized sequence.
                    continue

                uniq_trna_seq = UniqueMappedSequence(uniq_nontrna_seq.seq_string,
                                                     uniq_seq_represent_name,
                                                     uniq_nontrna_seq.read_count,
                                                     extra_fiveprime_length=len(uniq_nontrna_seq.seq_string) - norm_seq_length)
                uniq_trna_seq_dict[uniq_seq_represent_name] = uniq_trna_seq

                trimmed_seq = TrimmedMappedSequence(uniq_trna_seq)
                norm_seq.trimmed_seqs.append(trimmed_seq)
                trimmed_seq.norm_seq_represent_names.append(norm_seq.represent_name)
                norm_seq.start_positions.append(0)
                norm_seq.stop_positions.append(norm_seq_length)
                trimmed_trna_seq_dict[uniq_seq_represent_name] = trimmed_seq

            # Find all unique sequences in the cluster equal to in length or shorter than the
            # normalized tRNA sequence.
            for uniq_nontrna_seq in cluster.member_extras[cluster_norm_seq_index + 1: ]:
                uniq_seq_norm_seqs_dict[uniq_nontrna_seq.represent_name].append(norm_seq)

        for uniq_seq_represent_name, norm_seqs in uniq_seq_norm_seqs_dict.items():
            uniq_nontrna_seq = uniq_nontrna_seq_dict.pop(uniq_seq_represent_name)

            uniq_trna_seq = UniqueMappedSequence(uniq_nontrna_seq.seq_string,
                                                 uniq_seq_represent_name,
                                                 uniq_nontrna_seq.read_count)
            uniq_trna_seq_dict[uniq_seq_represent_name] = uniq_trna_seq

            trimmed_seq = TrimmedMappedSequence(uniq_trna_seq)
            uniq_seq_length = len(uniq_trna_seq.seq_string)
            # The same normalized sequence can be found in multiple clusters, so it can be
            # represented multiple times in `norm_seqs`.
            for norm_seq in set(norm_seqs):
                norm_seq.trimmed_seqs.append(trimmed_seq)
                trimmed_seq.norm_seq_represent_names.append(norm_seq.represent_name)
                norm_seq_length = len(norm_seq.seq_string)
                norm_seq.start_positions.append(norm_seq_length - uniq_seq_length)
                norm_seq.stop_positions.append(norm_seq_length)
            trimmed_trna_seq_dict[uniq_seq_represent_name] = trimmed_seq

        self.progress.end()


    def map_fragments(self):
        """Map unprofiled tRNA fragments to longer profiled tRNA sequences. Fragments only missing a
        3' terminus were already found by mapping with
        `threeprime_dereplicate_sequences_without_terminus` using an efficient shortcut or by profiling
        if '' was an accepted 3' terminus.

        EXAMPLE:
        Normalized tRNA:                 (GT)TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Mapped tRNA 1 (extra 5' bases) :   T TCCGTGATAGTTTAATGGTCAGAATGG
        Mapped tRNA 2 (interior)       :            TAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGG
        """
        start_time = time.time()

        pid = "Set up search of unprofiled reads to profiled tRNA"
        self.progress.new(pid)

        self.progress.update("Getting queries from unprofiled reads")
        temp_dir_path = filesnpaths.get_temp_directory_path()
        query_fasta_path = os.path.join(temp_dir_path, 'query.fa')

        uniq_nontrna_seq_dict = self.uniq_nontrna_seq_dict
        min_trna_frag_size = self.min_trna_frag_size
        with open(query_fasta_path, 'w') as query_fasta:
            for uniq_seq in [uniq_seq for uniq_seq in uniq_nontrna_seq_dict.values()
                             if len(uniq_seq.seq_string) >= min_trna_frag_size]:
                # Include the length of the sequence in the defline for the purposes of parsing
                # vmatch output.
                query_fasta.write(f">{uniq_seq.represent_name}-{len(uniq_seq.seq_string)}\n"
                                  f"{uniq_seq.seq_string}\n")


        self.progress.update_pid(pid)
        self.progress.update("Getting targets from profiled tRNAs")
        # Non-tRNA sequences are mapped to normalized tRNA sequences with extra 5' bases added when
        # present in underlying unique tRNA sequences. Multiple targets for each normalized sequence
        # are therefore produced for different 5' sequence extensions.
        target_fasta_path = os.path.join(temp_dir_path, 'target.fa')
        with open(target_fasta_path, 'w') as target_fasta:
            for norm_seq in self.norm_trna_seq_dict.values():
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
                        target_fasta.write(f">{norm_seq.represent_name}-{len(fiveprime_seq_string)}-{fiveprime_index}\n"
                                           f"{fiveprime_seq_string}{norm_seq_string}\n")
                else:
                    target_fasta.write(f">{norm_seq.represent_name}-0-0\n" # no extra 5' bases
                                       f"{norm_seq_string}\n")
        self.progress.end()


        # Use a 10x bigger query chunk size than the Vmatch default, as the default is rather
        # conservative for searches with mismatches/indels that take longer to process each chunk
        # and that may generate more alignments per chunk.
        match_df = Vmatch(argparse.Namespace(match_mode='exact_query_substring',
                                             fasta_db_file=target_fasta_path,
                                             fasta_query_file=query_fasta_path,
                                             num_threads=self.num_threads,
                                             query_chunk_size=1000000 // self.num_threads,
                                             temp_dir=temp_dir_path)).search_queries()

        pid = "Filtering matches"
        self.progress.new(pid)
        self.progress.update("...")

        self.restructure_fragment_match_table(match_df)

        # Process each unique sequence match. Each unique sequence can match more than one
        # normalized sequence.
        match_gb = match_df.groupby('query_name')
        del match_df
        gc.collect()

        fragment_filter_progress_interval = 25000
        total_matched_queries = len(match_gb)
        pp_total_matched_queries = pp(total_matched_queries)
        num_filtered_queries = -1
        uniq_trna_seq_dict = self.uniq_trna_seq_dict
        trimmed_trna_seq_dict = self.trimmed_trna_seq_dict
        norm_trna_seq_dict = self.norm_trna_seq_dict
        for uniq_seq_name, query_match_df in match_gb:
            num_filtered_queries += 1
            if num_filtered_queries % fragment_filter_progress_interval == 0:
                pp_progress_interval_end = pp(total_matched_queries if num_filtered_queries + fragment_filter_progress_interval > total_matched_queries else num_filtered_queries + fragment_filter_progress_interval)
                self.progress.update_pid(pid)
                self.progress.update(f"Queries {num_filtered_queries + 1}-{pp_progress_interval_end}/{pp_total_matched_queries}")

            # Each unique sequence with a validated match will yield a `UniqueMappedSequence` and
            # `TrimmedMappedSequence`.
            uniq_trna_seq = None
            trimmed_seq = None

            for norm_seq_name, target_fiveprime_length, query_start, uniq_seq_length in zip(query_match_df['target_name'],
                                                                                            query_match_df['fiveprime_length'],
                                                                                            query_match_df['query_start_in_target'],
                                                                                            query_match_df['query_length']):
                query_stop = query_start + uniq_seq_length
                trimmed_seq_stop_in_norm_seq = query_stop - target_fiveprime_length

                if trimmed_seq_stop_in_norm_seq <= 0:
                    # Ignore queries that align entirely to extra 5' bases. Sequences mapping
                    # exclusively to the 5' extension that are long enough to fulfill the minimum
                    # length requirement may be mapping to an artifactual chimeric sequence.
                    continue

                norm_seq = norm_trna_seq_dict[norm_seq_name]

                if not uniq_trna_seq:
                    # Enter this block the first time the unique sequence query validly matches a
                    # normalized sequence.
                    uniq_nontrna_seq = uniq_nontrna_seq_dict.pop(uniq_seq_name)

                    # Assume that 5' extensions are the same for the query regardless of the reference.
                    # This could be false when
                    # 1. tRNA profiling erroneously identified the end of the acceptor stem
                    # or 2. the query mapped to different places at the end of the acceptor stem in different tRNAs.
                    if target_fiveprime_length - query_start > 0:
                        query_fiveprime_length = target_fiveprime_length - query_start
                        trimmed_seq_start_in_norm_seq = 0
                    else:
                        query_fiveprime_length = 0
                        trimmed_seq_start_in_norm_seq = query_start - target_fiveprime_length

                    uniq_trna_seq = UniqueMappedSequence(uniq_nontrna_seq.seq_string,
                                                         uniq_seq_name,
                                                         uniq_nontrna_seq.read_count,
                                                         extra_fiveprime_length=query_fiveprime_length)
                    uniq_trna_seq_dict[uniq_seq_name] = uniq_trna_seq

                    trimmed_seq = TrimmedMappedSequence(uniq_trna_seq)
                    norm_seq.trimmed_seqs.append(trimmed_seq)
                    trimmed_seq.norm_seq_represent_names.append(norm_seq.represent_name)
                    norm_seq.start_positions.append(trimmed_seq_start_in_norm_seq)
                    norm_seq.stop_positions.append(trimmed_seq_stop_in_norm_seq)
                    trimmed_trna_seq_dict[uniq_seq_name] = trimmed_seq
                else:
                    for prev_trimmed_seq in norm_seq.trimmed_seqs[::-1]:
                        # Ensure that the trimmed sequence maps to the normalized sequence only
                        # once. Multiple targets can be created from the same normalized sequence
                        # for different 5' extensions. Trimmed mapped sequences are added after
                        # trimmed profiled sequences to the normalized sequence's list of trimmed
                        # sequences.
                        if trimmed_seq.represent_name == prev_trimmed_seq.represent_name:
                            break

                        if not isinstance(prev_trimmed_seq, TrimmedMappedSequence):
                            norm_seq.trimmed_seqs.append(trimmed_seq)
                            trimmed_seq.norm_seq_represent_names.append(norm_seq.represent_name)
                            if target_fiveprime_length - query_start > 0:
                                trimmed_seq_start_in_norm_seq = 0
                            else:
                                trimmed_seq_start_in_norm_seq = query_start - target_fiveprime_length
                            norm_seq.start_positions.append(trimmed_seq_start_in_norm_seq)
                            norm_seq.stop_positions.append(trimmed_seq_stop_in_norm_seq)
                            break
        self.progress.end()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed mapping tRNA fragments (min)", time.time() - start_time, is_time_value=True))


    def restructure_fragment_match_table(self, match_df):
        """Helper method for `map_fragments`."""
        uniq_seq_names = []
        for query_name in match_df['query_name']:
            # Sequence names in anvi'o cannot contain a hyphen.
            uniq_seq_name, uniq_seq_length = query_name.split('-')
            uniq_seq_names.append(uniq_seq_name)
        match_df.loc[:, 'query_name'] = uniq_seq_names

        norm_seq_names = []
        fiveprime_lengths = []
        for target_name in match_df['target_name']:
            norm_seq_name, fiveprime_length, fiveprime_index = target_name.split('-')
            norm_seq_names.append(norm_seq_name)
            fiveprime_lengths.append(int(fiveprime_length))
        match_df.loc[:, 'target_name'] = norm_seq_names
        match_df['fiveprime_length'] = fiveprime_lengths


    def find_substitutions(self):
        """Find sites of potential modification-induced substitutions."""
        start_time = time.time()
        pid = "Finding modification-induced substitutions"
        self.progress.new(pid)
        self.progress.update("...")

        # Cluster normalized tRNA sequences. Clusters agglomerate sequences that differ from at
        # least one other sequence in the cluster by no more than 3 nucleotides in 100 (by default)
        # in a gapless end-to-end alignment with no clipping.
        norm_trna_seq_dict = self.norm_trna_seq_dict
        norm_seq_represent_names = []
        norm_seq_strings = []
        norm_seq_feature_completeness_dict = {}
        for represent_name, norm_seq in norm_trna_seq_dict.items():
            norm_seq_represent_names.append(represent_name)
            norm_seq_strings.append(norm_seq.seq_string)
            norm_seq_feature_completeness_dict[represent_name] = norm_seq.has_complete_feature_set
        self.progress.end()

        agglomerator = Agglomerator(norm_seq_represent_names, norm_seq_strings, num_threads=self.num_threads)
        # Provide a priority function for seeding clusters that favors, in order:
        # 1. normalized sequences with a complete set of tRNA features,
        # 2. longer normalized sequences,
        # 3. normalized sequences with more alignments in the all-against-all search,
        # 4. alphanumeric order of the normalized sequence representative name.
        agglomerator.agglomerate(max_mismatch_freq=self.agglom_max_mismatch_freq,
                                 priority_function=lambda aligned_ref: (-norm_seq_feature_completeness_dict[aligned_ref.name],
                                                                        -len(aligned_ref.seq_string),
                                                                        -len(aligned_ref.alignments),
                                                                        aligned_ref.name))
        agglom_aligned_ref_dict = agglomerator.agglom_aligned_ref_dict

        pid = "Decomposing clusters"
        self.progress.new(pid)
        self.progress.update("...")

        excluded_norm_seq_names = [] # Used to exclude normalized sequences from being considered as aligned queries in clusters (see below)
        represent_norm_seq_names = [] # Used to prevent the same modified sequence from being created twice
        mod_trna_seq_dict = self.mod_trna_seq_dict
        num_processed_refs = -1
        total_ref_count = len(agglom_aligned_ref_dict)
        decomposition_progress_interval = 1000
        pp_total_ref_count = pp(total_ref_count)
        for ref_name, aligned_ref in agglom_aligned_ref_dict.items():
            num_processed_refs += 1
            if num_processed_refs % decomposition_progress_interval == 0:
                pp_progress_interval_end = pp(total_ref_count if num_processed_refs + decomposition_progress_interval > total_ref_count else num_processed_refs + decomposition_progress_interval)
                self.progress.update_pid(pid)
                self.progress.update(f"{pp(num_processed_refs + 1)}-{pp_progress_interval_end}/{pp_total_ref_count}")

            # A modification requires at least 3 different nucleotides to be detected, and each
            # normalized sequence differs by at least 1 nucleotide (substitution or gap), so for a
            # cluster to form a modified sequence, it must contain at least 3 normalized sequences.
            if len(aligned_ref.alignments) < 2:
                continue

            aligned_ref_length = len(aligned_ref.seq_string)
            valid_aligned_queries = []
            for alignment in aligned_ref.alignments:
                # Normalized tRNA sequences should only align at the 3' end. Alignments to the
                # interior of the sequence can theoretically occur when the reference is a tRNA-tRNA
                # chimera.
                if aligned_ref_length != alignment.target_start + alignment.alignment_length:
                    continue

                query_name = alignment.aligned_query.name
                # The normalized sequence query may have formed a modified sequence already. If the
                # normalized sequence had a complete feature profile, or if it was the same length
                # as such a sequence, then it should not be able to form a longer modified sequence
                # that would have 5' nucleotides beyond the end of a complete feature profile.
                if query_name in excluded_norm_seq_names:
                    continue

                valid_aligned_queries.append(norm_trna_seq_dict[query_name])

            # Confirm that 2 or more queries passed the filters, so at least 3 normalized sequences
            # are still in the cluster.
            if len(valid_aligned_queries) < 2:
                continue

            valid_aligned_queries.sort(key=lambda norm_seq: (-len(norm_seq.seq_string), -norm_seq.has_complete_feature_set, norm_seq.represent_name))

            seq_array = np.zeros((len(valid_aligned_queries) + 1, aligned_ref_length), dtype=int)
            # Rather than using the ASCII representation of each character, which saves some time in
            # converting the sequence string to a numpy array, constrain the integer representation
            # to the smallest possible range of integers to speed up the bincount method used to
            # determine the number of unique nucleotides at an alignment position.
            seq_array[0, :] += [NT_INT_DICT[nt] for nt in aligned_ref.seq_string]
            for i, aligned_query in enumerate(valid_aligned_queries, start=1):
                seq_array[i, aligned_ref_length - len(aligned_query.seq_string): ] += [NT_INT_DICT[nt] for nt in aligned_query.seq_string]

            norm_seqs = np.array([norm_trna_seq_dict[ref_name]] + valid_aligned_queries)

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
                next_clusters = deque() # Make a new object with each iteration rather than clearing the same one

                while clusters:
                    seq_array, norm_seqs, three_four_nt_alignment_positions = clusters.pop()

                    # A modification requires at least 3 different nucleotides to be detected, and
                    # each normalized sequence differs by at least 1 nucleotide, so for a cluster to
                    # form a modified sequence, it must contain at least 3 normalized sequences.
                    if norm_seqs.size < 3:
                        continue

                    aligned_nts = seq_array[:, alignment_pos]
                    nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]

                    if (nt_counts != 0).sum() < 2:
                        # There are now < 2 nucleotides at the alignment position in the cluster
                        # under consideration. 2 different nucleotides are needed to distinguish
                        # single nucleotide variants.
                        next_clusters.appendleft((seq_array, norm_seqs, three_four_nt_alignment_positions))
                        continue

                    # Add a new cluster for each single nucleotide variant to the stack of clusters
                    # to process if: 1. the new cluster contains at least 3 sequences and 2. the
                    # longest normalized sequence (with a complete feature profile, if applicable)
                    # in the new cluster has not yet formed a modified sequence. Agglomerative
                    # clustering ensures that the sequences agglomerated with the longest normalized
                    # sequence will be the same regardless of the original unsplit cluster.
                    represented_nts = nt_counts.nonzero()[0] + 1
                    for nt in represented_nts:
                        split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]

                        if split_cluster_seq_indices.size > 2:
                            split_cluster_norm_seqs = norm_seqs[split_cluster_seq_indices]

                            if split_cluster_norm_seqs[0].represent_name in represent_norm_seq_names:
                                continue

                            next_clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                      split_cluster_norm_seqs,
                                                      three_four_nt_alignment_positions.copy()))
                if next_clusters:
                    clusters = next_clusters
                else:
                    break
            if not clusters:
                continue

            # Check alignment positions previously found to have 3-4 nucleotides. Further split
            # clusters when positions now have 2 nucleotides.
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

                            # At least 3 normalized sequences are needed, and the split cluster
                            # cannot have already formed a modified sequence.
                            if split_cluster_seq_indices.size > 2:
                                split_cluster_norm_seqs = norm_seqs[split_cluster_seq_indices]

                                if split_cluster_norm_seqs[0].represent_name in represent_norm_seq_names:
                                    continue

                                clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                     split_cluster_norm_seqs,
                                                     np.delete(three_four_nt_alignment_positions, candidates_to_remove)))
                        # Reevaluate previous alignment positions in the split clusters.
                        break
                else:
                    # At least 1 position was discounted as no longer having 3-4 different
                    # nucleotides, but these positions had fewer than 2 nucleotides, and so did not
                    # cause the cluster to be split into new clusters. Therefore, do not cycle
                    # through the remaining positions again to find any more with fewer than 3
                    # nucleotides.
                    if candidates_to_remove:
                        next_clusters.appendleft((norm_seqs, np.delete(three_four_nt_alignment_positions, candidates_to_remove)))
                    else:
                        next_clusters.appendleft((norm_seqs, three_four_nt_alignment_positions))

            if not next_clusters:
                continue
            clusters = next_clusters

            while clusters:
                norm_seqs, mod_positions = clusters.pop() # Normalized sequences should have retained their order
                norm_seqs = list(norm_seqs) # Turn the array into a list
                represent_norm_seq = norm_seqs[0]

                represent_norm_seq_length = len(represent_norm_seq.seq_string)
                represent_norm_seq_start_in_array = aligned_ref_length - represent_norm_seq_length
                mod_positions -= represent_norm_seq_start_in_array
                mod_seq = ModifiedSequence(norm_seqs, mod_positions.tolist())

                if represent_norm_seq.has_complete_feature_set:
                    for norm_seq in norm_seqs:
                        if len(norm_seq.seq_string) < represent_norm_seq_length:
                            break
                        excluded_norm_seq_names.append(norm_seq.represent_name)

                mod_trna_seq_dict[mod_seq.represent_name] = mod_seq

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed finding modification-induced substitutions (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def find_deletions(self):
        """Find potential modification-induced deletions. First, in silico "test" deletions are
        introduced at and around substitution sites in sequences with potential modification-induced
        substitutions. Notate such modified sequences with in silico deletions as *M'*. M' are
        searched against two pools of sequence targets:
        1. normalized tRNA sequences not assigned to modified sequences, notated *Nf*, for
           normalized sequences with a full feature profile, and
        2. normalized "non-tRNA" sequences with truncated tRNA profiles, notated *Nt*.

        Why these two pools?
        1. Why are Nf considered at all, when they have been successfully profiled, and therefore,
           presumably, do not have deletions interrupting the profile? Deletions can be erroneously
           accommodated by flexibility in feature lengths. For example, a deletion associated with a
           modification in the D loop can cause the variable-length alpha or beta sections of the D
           loop to be assigned one fewer nucleotide than is correct in order to optimize the
           profile.
        2. Why not all "non-tRNAs," why just those with a truncated feature profile (Nt)? Deletions
           can cause truncation of the profile. "Non-tRNAs" without even a truncated profile
           (sequences that were not profiled past the minimum length threshold of the T arm) also
           have fewer opportunities for modification-induced mutations.

        Normalized sequences rather than trimmed or unique sequences are searched for the sake of
        speed and simplicity. Ideally, normalized sequences found to have deletions would be further
        processed, finding which of their constituent trimmed and unique sequences actually contain
        the deletions. However, *nonspecific* trimmed and unique sequences are, by definition, in
        other normalized sequences, creating the possible ambiguity that the trimmed sequence would
        be marked as having a deletion in one but not another normalized sequence. This would not
        necessarily be an error, as identical underlying reads could theoretically originate from
        different cDNA sequences, with some containing a deletion, and others, representing a
        different tRNA, not containing a deletion."""
        start_time = time.time()
        pid = "Finding seqs with mod-induced dels"
        self.progress.new(pid)

        self.progress.update("Generating mod seqs with in silico dels")
        # A "child" here is a modified sequence with in silico deletions, M'.
        mod_seq_child_names = []
        mod_seq_child_reversed_seq_strings = []
        mod_seq_child_extras = []
        for mod_seq in self.mod_trna_seq_dict.values():
            mod_seq_child_index = 0
            for seq_string_with_del, del_config in self.get_sequences_with_deletions(mod_seq):
                mod_seq_child_names.append(mod_seq.represent_name + '_' + str(mod_seq_child_index))
                mod_seq_child_reversed_seq_strings.append(seq_string_with_del[::-1])
                mod_seq_child_extras.append((del_config, mod_seq))
                mod_seq_child_index += 1


        self.progress.update_pid(pid)
        self.progress.update("Gathering norm seqs with trunc profiles")
        norm_trunc_seq_names = []
        norm_trunc_seq_reversed_seq_strings = []
        norm_trunc_seq_extras = []
        for norm_trunc_seq_name, norm_trunc_seq in self.norm_trunc_seq_dict.items():
            norm_trunc_seq_names.append(norm_trunc_seq_name)
            norm_trunc_seq_reversed_seq_strings.append(norm_trunc_seq.seq_string[::-1])
            norm_trunc_seq_extras.append(norm_trunc_seq)

        self.progress.update_pid(pid)
        self.progress.update("3'-derep seqs with in silico dels & seqs with trunc profiles")
        clusters = self.prefix_dereplicate_deletion_candidates(norm_trunc_seq_names,
                                                               norm_trunc_seq_reversed_seq_strings,
                                                               norm_trunc_seq_extras,
                                                               mod_seq_child_names,
                                                               mod_seq_child_reversed_seq_strings,
                                                               mod_seq_child_extras)

        if clusters:
            self.progress.update_pid(pid)
            self.progress.update("Searching norm seqs with trunc profiles")
            self.process_deletion_clusters(clusters, 'trunc')


        self.progress.update_pid(pid)
        self.progress.update("Gathering norm seqs with full profiles")
        norm_trna_seq_names = []
        norm_trna_seq_reversed_seq_strings = []
        norm_trna_seq_extras = []
        for norm_trna_seq_name, norm_trna_seq in self.norm_trna_seq_dict.items():
            if not norm_trna_seq.mod_seqs:
                norm_trna_seq_names.append(norm_trna_seq_name)
                norm_trna_seq_reversed_seq_strings.append(norm_trna_seq.seq_string[::-1])
                norm_trna_seq_extras.append(norm_trna_seq)

        self.progress.update_pid(pid)
        self.progress.update("3'-derep seqs with in silico dels and seqs with full profiles")
        clusters = self.prefix_dereplicate_deletion_candidates(norm_trna_seq_names,
                                                               norm_trna_seq_reversed_seq_strings,
                                                               norm_trna_seq_extras,
                                                               mod_seq_child_names,
                                                               mod_seq_child_reversed_seq_strings,
                                                               mod_seq_child_extras)

        if clusters:
            self.progress.update_pid(pid)
            self.progress.update("Searching norm seqs with full profiles")
            self.process_deletion_clusters(clusters, 'trna')


        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed finding mod-induced dels (min)",
                                          time.time() - start_time,
                                          is_time_value=True))

        self.progress.end()


    def get_sequences_with_deletions(self, mod_seq):
        """Generate in silico modified sequences with deletions, *M'*, at and/or around substitution
        sites. This method is called by `find_deletions`.

        This method first introduces in silico nucleotide *substitutions* at possible substitution
        sites to create template sequences. Substitutions are only allowed to the nucleotides
        observed in the normalized sequences comprising the modified sequence, i.e., if 3
        nucleotides are observed at a substitution site, then 3 are introduced in silico, rather
        than the maximum 4. For each template sequence produced, `introduce_deletions` is called to
        introduce deletions of different lengths at the possible configurations of substitution
        sites. This may produce redundant sequences from the different templates, which are removed
        in the present method.

        Returns
        =======
        del_set : set
            A set of tuples, each with 2 elements. The first element is a sequence string containing
            deletions. The second element is a tuple of the indices of these deletions in the input
            sequence.
        """
        longest_norm_seq_string = mod_seq.norm_seqs_without_dels[0].seq_string
        mod_seq_length = len(longest_norm_seq_string)
        sub_positions = mod_seq.sub_positions

        # Record the configurations of nucleotides at the substitution positions found in the
        # normalized sequences comprising the modified sequence. For normalized sequences shorter
        # than the modified sequence (shorter than its longest normalized sequence) that do not
        # include certain substitution positions nearer the 5' end of the modified sequence, use the
        # nucleotides at the missing substitution positions from the representative longest
        # normalized sequence in the substitution configuration.
        longest_norm_seq_sub_nts = tuple([longest_norm_seq_string[sub_pos] for sub_pos in sub_positions])
        sub_nt_config_set = set([longest_norm_seq_sub_nts])
        for norm_seq in mod_seq.norm_seqs_without_dels[1: ]:
            norm_seq_string = norm_seq.seq_string
            norm_seq_start_in_mod_seq = mod_seq_length - len(norm_seq_string)
            norm_seq_sub_nts = []
            for sub_num, sub_pos in enumerate(sub_positions):
                if sub_pos < norm_seq_start_in_mod_seq:
                    norm_seq_sub_nts.append(longest_norm_seq_sub_nts[sub_num])
                else:
                    norm_seq_sub_nts.append(norm_seq_string[sub_pos - norm_seq_start_in_mod_seq])
            sub_nt_config_set.add(tuple(norm_seq_sub_nts))

        # Make template sequences (without deletions) from the substitution configurations.
        seq_strings_without_dels = []
        for sub_nt_config in sub_nt_config_set:
            altered_seq_string = longest_norm_seq_string
            for sub_nt, sub_pos in zip(sub_nt_config, sub_positions):
                altered_seq_string = altered_seq_string[: sub_pos] + sub_nt + altered_seq_string[sub_pos + 1: ]
            seq_strings_without_dels.append(altered_seq_string)

        del_set = set()
        for max_distinct_dels in range(self.max_distinct_dels, 0, -1):
            # Introduce deletions into each template sequence, potentially producing a number of new
            # sequences.
            del_dict = {}
            for seq_string in seq_strings_without_dels:
                del_dict_for_seq = self.introduce_deletions(seq_string, sub_positions, max_distinct_dels=max_distinct_dels)
                # The same sequence with deletions may sometimes be generated from different
                # template normalized sequences. In case of redundancy, favor the deletion
                # configuration closest to the substitution site around which deletions were
                # introduced, with equally close configurations resolved by favoring more 5' over
                # more 3' configurations.
                for seq_string_with_del, del_pos_info in del_dict_for_seq.items():
                    if seq_string_with_del in del_dict:
                        if del_pos_info[1] >= del_dict[seq_string_with_del][1]:
                            continue
                    del_dict[seq_string_with_del] = del_pos_info
            # Nota bene: Generated sequences may be 3' subsequences of each other. These are
            # 3'-dereplicated in `find_deletions`, which calls the present method.
            del_set = set([(seq_string_with_del, del_pos_info[0]) for seq_string_with_del, del_pos_info in del_dict.items()])
            if len(del_set) <= self.max_del_configs:
                # A very large number of sequences with in silico deletions was generated,
                # potentially preventing processing in a reasonable time, so decrement the number of
                # distinct sites at which deletions can be introduced in a single template sequence.
                break
        return del_set


    def introduce_deletions(self, seq_string, sub_positions, max_distinct_dels):
        """Generate in silico sequences with deletions at and/or around substitution sites in the
        input sequence. This method is called by `get_sequences_with_deletions`.

        Parameters
        ==========
        seq_string : str
            The sequence in which deletions will be introduced

        sub_positions : list-like
            Where substitutions are located in the input sequence

        Returns
        =======
        del_dict : dict
            Each dict key is a sequence string resulting from the introduction of deletions in the
            input sequence. Each dict value is a length-two tuple. The first element of the tuple is
            another tuple of the deletion positions in the input. The second element is a score used
            to measure the distance of the deletions from the substitution loci of the input
            sequence.
        """
        del_ranges = self.del_ranges
        min_fiveprime_del_pos = self.min_length_of_long_fiveprime_extension - 1
        min_dist_between_dels = self.min_dist_between_dels
        # Find all the ways deletions can be introduced into the sequence given the
        # parameterization.
        del_pos_configs = []
        seq_string_length = len(seq_string)
        # Deletions of different sizes can be situated at each substitution site. Deletions may be
        # found at one or multiple sites, if multiple substititutions are present. Call each
        # deletion site a "locus".
        for num_del_loci in range(1, max_distinct_dels + 1):
            # Example: `num_del_loci` of 2 means deletions will be introduced in 2 places.
            del_range_configs = product(*[del_ranges for _ in range(num_del_loci)])
            for del_locus_config in combinations(sub_positions, num_del_loci):
                # Example: `del_locus_config` of (10, 20) means deletions will be introduced around
                # the substitutions at positions 10 and 20 in the input sequence.
                for del_range_config in del_range_configs:
                    # Example: `del_range_config` of (range(-2, 0), range(0, 1)) means the 2
                    # nucleotides 5' of the first substitution position and the nucleotide at the
                    # second substitution position will be deleted.
                    all_del_positions = []
                    all_del_pos_scores = []
                    for del_locus_index, del_range in enumerate(del_range_config):
                        sub_pos = del_locus_config[del_locus_index]
                        locus_del_positions = []
                        locus_del_pos_scores = []
                        for del_pos_relative_to_sub in del_range:
                            del_pos = sub_pos + del_pos_relative_to_sub

                            # Deletion positions must be within the template sequence. There must be
                            # sufficient sequence length on the 5' end to confirm the deletion. This
                            # is related to the issue of nontemplated nucleotides. For example,
                            # there may be a *biological* 3' tRNA fragment that ends at a
                            # modification site, but a nontemplated nucleotide is added to the 5'
                            # end of the read and is the same as the nucleotide that occurs 3
                            # nucleotides 5' of the modification site, so without an "anchoring"
                            # stretch of nucleotides on the 5' end, a deletion of length 2 could be
                            # mistakenly identified as lying between the modification and the
                            # nontemplated nucleotide. Therefore, only allow silico deletions that
                            # occur some distance from the 5' end.
                            if not seq_string_length > del_pos >= min_fiveprime_del_pos:
                                continue

                            if not locus_del_positions:
                                if all_del_positions:
                                    if del_pos - all_del_positions[-1] < min_dist_between_dels:
                                        continue

                            locus_del_positions.append(del_pos)

                            # Lower deletion position scores are better (see below for the purpose
                            # of the score). Distinguish 3' and 5' deletions that are the same
                            # distance from the substitution locus by adding 0.5 to the magnitude of
                            # the 3' distance.
                            if del_pos_relative_to_sub <= 0:
                                locus_del_pos_scores.append(abs(del_pos_relative_to_sub))
                            else:
                                locus_del_pos_scores.append(abs(del_pos_relative_to_sub) + 0.5)
                        all_del_positions.extend(locus_del_positions)
                        all_del_pos_scores.extend(locus_del_pos_scores)
                    if all_del_positions:
                        # It is possible that in silico deletions originating from one substitution
                        # locus may coincide with those originating from another, so these
                        # redundancies must be resolved.
                        del_pos_config = []
                        prev_del_pos = -1
                        for del_pos_info in sorted(zip(all_del_positions, all_del_pos_scores), key=lambda del_pos_info: (del_pos_info[0], del_pos_info[1])):
                            if del_pos_info[0] > prev_del_pos:
                                del_pos_config.append(del_pos_info)
                        del_pos_configs.append(tuple(del_pos_config))
        # It is possible to generate the same sequence with deletions given different deletion
        # configurations. For example, ACCG can become ACG by deleting either C. If the second C is
        # the sole substitution locus in the sequence, then only record the in silico deletion as
        # occurring at that nucleotide rather one nucleotide 5' of it. This is resolved by choosing
        # the lower deletion position score. The deletion at the first C has a score of 1, whereas
        # the deletion at the second has a score of 0.
        del_dict = {}
        for del_pos_config in del_pos_configs:
            seq_string_with_dels = seq_string
            del_positions = []
            total_del_pos_score = 0
            for del_pos, del_pos_score in del_pos_config[::-1]:
                seq_string_with_dels = seq_string_with_dels[: del_pos] + seq_string_with_dels[del_pos + 1: ]
                del_positions.append(del_pos)
                total_del_pos_score += del_pos_score

            if seq_string_with_dels in del_dict:
                if total_del_pos_score >= del_dict[seq_string_with_dels][1]:
                    continue
            del_dict[seq_string_with_dels] = (tuple(sorted(del_positions)), total_del_pos_score)
        return del_dict


    def prefix_dereplicate_deletion_candidates(self,
                                               norm_seq_names,
                                               norm_seq_reversed_seq_strings,
                                               norm_seq_extras,
                                               mod_seq_child_names,
                                               mod_seq_child_reversed_seq_strings,
                                               mod_seq_child_extras):
        """3'-dereplicate modified sequence "children" with speculative in silico deletions (M') and
        normalized sequences (N). Clusters seeded by M' are produced when an N is found that is a 3'
        subsequence of M', potentially of equal length. Other M' that are 3' subsequences of the
        seed M' may also be part of these clusters. Clusters seeded by N are produced. Only M' are
        3' subsequences in these clusters, as N have already been dereplicated. No cluster, seeded
        either by M' or N, can contain more than one N."""
        pid = "Finding seqs with mod-induced dels"
        hashed_norm_seq_strings = [sha1(seq_string.encode('utf-8')).hexdigest() for seq_string in norm_seq_reversed_seq_strings]
        hashed_mod_seq_child_strings = [sha1(seq_string.encode('utf-8')).hexdigest() for seq_string in mod_seq_child_reversed_seq_strings]

        # FIND N and other M' in M'.
        ############################
        self.progress.update_pid(pid)
        self.progress.update("Finding seqs without in silico dels in seqs with them")
        # Make an `AlignedTarget` object for each M' target that is hit, and add info on matching N
        # or M' queries to its alignment attribute. Store these objects in the following dict.
        mod_seq_child_aligned_target_dict = {}
        # If an N is a 3' subsequence of an M', then do not search for it later as a target to form
        # a cluster seeded by N.
        matched_norm_seq_names = set()
        # Generate k-mers of the size of each N or M' query from the 3' end of each M' target.
        kmer_sizes = sorted(set([len(seq_string) for seq_string in mod_seq_child_reversed_seq_strings + norm_seq_reversed_seq_strings]))
        chunk_start = 0
        target_chunk_size = self.alignment_target_chunk_size
        chunk_stop = target_chunk_size
        while chunk_start < len(mod_seq_child_names):

            # GENERATE A DICT mapping M' target 3' prefix k-mers to target info.
            ####################################################################
            # This is like a `sequence.Kmerizer` method that uses multiple k-mer sizes and doesn't
            # generate `AlignedTarget` objects for each target. To speed up lookup, create a nested
            # dict, with the outer dict keyed by k-mer length.
            kmer_dict = {kmer_size: {} for kmer_size in kmer_sizes}
            for name, seq_string, target_extra in zip(mod_seq_child_names[chunk_start: chunk_stop],
                                                      mod_seq_child_reversed_seq_strings[chunk_start: chunk_stop],
                                                      mod_seq_child_extras[chunk_start: chunk_stop]):
                seq_length = len(seq_string)
                # The target extra here is a tuple with the first element being the in silico
                # deletion configuration of M' and the second being the modified sequence object, M.
                target_info = (name, seq_string, target_extra)
                for kmer_size in kmer_sizes:
                    if kmer_size > seq_length:
                        break

                    hashed_kmer = sha1(seq_string[: kmer_size].encode('utf-8')).hexdigest()
                    try:
                        kmer_dict[kmer_size][hashed_kmer].append(target_info)
                    except KeyError:
                        kmer_dict[kmer_size][hashed_kmer] = [target_info]


            # SEARCH for N queries in M' targets.
            #####################################
            for hashed_norm_seq_string, norm_seq_name, norm_seq_reversed_seq_string, norm_seq_extra in zip(hashed_norm_seq_strings,
                                                                                                           norm_seq_names,
                                                                                                           norm_seq_reversed_seq_strings,
                                                                                                           norm_seq_extras):
                try:
                    # N may be found in multiple M'.
                    target_info_list = kmer_dict[len(norm_seq_reversed_seq_string)][hashed_norm_seq_string]
                except KeyError:
                    # N is not a 3' subsequence of any M'.
                    continue

                for target_name, target_seq_string, target_extra in target_info_list:
                    try:
                        # The M' target has already been hit by another sequence -- presumably an M'
                        # query, since only one N can hit any target, as no N can be a 3' subsequence of
                        # another N.
                        target = mod_seq_child_aligned_target_dict[target_name]
                    except KeyError:
                        target = AlignedTarget(target_seq_string, name=target_name)
                        mod_seq_child_aligned_target_dict[target_name] = target
                        # Record information on the target as the target object's first alignment item.
                        target.alignments.append((target_name, target_seq_string, target_extra))
                    target.alignments.append((norm_seq_name, norm_seq_reversed_seq_string, norm_seq_extra))
                    matched_norm_seq_names.add(norm_seq_name)


            # SEARCH for M' queries in M' targets.
            ######################################
            for hashed_mod_seq_child_string, mod_seq_child_name, mod_seq_child_reversed_seq_string, mod_seq_child_extra in zip(hashed_mod_seq_child_strings,
                                                                                                                               mod_seq_child_names,
                                                                                                                               mod_seq_child_reversed_seq_strings,
                                                                                                                               mod_seq_child_extras):
                try:
                    target_info_list = kmer_dict[len(mod_seq_child_reversed_seq_string)][hashed_mod_seq_child_string]
                except KeyError:
                    continue

                self_match = False
                for target_name, target_seq_string, target_extra in target_info_list:
                    if not self_match:
                        # Use the Boolean, `self_match`, to limit string comparisons.
                        if mod_seq_child_name == target_name:
                            # Ignore self matches
                            self_match = True
                            continue
                    try:
                        target = mod_seq_child_aligned_target_dict[target_name]
                    except KeyError:
                        target = AlignedTarget(target_seq_string, name=target_name)
                        mod_seq_child_aligned_target_dict[target_name] = target
                        target.alignments.append((target_name, target_seq_string, target_extra))
                    target.alignments.append((mod_seq_child_name, mod_seq_child_reversed_seq_string, mod_seq_child_extra))

            # Delete the k-mer dict for the current chunk.
            del(kmer_dict)
            gc.collect()
            chunk_start = chunk_stop
            chunk_stop += target_chunk_size

        # Find N that DID NOT MATCH any M' to use as targets in the next search.
        matched_norm_seq_names = list(matched_norm_seq_names)
        unmatched_norm_seq_names = []
        unmatched_norm_seq_reversed_seq_strings = []
        unmatched_norm_seq_extras = []
        for name, seq_string, extra in zip(norm_seq_names, norm_seq_reversed_seq_strings, norm_seq_extras):
            if name not in matched_norm_seq_names:
                unmatched_norm_seq_names.append(name)
                unmatched_norm_seq_reversed_seq_strings.append(seq_string)
                unmatched_norm_seq_extras.append(extra)

        del(matched_norm_seq_names)
        # Delete lists for all N rather than just unmatched N.
        del(hashed_norm_seq_strings)
        del(norm_seq_names)
        del(norm_seq_reversed_seq_strings)
        del(norm_seq_extras)
        gc.collect()

        # FIND M' in unmatched N.
        #########################
        self.progress.update_pid(pid)
        self.progress.update("Find seqs with in silico dels in seqs without them")
        # Make an `AlignedTarget` object for each N target that is hit, and add info on matching M'
        # queries to its alignment attribute. Store these objects in the following dict.
        norm_seq_aligned_target_dict = {}
        # Generate k-mers of the size of each M' query from the 3' end of each N target.
        kmer_sizes = sorted(set([len(seq_string) for seq_string in mod_seq_child_reversed_seq_strings]))
        chunk_start = 0
        chunk_stop = target_chunk_size
        while chunk_start < len(unmatched_norm_seq_names):
            chunk_norm_seq_lengths = sorted(set([len(seq_string) for seq_string in unmatched_norm_seq_reversed_seq_strings[chunk_start: chunk_stop]]))

            # GENERATE A DICT mapping N target 3' prefix k-mers to target info.
            ###################################################################
            kmer_dict = {kmer_size: {} for kmer_size in kmer_sizes}
            for name, seq_string, target_extra in zip(unmatched_norm_seq_names[chunk_start: chunk_stop],
                                                      unmatched_norm_seq_reversed_seq_strings[chunk_start: chunk_stop],
                                                      unmatched_norm_seq_extras[chunk_start: chunk_stop]):
                seq_length = len(seq_string)
                target_info = (name, seq_string, target_extra)
                for kmer_size in kmer_sizes:
                    if kmer_size > seq_length:
                        break

                    hashed_kmer = sha1(seq_string[: kmer_size].encode('utf-8')).hexdigest()
                    try:
                        kmer_dict[kmer_size][hashed_kmer].append(target_info)
                    except KeyError:
                        kmer_dict[kmer_size][hashed_kmer] = [target_info]


            # SEARCH for M' queries in N targets.
            #####################################
            for hashed_mod_seq_child_string, mod_seq_child_name, mod_seq_child_reversed_seq_string, mod_seq_child_extra in zip(hashed_mod_seq_child_strings,
                                                                                                                               mod_seq_child_names,
                                                                                                                               mod_seq_child_reversed_seq_strings,
                                                                                                                               mod_seq_child_extras):
                try:
                    target_info_list = kmer_dict[len(mod_seq_child_reversed_seq_string)][hashed_mod_seq_child_string]
                except KeyError:
                    continue

                for target_name, target_seq_string, target_extra in target_info_list:
                    try:
                        target = norm_seq_aligned_target_dict[target_name]
                    except KeyError:
                        target = AlignedTarget(target_seq_string, name=target_name)
                        norm_seq_aligned_target_dict[target_name] = target
                        target.alignments.append((target_name, target_seq_string, target_extra))
                    target.alignments.append((mod_seq_child_name, mod_seq_child_reversed_seq_string, mod_seq_child_extra))

            del(kmer_dict)
            gc.collect()
            chunk_start = chunk_stop
            chunk_stop += target_chunk_size

        del(hashed_mod_seq_child_strings)
        gc.collect()

        clusters = []
        # CREATE CLUSTER objects for each M' containing an N.
        #####################################################
        self.progress.update_pid(pid)
        self.progress.update("Creating clusters seeded by seqs with in silico dels")
        for mod_seq_child_target in mod_seq_child_aligned_target_dict.values():
            norm_seq_added = False
            if mod_seq_child_target.alignments:
                # Start filling out a cluster object, but only add meaningful clusters containing an
                # N to the final list.
                cluster = Cluster()
                # The first item in the alignment is info on the M' seed. Do not add info on the
                # seed to the cluster until it is established that the cluster contains an N. For
                # full reproducibility, sort matching queries in descending order of length, and in
                # case of ties, by name.
                for name, reversed_seq_string, extra in sorted(mod_seq_child_target.alignments[1: ],
                                                               key=lambda alignment_with_mod_seq_child: (-len(alignment_with_mod_seq_child[1]), alignment_with_mod_seq_child[0])):
                    cluster.member_names.append(name)
                    cluster.member_seqs.append(reversed_seq_string)
                    cluster.member_extras.append(extra)
                    if not norm_seq_added:
                        if not isinstance(extra, tuple):
                            norm_seq_added = True
                if norm_seq_added:
                    name, reversed_seq_string, extra = mod_seq_child_target.alignments[0]
                    cluster.member_names.insert(0, name)
                    cluster.member_seqs.insert(0, reversed_seq_string)
                    cluster.member_extras.insert(0, extra)
                    clusters.append(cluster)

        # CREATE CLUSTER objects for each N containing but not contained by an M'.
        ##########################################################################
        self.progress.update_pid(pid)
        self.progress.update("Creating clusters seeded by seqs without in silico dels")
        for norm_seq_target in norm_seq_aligned_target_dict.values():
            if norm_seq_target.alignments:
                cluster = Cluster()
                # The first item is info on the N seed.
                name, reversed_seq_string, extra = norm_seq_target.alignments[0]
                cluster.member_names.append(name)
                cluster.member_seqs.append(reversed_seq_string)
                cluster.member_extras.append(extra)
                # Add info on matching M' queries.
                for name, reversed_seq_string, extra in sorted(norm_seq_target.alignments[1: ],
                                                               key=lambda mod_seq_child_alignment_with_norm_seq: (-len(mod_seq_child_alignment_with_norm_seq[1]), mod_seq_child_alignment_with_norm_seq[0])):
                    cluster.member_names.append(name)
                    cluster.member_seqs.append(reversed_seq_string)
                    cluster.member_extras.append(extra)
                clusters.append(cluster)

        # For full reproducibility, sort clusters in descending order of size, and in case of ties,
        # by name.
        clusters.sort(key=lambda cluster: (-len(cluster.member_names), cluster.member_names[0]))
        return clusters


    def process_deletion_clusters(self, clusters, norm_seq_type):
        """Process 3'-dereplicated clusters comprised of modified sequences with in silico
        deletions, *M'*, and normalized sequences. The normalized sequences either all have a
        truncated feature profile, *Nt*, or all have a full feature profile, *Nf*.

        We are interested in clusters with both M' and a normalized sequence, verifying the
        speculative deletions in M' by the existence of a corresponding normalized sequence. There
        can only be one N per cluster, as normalized sequences have previously been 3'-dereplicated.
        There can be multiple M'.

        N with supported deletions are removed and reconstituted as `NormalizedDeletionSequence`
        objects.

        Normalized sequences can often be sequences with deletions that also have extra 5'
        nucleotides. When this type of normalized sequence is present in the cluster, M' will be
        shorter. In the case of Nf, the extra nucleotides went undetected in the initial assignment
        of the feature profile, as they are within a normalized sequence, which, by definition, is
        meant to have trimmed 3' and 5' ends. For the extra 5' nucleotides to be detected, it must
        be confirmed that the longest modified sequence is a full-length tRNA.

        N can also arise from tRNA fragments. In a cluster, these would be 3' subsequences of one or
        more M'.
        """
        norm_seq_mod_seqs_dict = self.get_normalized_sequences_containing_modified_sequences_with_deletions(clusters)

        # Process the matches between N and one or more M'.
        if norm_seq_type == 'trna':
            norm_seq_dict = self.norm_trna_seq_dict
        elif norm_seq_type == 'trunc':
            norm_seq_dict = self.norm_trunc_seq_dict
        new_norm_del_seq_dict = {}
        # It is useful to know which of a trimmed sequence's normalized sequences are found to
        # contain deletions.
        trimmed_seq_norm_del_seq_dict = defaultdict(list)
        # Winnow the matches down to one-to-one matches between a normalized sequence and modified
        # sequence, ignoring N that can be formed from the introduction of deletions in different
        # modified sequences, as indicated by the following variable, which is set to `False` and
        # cannot currently be changed by the user.
        allow_norm_seq_with_dels_from_multiple_mod_seqs = self.allow_norm_seq_with_dels_from_multiple_mod_seqs
        for norm_seq_name, match_info in norm_seq_mod_seqs_dict.items():
            if len(match_info) == 1:
                norm_seq, mod_seq, del_config, extra_fiveprime_length = match_info[0]
            else:
                # The normalized sequence was found in multiple M', which may be from the same or
                # different modified sequences.
                uniq_mod_seq_info_dict = {}
                # In the following loop, `norm_seq` is the same in every iteration. This same
                # variable is referenced after the loop.
                for norm_seq, mod_seq, del_config, extra_fiveprime_length in match_info:
                    if mod_seq.represent_name in uniq_mod_seq_info_dict:
                        # Multiple deletion configurations in the same modified sequence are
                        # apparently able to produce the normalized sequence.
                        if len(del_config) < len(uniq_mod_seq_info_dict[mod_seq.represent_name]):
                            # Favor the most parsimonious configuration of deletions producing the
                            # normalized sequence.
                            uniq_mod_seq_info_dict[mod_seq.represent_name] = del_config
                    else:
                        uniq_mod_seq_info_dict[mod_seq.represent_name] = del_config

                if len(uniq_mod_seq_info_dict) > 1:
                    if not allow_norm_seq_with_dels_from_multiple_mod_seqs:
                        # The normalized sequence with deletions can arise from multiple modified
                        # sequences, so ignore it.
                        continue

                mod_seq = match_info[0][1]
                del_config = uniq_mod_seq_info_dict[mod_seq.represent_name]
                extra_fiveprime_length = match_info[0][3]
            for trimmed_seq in norm_seq.trimmed_seqs:
                trimmed_seq_norm_del_seq_dict[trimmed_seq.represent_name].append(norm_seq.represent_name)

            # Transfer the contents of a normalized sequence with supported deletions to a
            # `NormalizedDeletionSequence` object. Some if not all of the feature profiles of the
            # trimmed/unique sequences underlying the normalized sequence are invalidated by the
            # deletions. However, do not alter the trimmed and unique sequences. A certain amount of
            # information contradicting the trimmed and unique sequence profiles is stored in the
            # `NormalizedDeletionSequence`, such as any additional 5' extension contained in the
            # normalized sequence.
            if not extra_fiveprime_length:
                # N was the same length as M'. No other normalized sequence will contain M', as
                # normalized sequences were 3' dereplicated earlier in the workflow.
                norm_seq_dict.pop(norm_seq.represent_name)
                norm_del_seq_string = norm_seq.seq_string
                # Offload the work of finding the position of N in M and whether N contains the
                # anticodon to `TRNASeqDataset` rather than having methods for this in each instance
                # of `NormalizedDeletionSequence`.
                norm_del_seq_start_in_mod_seq, norm_seq_del_config = self.find_normalized_deletion_sequence_in_modified_sequence(len(norm_del_seq_string), mod_seq, del_config)
                norm_del_seq_contains_anticodon = self.check_normalized_deletion_sequence_for_anticodon(mod_seq, norm_del_seq_start_in_mod_seq)
                norm_del_seq = NormalizedDeletionSequence(norm_del_seq_string, norm_seq, mod_seq, del_config, norm_del_seq_start_in_mod_seq, norm_seq_del_config, norm_del_seq_contains_anticodon)
                new_norm_del_seq_dict[norm_del_seq_string] = norm_del_seq
                mod_seq.norm_seqs_with_dels.append(norm_del_seq)
                mod_seq.del_configs.append(del_config)

                continue

            # By reaching this point, N was found to be longer than M'. Multiple N can contain a
            # given M'. N is longer than M' when it contains previously unidentified extra 5'
            # nucleotides. The new 5' extension is recorded in the `NormalizedDeletionSequence`
            # object. When multiple N are the same as M' except for the new 5' extension, they
            # are consolidated into the same object.
            norm_del_seq_string = norm_seq.seq_string[extra_fiveprime_length: ]
            if norm_del_seq_string in new_norm_del_seq_dict:
                norm_del_seq = new_norm_del_seq_dict[norm_del_seq_string]
                # Avoid adding duplicate trimmed sequences to the object (those shorter than N).
                norm_del_seq_length = len(norm_del_seq_string)
                trimmed_seqs = [trimmed_seq for trimmed_seq in norm_seq.trimmed_seqs if len(trimmed_seq.seq_string) > norm_del_seq_length]
                norm_del_seq.trimmed_seqs.extend(trimmed_seqs)
                norm_del_seq.defunct_norm_seqs.extend([norm_seq for _ in trimmed_seqs])
            else:
                # Prevent N from actually being a slightly shorter, deletion-free 3' subsequence
                # fragment of the modified sequence, M.
                unsupported_dels = False
                for norm_seq_without_dels in mod_seq.norm_seqs_without_dels:
                    if norm_del_seq_string == norm_seq_without_dels.seq_string[len(del_config): ]:
                        unsupported_dels = True
                        break
                if unsupported_dels:
                    continue

                # The following code is the same as above, when dealing with normalized sequences
                # without a newly discovered 5' extension; perhaps it should be a separate method.
                norm_seq_dict.pop(norm_seq.represent_name)
                norm_del_seq_start_in_mod_seq, norm_seq_del_config = self.find_normalized_deletion_sequence_in_modified_sequence(len(norm_del_seq_string), mod_seq, del_config)
                norm_del_seq_contains_anticodon = self.check_normalized_deletion_sequence_for_anticodon(mod_seq, norm_del_seq_start_in_mod_seq)
                norm_del_seq = NormalizedDeletionSequence(norm_del_seq_string, norm_seq, mod_seq, del_config, norm_del_seq_start_in_mod_seq, norm_seq_del_config, norm_del_seq_contains_anticodon)
                new_norm_del_seq_dict[norm_del_seq_string] = norm_del_seq
                mod_seq.norm_seqs_with_dels.append(norm_del_seq)
                mod_seq.del_configs.append(del_config)

        norm_del_seq_dict = self.norm_del_seq_dict
        trimmed_del_seq_dict = self.trimmed_del_seq_dict
        uniq_del_seq_dict = self.uniq_del_seq_dict
        for norm_seq in new_norm_del_seq_dict.values():
            norm_seq.init()
            norm_del_seq_dict[norm_seq.represent_name] = norm_seq
            # Record the constituent sequences of normalized sequences with deletions.
            for trimmed_seq in norm_seq.trimmed_seqs:
                trimmed_del_seq_dict[trimmed_seq.represent_name] = (trimmed_seq, set(trimmed_seq_norm_del_seq_dict[trimmed_seq.represent_name]))
                for uniq_seq in trimmed_seq.uniq_seqs:
                    uniq_del_seq_dict[uniq_seq.represent_name] = uniq_seq
        if norm_seq_type == 'trunc':
            # If a sequence with a truncated feature profile is identified as tRNA, albeit with a
            # deletion, then record it as such. Trimmed truncated sequences may be nonspecific,
            # found in other normalized truncated sequences not recovered by this method. This
            # nuance should not affect anything downstream.
            trimmed_trna_seq_dict = self.trimmed_trna_seq_dict
            trimmed_trunc_seq_dict = self.trimmed_trunc_seq_dict
            uniq_trna_seq_dict = self.uniq_trna_seq_dict
            uniq_trunc_seq_dict = self.uniq_trunc_seq_dict
            for norm_seq in new_norm_del_seq_dict.values():
                for trimmed_seq in norm_seq.trimmed_seqs:
                    if isinstance(trimmed_seq, TrimmedTruncatedProfileSequence):
                        try:
                            trimmed_trunc_seq_dict.pop(trimmed_seq.represent_name)
                        except KeyError:
                            # The trimmed sequence was part of another normalized sequence with
                            # deletions that was already processed.
                            continue
                        trimmed_trna_seq_dict[trimmed_seq.represent_name] = trimmed_seq
                        for uniq_seq in trimmed_seq.uniq_seqs:
                            uniq_trunc_seq_dict.pop(uniq_seq.represent_name)
                            uniq_trna_seq_dict[uniq_seq.represent_name] = uniq_seq


    def get_normalized_sequences_containing_modified_sequences_with_deletions(self, clusters):
        """The first step of `process_deletion_clusters`. Clusters are inspected to find normalized
        sequences that are the same as or contain as 3' subsequences modified sequences with in
        silico deletions."""
        norm_seq_mod_seqs_dict = defaultdict(list)
        min_fiveprime_del_pos = self.min_length_of_long_fiveprime_extension - 1
        for cluster in clusters:
            norm_seq_name = None
            # Track the previous (the only) normalized sequence found in the cluster.
            norm_seq = None
            norm_seq_length = None
            # Track the previous modified sequences found in the cluster before a normalized
            # sequence is found.
            mod_seqs = []
            del_configs = []
            for seq_name, member_extra in zip(cluster.member_names, cluster.member_extras):
                if isinstance(member_extra, tuple):
                    del_config = member_extra[0]
                    mod_seq = member_extra[1]
                    del_configs.append(del_config)
                    mod_seqs.append(mod_seq)
                    if norm_seq:
                        # The normalized sequence came before the modified sequence in the cluster.
                        extra_fiveprime_length = norm_seq_length + len(del_config) - len(mod_seq.norm_seqs_without_dels[0].seq_string)
                        if extra_fiveprime_length == 0:
                            # The normalized sequence is equal in length to the longest M' in the
                            # cluster, and happened to come before it in the cluster. It seems that
                            # this does not actually happen when this method is called from
                            # `find_deletions` due to the order in which sequences are processed
                            # during clustering, but entertain the possibility for the sake of
                            # completeness were the clusters to be formed slightly differently.
                            for mod_seq, del_config in zip(mod_seqs, del_configs):
                                norm_seq_mod_seqs_dict[norm_seq_name].append((norm_seq, mod_seq, del_config, 0))
                            mod_seqs = []
                            del_configs = []
                        elif extra_fiveprime_length > 0:
                            # The normalized sequence is longer than the longest M' in the cluster.
                            if mod_seq.norm_seqs_without_dels[0].has_complete_feature_set:
                                # The longest modified sequence in the cluster is a full-length
                                # tRNA. Therefore, the normalized sequence contains extra 5'
                                # nucleotides that should be trimmed.
                                for mod_seq, del_config in zip(mod_seqs, del_configs):
                                    norm_seq_mod_seqs_dict[norm_seq_name].append((norm_seq, mod_seq, del_config, extra_fiveprime_length))
                                mod_seqs = []
                                del_configs = []
                else: # considering normalized sequence
                    norm_seq_name = seq_name
                    norm_seq = member_extra
                    norm_seq_length = len(norm_seq.seq_string)
                    if not mod_seqs:
                        continue
                    for mod_seq, del_config in zip(mod_seqs, del_configs):
                        if del_config[0] - len(mod_seq.norm_seqs_without_dels[0].seq_string) + norm_seq_length + len(del_config) < min_fiveprime_del_pos:
                            # There must be sufficient matching sequence length on the 5' end of the
                            # normalized sequence to confirm the 5'-most deletion.
                            continue
                        norm_seq_mod_seqs_dict[norm_seq_name].append((norm_seq, mod_seq, del_config, 0))
                    mod_seqs = []
                    del_configs = []
        return norm_seq_mod_seqs_dict


    def find_normalized_deletion_sequence_in_modified_sequence(self, norm_del_seq_length, mod_seq, del_config):
        """Normalized sequences with deletions and modified sequences are aligned at the 3' end, so
        find the start position of the normalized sequence in the modified sequence by working
        backward from the 3' end. Also record the normalized sequence positions immediately 5' of
        the deletions."""
        # Example:
        # ACGTAAC (mod seq)
        #    T  C (norm seq)
        # del_pos = 4
        # norm_seq_pos == 1
        # mod_seq_pos == 6
        # mod_seq_pos -> 5
        # norm_seq_pos -> 0
        # mod_seq_pos -> 4 (continue)
        # mod_seq_pos -> 3 (del_pos = -1)
        # mod_seq_pos -> 2
        # norm_seq_pos -> -1 (break)
        # norm_seq_start = mod_seq_pos + 1 == 3

        del_config_iterator = iter(del_config[::-1])
        del_pos = next(del_config_iterator)
        norm_seq_pos = norm_del_seq_length - 1
        mod_seq_pos = len(mod_seq.norm_seqs_without_dels[0].seq_string) - 1
        norm_seq_del_config = []
        while norm_seq_pos > -1:
            if mod_seq_pos == del_pos:
                mod_seq_pos -= 1
                norm_seq_del_config.append(norm_seq_pos)
                try:
                    del_pos = next(del_config_iterator)
                    continue
                except StopIteration:
                    del_pos = -1
            mod_seq_pos -= 1
            norm_seq_pos -= 1
        norm_seq_start_in_mod_seq = mod_seq_pos + 1
        return norm_seq_start_in_mod_seq, norm_seq_del_config


    def check_normalized_deletion_sequence_for_anticodon(self, mod_seq, norm_del_seq_start_in_mod_seq):
        # Determine whether the normalized sequence with deletions contains the anticodon.
        try:
            anticodon_loop_start = mod_seq.norm_seqs_without_dels[0].trimmed_seqs[0].feature_start_indices[self.RELATIVE_ANTICODON_LOOP_INDEX]
        except IndexError:
            # The anticodon loop was not reached in the profile.
            return False
        if anticodon_loop_start >= norm_del_seq_start_in_mod_seq:
            return True
        else:
            return False


    def report_statistics(self):
        """Add run statistics to the database and write them to the summary file."""
        profiled_trna_reads = 0
        trna_reads_with_threeprime_cca = 0
        trna_reads_with_threeprime_cc = 0
        trna_reads_with_threeprime_c = 0
        trna_reads_with_other_threeprime_termini = 0
        profiled_trna_reads_containing_anticodon = 0
        full_length_trna_reads = 0
        trna_reads_with_extrapolated_fiveprime_feature = 0
        min_length_of_long_fiveprime_extension = self.min_length_of_long_fiveprime_extension
        profiled_trna_reads_with_short_fiveprime_extension = 0
        profiled_trna_reads_with_long_fiveprime_extension = 0
        trna_reads_with_recovered_trunc_profile = 0
        uniq_profiled_trna_seqs = 0
        uniq_trna_seqs_with_threeprime_cca = 0
        uniq_trna_seqs_with_threeprime_cc = 0
        uniq_trna_seqs_with_threeprime_c = 0
        uniq_trna_seqs_with_other_threeprime_termini = 0
        uniq_trna_seqs_containing_anticodon = 0
        full_length_uniq_seqs = 0
        uniq_seqs_with_extrapolated_fiveprime_feature = 0
        uniq_seqs_with_short_fiveprime_extension = 0
        uniq_seqs_with_long_fiveprime_extension = 0
        uniq_seqs_with_recovered_trunc_profile = 0
        interior_mapped_reads = 0
        fiveprime_mapped_reads = 0
        interior_mapped_uniq_seqs = 0
        fiveprime_mapped_uniq_seqs = 0
        for uniq_trna_seq in self.uniq_trna_seq_dict.values():
            read_count = uniq_trna_seq.read_count

            if isinstance(uniq_trna_seq, UniqueMappedSequence):
                if uniq_trna_seq.extra_fiveprime_length:
                    fiveprime_mapped_reads += read_count
                    fiveprime_mapped_uniq_seqs += 1
                else:
                    interior_mapped_reads += read_count
                    interior_mapped_uniq_seqs += 1
                continue

            # To reach this point, the unique sequence has a full or truncated feature profile.
            profiled_trna_reads += read_count
            uniq_profiled_trna_seqs += 1

            has_full_profile = True
            if isinstance(uniq_trna_seq, UniqueTruncatedProfileSequence):
                has_full_profile = False

            if uniq_trna_seq.threeprime_terminus_length:
                threeprime_terminus = uniq_trna_seq.seq_string[-uniq_trna_seq.threeprime_terminus_length: ]
                if threeprime_terminus == 'CCA':
                    trna_reads_with_threeprime_cca += read_count
                    uniq_trna_seqs_with_threeprime_cca += 1
                elif threeprime_terminus == 'CC':
                    trna_reads_with_threeprime_cc += read_count
                    uniq_trna_seqs_with_threeprime_cc += 1
                elif threeprime_terminus == 'C':
                    trna_reads_with_threeprime_c += read_count
                    uniq_trna_seqs_with_threeprime_c += 1
                else:
                    trna_reads_with_other_threeprime_termini += read_count
                    uniq_trna_seqs_with_other_threeprime_termini += 1
            else:
                trna_reads_with_other_threeprime_termini += read_count
                uniq_trna_seqs_with_other_threeprime_termini += 1

            if uniq_trna_seq.contains_anticodon:
                profiled_trna_reads_containing_anticodon += read_count
                uniq_trna_seqs_containing_anticodon += 1

            if has_full_profile:
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
            else:
                trna_reads_with_recovered_trunc_profile += read_count
                uniq_seqs_with_recovered_trunc_profile += 1

        total_trna_reads = profiled_trna_reads + interior_mapped_reads + fiveprime_mapped_reads
        total_uniq_trna_seqs = uniq_profiled_trna_seqs + interior_mapped_uniq_seqs + fiveprime_mapped_uniq_seqs
        mean_profiled_reads_per_uniq_seq = profiled_trna_reads / uniq_profiled_trna_seqs
        mean_mapped_reads_per_uniq_seq = (fiveprime_mapped_reads + interior_mapped_reads) / (fiveprime_mapped_uniq_seqs + interior_mapped_uniq_seqs)


        trimmed_del_seq_dict = self.trimmed_del_seq_dict
        trna_reads_specific_to_norm_seqs_with_del_signature = 0
        trna_reads_nonspecific_to_norm_seqs_with_del_signature = 0
        uniq_seqs_specific_to_norm_seqs_with_del_signature = 0
        uniq_seqs_nonspecific_to_norm_seqs_with_del_signature = 0
        # Count unique sequences that are part of normalized sequences with deletions. For unique
        # sequences that are only part of normalized sequences with deletions, subtract their counts
        # from the previous counts.
        for uniq_trna_seq in self.uniq_del_seq_dict.values():
            read_count = uniq_trna_seq.read_count

            trimmed_seq, defunct_norm_seq_names = trimmed_del_seq_dict[uniq_trna_seq.trimmed_seq_represent_name]
            if len(trimmed_seq.norm_seq_represent_names) == len(defunct_norm_seq_names):
                specific_to_norm_seqs_with_del_signature = True
            else:
                specific_to_norm_seqs_with_del_signature = False

            if specific_to_norm_seqs_with_del_signature:
                trna_reads_specific_to_norm_seqs_with_del_signature += read_count
                uniq_seqs_specific_to_norm_seqs_with_del_signature += 1
            else:
                trna_reads_nonspecific_to_norm_seqs_with_del_signature += read_count
                uniq_seqs_nonspecific_to_norm_seqs_with_del_signature += 1
                continue

            # To reach this point, the unique sequence must be specific to a normalized sequence
            # with deletions.
            if isinstance(uniq_trna_seq, UniqueMappedSequence):
                if uniq_trna_seq.extra_fiveprime_length:
                    fiveprime_mapped_reads -= read_count
                    fiveprime_mapped_uniq_seqs -= 1
                else:
                    interior_mapped_reads -= read_count
                    interior_mapped_uniq_seqs -= 1
                continue

            # To reach this point, the unique sequence has a full or truncated profile.
            profiled_trna_reads -= read_count
            uniq_profiled_trna_seqs -= 1

            has_full_profile = True
            if isinstance(uniq_trna_seq, UniqueTruncatedProfileSequence):
                has_full_profile = False

            if uniq_trna_seq.threeprime_terminus_length:
                threeprime_terminus = uniq_trna_seq.seq_string[-uniq_trna_seq.threeprime_terminus_length: ]
                if threeprime_terminus == 'CCA':
                    trna_reads_with_threeprime_cca -= read_count
                    uniq_trna_seqs_with_threeprime_cca -= 1
                elif threeprime_terminus == 'CC':
                    trna_reads_with_threeprime_cc -= read_count
                    uniq_trna_seqs_with_threeprime_cc -= 1
                elif threeprime_terminus == 'C':
                    trna_reads_with_threeprime_c -= read_count
                    uniq_trna_seqs_with_threeprime_c -= 1
                else:
                    trna_reads_with_other_threeprime_termini -= read_count
                    uniq_trna_seqs_with_other_threeprime_termini -= 1
            else:
                trna_reads_with_other_threeprime_termini -= read_count
                uniq_trna_seqs_with_other_threeprime_termini -= 1

            if uniq_trna_seq.contains_anticodon:
                profiled_trna_reads_containing_anticodon -= read_count
                uniq_trna_seqs_containing_anticodon -= 1

            if has_full_profile:
                if uniq_trna_seq.has_complete_feature_set:
                    full_length_trna_reads -= read_count
                    full_length_uniq_seqs -= 1

                if uniq_trna_seq.num_extrapolated_fiveprime_nts:
                    trna_reads_with_extrapolated_fiveprime_feature -= read_count
                    uniq_seqs_with_extrapolated_fiveprime_feature -= 1

                if uniq_trna_seq.extra_fiveprime_length:
                    if uniq_trna_seq.extra_fiveprime_length >= min_length_of_long_fiveprime_extension:
                        profiled_trna_reads_with_long_fiveprime_extension -= read_count
                        uniq_seqs_with_long_fiveprime_extension -= 1
                    else:
                        profiled_trna_reads_with_short_fiveprime_extension -= read_count
                        uniq_seqs_with_short_fiveprime_extension -= 1
            else:
                trna_reads_with_recovered_trunc_profile -= read_count
                uniq_seqs_with_recovered_trunc_profile -= 1


        # The TOTAL count of trimmed sequence objects combines apples and oranges (profiled and
        # mapped sequences) and so should not be reported. Separate trimmed sequences are formed for
        # each unique mapped sequence, because the 5' extension of a mapped sequence can be most of
        # its length, so consolidation of unique mapped sequences with different 5' extensions would
        # result in very short trimmed sequences of only a handful of nucleotides.
        trimmed_profiled_seqs_containing_anticodon = 0
        full_length_trimmed_seqs = 0
        trimmed_seqs_with_extrapolated_fiveprime_feature = 0
        trimmed_seqs_with_recovered_trunc_profile = 0
        for trimmed_trna_seq in self.trimmed_trna_seq_dict.values():
            if isinstance(trimmed_trna_seq, TrimmedMappedSequence):
                continue

            represent_uniq_seq = trimmed_trna_seq.uniq_seqs[0]

            if trimmed_trna_seq.contains_anticodon:
                trimmed_profiled_seqs_containing_anticodon += 1

            has_full_profile = True
            if isinstance(trimmed_trna_seq, TrimmedTruncatedProfileSequence):
                has_full_profile = False

            if has_full_profile:
                if trimmed_trna_seq.has_complete_feature_set:
                    full_length_trimmed_seqs += 1

                if represent_uniq_seq.num_extrapolated_fiveprime_nts:
                    trimmed_seqs_with_extrapolated_fiveprime_feature += 1
            else:
                trimmed_seqs_with_recovered_trunc_profile += 1


        # Count trimmed sequences that are part of normalized sequences with deletions. For trimmed
        # sequences that are only part of normalized sequences with deletions, subtract their counts
        # from the previous counts.
        trimmed_seqs_specific_to_norm_seqs_with_del_signature = 0
        trimmed_seqs_nonspecific_to_norm_seqs_with_del_signature = 0
        for trimmed_trna_seq, defunct_norm_seq_names in self.trimmed_del_seq_dict.values():
            if len(trimmed_trna_seq.norm_seq_represent_names) == len(defunct_norm_seq_names):
                specific_to_norm_seqs_with_del_signature = True
            else:
                specific_to_norm_seqs_with_del_signature = False

            if specific_to_norm_seqs_with_del_signature:
                trimmed_seqs_specific_to_norm_seqs_with_del_signature += 1
            else:
                trimmed_seqs_nonspecific_to_norm_seqs_with_del_signature += 1
                continue

            if isinstance(trimmed_trna_seq, TrimmedMappedSequence):
                continue

            # To reach this point, the trimmed sequence must be specific to a normalized sequence
            # with deletions.
            has_full_profile = True
            if isinstance(trimmed_trna_seq, TrimmedTruncatedProfileSequence):
                has_full_profile = False

            if has_full_profile:
                if trimmed_trna_seq.has_complete_feature_set:
                    full_length_trimmed_seqs -= 1

                if trimmed_trna_seq.uniq_seqs[0].num_extrapolated_fiveprime_nts:
                    trimmed_seqs_with_extrapolated_fiveprime_feature -= 1
            else:
                trimmed_seqs_with_recovered_trunc_profile -= 1


        norm_trna_seqs = len(self.norm_trna_seq_dict)
        norm_trna_seqs_without_mods = 0
        norm_trna_seqs_containing_anticodon = 0
        full_length_norm_seqs = 0
        for norm_trna_seq in self.norm_trna_seq_dict.values():
            if not norm_trna_seq.mod_seqs:
                norm_trna_seqs_without_mods += 1

            if norm_trna_seq.trimmed_seqs[0].contains_anticodon:
                norm_trna_seqs_containing_anticodon += 1

            if norm_trna_seq.has_complete_feature_set:
                full_length_norm_seqs += 1

        norm_seqs_with_del_signature = len(self.norm_del_seq_dict)

        mod_trna_seqs = len(self.mod_trna_seq_dict)
        mod_seqs_with_del_signature = 0
        for mod_trna_seq in self.mod_trna_seq_dict.values():
            if mod_trna_seq.del_configs:
                mod_seqs_with_del_signature += 1

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        set_meta_value = trnaseq_db.db.set_meta_value

        set_meta_value('total_trna_reads', total_trna_reads)
        set_meta_value('profiled_trna_reads', profiled_trna_reads)
        set_meta_value('trna_reads_with_threeprime_cca', trna_reads_with_threeprime_cca)
        set_meta_value('trna_reads_with_threeprime_cc', trna_reads_with_threeprime_cc)
        set_meta_value('trna_reads_with_threeprime_c', trna_reads_with_threeprime_c)
        set_meta_value('trna_reads_with_other_threeprime_termini', trna_reads_with_other_threeprime_termini)
        set_meta_value('profiled_trna_reads_containing_anticodon', profiled_trna_reads_containing_anticodon)
        set_meta_value('full_length_trna_reads', full_length_trna_reads)
        set_meta_value('trna_reads_with_extrapolated_fiveprime_feature', trna_reads_with_extrapolated_fiveprime_feature)
        set_meta_value('profiled_trna_reads_with_short_fiveprime_extension', profiled_trna_reads_with_short_fiveprime_extension)
        set_meta_value('profiled_trna_reads_with_long_fiveprime_extension', profiled_trna_reads_with_long_fiveprime_extension)
        set_meta_value('trna_reads_with_recovered_truncated_profile', trna_reads_with_recovered_trunc_profile)
        set_meta_value('interior_mapped_reads', interior_mapped_reads)
        set_meta_value('fiveprime_mapped_reads', fiveprime_mapped_reads)
        set_meta_value('trna_reads_specific_to_norm_seqs_with_del_signature', trna_reads_specific_to_norm_seqs_with_del_signature)
        set_meta_value('trna_reads_nonspecific_to_norm_seqs_with_del_signature', trna_reads_nonspecific_to_norm_seqs_with_del_signature)

        set_meta_value('total_unique_trna_seqs', total_uniq_trna_seqs)
        set_meta_value('unique_profiled_trna_seqs', uniq_profiled_trna_seqs)
        set_meta_value('unique_trna_seqs_with_threeprime_cca', uniq_trna_seqs_with_threeprime_cca)
        set_meta_value('unique_trna_seqs_with_threeprime_cc', uniq_trna_seqs_with_threeprime_cc)
        set_meta_value('unique_trna_seqs_with_threeprime_c', uniq_trna_seqs_with_threeprime_c)
        set_meta_value('unique_trna_seqs_with_other_threeprime_termini', uniq_trna_seqs_with_other_threeprime_termini)
        set_meta_value('unique_trna_seqs_containing_anticodon', uniq_trna_seqs_containing_anticodon)
        set_meta_value('full_length_unique_seqs', full_length_uniq_seqs)
        set_meta_value('unique_seqs_with_extrapolated_fiveprime_feature', uniq_seqs_with_extrapolated_fiveprime_feature)
        set_meta_value('unique_seqs_with_short_fiveprime_extension', uniq_seqs_with_short_fiveprime_extension)
        set_meta_value('unique_seqs_with_long_fiveprime_extension', uniq_seqs_with_long_fiveprime_extension)
        set_meta_value('unique_seqs_with_recovered_truncated_profile', uniq_seqs_with_recovered_trunc_profile)
        set_meta_value('interior_mapped_uniq_seqs', interior_mapped_uniq_seqs)
        set_meta_value('fiveprime_mapped_uniq_seqs', fiveprime_mapped_uniq_seqs)
        set_meta_value('uniq_seqs_specific_to_norm_seqs_with_del_signature', uniq_seqs_specific_to_norm_seqs_with_del_signature)
        set_meta_value('uniq_seqs_nonspecific_to_norm_seqs_with_del_signature', uniq_seqs_nonspecific_to_norm_seqs_with_del_signature)

        set_meta_value('trimmed_profiled_seqs_containing_anticodon', trimmed_profiled_seqs_containing_anticodon)
        set_meta_value('full_length_trimmed_seqs', full_length_trimmed_seqs)
        set_meta_value('trimmed_seqs_with_extrapolated_fiveprime_feature', trimmed_seqs_with_extrapolated_fiveprime_feature)
        set_meta_value('trimmed_seqs_with_recovered_trunc_profile', trimmed_seqs_with_recovered_trunc_profile)
        set_meta_value('trimmed_seqs_specific_to_norm_seqs_with_del_signature', trimmed_seqs_specific_to_norm_seqs_with_del_signature)
        set_meta_value('trimmed_seqs_nonspecific_to_norm_seqs_with_del_signature', trimmed_seqs_nonspecific_to_norm_seqs_with_del_signature)

        set_meta_value('normalized_trna_seqs', norm_trna_seqs)
        set_meta_value('normalized_trna_seqs_without_potential_modifications', norm_trna_seqs_without_mods)
        set_meta_value('normalized_trna_seqs_containing_anticodon', norm_trna_seqs_containing_anticodon)
        set_meta_value('full_length_normalized_seqs', full_length_norm_seqs)
        set_meta_value('normalized_seqs_with_del_signature', norm_seqs_with_del_signature)

        set_meta_value('potentially_modified_seqs', mod_trna_seqs)
        set_meta_value('potentially_modified_seqs_with_del_signature', mod_seqs_with_del_signature)

        trnaseq_db.disconnect()

        # Use debug flag to write additional statistics of technical but not general interest to the
        # summary text file.
        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("Total tRNA reads (profiled and mapped)", total_trna_reads))
            f.write(get_summary_line("Reads profiled as tRNA", profiled_trna_reads))
            f.write(get_summary_line("Profiled reads ending 3'-CCA", trna_reads_with_threeprime_cca))
            f.write(get_summary_line("Profiled reads ending 3'-CC", trna_reads_with_threeprime_cc))
            f.write(get_summary_line("Profiled reads ending 3'-C", trna_reads_with_threeprime_c))
            f.write(get_summary_line("Profiled reads ending in other 3' termini", trna_reads_with_other_threeprime_termini))
            f.write(get_summary_line("Profiled reads containing anticodon", profiled_trna_reads_containing_anticodon))
            f.write(get_summary_line("Profiled reads spanning acceptor stem", full_length_trna_reads))
            if anvio.DEBUG:
                f.write(get_summary_line("Profiled reads with extrapolated 5' feature", trna_reads_with_extrapolated_fiveprime_feature))
            f.write(get_summary_line("Profiled reads with 1-%d extra 5' nts" % (self.min_length_of_long_fiveprime_extension - 1), profiled_trna_reads_with_short_fiveprime_extension))
            f.write(get_summary_line("Profiled reads with >%d extra 5' nts" % (self.min_length_of_long_fiveprime_extension - 1), profiled_trna_reads_with_long_fiveprime_extension))
            if anvio.DEBUG:
                f.write(get_summary_line("Profiled reads with recovered truncated profile", trna_reads_with_recovered_trunc_profile))
            f.write(get_summary_line("Reads mapped to tRNA interior", interior_mapped_reads))
            f.write(get_summary_line("Reads mapped to tRNA with extra 5' nts", fiveprime_mapped_reads))
            f.write(get_summary_line("Profiled reads specific to normalized seqs with del signature", trna_reads_specific_to_norm_seqs_with_del_signature))
            f.write(get_summary_line("Profiled reads nonspecific to normalized seqs with del signature", trna_reads_nonspecific_to_norm_seqs_with_del_signature))

            f.write(get_summary_line("Total unique tRNA seqs", total_uniq_trna_seqs))
            f.write(get_summary_line("Unique profiled tRNA seqs", uniq_profiled_trna_seqs))
            f.write(get_summary_line("Unique tRNA seqs ending 3'-CCA", uniq_trna_seqs_with_threeprime_cca))
            f.write(get_summary_line("Unique tRNA seqs ending 3'-CC", uniq_trna_seqs_with_threeprime_cc))
            f.write(get_summary_line("Unique tRNA seqs ending 3'-C", uniq_trna_seqs_with_threeprime_c))
            f.write(get_summary_line("Unique tRNA seqs ending in other 3' termini", uniq_trna_seqs_with_other_threeprime_termini))
            f.write(get_summary_line("Unique tRNA seqs containing anticodon", uniq_trna_seqs_containing_anticodon))
            f.write(get_summary_line("Unique tRNA seqs spanning acceptor stem", full_length_uniq_seqs))
            if anvio.DEBUG:
                f.write(get_summary_line("Unique tRNA seqs with extrapolated 5' feature", uniq_seqs_with_extrapolated_fiveprime_feature))
            f.write(get_summary_line("Unique tRNA seqs with 1-%d extra 5' nts" % (self.min_length_of_long_fiveprime_extension - 1), uniq_seqs_with_short_fiveprime_extension))
            f.write(get_summary_line("Unique tRNA seqs with >%d extra 5' nts" % (self.min_length_of_long_fiveprime_extension - 1), uniq_seqs_with_long_fiveprime_extension))
            if anvio.DEBUG:
                f.write(get_summary_line("Unique tRNA seqs with recovered truncated profile", uniq_seqs_with_recovered_trunc_profile))
            f.write(get_summary_line("Mean profiled reads per unique seq", mean_profiled_reads_per_uniq_seq))
            f.write(get_summary_line("Unique seqs mapped to tRNA interior", interior_mapped_uniq_seqs))
            f.write(get_summary_line("Unique seqs mapped to tRNA with extra 5' nts", fiveprime_mapped_uniq_seqs))
            f.write(get_summary_line("Mean mapped reads per unique seq", mean_mapped_reads_per_uniq_seq))

            f.write(get_summary_line("Trimmed tRNA seqs containing anticodon", trimmed_profiled_seqs_containing_anticodon))
            f.write(get_summary_line("Trimmed tRNA seqs spanning acceptor stem", full_length_trimmed_seqs))
            if anvio.DEBUG:
                f.write(get_summary_line("Trimmed tRNA seqs with extrapolated 5' feature", trimmed_seqs_with_extrapolated_fiveprime_feature))
                f.write(get_summary_line("Trimmed tRNA seqs specific to normalized seqs with del signature", trimmed_seqs_specific_to_norm_seqs_with_del_signature))
                f.write(get_summary_line("Trimmed tRNA seqs nonspecific to normalized seqs with del signature", trimmed_seqs_nonspecific_to_norm_seqs_with_del_signature))

            f.write(get_summary_line("Normalized tRNA seqs", norm_trna_seqs))
            f.write(get_summary_line("Normalized tRNA seqs without potential modifications", norm_trna_seqs_without_mods))
            f.write(get_summary_line("Normalized tRNA seqs containing anticodon", norm_trna_seqs_containing_anticodon))
            f.write(get_summary_line("Normalized tRNA seqs spanning acceptor stem", full_length_norm_seqs))
            f.write(get_summary_line("Normalized seqs with del signature", norm_seqs_with_del_signature))

            f.write(get_summary_line("Potentially modified seqs", mod_trna_seqs))
            f.write(get_summary_line("Potentially modified seqs with del signature", mod_seqs_with_del_signature))

        self.run.info("Summary", self.analysis_summary_path)


    def write_feature_table(self):
        self.progress.new("Writing tRNA-seq db table of profiled tRNA features")
        self.progress.update("...")

        feature_table_entries = []
        for uniq_seq in self.uniq_trna_seq_dict.values():
            class_name = type(uniq_seq).__name__
            if class_name == 'UniqueMappedSequence':
                continue
            if class_name == 'UniqueFullProfileSequence' or class_name == 'UniqueTransferredProfileSequence':
                has_complete_feature_set = uniq_seq.has_complete_feature_set
                num_extrapolated_fiveprime_nts = uniq_seq.num_extrapolated_fiveprime_nts
                extra_fiveprime_length = uniq_seq.extra_fiveprime_length
            elif class_name == 'UniqueTruncatedProfileSequence':
                has_complete_feature_set = False
                num_extrapolated_fiveprime_nts = 0
                extra_fiveprime_length = 0
            else:
                raise ConfigError(f"An object of class, {class_name}, is not recognized as a profiled tRNA sequence.")
            feature_table_entries.append(
                (uniq_seq.represent_name,
                 has_complete_feature_set,
                 uniq_seq.anticodon_string,
                 uniq_seq.anticodon_aa,
                 len(uniq_seq.seq_string),
                 # Zero-based start position of identified tRNA features within the read.
                 len(uniq_seq.seq_string) - uniq_seq.profiled_seq_length,
                 uniq_seq.num_conserved,
                 uniq_seq.num_unconserved,
                 uniq_seq.num_paired,
                 uniq_seq.num_unpaired,
                 num_extrapolated_fiveprime_nts,
                 extra_fiveprime_length,
                 uniq_seq.threeprime_terminus_length)
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
        self.progress.new("Writing tRNA-seq db table of unconserved nts in fully profiled tRNA")
        self.progress.update("...")

        unconserved_table_entries = []
        for uniq_seq in self.uniq_trna_seq_dict.values():
            if not isinstance(uniq_seq, UniqueFullProfileSequence):
                continue

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
        self.progress.new("Writing tRNA-seq db table of unpaired nts in profiled tRNA")
        self.progress.update("...")

        unpaired_table_entries = []
        for uniq_seq in self.uniq_trna_seq_dict.values():
            if not isinstance(uniq_seq, UniqueFullProfileSequence):
                continue

            for unpaired_tuple in uniq_seq.unpaired_info:
                unpaired_table_entries.append((uniq_seq.represent_name, ) + unpaired_tuple)

        if unpaired_table_entries:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('unpaired', ','.join('?' * len(tables.trnaseq_unpaired_table_structure))),
                                     unpaired_table_entries)
            trnaseq_db.disconnect()

        self.progress.end()


    def write_sequences_table(self):
        self.progress.new("Writing tRNA-seq db table of unique tRNA seqs")
        self.progress.update("...")

        sequences_table_entries = []
        class_name_id_info_dict = {'UniqueFullProfileSequence': 'full_profile',
                                   'UniqueTruncatedProfileSequence': 'truncated_profile',
                                   'UniqueTransferredProfileSequence': 'transferred_profile',
                                   'UniqueMappedSequence': 'mapped'}
        for uniq_seq in self.uniq_trna_seq_dict.values():
            class_name = type(uniq_seq).__name__
            try:
                id_info = class_name_id_info_dict[class_name]
            except KeyError:
                raise ConfigError(f"An object of class, {class_name}, is not recognized as a unique tRNA sequence. "
                                  f"The following classes are recognized: f{', '.join(class_name_id_info_dict)}.")

            sequences_table_entries.append(
                (uniq_seq.represent_name,
                 uniq_seq.read_count,
                 id_info,
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
        self.progress.new("Writing tRNA-seq db table of trimmed tRNA seqs")
        self.progress.update("...")

        trimmed_table_entries = []
        class_name_id_info_dict = {'TrimmedFullProfileSequence': 'full_profile',
                                   'TrimmedTruncatedProfileSequence': 'truncated_profile',
                                   'TrimmedMappedSequence': 'mapped'}
        no_threeprime_variants = tuple([0 for _ in THREEPRIME_VARIANTS])
        for trimmed_seq in self.trimmed_trna_seq_dict.values():
            class_name = type(trimmed_seq).__name__
            try:
                id_info = class_name_id_info_dict[class_name]
            except KeyError:
                raise ConfigError(f"An object of class, {class_name}, is not recognized as a trimmed tRNA sequence. "
                                  f"The following classes are recognized: f{', '.join(class_name_id_info_dict)}.")

            threeprime_termini = ''
            threeprime_terminus_read_counts = ''
            if id_info != 'mapped':
                for threeprime_terminus_string, read_count in trimmed_seq.read_threeprime_terminus_count_dict.items():
                    threeprime_termini += threeprime_terminus_string + ','
                    threeprime_terminus_read_counts += str(read_count) + ','

            trimmed_table_entries.append(
                (trimmed_seq.represent_name,
                 len(trimmed_seq.uniq_seqs),
                 trimmed_seq.read_count,
                 id_info,
                 trimmed_seq.seq_string,
                 len(trimmed_seq.norm_seq_represent_names),
                 trimmed_seq.uniq_with_extra_fiveprime_count if id_info != 'truncated_profile' else 0,
                 trimmed_seq.read_with_extra_fiveprime_count if id_info != 'truncated_profile' else 0,
                 threeprime_termini,
                 threeprime_terminus_read_counts)
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
        self.progress.new("Writing tRNA-seq db table of fragment-dereplicated tRNA seqs")
        self.progress.update("...")

        norm_table_entries = []
        for norm_seqs, id_info in ((self.norm_trna_seq_dict.values(), 'full_profile'),
                                   (self.norm_del_seq_dict.values(), 'deletion')):
            for norm_seq in norm_seqs:
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

                specific_threeprime_termini = ''
                specific_threeprime_terminus_read_counts = ''
                for threeprime_terminus_string, read_count in norm_seq.specific_read_threeprime_terminus_count_dict.items():
                    specific_threeprime_termini += threeprime_terminus_string + ','
                    specific_threeprime_terminus_read_counts += str(read_count) + ','

                nonspecific_threeprime_termini = ''
                nonspecific_threeprime_terminus_read_counts = ''
                for threeprime_terminus_string, read_count in norm_seq.nonspecific_read_threeprime_terminus_count_dict.items():
                    nonspecific_threeprime_termini += threeprime_terminus_string + ','
                    nonspecific_threeprime_terminus_read_counts += str(read_count) + ','

                norm_table_entries.append(
                    (norm_seq.represent_name,
                    len(norm_seq.trimmed_seqs),
                    id_info,
                    norm_seq.mean_specific_cov,
                    norm_seq.mean_nonspecific_cov,
                    ','.join(map(str, norm_seq.specific_covs)) + ',',
                    ','.join(map(str, norm_seq.nonspecific_covs)) + ',',
                    len(norm_seq.mod_seqs),
                    norm_seq.specific_read_count,
                    norm_seq.nonspecific_read_count,
                    norm_seq.specific_read_with_extra_fiveprime_count,
                    norm_seq.nonspecific_read_with_extra_fiveprime_count,
                    norm_seq.specific_mapped_read_count,
                    norm_seq.nonspecific_mapped_read_count,
                    specific_long_fiveprime_extensions,
                    specific_long_fiveprime_extension_read_counts,
                    nonspecific_long_fiveprime_extensions,
                    nonspecific_long_fiveprime_extension_read_counts,
                    specific_threeprime_termini,
                    specific_threeprime_terminus_read_counts,
                    nonspecific_threeprime_termini,
                    nonspecific_threeprime_terminus_read_counts)
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
        self.progress.new("Writing tRNA-seq db table of modified tRNA seqs")
        self.progress.update("...")

        mod_table_entries = []
        for mod_seq in self.mod_trna_seq_dict.values():
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

            specific_threeprime_termini = ''
            specific_threeprime_terminus_read_counts = ''
            for threeprime_terminus_string, read_count in mod_seq.specific_read_threeprime_terminus_count_dict.items():
                specific_threeprime_termini += threeprime_terminus_string + ','
                specific_threeprime_terminus_read_counts += str(read_count) + ','

            nonspecific_threeprime_termini = ''
            nonspecific_threeprime_terminus_read_counts = ''
            for threeprime_terminus_string, read_count in mod_seq.nonspecific_read_threeprime_terminus_count_dict.items():
                nonspecific_threeprime_termini += threeprime_terminus_string + ','
                nonspecific_threeprime_terminus_read_counts += str(read_count) + ','

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
                   mod_seq.specific_read_with_extra_fiveprime_count,
                   mod_seq.nonspecific_read_with_extra_fiveprime_count,
                   mod_seq.specific_mapped_read_count,
                   mod_seq.nonspecific_mapped_read_count,
                   specific_long_fiveprime_extensions,
                   specific_long_fiveprime_extension_read_counts,
                   nonspecific_long_fiveprime_extensions,
                   nonspecific_long_fiveprime_extension_read_counts,
                   specific_threeprime_termini,
                   specific_threeprime_terminus_read_counts,
                   nonspecific_threeprime_termini,
                   nonspecific_threeprime_terminus_read_counts)
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


    def write_unique_nontrna_supplement(self):
        self.progress.new("Writing file of unique seqs not identified as tRNA")
        self.progress.update("...")

        with open(self.uniq_nontrna_path, 'w') as nontrna_file:
            nontrna_file.write("\t".join(self.UNIQ_NONTRNA_HEADER) + "\n")
            for uniq_seq in self.uniq_nontrna_seq_dict.values():
                nontrna_file.write(uniq_seq.represent_name + "\t"
                                   + str(uniq_seq.read_count) + "\t"
                                   + "\t"
                                   + uniq_seq.seq_string + "\n")
            for uniq_seq in self.uniq_trunc_seq_dict.values():
                nontrna_file.write(uniq_seq.represent_name + "\t"
                                   + str(uniq_seq.read_count) + "\t"
                                   + str(uniq_seq.trunc_profile_index) + "\t"
                                   + uniq_seq.seq_string + "\n")

        self.progress.end()

        self.run.info("Unique non-tRNA supplement", self.uniq_nontrna_path)


    def write_trimmed_supplement(self):
        """Write a supplementary file showing the spectrum of 5'/3' extensions of trimmed, fully
        profiled tRNA sequences."""
        self.progress.new("Writing file showing 5'/3' ends of trimmed, fully profiled tRNA seqs")
        self.progress.update("...")

        with open(self.trimmed_ends_path, 'w') as trimmed_file:
            trimmed_file.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for trimmed_seq in sorted(self.trimmed_trna_seq_dict.values(), key=lambda trimmed_seq: -trimmed_seq.read_count):
                if not isinstance(trimmed_seq, TrimmedFullProfileSequence):
                    continue

                represent_name = trimmed_seq.represent_name
                for uniq_seq in sorted(trimmed_seq.uniq_seqs,
                                       key=lambda uniq_seq: (-uniq_seq.extra_fiveprime_length, -uniq_seq.threeprime_terminus_length)):
                    trimmed_file.write(represent_name + "\t"
                                       + uniq_seq.represent_name + "\t"
                                       + uniq_seq.seq_string[: uniq_seq.extra_fiveprime_length] + "\t"
                                       + uniq_seq.seq_string[len(uniq_seq.seq_string) - uniq_seq.threeprime_terminus_length: ] + "\t"
                                       + str(uniq_seq.read_count) + "\n")

        self.progress.end()

        self.run.info("Trimmed tRNA supplement", self.trimmed_ends_path)


def profile_worker(input_queue, output_queue, profiler):
    """This client for `trnaidentifier.Profiler.profile` is located outside the `TRNASeqDataset`
    class to allow multiprocessing."""
    while True:
        seq_string, represent_name, read_count = input_queue.get()
        output_queue.put((profiler.profile(seq_string, name=represent_name), read_count))


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
        'id_info',
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
        self.min_del_fraction = A('min_deletion_fraction')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default

        if not self.project_name:
            raise ConfigError("Please specify a name for the collection of input tRNA-seq dbs "
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
            raise ConfigError("Sample IDs in each input tRNA-seq db must be unique. This is not "
                              "the case with your input. Here are the sample names so you can see "
                              "which ones occur more than once: '%s'" % (", ".join(self.trnaseq_db_sample_ids)))
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
                              "To remove the limit on tRNA seeds reported to the contigs db, "
                              "provide a value of -1. Otherwise provide an integer greater than 0.")

        self.run.info("Input tRNA-seq dbs", ", ".join(self.trnaseq_db_paths))
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
                self.run.warning("Not all input tRNA-seq dbs have the same version number, "
                                 "but since you have used the `--force` flag, `anvi-convert-trnaseq-database` "
                                 "will proceed though this is dangerous and may lead to errors. "
                                 f"Here is the version number of each database:\n{trnaseq_db_version_report}")
            else:
                raise ConfigError("Not all input tRNA-seq dbs have the same version number. "
                                  f"Here is the version number of each db:\n{trnaseq_db_version_report}")


    def set_treatment_preference(self):
        if not self.preferred_treatment:
            return

        input_treatments = [inner_dict['treatment'] for inner_dict in self.trnaseq_dbs_info_dict.values()]
        self.preferred_trnaseq_db_sample_ids = []
        self.preferred_trnaseq_db_nums = []
        if self.preferred_treatment not in input_treatments:
            raise ConfigError("You provided a preferred treatment type, %s, "
                              "but it was not found in any of the input dbs, "
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
                raise ConfigError("The nonspecific profile db types provided by `--nonspecific-output` are not recognized. "
                                  "The db types must be comma separated without spaces, "
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
        self.progress.new("Loading seq info from tRNA-seq dbs")
        self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} dbs loaded")

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
            self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} dbs loaded")

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
             id_info,
             mean_specific_cov,
             specific_covs_string,
             nonspecific_covs_string) = norm_seq_info

            if id_info == 'deletion':
                # Ignore normalized sequences with deletions. The coverage of deletions themselves
                # is recorded in the parent modified sequence, but the contribution of these
                # sequences to nucleotide coverage is ignored. Inclusion of these sequences would
                # produce numerous complications (e.g., they don't have feature profiles).
                continue

            norm_seq_summary = NormalizedSeqSummary()
            norm_seq_summary.name = name
            norm_seq_summary.sample_id = sample_id
            norm_seq_summary.mean_specific_cov = mean_specific_cov
            # There is always a trailing comma in the coverage strings.
            norm_seq_summary.specific_covs = np.fromiter(map(int, specific_covs_string.split(',')[: -1]), int)
            norm_seq_summary.nonspecific_covs = np.fromiter(map(int, nonspecific_covs_string.split(',')[: -1]), int)
            (norm_seq_summary.seq_string,
             norm_seq_summary.anticodon_seq_string,
             norm_seq_summary.feature_threshold_start) = seq_string_and_feature_df.loc[norm_seq_summary.name, ['sequence', 'anticodon_sequence', self.feature_threshold + '_start']]

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
                     nonspecific_del_cov) in zip(sorted(set(map(int, del_positions.replace(';', ',').split(',')[: -1]))),
                                                 map(int, del_specific_covs.split(',')[: -1]),
                                                 map(int, del_nonspecific_covs.split(',')[: -1])):
                    specific_del_covs[del_pos] = specific_del_cov
                    nonspecific_del_covs[del_pos] = nonspecific_del_cov
            mod_seq_summary.specific_del_covs = specific_del_covs
            mod_seq_summary.nonspecific_del_covs = nonspecific_del_covs

            mod_seq_summary.norm_seq_summaries = []
            for norm_seq_name in names_of_norm_seqs_without_dels.split(','):
                norm_seq_summary = norm_seq_summary_dict[norm_seq_name]

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
        self.progress.new("Forming seed seqs from input samples")

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
        self.progress.new("Generating a contigs db of tRNA seeds")
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
                sha1(seed_seq.seq_string.encode('utf-8')).hexdigest(), # "Gene unique identifier"
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
        self.progress.new(f"Generating {db_cov_type} profile db")
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
