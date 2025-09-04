# -*- coding: utf-8
# pylint: disable=line-too-long
"""Library for tRNA-seq dataset operations

`bin/anvi-trnaseq` and `bin/anvi-merge-trnaseq` are the default clients using this module.
`anvi-trnaseq` instantiates a `TRNASeqDataset` object. `anvi-merge-trnaseq` instantiates a
`DatabaseMerger` object. The clients call the objects' `process` methods to start the analytic
workflows.

Each sequence library in an experiment is processed separately as a `TRNASeqDataset`, storing an
information-rich anvi'o tRNA-seq database. `DatabaseMerger` finds reference seed sequences from a
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
import math
import time
import queue
import random
import shutil
import argparse
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

from hashlib import sha1
from itertools import chain
from functools import partial
from bisect import bisect_left
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from collections import defaultdict, deque, OrderedDict

import anvio
import anvio.dbops as dbops
import anvio.tables as tables
import anvio.fastalib as fastalib
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.tables.miscdata as miscdata
import anvio.trnaidentifier as trnaidentifier
import anvio.auxiliarydataops as auxiliarydataops

from anvio.dbinfo import DBInfo
from anvio.errors import ConfigError
from anvio.sequence import Dereplicator
from anvio.drivers.vmatch import Vmatch
from anvio.agglomeration import Agglomerator
from anvio.tables.views import TablesForViews
from anvio.tables.miscdata import TableForLayerOrders
from anvio.dbinfo import is_trnaseq_db
from anvio.utils.algorithms import convert_binary_blob_to_numpy_array
from anvio.utils.fasta import check_fasta_id_formatting
from anvio.utils.validation import check_sample_id


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


pp = terminal.pretty_print

MAXSIZE = sys.maxsize

ALL_NTS = tuple(constants.nucleotides)
UNAMBIG_NTS = ('A', 'C', 'G', 'T')
UNAMBIG_NTS_LIST = list(UNAMBIG_NTS)
NT_INT_DICT = {nt: i for i, nt in enumerate(UNAMBIG_NTS, start=1)}
INT_NT_DICT = {i: nt for i, nt in enumerate(UNAMBIG_NTS, start=1)}
# NUM_NT_BINS is used in counting the number of distinct nucleotides (1, 2, 3, or 4) at positions in
# internally ungapped alignments: there is one bin (value 0) for end gaps in the alignment.
NUM_NT_BINS = len(UNAMBIG_NTS) + 1

AMINO_ACIDS = constants.amino_acids
ANTICODON_AA_DICT = constants.anticodon_to_AA
ANTICODONS = list(ANTICODON_AA_DICT)
AA_ANTICODON_DICT = constants.AA_to_anticodons
TRNA_FEATURE_NAMES = constants.TRNA_FEATURE_NAMES

# The user can specify in `anvi-trnaseq` what defines a long (biological vs. non-templated) 5'
# extension. The variable is set by `TRNASeqDataset.__init__`.
MIN_LENGTH_LONG_5PRIME_EXTENSION = 4
# The user can profile tRNA that does not end in any 3' terminus, such as CCA. Normally, such
# sequences are not profiled but may be mapped.
PROFILE_ABSENT_3PRIME_TERMINUS = False


class UniqueSequence(object):
    """A dereplicated tRNA-seq read."""
    # Instances of this class are called `seq_U`.

    __slots__ = ('string', 'name', 'read_count')

    def __init__(self, string, name, read_count):
        self.string = string
        self.name = name
        self.read_count = read_count


class UniqueProfileSequence(UniqueSequence):
    """A tRNA feature profile, either full or truncated, was assigned to the sequence."""
    # Instances of this class are called `seq_Up`.

    __slots__ = (
        'feature_starts',
        'feature_stops',
        'has_complete_feature_set',
        'has_His_G',
        'alpha_start',
        'alpha_stop',
        'beta_start',
        'beta_stop',
        'anticodon_string',
        'anticodon_aa',
        'contains_anticodon',
        'length_3prime_terminus',
        'num_conserved',
        'num_unconserved',
        'num_paired',
        'num_unpaired',
        'unconserved_info',
        'unpaired_info',
        'profiled_seq_length',
        'name_T'
    )

    def __init__(self, string, name, read_count, profile):
        super().__init__(string, name, read_count)

        self.feature_starts = tuple([f.start_pos if hasattr(f, 'start_pos') else f.start_positions for f in profile.features])
        self.feature_stops = tuple([f.stop_pos if hasattr(f, 'stop_pos') else f.stop_positions for f in profile.features])
        self.alpha_start = profile.alpha_start
        self.alpha_stop = profile.alpha_stop
        self.beta_start = profile.beta_start
        self.beta_stop = profile.beta_stop
        self.anticodon_string = anticodon = profile.anticodon_seq
        self.anticodon_aa = profile.anticodon_aa if profile.anticodon_aa else None
        self.contains_anticodon = True if anticodon else False
        self.length_3prime_terminus = len(profile.threeprime_terminus_seq)
        self.num_conserved = profile.num_conserved
        self.num_unconserved = profile.num_unconserved
        self.num_paired = profile.num_paired
        self.num_unpaired = profile.num_unpaired
        self.unconserved_info = tuple(profile.unconserved_info) if profile.unconserved_info else tuple()
        self.unpaired_info = tuple(profile.unpaired_info) if profile.unpaired_info else tuple()
        self.profiled_seq_length = len(profile.profiled_seq)
        self.name_T = None


class UniqueFullProfileSequence(UniqueProfileSequence):
    """A full tRNA feature profile was assigned to the sequence."""
    # Instances of this class are called `seq_Uf`.

    __slots__ = (
        'has_complete_feature_set',
        'has_His_G',
        'num_extrap_5prime_nts',
        'xtra_5prime_length'
    )

    def __init__(self, string, name, read_count, profile):
        super().__init__(string, name, read_count, profile)

        self.has_complete_feature_set = profile.has_complete_feature_set
        self.has_His_G = True if profile.features[0].name == 'tRNA-His position 0' else False
        self.num_extrap_5prime_nts = profile.num_in_extrapolated_fiveprime_feature
        self.xtra_5prime_length = 0 if profile.num_extra_fiveprime is None else profile.num_extra_fiveprime


class UniqueTruncatedProfileSequence(UniqueProfileSequence):
    """A truncated tRNA feature profile was assigned to the sequence."""
    # Instances of this class are called `seq_Uc`.

    __slots__ = ('trunc_profile_index', )

    def __init__(self, string, name, read_count, profile):
        super().__init__(string, name, read_count, profile)

        self.trunc_profile_index = profile.trunc_profile_index


class UniqueTransferredProfileSequence(UniqueFullProfileSequence):
    """This object is generated as part of the determination of Nf from Tf. This type of seq is
    produced in the special circumstance that the profile of a shorter Tf is transferred to a longer
    Tf, because the longer Tf was originally found to have a complete profile, but a shorter 3'
    subseq also had a complete profile; so, parsimoniously, the profile of the shorter Tf was
    transferred to the longer, and the additional 5' nts of the longer Tf reclassified as extra nts
    beyond the 5' end of a mature tRNA sequence."""
    # Instances of this class are called `seq_Us`.

    __slots__ = ('defunct_Uf', )

    def __init__(self, defunct_Uf, replacement_dict):
        UniqueSequence.__init__(self, defunct_Uf.string, defunct_Uf.name, defunct_Uf.read_count)

        string_U = defunct_Uf.string
        length_U = len(string_U)
        stop_T_in_U = length_U - string_U[::-1].index(replacement_dict['string_T'][::-1])

        feature_starts = []
        for feature_start_from_T_3prime in replacement_dict['feature_starts_from_T_3prime']:
            if isinstance(feature_start_from_T_3prime, int):
                feature_starts.append(stop_T_in_U + feature_start_from_T_3prime)
            else:
                feature_starts.append(tuple([stop_T_in_U + strand_start_from_T_3prime for strand_start_from_T_3prime in feature_start_from_T_3prime]))
        self.feature_starts = tuple(feature_starts)

        feature_stops = []
        for feature_stop_from_T_3prime in replacement_dict['feature_stops_from_T_3prime']:
            if isinstance(feature_stop_from_T_3prime, int):
                feature_stops.append(stop_T_in_U + feature_stop_from_T_3prime)
            else:
                feature_stops.append(tuple([stop_T_in_U + strand_stop_from_T_3prime for strand_stop_from_T_3prime in feature_stop_from_T_3prime]))
        self.feature_stops = tuple(feature_stops)

        self.has_complete_feature_set = True
        self.num_extrap_5prime_nts = 0
        self.has_His_G = replacement_dict['has_His_G']
        self.alpha_start = None if replacement_dict['alpha_start_from_T_3prime'] is None else stop_T_in_U + replacement_dict['alpha_start_from_T_3prime']
        self.alpha_stop = None if replacement_dict['alpha_stop_from_T_3prime'] is None else stop_T_in_U + replacement_dict['alpha_stop_from_T_3prime']
        self.beta_start = None if replacement_dict['beta_start_from_T_3prime'] is None else stop_T_in_U + replacement_dict['beta_start_from_T_3prime']
        self.beta_stop = None if replacement_dict['beta_stop_from_T_3prime'] is None else stop_T_in_U + replacement_dict['beta_stop_from_T_3prime']
        self.anticodon_string = replacement_dict['anticodon_string']
        self.anticodon_aa = replacement_dict['anticodon_aa']
        self.contains_anticodon = replacement_dict['contains_anticodon']
        self.length_3prime_terminus = length_U - stop_T_in_U
        self.num_conserved = replacement_dict['num_conserved']
        self.num_unconserved = replacement_dict['num_unconserved']
        self.num_paired = replacement_dict['num_paired']
        self.num_unpaired = replacement_dict['num_unpaired']

        unconserved_info = []
        for unconserved_tuple in replacement_dict['unconserved_info_from_T_3prime']:
            unconserved_info.append((stop_T_in_U + unconserved_tuple[0],
                                     unconserved_tuple[1],
                                     unconserved_tuple[2]))
        self.unconserved_info = tuple(unconserved_info)

        unpaired_info = []
        for unpaired_tuple in replacement_dict['unpaired_info_from_T_3prime']:
            unpaired_info.append((stop_T_in_U + unpaired_tuple[0],
                                  stop_T_in_U + unpaired_tuple[1],
                                  unpaired_tuple[2],
                                  unpaired_tuple[3]))
        self.unpaired_info = tuple(unpaired_info)

        self.profiled_seq_length = replacement_dict['profiled_seq_without_terminus_length'] + self.length_3prime_terminus
        self.xtra_5prime_length = length_U - self.profiled_seq_length
        self.name_T = None

        # Store the defunct profile information for posterity.
        self.defunct_Uf = defunct_Uf


class UniqueMappedSequence(UniqueSequence):
    """This object is generated in the identification of tRNA fragments by mapping."""
    # Instances of this class are called `seq_Um`.

    __slots__ = ('xtra_5prime_length', 'name_T')

    def __init__(self, string, name, read_count, xtra_5prime_length=0):
        super().__init__(string, name, read_count)

        self.xtra_5prime_length = xtra_5prime_length
        self.name_T = None


class UniqueIndelSequence(UniqueSequence):
    """This object is generated in the identification of tRNA sequences with indels. Nq are found to
    have indels, which can contradict feature profiles and the lengths of 5' extensions and 3'
    termini in existing U. No profile is assigned to this object."""
    # Instances of this class are called `seq_Ui`.

    __slots__ = ('orig_U', 'xtra_5prime_length', 'length_3prime_terminus', 'name_T')

    def __init__(self, seq_U, length_3prime_terminus=0, xtra_5prime_length=0):
        super().__init__(seq_U.string, seq_U.name, seq_U.read_count)

        self.orig_U = seq_U
        self.length_3prime_terminus = length_3prime_terminus
        self.xtra_5prime_length = xtra_5prime_length
        self.name_T = None


class TrimmedSequence(object):
    """A tRNA sequence with bases trimmed 5' of the acceptor stem (or 5'-G in the case of tRNA-His)
    and 3' of the discriminator.

    The purpose of trimming is to collapse non-biological variability prevalent at the ends of
    reads.

    EXAMPLE 1:
    E. coli tRNA-Ala-GGC-1-1
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA

    This collapses to the following T, removing the 3' terminus (the acceptor happens to be genomic
    rather than post-transcriptionally added in this example, but it doesn't matter):
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA

    Examples of possible Up that collapse to T:
    AGGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA
     GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCAC

    EXAMPLE 2:
    3' fragment of the same tRNA, ending in 3'-CC rather than canonical 3'-CCA
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACC

    This collapses to the following T, removing 3'-CC:
                                TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA
    """
    # Instances of this class are called T.

    __slots__ = ('string', 'read_count', 'names_U', 'names_N')

    def __init__(self, string, seqs_U):
        self.string = string
        self.names_U = tuple([seq_U.name for seq_U in seqs_U])
        for seq_U in seqs_U:
            if seq_U.name_T is not None:
                raise ConfigError(f"The unique sequence with the name {seq_U.name} "
                                  f"was already assigned to a trimmed sequence with the name {seq_U.name_T} "
                                  "and so cannot be assigned to a new trimmed sequence.")
        self.read_count = sum([seq_U.read_count for seq_U in seqs_U])
        self.names_N = []


def get_representative_unique_sequence(seqs_U):
    """The representative U in Tf and Ti is chosen as follows:
    1. Most abundant full-length tRNA without extra 5' nts, ignoring the 3' terminus nts, OR
    2. Most abundant full-length tRNA with extra 5' nts, OR
    3. Most abundant 3' tRNA fragment
    Sort U such that the first is the most abundant+longest and the last is the least
    abundant+shortest."""
    seqs_U.sort(key=lambda seq_U: (-seq_U.xtra_5prime_length, -seq_U.read_count, seq_U.name))
    if seqs_U[0].xtra_5prime_length > 0:
        if seqs_U[-1].xtra_5prime_length == 0:
            # If the first U has extra 5' nts and the last has none, then the last U and others
            # without extra 5' nts must be a full-length tRNA (ignoring the 3' terminus). Therefore,
            # select the most abundant of these full-length tRNAs with extra 5' nts as the
            # representative seq.
            represent_seq_U = sorted(seqs_U, key=lambda seq_U: (-seq_U.xtra_5prime_length, seq_U.read_count, seq_U.name))[-1]
        else:
            represent_seq_U = seqs_U[0]
    else:
        # In this case, ALL U are EITHER full-length tRNA OR a 3' tRNA fragment.
        represent_seq_U = seqs_U[0]

    return represent_seq_U


def get_read_threeprime_terminus_count_dict(seqs_U):
    """Count the number of reads with each 3' terminus in U."""
    read_3prime_terminus_count_dict = defaultdict(int)
    for seq_U in seqs_U:
        if seq_U.length_3prime_terminus:
            read_3prime_terminus_count_dict[seq_U.string[-seq_U.length_3prime_terminus: ]] += seq_U.read_count
        else:
            read_3prime_terminus_count_dict[''] += seq_U.read_count

    return read_3prime_terminus_count_dict


def get_extra_fiveprime_info(seqs_U):
    """Get information on the extra 5' nucleotides in unique sequences."""
    uniq_with_xtra_5prime_count = 0
    read_with_xtra_5prime_count = 0
    # Find the number of reads containing each unique 5' extension.
    long_5prime_extension_dict = {}
    for seq_U in seqs_U:
        if seq_U.xtra_5prime_length:
            uniq_with_xtra_5prime_count += 1
            read_with_xtra_5prime_count += seq_U.read_count
            if seq_U.xtra_5prime_length >= MIN_LENGTH_LONG_5PRIME_EXTENSION:
                long_5prime_extension_dict[seq_U.string[: seq_U.xtra_5prime_length]] = seq_U.read_count

    return uniq_with_xtra_5prime_count, read_with_xtra_5prime_count, long_5prime_extension_dict


class TrimmedFullProfileSequence(TrimmedSequence):
    """This object is formed from sequences with a full tRNA feature profile."""
    # Instances of this class are called `Tf` and can derive from `Uf` and optionally `Us` objects.

    __slots__ = (
        'name',
        'categories_U',
        'feature_starts',
        'feature_stops',
        'contains_anticodon',
        'read_3prime_terminus_count_dict',
        'has_complete_feature_set',
        'has_His_G',
        'num_extrap_5prime_nts',
        'uniq_with_xtra_5prime_count',
        'read_with_xtra_5prime_count',
        'long_5prime_extension_dict'
    )

    def __init__(self, string, seqs_U):
        # U are sorted in place when finding the representative unique sequence. Call the parent
        # class constructor after this to ensure recorded U names are in the correct order, with the
        # representative name coming first.
        represent_U = get_representative_unique_sequence(seqs_U)
        super().__init__(string, seqs_U)
        self.name = name = represent_U.name
        categories_U = []
        for seq_U in seqs_U:
            if isinstance(seq_U, UniqueTransferredProfileSequence):
                categories_U.append('Us')
            elif isinstance(seq_U, UniqueFullProfileSequence):
                categories_U.append('Uf')
            else:
                raise Exception(f"A unique sequence with name `{seq_U.name}` of class `{type(seq_U)}` was encountered.")
            seq_U.name_T = name
        self.categories_U = tuple(categories_U)

        # Assume that the feature profile indices of the representative U are the same as the other
        # U. The 3' terminus is the last feature in the profile and not part of T.
        self.feature_starts = represent_U.feature_starts[: -1] if represent_U.feature_starts else None
        self.feature_stops = represent_U.feature_stops[: -1] if represent_U.feature_stops else None
        self.contains_anticodon = represent_U.contains_anticodon
        self.has_complete_feature_set = represent_U.has_complete_feature_set
        self.has_His_G = represent_U.has_His_G
        self.num_extrap_5prime_nts = represent_U.num_extrap_5prime_nts

        self.read_3prime_terminus_count_dict = get_read_threeprime_terminus_count_dict(seqs_U)

        (self.uniq_with_xtra_5prime_count,
         self.read_with_xtra_5prime_count,
         self.long_5prime_extension_dict) = get_extra_fiveprime_info(seqs_U)


class TrimmedTruncatedProfileSequence(TrimmedSequence):
    """This object is formed from sequences with a truncated tRNA feature profile."""
    # Instances of this class are called `Tc` and derive from `Uc`. All instances are initially
    # categorized as `nontrna`, which can later change to `trna`. The category encompasses all
    # component `Uc`.

    __slots__ = (
        'name',
        'category',
        'feature_starts',
        'feature_stops',
        'contains_anticodon',
        'read_3prime_terminus_count_dict',
        'trunc_profile_index'
    )

    def __init__(self, string, seqs_U):
        # Make the most abundant U the representative sequence.
        seqs_U.sort(key=lambda seq_U: (-seq_U.read_count, seq_U.name))
        super().__init__(string, seqs_U)
        represent_U = seqs_U[0]
        self.name = name = represent_U.name
        self.category = 'nontrna'
        for seq_U in seqs_U:
            seq_U.name_T = name

        # Assume that the feature profile indices of the representative U are the same as the other
        # U. The 3' terminus is the last feature in the profile and not part of T.
        self.feature_starts = represent_U.feature_starts[: -1] if represent_U.feature_starts else None
        self.feature_stops = represent_U.feature_stops[: -1] if represent_U.feature_stops else None
        self.contains_anticodon = represent_U.contains_anticodon

        self.read_3prime_terminus_count_dict = get_read_threeprime_terminus_count_dict(seqs_U)

        self.trunc_profile_index = represent_U.trunc_profile_index


class TrimmedMappedSequence(TrimmedSequence):
    """This object is formed from a single Um in the process of mapping unique unprofiled seqs to N.
    It is not like the other T objects. Its purpose is to be one of the T objects added to an N,
    however, unliked Tp, no 5' nts are trimmed from the U string in creating the T string (and
    mapped seqs do not have extra 3' nts). The reason for this is that the 5' extension may
    represent all but a small number of nts in the seq, so it is best not to dereplicate mapped seqs
    identical in the non-5' section by lumping them together as the same T."""
    # Instances of this class are called `Tm` and derive from a single `Um`.

    __slots__ = (
        'name',
        'uniq_with_xtra_5prime_count',
        'read_with_xtra_5prime_count',
        'long_5prime_extension_dict'
    )

    def __init__(self, seq_U):
        super().__init__(seq_U.string, [seq_U])

        self.name = seq_U.name
        seq_U.name_T = self.name

        xtra_5prime_length = seq_U.xtra_5prime_length
        self.uniq_with_xtra_5prime_count = 1 if xtra_5prime_length else 0
        self.read_with_xtra_5prime_count = seq_U.read_count if xtra_5prime_length else 0
        self.long_5prime_extension_dict = {}
        if xtra_5prime_length >= MIN_LENGTH_LONG_5PRIME_EXTENSION:
            self.long_5prime_extension_dict[seq_U.string[: xtra_5prime_length]] = seq_U.read_count


class TrimmedIndelSequence(TrimmedSequence):
    """This object is formed in the identification of tRNA sequences with indels."""
    # Instances of this class are called `Ti` and derive from `Ui`.

    __slots__ = (
        'name',
        'read_3prime_terminus_count_dict',
        'uniq_with_xtra_5prime_count',
        'read_with_xtra_5prime_count',
        'long_5prime_extension_dict'
    )

    def __init__(self, string, seqs_U):
        represent_U = get_representative_unique_sequence(seqs_U)
        super().__init__(string, seqs_U)
        self.name = name = represent_U.name
        for seq_U in seqs_U:
            seq_U.name_T = name

        self.read_3prime_terminus_count_dict = get_read_threeprime_terminus_count_dict(seqs_U)

        (self.uniq_with_xtra_5prime_count,
         self.read_with_xtra_5prime_count,
         self.long_5prime_extension_dict) = get_extra_fiveprime_info(seqs_U)


class NormalizedSequence(object):
    """A tRNA sequence that can contain shorter tRNA fragment subsequences.

    N are derived from T. Tp are first prefix-dereplicated from the 3' end by the method,
    `TRNASeqDataset.dereplicate_threeprime`. The longest Tp in a cluster of dereplicated Tp becomes
    the representative seq of N. `TRNASeqDataset.map_fragments` subsequently maps unprofiled reads
    to the set of N. Mapped tRNA fragments are added as T, and the `init` method is called to
    finalize attributes of each N.

    EXAMPLE:
    Consider the full-length and fragmentary E. coli tRNA-Ala-GGC-1-1 used in T examples.
    GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA
                               TTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA

    The seqs collapse into a single N when dereplicated from the 3' end. N is represented by the
    longer seq, which should be the first in the list of T added to N.
    """
    # Instances of this class are called `N`.

    __slots__ = (
        'name',
        'names_T',
        'string',
        'starts_T_in_N',
        'stops_T_in_N',
        'spec_read_count',
        'nonspec_read_count',
        'spec_covs',
        'nonspec_covs',
        'mean_spec_cov',
        'mean_nonspec_cov',
        'names_M'
    )

    def __init__(self, seqs_T, starts_T_in_N=None, stops_T_in_N=None, skip_init=False):
        represent_T = seqs_T[0]
        self.name = name = represent_T.name
        for seq_T in seqs_T:
            seq_T.names_N.append(name)
        self.names_T = [seq_T.name for seq_T in seqs_T]
        self.string = represent_T.string

        if starts_T_in_N and stops_T_in_N:
            self.starts_T_in_N = starts_T_in_N
            self.stops_T_in_N = stops_T_in_N
        elif not starts_T_in_N and not stops_T_in_N:
            # This is what occurs in `anvi-trnaseq` for the instantiation of
            # `NormalizedFullProfileSequence` and `NormalizedTruncatedProfileSequence` objects. All
            # Tp provided as input are aligned from the 3' end. Tm that can be aligned to other
            # places in N are added later.
            length_N = len(self.string)
            self.starts_T_in_N = [length_N - len(seq_T.string) for seq_T in seqs_T]
            self.stops_T_in_N = [length_N] * len(seqs_T)
        else:
            self.starts_T_in_N = None
            self.stops_T_in_N = None

        # It is useful to know which M, if any, contain this N. Though Ni can theoretically be
        # assigned to multiple M, this is not allowed for simplicity's sake.
        self.names_M = []

        if skip_init:
            self.spec_read_count = None
            self.nonspec_read_count = None
            self.spec_covs = None
            self.nonspec_covs = None
            self.mean_spec_cov = None
            self.mean_nonspec_cov = None
        else:
            self.init(seqs_T)


    def init(self, seqs_T):
        """Set attributes for a "finalized" set of input T, potentially including those that were
        added after object instantiation, such as Tm."""
        self.names_T = tuple(self.names_T)

        # Specific reads are those that are only contained in this N.
        spec_read_count = 0
        nonspec_read_count = 0
        spec_covs = np.zeros(len(self.string), dtype=int)
        nonspec_covs = np.zeros(len(self.string), dtype=int)

        for seq_T, start_T_in_N, stop_T_in_N in zip(seqs_T, self.starts_T_in_N, self.stops_T_in_N):
            if len(seq_T.names_N) == 1:
                spec_read_count += seq_T.read_count
                spec_covs[start_T_in_N: stop_T_in_N] += seq_T.read_count
            else:
                nonspec_read_count += seq_T.read_count
                nonspec_covs[start_T_in_N: stop_T_in_N] += seq_T.read_count

        self.spec_read_count = spec_read_count
        self.nonspec_read_count = nonspec_read_count
        self.spec_covs = spec_covs
        self.nonspec_covs = nonspec_covs
        self.mean_spec_cov = spec_covs.mean()
        self.mean_nonspec_cov = nonspec_covs.mean()


class NormalizedFullProfileSequence(NormalizedSequence):
    """This object is instantiated with `TrimmedFullProfileSequence` objects in the workflow.
    `TrimmedTruncatedProfileSequence` and `TrimmedMappedSequence` objects are potentially added
    after instantiation."""
    # Instances of this class are called `Nf`.

    __slots__ = (
        'categories_T',
        'feature_starts',
        'feature_stops',
        'contains_anticodon',
        'has_complete_feature_set',
        'has_His_G',
        'spec_read_xtra_5prime_count',
        'nonspec_read_xtra_5prime_count',
        'spec_map_seq_count',
        'nonspec_map_seq_count',
        'spec_map_read_count',
        'nonspec_map_read_count',
        'absent_3prime_terminus_seq_count',
        'absent_3prime_terminus_read_count',
        'spec_long_5prime_extension_dict',
        'nonspec_long_5prime_extension_dict',
        'spec_read_3prime_terminus_count_dict',
        'nonspec_read_3prime_terminus_count_dict'
    )

    def __init__(self, seqs_T, starts_T_in_N=None, stops_T_in_N=None):
        super().__init__(seqs_T, starts_T_in_N=None, stops_T_in_N=None, skip_init=True)

        categories_T = []
        for seq_T in seqs_T:
            if isinstance(seq_T, TrimmedFullProfileSequence):
                categories_T.append('Tf')
            elif isinstance(seq_T, TrimmedTruncatedProfileSequence):
                categories_T.append('Tc_trna')
            else:
                raise Exception(f"A trimmed seq ({seq_T.name}) of the unexpected class `{type(seq_T)}` was used to instantiate `Nf`.")
        self.categories_T = categories_T

        represent_T = seqs_T[0]
        self.feature_starts = represent_T.feature_starts
        self.feature_stops = represent_T.feature_stops
        self.contains_anticodon = represent_T.contains_anticodon
        self.has_complete_feature_set = represent_T.has_complete_feature_set
        self.has_His_G = represent_T.has_His_G

        self.spec_read_xtra_5prime_count = None
        self.nonspec_read_xtra_5prime_count = None
        self.spec_map_read_count = None
        self.nonspec_map_read_count = None
        self.spec_long_5prime_extension_dict = None
        self.nonspec_long_5prime_extension_dict = None
        self.spec_read_3prime_terminus_count_dict = None
        self.nonspec_read_3prime_terminus_count_dict = None


    def init(self, seqs_T):
        """Set attributes for a "finalized" set of input trimmed sequences."""
        super().init(seqs_T)
        self.categories_T = tuple(self.categories_T)

        spec_map_seq_count = 0
        nonspec_map_seq_count = 0
        spec_map_read_count = 0
        nonspec_map_read_count = 0
        absent_3prime_terminus_seq_count = 0
        absent_3prime_terminus_read_count = 0
        spec_read_xtra_5prime_count = 0
        nonspec_read_xtra_5prime_count = 0
        spec_long_5prime_extension_dict = defaultdict(int)
        nonspec_long_5prime_extension_dict = defaultdict(int)
        spec_read_3prime_terminus_count_dict = defaultdict(int)
        nonspec_read_3prime_terminus_count_dict = defaultdict(int)

        length_N = len(self.string)
        for seq_T, start_T_in_N, stop_T_in_N in zip(seqs_T, self.starts_T_in_N, self.stops_T_in_N):
            if len(seq_T.names_N) == 1:
                if isinstance(seq_T, TrimmedMappedSequence):
                    spec_map_seq_count += 1
                    spec_map_read_count += seq_T.read_count
                    if not PROFILE_ABSENT_3PRIME_TERMINUS:
                        if start_T_in_N == 0 and stop_T_in_N == length_N:
                            absent_3prime_terminus_seq_count += 1
                            absent_3prime_terminus_read_count += seq_T.read_count
                    spec_read_xtra_5prime_count += seq_T.read_with_xtra_5prime_count
                    for string_5prime_extension, read_count in seq_T.long_5prime_extension_dict.items():
                        spec_long_5prime_extension_dict[string_5prime_extension] += read_count
                else:
                    for string_3prime_terminus, read_count in seq_T.read_3prime_terminus_count_dict.items():
                        spec_read_3prime_terminus_count_dict[string_3prime_terminus] += read_count

                    if not isinstance(seq_T, TrimmedTruncatedProfileSequence):
                        spec_read_xtra_5prime_count += seq_T.read_with_xtra_5prime_count
                        for string_5prime_extension, read_count in seq_T.long_5prime_extension_dict.items():
                            spec_long_5prime_extension_dict[string_5prime_extension] += read_count
            else:
                if isinstance(seq_T, TrimmedMappedSequence):
                    nonspec_map_seq_count += 1
                    nonspec_map_read_count += seq_T.read_count
                    nonspec_read_xtra_5prime_count += seq_T.read_with_xtra_5prime_count
                else:
                    for string_3prime_terminus, read_count in seq_T.read_3prime_terminus_count_dict.items():
                        nonspec_read_3prime_terminus_count_dict[string_3prime_terminus] += read_count

                    if not isinstance(seq_T, TrimmedTruncatedProfileSequence):
                        nonspec_read_xtra_5prime_count += seq_T.read_with_xtra_5prime_count
                        for string_5prime_extension, read_count in seq_T.long_5prime_extension_dict.items():
                            nonspec_long_5prime_extension_dict[string_5prime_extension] += read_count

        self.spec_map_seq_count = spec_map_seq_count
        self.nonspec_map_seq_count = nonspec_map_seq_count
        self.spec_map_read_count = spec_map_read_count
        self.nonspec_map_read_count = nonspec_map_read_count
        self.absent_3prime_terminus_seq_count = absent_3prime_terminus_seq_count
        self.absent_3prime_terminus_read_count = absent_3prime_terminus_read_count
        self.spec_read_xtra_5prime_count = spec_read_xtra_5prime_count
        self.nonspec_read_xtra_5prime_count = nonspec_read_xtra_5prime_count
        self.spec_long_5prime_extension_dict = spec_long_5prime_extension_dict
        self.nonspec_long_5prime_extension_dict = nonspec_long_5prime_extension_dict
        self.spec_read_3prime_terminus_count_dict = spec_read_3prime_terminus_count_dict
        self.nonspec_read_3prime_terminus_count_dict = nonspec_read_3prime_terminus_count_dict


class NormalizedTruncatedProfileSequence(NormalizedSequence):
    """This object is formed exclusively from `TrimmedTruncatedProfileSequence` objects."""
    # Instances of this class are called `Nc`.

    __slots__ = (
        'feature_starts',
        'feature_stops',
        'contains_anticodon',
        'trunc_profile_index',
        'spec_read_3prime_terminus_count_dict',
        'nonspec_read_3prime_terminus_count_dict'
    )

    def __init__(self, seqs_T, starts_T_in_N=None, stops_T_in_N=None):
        super().__init__(seqs_T, starts_T_in_N=None, stops_T_in_N=None, skip_init=False)

        represent_T = seqs_T[0]
        self.feature_starts = represent_T.feature_starts
        self.feature_stops = represent_T.feature_stops
        self.contains_anticodon = represent_T.contains_anticodon
        self.trunc_profile_index = represent_T.trunc_profile_index

        spec_read_3prime_terminus_count_dict = defaultdict(int)
        nonspec_read_3prime_terminus_count_dict = defaultdict(int)
        for seq_T in seqs_T:
            if len(seq_T.names_N) == 1:
                for string_3prime_terminus, read_count in seq_T.read_3prime_terminus_count_dict.items():
                    spec_read_3prime_terminus_count_dict[string_3prime_terminus] += read_count
            else:
                for string_3prime_terminus, read_count in seq_T.read_3prime_terminus_count_dict.items():
                    nonspec_read_3prime_terminus_count_dict[string_3prime_terminus] += read_count
        self.spec_read_3prime_terminus_count_dict = spec_read_3prime_terminus_count_dict
        self.nonspec_read_3prime_terminus_count_dict = nonspec_read_3prime_terminus_count_dict


class NormalizedIndelSequence(NormalizedSequence):
    """This object is formed exclusively from `TrimmedIndelSequence` objects."""
    # Instances of this class are called `Ni`.

    __slots__ = (
        'insert_starts_Ni',
        'insert_starts_M',
        'insert_lengths',
        'del_starts_Ni',
        'del_starts_M',
        'del_lengths',
        'contains_anticodon',
        'spec_insert_covs',
        'nonspec_insert_covs',
        'spec_del_covs',
        'nonspec_del_covs',
        'spec_read_xtra_5prime_count',
        'nonspec_read_xtra_5prime_count',
        'spec_map_read_count',
        'nonspec_map_read_count',
        'spec_long_5prime_extension_dict',
        'nonspec_long_5prime_extension_dict',
        'spec_read_3prime_terminus_count_dict',
        'nonspec_read_3prime_terminus_count_dict'
    )

    def __init__(self,
                 string,
                 seqs_T,
                 starts_T_in_N,
                 name_M,
                 insert_starts_Ni,
                 insert_starts_M,
                 insert_lengths,
                 del_starts_Ni,
                 del_starts_M,
                 del_lengths,
                 contains_anticodon):
        self.name = name = seqs_T[0].name
        for seq_Ti in seqs_T:
            seq_Ti.names_N.append(name)
        self.names_T = [seq_T.name for seq_T in seqs_T]
        self.string = string
        self.starts_T_in_N = starts_T_in_N
        self.names_M = [name_M]
        self.insert_starts_Ni = insert_starts_Ni
        self.insert_starts_M = insert_starts_M
        self.insert_lengths = insert_lengths
        self.del_starts_Ni = del_starts_Ni
        self.del_starts_M = del_starts_M
        self.del_lengths = del_lengths
        self.contains_anticodon = contains_anticodon


    def init(self, seqs_T, dict_Ui):
        """Set attributes for a "finalized" set of input Ti."""
        self.names_T = tuple(self.names_T)
        self.starts_T_in_N = tuple(self.starts_T_in_N)
        spec_read_count = 0
        nonspec_read_count = 0
        length_Ni = len(self.string)
        spec_covs = np.zeros(length_Ni, dtype=int)
        nonspec_covs = np.zeros(length_Ni, dtype=int)
        insert_starts_Ni = self.insert_starts_Ni
        insert_stops_Ni = [insert_start + insert_length for insert_start, insert_length in zip(insert_starts_Ni, self.insert_lengths)]
        spec_insert_covs = np.zeros(len(insert_starts_Ni), dtype=int)
        nonspec_insert_covs = np.zeros(len(insert_starts_Ni), dtype=int)
        del_starts_Ni = self.del_starts_Ni
        spec_del_covs = np.zeros(len(del_starts_Ni), dtype=int)
        nonspec_del_covs = np.zeros(len(del_starts_Ni), dtype=int)
        spec_read_xtra_5prime_count = 0
        nonspec_read_xtra_5prime_count = 0
        spec_map_read_count = 0
        nonspec_map_read_count = 0
        spec_long_5prime_extension_dict = defaultdict(int)
        nonspec_long_5prime_extension_dict = defaultdict(int)
        spec_read_3prime_terminus_count_dict = defaultdict(int)
        nonspec_read_3prime_terminus_count_dict = defaultdict(int)
        for seq_Ti, start_Ti_in_Ni in zip(seqs_T, self.starts_T_in_N):
            read_count = seq_Ti.read_count
            orig_U = dict_Ui[seq_Ti.names_U[0]].orig_U
            if len(seq_Ti.names_N) == 1:
                spec_read_count += read_count
                if seq_Ti.read_with_xtra_5prime_count:
                    spec_read_xtra_5prime_count += seq_Ti.read_with_xtra_5prime_count
                if isinstance(orig_U, UniqueMappedSequence):
                    spec_map_read_count += 1
                covs = spec_covs
                insert_covs = spec_insert_covs
                del_covs = spec_del_covs
                long_5prime_extension_dict = spec_long_5prime_extension_dict
                read_3prime_terminus_count_dict = spec_read_3prime_terminus_count_dict
            else:
                nonspec_read_count += read_count
                if seq_Ti.read_with_xtra_5prime_count:
                    nonspec_read_xtra_5prime_count += seq_Ti.read_with_xtra_5prime_count
                if isinstance(orig_U, UniqueMappedSequence):
                    nonspec_map_read_count += 1
                covs = nonspec_covs
                insert_covs = nonspec_insert_covs
                del_covs = nonspec_del_covs
                long_5prime_extension_dict = nonspec_long_5prime_extension_dict
                read_3prime_terminus_count_dict = nonspec_read_3prime_terminus_count_dict

            stop_Ti_in_Ni = start_Ti_in_Ni + len(seq_Ti.string)
            covs[start_Ti_in_Ni: stop_Ti_in_Ni] += read_count

            insert_index = 0
            for insert_start, insert_stop in zip(insert_starts_Ni, insert_stops_Ni):
                if start_Ti_in_Ni <= insert_start and stop_Ti_in_Ni >= insert_stop:
                    insert_covs[insert_index] += read_count
                insert_index += 1

            for del_index, del_start in enumerate(del_starts_Ni):
                if start_Ti_in_Ni <= del_start and stop_Ti_in_Ni > del_start + 1:
                    del_covs[del_index] += read_count

            for string_5prime, extension_read_count in seq_Ti.long_5prime_extension_dict.items():
                long_5prime_extension_dict[string_5prime] += extension_read_count

            for string_3prime, terminus_read_count in seq_Ti.read_3prime_terminus_count_dict.items():
                read_3prime_terminus_count_dict[string_3prime] += terminus_read_count
        self.spec_read_count = spec_read_count
        self.nonspec_read_count = nonspec_read_count
        self.spec_covs = spec_covs
        self.nonspec_covs = nonspec_covs
        self.mean_spec_cov = spec_covs.mean()
        self.mean_nonspec_cov = nonspec_covs.mean()
        self.spec_insert_covs = spec_insert_covs
        self.nonspec_insert_covs = nonspec_insert_covs
        self.spec_del_covs = spec_del_covs
        self.nonspec_del_covs = nonspec_del_covs
        self.spec_read_xtra_5prime_count = spec_read_xtra_5prime_count
        self.nonspec_read_xtra_5prime_count = nonspec_read_xtra_5prime_count
        self.spec_map_read_count = spec_map_read_count
        self.nonspec_map_read_count = nonspec_map_read_count
        self.spec_long_5prime_extension_dict = spec_long_5prime_extension_dict
        self.nonspec_long_5prime_extension_dict = nonspec_long_5prime_extension_dict
        self.spec_read_3prime_terminus_count_dict = spec_read_3prime_terminus_count_dict
        self.nonspec_read_3prime_terminus_count_dict = nonspec_read_3prime_terminus_count_dict


class ModifiedSequence(object):
    """This object represents a tRNA sequence with sites of predicted potential mod-induced subs
    and, optionally, indels.

    The `anvi-trnaseq` workflow aggregates similar Nf. The aggregations are decomposed into clusters
    of Nf distinguished by potential mod-induced subs (3-4 different nts at ≥1 aligned positions). M
    is instantiated with the list of Nf, with the first Nf being longest or tied for longest.
    Corresponding lists of Nf sub positions in M are required. The workflow later finds Ni and adds
    them to M. Lists of indel positions in M and indel lengths are needed. The `init` method of M is
    called to calculate coverages and other information from nts, subs, and indels, if present.

    The workflow requires that Ni be assigned to 1 M for simplicity's sake. If the same Ni can arise
    from indels in multiple M, then Ni is disregarded.

    M EXAMPLE:
    Consider E. coli tRNA-Ala-GGC-1-1, with detected mods at positions 17 and 46. As seen in the N
    example, the first sequence is the N with unmutated nucleotides. The next set of N include
    possible mod-induced subs. The next set of N are Ni with insertions, and the last set are Ni
    with deletions.

                     |                              |
    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA

    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG A TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAG C TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC C GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA
              CTCAGC G GGGAGAGCGCTTGCATGGCATGCAAGAG G TCAGCGGTTCGATCCCGCTTAGCTCCA
                                    CATGGCATGCAAGAG T TCAGCGGTTCGATCCCGCTTAGCTCCA

    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGAGAG TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGAGGG TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC A GGGAGAGCGCTTGCATGGCATGCAAGAGGAATCAGCGGTTCGATCCCGCTTAGCTCCA

    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGAG - TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGA- - TCAGCGGTTCGATCCCGCTTAGCTCCA
    GGGGCTATAGCTCAGC T GGGAGAGCGCTTGCATGGCATGCAAGA- G TCAGCGGTTCGATCCCGCTTAGCTCCA
    """

    __slots__ = (
        'names_Nb',
        'sub_positions',
        'name',
        'length',
        'names_Ni',
        'starts_Ni_in_M',
        'insert_starts',
        'insert_strings',
        'del_starts',
        'del_lengths',
        'spec_read_count',
        'nonspec_read_count',
        'spec_map_read_count',
        'nonspec_map_read_count',
        'spec_read_xtra_5prime_count',
        'nonspec_read_xtra_5prime_count',
        'spec_long_5prime_extension_dict',
        'nonspec_long_5prime_extension_dict',
        'spec_read_3prime_terminus_count_dict',
        'nonspec_read_3prime_terminus_count_dict',
        'spec_covs',
        'nonspec_covs',
        'mean_spec_cov',
        'mean_nonspec_cov',
        'spec_sub_covs',
        'nonspec_sub_covs',
        'spec_insert_covs',
        'nonspec_insert_covs',
        'spec_del_covs',
        'nonspec_del_covs',
        'consensus_string'
    )

    def __init__(self, seqs_Nb, sub_positions):
        self.names_Nb = [seq_Nb.name for seq_Nb in seqs_Nb]
        self.sub_positions = sub_positions
        self.name = name = self.names_Nb[0]
        for seq_Nb in seqs_Nb:
            seq_Nb.names_M.append(name)
        self.length = len(seqs_Nb[0].string)
        self.names_Ni = tuple()
        self.starts_Ni_in_M = tuple()
        self.insert_starts = None
        self.insert_strings = None
        self.del_starts = None
        self.del_lengths = None
        self.spec_read_count = None
        self.nonspec_read_count = None
        self.spec_map_read_count = None
        self.nonspec_map_read_count = None
        self.spec_read_xtra_5prime_count = None
        self.nonspec_read_xtra_5prime_count = None
        self.spec_long_5prime_extension_dict = None
        self.nonspec_long_5prime_extension_dict = None
        self.spec_read_3prime_terminus_count_dict = None
        self.nonspec_read_3prime_terminus_count_dict = None
        self.spec_covs = None
        self.nonspec_covs = None
        self.mean_spec_cov = None
        self.mean_nonspec_cov = None
        self.spec_sub_covs = None
        self.nonspec_sub_covs = None
        self.spec_insert_covs = None
        self.nonspec_insert_covs = None
        self.spec_del_covs = None
        self.nonspec_del_covs = None
        self.consensus_string = None


    def init(self, seqs_Nb, seqs_Ni):
        """Set attributes for a "finalized" set of input Ns and Ni objects."""
        self.names_Nb = tuple(self.names_Nb)

        spec_read_count = 0
        nonspec_read_count = 0
        spec_map_read_count = 0
        nonspec_map_read_count = 0
        spec_read_xtra_5prime_count = 0
        nonspec_read_xtra_5prime_count = 0
        spec_long_5prime_extension_dict = defaultdict(int)
        nonspec_long_5prime_extension_dict = defaultdict(int)
        spec_read_3prime_terminus_count_dict = defaultdict(int)
        nonspec_read_3prime_terminus_count_dict = defaultdict(int)
        length_M = self.length
        spec_covs_M = np.zeros(length_M, dtype=int)
        nonspec_covs_M = np.zeros(length_M, dtype=int)


        # Find the read counts of different types of seqs composing M.
        for seq_N in seqs_Nb + seqs_Ni:
            spec_read_count += seq_N.spec_read_count
            nonspec_read_count += seq_N.nonspec_read_count
            spec_map_read_count += seq_N.spec_map_read_count
            nonspec_map_read_count += seq_N.nonspec_map_read_count
            spec_read_xtra_5prime_count += seq_N.spec_read_xtra_5prime_count
            nonspec_read_xtra_5prime_count += seq_N.nonspec_read_xtra_5prime_count
            for string_5prime, read_count in seq_N.spec_long_5prime_extension_dict.items():
                spec_long_5prime_extension_dict[string_5prime] += read_count
            for string_5prime, read_count in seq_N.nonspec_long_5prime_extension_dict.items():
                nonspec_long_5prime_extension_dict[string_5prime] += read_count
            for string_3prime, read_count in seq_N.spec_read_3prime_terminus_count_dict.items():
                spec_read_3prime_terminus_count_dict[string_3prime] += read_count
            for string_3prime, read_count in seq_N.nonspec_read_3prime_terminus_count_dict.items():
                nonspec_read_3prime_terminus_count_dict[string_3prime] += read_count


        # Find the covs of subs in Nb.
        sub_positions = self.sub_positions
        spec_sub_covs = np.zeros((len(sub_positions), len(UNAMBIG_NTS)), dtype=int)
        nonspec_sub_covs = np.zeros((len(sub_positions), len(UNAMBIG_NTS)), dtype=int)

        reverse_sub_positions = sub_positions[::-1]
        for seq_Nb in seqs_Nb:
            # N are aligned with the 3' end of M.
            start_Nb_in_M = length_M - len(seq_Nb.string)
            spec_covs_M[start_Nb_in_M: ] += seq_Nb.spec_covs
            nonspec_covs_M[start_Nb_in_M: ] += seq_Nb.nonspec_covs

            # Find the covs of subs in Nb. Loop through subs from the 3' end of Nb and M.
            string_Nb = seq_Nb.string
            spec_covs_Nb = seq_Nb.spec_covs
            nonspec_covs_Nb = seq_Nb.nonspec_covs
            for sub_index, sub_pos_M in enumerate(reverse_sub_positions, 1):
                if sub_pos_M < start_Nb_in_M:
                    # Remaining subs are 5' of Nb.
                    break
                sub_pos_Nb = sub_pos_M - start_Nb_in_M
                nt_col = NT_INT_DICT[string_Nb[sub_pos_Nb]] - 1
                spec_sub_covs[-sub_index, nt_col] += spec_covs_Nb[sub_pos_Nb]
                nonspec_sub_covs[-sub_index, nt_col] += nonspec_covs_Nb[sub_pos_Nb]


        # Find the covs of M positions, subs, and indels.
        insert_covs_dict = {}
        del_covs_dict = {}
        for seq_Ni, start_Ni_in_M in zip(seqs_Ni, self.starts_Ni_in_M):
            # Find insertion covs in Ni.
            for insert_Ni_start, insert_M_start, insert_length, spec_insert_cov, nonspec_insert_cov in zip(seq_Ni.insert_starts_Ni,
                                                                                                           seq_Ni.insert_starts_M,
                                                                                                           seq_Ni.insert_lengths,
                                                                                                           seq_Ni.spec_insert_covs,
                                                                                                           seq_Ni.nonspec_insert_covs):
                try:
                    insert_covs = insert_covs_dict[(insert_M_start, seq_Ni.string[insert_Ni_start: insert_Ni_start + insert_length])]
                    insert_covs[0] += spec_insert_cov
                    insert_covs[1] += nonspec_insert_cov
                except KeyError:
                    insert_covs_dict[(insert_M_start, seq_Ni.string[insert_Ni_start: insert_Ni_start + insert_length])] = np.array((spec_insert_cov, nonspec_insert_cov))

            # Find deletion covs in Ni.
            for del_M_start, del_length, spec_del_cov, nonspec_del_cov in zip(seq_Ni.del_starts_M,
                                                                              seq_Ni.del_lengths,
                                                                              seq_Ni.spec_del_covs,
                                                                              seq_Ni.nonspec_del_covs):
                try:
                    del_covs = del_covs_dict[(del_M_start, del_length)]
                    del_covs[0] += spec_del_cov
                    del_covs[1] += nonspec_del_cov
                except KeyError:
                    del_covs_dict[(del_M_start, del_length)] = np.array((spec_del_cov, nonspec_del_cov))

            # Find nt and sub covs in Ni. Loop through each nt of Ni to account for the positions of
            # indels, a more complex process than that used for Nb.
            string_Ni = seq_Ni.string
            spec_covs_Ni = seq_Ni.spec_covs
            nonspec_covs_Ni = seq_Ni.nonspec_covs

            # Make an iterator of the positions of inserted nts in Ni.
            iter_Ni_insert_positions = iter(chain.from_iterable(
                [range(insert_start, insert_start + insert_length)
                 for insert_start, insert_length in zip(seq_Ni.insert_starts_Ni, seq_Ni.insert_lengths)]))
            try:
                next_Ni_insert_pos = next(iter_Ni_insert_positions)
            except StopIteration:
                next_Ni_insert_pos = -1

            # Make iterators of the positions and lengths of deletions in Ni.
            iter_Ni_del_positions = iter([del_Ni_start + 1 for del_Ni_start in seq_Ni.del_starts_Ni])
            iter_Ni_del_lengths = iter(seq_Ni.del_lengths)
            try:
                next_Ni_del_pos = next(iter_Ni_del_positions)
                next_Ni_del_length = next(iter_Ni_del_lengths)
            except StopIteration:
                next_Ni_del_pos = -1
                next_Ni_del_length = -1

            nt_pos_M = start_Ni_in_M
            sub_index = bisect_left(sub_positions, start_Ni_in_M)
            iter_M_sub_positions = iter(sub_positions[sub_index: ])
            try:
                next_M_sub_pos = next(iter_M_sub_positions)
            except StopIteration:
                next_M_sub_pos = MAXSIZE
            for nt_pos_Ni, nt in enumerate(string_Ni):
                if nt_pos_Ni == next_Ni_insert_pos:
                    try:
                        next_Ni_insert_pos = next(iter_Ni_insert_positions)
                    except StopIteration:
                        next_Ni_insert_pos = -1
                    # Since insertions are not nts in M, do not increment the nt position in M.
                    continue

                if nt_pos_Ni == next_Ni_del_pos:
                    # Increment the nt position in M by the size of the del.
                    nt_pos_M += next_Ni_del_length
                    try:
                        next_Ni_del_pos = next(iter_Ni_del_positions)
                        next_Ni_del_length = next(iter_Ni_del_lengths)
                    except StopIteration:
                        next_Ni_del_pos = -1
                        next_Ni_del_length = -1
                    while nt_pos_M > next_M_sub_pos:
                        # The del contained a sub. Find the position of the next sub after the del.
                        try:
                            next_M_sub_pos = next(iter_M_sub_positions)
                            sub_index += 1
                        except StopIteration:
                            break
                    nt_pos_M += 1
                    continue

                # To make it to this point, the position is not an indel.
                if nt_pos_M == next_M_sub_pos:
                    spec_cov_Ni = spec_covs_Ni[nt_pos_Ni]
                    nonspec_cov_Ni = nonspec_covs_Ni[nt_pos_Ni]
                    spec_covs_M[nt_pos_M] += spec_cov_Ni
                    nonspec_covs_M[nt_pos_M] += nonspec_cov_Ni
                    nt_col = NT_INT_DICT[string_Ni[nt_pos_Ni]] - 1
                    spec_sub_covs[sub_index, nt_col] += spec_cov_Ni
                    nonspec_sub_covs[sub_index, nt_col] += nonspec_cov_Ni
                    try:
                        next_M_sub_pos = next(iter_M_sub_positions)
                        sub_index += 1
                    except StopIteration:
                        next_M_sub_pos = -1
                else:
                    spec_covs_M[nt_pos_M] += seq_Ni.spec_covs[nt_pos_Ni]
                    nonspec_covs_M[nt_pos_M] += seq_Ni.nonspec_covs[nt_pos_Ni]

                nt_pos_M += 1

        # Record the positions, lengths, and coverages of insertions in M.
        insert_starts = []
        insert_strings = []
        spec_insert_covs = []
        nonspec_insert_covs = []
        for insert_config, insert_covs in sorted(insert_covs_dict.items()):
            insert_starts.append(insert_config[0])
            insert_strings.append(insert_config[1])
            spec_insert_covs.append(insert_covs[0])
            nonspec_insert_covs.append(insert_covs[1])

        # Record the positions, lengths, and coverages of deletions in M.
        del_starts = []
        del_lengths = []
        spec_del_covs = []
        nonspec_del_covs = []
        for del_config, del_covs in sorted(del_covs_dict.items()):
            del_starts.append(del_config[0])
            del_lengths.append(del_config[1])
            spec_del_covs.append(del_covs[0])
            nonspec_del_covs.append(del_covs[1])

        # Set a consensus seq using the nts with the highest specific cov at each sub position.
        consensus_string = seqs_Nb[0].string
        for sub_pos, sub_nt_covs in zip(sub_positions, spec_sub_covs):
            nt_int = sub_nt_covs.argmax() + 1
            consensus_string = consensus_string[: sub_pos] + INT_NT_DICT[nt_int] + consensus_string[sub_pos + 1: ]

        self.insert_starts = tuple(insert_starts)
        self.insert_strings = tuple(insert_strings)
        self.del_starts = tuple(del_starts)
        self.del_lengths = tuple(del_lengths)
        self.spec_read_count = spec_read_count
        self.nonspec_read_count = nonspec_read_count
        self.spec_map_read_count = spec_map_read_count
        self.nonspec_map_read_count = nonspec_map_read_count
        self.spec_read_xtra_5prime_count = spec_read_xtra_5prime_count
        self.nonspec_read_xtra_5prime_count = nonspec_read_xtra_5prime_count
        self.spec_long_5prime_extension_dict = spec_long_5prime_extension_dict
        self.nonspec_long_5prime_extension_dict = nonspec_long_5prime_extension_dict
        self.spec_read_3prime_terminus_count_dict = spec_read_3prime_terminus_count_dict
        self.nonspec_read_3prime_terminus_count_dict = nonspec_read_3prime_terminus_count_dict
        self.spec_covs = tuple(spec_covs_M)
        self.nonspec_covs = tuple(nonspec_covs_M)
        self.mean_spec_cov = spec_covs_M.mean()
        self.mean_nonspec_cov = nonspec_covs_M.mean()
        self.spec_sub_covs = spec_sub_covs
        self.nonspec_sub_covs = nonspec_sub_covs
        self.spec_insert_covs = tuple(spec_insert_covs)
        self.nonspec_insert_covs = tuple(nonspec_insert_covs)
        self.spec_del_covs = tuple(spec_del_covs)
        self.nonspec_del_covs = tuple(nonspec_del_covs)
        self.consensus_string = consensus_string


class TRNASeqDataset(object):
    """Processes reads from a tRNA-seq library. `bin/anvi-trnaseq` is the client."""

    RELATIVE_ANTICODON_LOOP_INDEX = TRNA_FEATURE_NAMES.index('anticodon_loop') - len(TRNA_FEATURE_NAMES) + 1

    # Column headers for supplementary tables written to text files
    UNIQ_NONTRNA_HEADER = [
        "name",
        "read_count",
        "truncated_profile_index",
        "sequence"
    ]
    TRIMMED_ENDS_HEADER = [
        "name",
        "unique_name",
        "fiveprime_sequence",
        "threeprime_sequence",
        "read_count"
    ]

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # Argument group 1A: MANDATORY
        self.input_fasta_path = A('trnaseq_fasta')
        self.sample_id = A('sample_name')
        self.out_dir = os.path.abspath(A('output_dir')) if A('output_dir') else None
        self.checkpoint_dir = os.path.join(self.out_dir, "CHECKPOINT")
        get_checkpoint_subdir = partial(os.path.join, self.checkpoint_dir)
        self.checkpoint_subdir_dict = {checkpoint: get_checkpoint_subdir(checkpoint.upper()) for checkpoint in constants.TRNASEQ_CHECKPOINTS}

        # Argument group 1B: EXTRAS
        self.treatment = A('treatment')
        self.overwrite_out_dest = A('overwrite_output_destinations')
        self.descrip_path = os.path.abspath(A('description')) if A('description') else None

        # Argument group 1C: ADVANCED
        self.write_checkpoints = A('write_checkpoints')
        self.load_checkpoint = A('load_checkpoint')
        self.feature_param_path = os.path.abspath(A('feature_param_file')) if A('feature_param_file') else None
        self.param_3prime_termini = A('threeprime_termini')
        global MIN_LENGTH_LONG_5PRIME_EXTENSION
        MIN_LENGTH_LONG_5PRIME_EXTENSION = A('min_length_long_fiveprime')
        self.min_trna_frag_size = A('min_trna_fragment_size')
        agglom_max_mismatch_freq = A('agglomeration_max_mismatch_freq')
        self.agglom_max_mismatch_freq = round(agglom_max_mismatch_freq * 100) / 100
        self.skip_indel_profiling = A('skip_INDEL_profiling')
        max_indel_freq = A('max_indel_freq')
        self.max_indel_freq = round(max_indel_freq * 100) / 100
        self.left_indel_buffer = A('left_indel_buffer')
        self.right_indel_buffer = A('right_indel_buffer')

        # Argument group 1D: PERFORMANCE
        self.num_threads = A('num_threads')
        self.skip_fasta_check = A('skip_fasta_check')
        self.profiling_chunk_size = A('profiling_chunk_size')
        self.alignment_target_chunk_size = A('alignment_target_chunk_size')

        if not self.input_fasta_path:
            raise ConfigError("Please specify the path to a FASTA file of tRNA-seq reads using `--fasta-file` or `-f`.")
        if not self.sample_id:
            raise ConfigError("Please provide a sample name using `--sample-name` or `-S`.")
        if not self.out_dir:
            raise ConfigError("Please provide an output directory using `--output-dir` or `-o`.")

        self.descrip = None

        get_out_file_path = partial(os.path.join, self.out_dir)

        self.trnaseq_db_path = get_out_file_path(self.sample_id + "-TRNASEQ.db")

        self.analysis_summary_path = get_out_file_path(self.sample_id + "-ANALYSIS_SUMMARY.txt")

        # Supplementary text file paths produced by DEBUG flag
        self.path_Un_supplement = get_out_file_path(self.sample_id + "-UNIQUED_NONTRNA.txt")
        self.path_Tf_ends = get_out_file_path(self.sample_id + "-TRIMMED_ENDS.txt")
        self.consol_seqs_with_inconsis_profiles_path = get_out_file_path(self.sample_id + "-CONSOLIDATED_SEQS_WITH_INCONSISTENT_PROFILES.txt")
        self.count_consol_Tf = None

        # The identification of sequences as tRNA occurs through different means. By the time of the
        # first "profile" checkpoint, only Uf are recognized as tRNA. Further processing before the
        # second "normalize" checkpoint can show some Ut to be tRNA; mapping can reveal some Un to
        # be tRNA by the third "mapping" checkpoint. The changing classification of sequences over
        # the workflow means that the contents of the dictionaries and intermediate files storing U,
        # T, and N also change. It is important to note that object names do not correspond to these
        # classifications. For example, all Ut will be in `uniq_trunc_dict` before the "profile"
        # checkpoint, but some Ut can move to `uniq_trna_dict` after being confirmed as tRNA by the
        # "normalize" checkpoint. The underlying sequences encapsulated in Ut objects stay in these
        # objects despite being recognized as tRNA.

        # Not every dict of seq objects changes between checkpoints, yet all existing dicts are
        # written at every checkpoint. The alternative would be to reduce the number of intermediate
        # files stored at later checkpoints by relying on any unchanged intermediate files from
        # earlier checkpoints.

        self.intermed_file_path_dict = {}
        get_intermed_file_path = partial(os.path.join, self.checkpoint_subdir_dict['profile'])
        self.intermed_file_path_dict['profile'] = {
            'dict_Uf': get_intermed_file_path("Uf.pkl"),
            'dict_Uc_nontrna': get_intermed_file_path("Uc_nontRNA.pkl"),
            'dict_Un': get_intermed_file_path("Un.pkl")
        }
        get_intermed_file_path = partial(os.path.join, self.checkpoint_subdir_dict['normalize'])
        self.intermed_file_path_dict['normalize'] = {
            'dict_Uf': get_intermed_file_path("Uf.pkl"),
            'dict_Us': get_intermed_file_path("Us.pkl"),
            'dict_Uc_trna': get_intermed_file_path("Uc_tRNA.pkl"),
            'dict_Uc_nontrna': get_intermed_file_path("Uc_nontRNA.pkl"),
            'dict_Un': get_intermed_file_path("Un.pkl"),
            'dict_Tf': get_intermed_file_path("Tf.pkl"),
            'dict_Tc_trna': get_intermed_file_path("Tc_tRNA.pkl"),
            'dict_Tc_nontrna': get_intermed_file_path("Tc_nontrna.pkl"),
            'dict_Nf': get_intermed_file_path("Nf.pkl"),
            'dict_Nc': get_intermed_file_path("Nc.pkl")
        }
        get_intermed_file_path = partial(os.path.join, self.checkpoint_subdir_dict['map_fragments'])
        self.intermed_file_path_dict['map_fragments'] = {
            'dict_Uf': get_intermed_file_path("Uf.pkl"),
            'dict_Us': get_intermed_file_path("Us.pkl"),
            'dict_Uc_trna': get_intermed_file_path("Uc_tRNA.pkl"),
            'dict_Um': get_intermed_file_path("Um.pkl"),
            'dict_Uc_nontrna': get_intermed_file_path("Uc_nontRNA.pkl"),
            'dict_Un': get_intermed_file_path("Un.pkl"),
            'dict_Tf': get_intermed_file_path("Tf.pkl"),
            'dict_Tc_trna': get_intermed_file_path("Tc_tRNA.pkl"),
            'dict_Tc_nontrna': get_intermed_file_path("Tc_nontrna.pkl"),
            'dict_Tm': get_intermed_file_path("Tm.pkl"),
            'dict_Nf': get_intermed_file_path("Nf.pkl"),
            'dict_Nc': get_intermed_file_path("Nc.pkl")
        }
        get_intermed_file_path = partial(os.path.join, self.checkpoint_subdir_dict['substitutions'])
        self.intermed_file_path_dict['substitutions'] = {
            'dict_Uf': get_intermed_file_path("Uf.pkl"),
            'dict_Us': get_intermed_file_path("Us.pkl"),
            'dict_Uc_trna': get_intermed_file_path("Uc_tRNA.pkl"),
            'dict_Um': get_intermed_file_path("Um.pkl"),
            'dict_Uc_nontrna': get_intermed_file_path("Uc_nontRNA.pkl"),
            'dict_Un': get_intermed_file_path("Un.pkl"),
            'dict_Tf': get_intermed_file_path("Tf.pkl"),
            'dict_Tc_trna': get_intermed_file_path("Tc_tRNA.pkl"),
            'dict_Tc_nontrna': get_intermed_file_path("Tc_nontrna.pkl"),
            'dict_Tm': get_intermed_file_path("Tm.pkl"),
            'dict_Nf': get_intermed_file_path("Nf.pkl"),
            'dict_Nc': get_intermed_file_path("Nc.pkl"),
            'dict_M': get_intermed_file_path("M.pkl")
        }
        get_intermed_file_path = partial(os.path.join, self.checkpoint_subdir_dict['indels'])
        self.intermed_file_path_dict['indels'] = {
            'dict_Uf': get_intermed_file_path("Uf.pkl"),
            'dict_Us': get_intermed_file_path("Us.pkl"),
            'dict_Uc_trna': get_intermed_file_path("Uc_tRNA.pkl"),
            'dict_Um': get_intermed_file_path("Um.pkl"),
            'dict_Ui': get_intermed_file_path("Ui.pkl"),
            'dict_Uc_nontrna': get_intermed_file_path("Uc_nontRNA.pkl"),
            'dict_Un': get_intermed_file_path("Un.pkl"),
            'dict_Tf': get_intermed_file_path("Tf.pkl"),
            'dict_Tc_trna': get_intermed_file_path("Tc_tRNA.pkl"),
            'dict_Tc_nontrna': get_intermed_file_path("Tc_nontrna.pkl"),
            'dict_Tm': get_intermed_file_path("Tm.pkl"),
            'dict_Ti': get_intermed_file_path("Ti.pkl"),
            'dict_Nf': get_intermed_file_path("Nf.pkl"),
            'dict_Nc': get_intermed_file_path("Nc.pkl"),
            'dict_Ni': get_intermed_file_path("Ni.pkl"),
            'dict_M': get_intermed_file_path("M.pkl")
        }

        self.dict_Uf = {}
        self.dict_Uc_trna = {}
        self.dict_Uc_nontrna = {}
        self.dict_Un = {}
        self.dict_Us = {}
        self.dict_Um = {}
        self.dict_Ui = {}
        self.dict_Tf = {}
        self.dict_Tc_trna = {}
        self.dict_Tc_nontrna = {}
        self.dict_Tm = {}
        self.dict_Ti = {}
        self.dict_Nf = {}
        self.dict_Nc = {}
        self.dict_Ni = {}
        self.dict_M = {}


    def process(self):
        """The entry method of TRNASeqDataset, called by `anvi-trnaseq`."""
        total_time_start = time.time()

        self.sanity_check()
        if '' in self.parsed_3prime_termini:
            global PROFILE_ABSENT_3PRIME_TERMINUS
            PROFILE_ABSENT_3PRIME_TERMINUS = True

        load_checkpoint = self.load_checkpoint
        if not load_checkpoint:
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
            self.report_indel_analysis_parameters()

            # Profile each (unique) read for tRNA features.
            self.profile_trna()
            self.report_profile_stats()

            if self.write_checkpoints:
                self.write_checkpoint_files('profile')
        elif load_checkpoint == 'profile':
            self.load_checkpoint_files('profile')
            self.report_fragment_mapping_parameters()
            self.report_substitution_analysis_parameters()
            self.report_indel_analysis_parameters()


        if (load_checkpoint == 'profile'
            or not load_checkpoint):
            # Do the steps between the "profile" and "normalize" checkpoints.

            # Trim 5' and 3' ends of Uf, forming Tf.
            self.trim_trna_ends()
            # Trim 3' ends of Uc, forming Tc.
            self.trim_truncated_profile_ends()
            self.report_trim_stats()

            # Consolidate 3' fragments of longer Tf, forming Nf.
            self.threeprime_dereplicate_profiled_trna()

            # Recover Tc as tRNA by comparing to Nf.
            self.threeprime_dereplicate_truncated_sequences()
            self.report_threeprime_dereplication_statistics()

            if self.write_checkpoints:
                self.write_checkpoint_files('normalize')
        elif load_checkpoint == 'normalize':
            self.load_checkpoint_files('normalize')
            self.report_fragment_mapping_parameters()
            self.report_substitution_analysis_parameters()
            self.report_indel_analysis_parameters()


        if (load_checkpoint == 'normalize'
            or load_checkpoint == 'profile'
            or not load_checkpoint):
            # Do the steps between the "normalize" and "map_fragments" checkpoints.

            if not PROFILE_ABSENT_3PRIME_TERMINUS:
                # Recover 3' tRNA sequences lacking a 3' terminus.
                self.threeprime_dereplicate_sequences_without_terminus()

            # Map fragments derived from the interior and 5' end of tRNA.
            self.map_fragments()

            # Finalize Nf now that all T found through various means have been added to them.
            self.progress.new("Finalizing normalized tRNA sequences")
            self.progress.update("...")
            for seq_Nf in self.dict_Nf.values():
                seq_Nf.init([getattr(self, 'dict_' + category_T)[name_T]
                             for category_T, name_T in zip(seq_Nf.categories_T, seq_Nf.names_T)])
            self.progress.end()
            self.report_mapping_statistics()
            self.report_initialized_normalized_sequence_coverage_statistics()

            if self.write_checkpoints:
                self.write_checkpoint_files('map_fragments')
        elif load_checkpoint == 'map_fragments':
            self.load_checkpoint_files('map_fragments')
            self.report_substitution_analysis_parameters()
            self.report_indel_analysis_parameters()


        if (load_checkpoint == 'map_fragments'
            or load_checkpoint == 'normalize'
            or load_checkpoint == 'profile'
            or not load_checkpoint):
            # Do the steps between the "map_fragments" and "substitutions" checkpoints.

            # Find modified nucleotides, grouping normalized sequences into modified sequences.
            self.find_substitutions()

            self.report_sub_stats()

            if self.write_checkpoints:
                self.write_checkpoint_files('substitutions')
        elif load_checkpoint == 'substitutions':
            self.load_checkpoint_files('substitutions')
            self.report_indel_analysis_parameters()


        if (load_checkpoint == 'substitutions'
            or load_checkpoint == 'map_fragments'
            or load_checkpoint == 'normalize'
            or load_checkpoint == 'profile'
            or not load_checkpoint):
            # Do the steps between the "substitutions" and "indels" checkpoint.

            if self.skip_indel_profiling:
                trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
                set_meta_value = trnaseq_db.db.set_meta_value
                set_meta_value('count_Nqf_with_indels', 0)
                set_meta_value('count_Nc_with_indels', 0)
                trnaseq_db.disconnect()
            else:
                self.find_indels()

            # "Finalize" M now that all Nf and Ni have been added to them.
            self.progress.new("Finalizing modified tRNA sequences")
            self.progress.update("...")
            dict_Nf = self.dict_Nf
            dict_Ni = self.dict_Ni
            for seq_M in self.dict_M.values():
                seq_M.init([dict_Nf[name_Nb] for name_Nb in seq_M.names_Nb],
                           [dict_Ni[name_Ni] for name_Ni in seq_M.names_Ni])
            self.progress.end()

            self.report_M_stats()

            if self.write_checkpoints:
                self.write_checkpoint_files('indels')
        elif load_checkpoint == 'indels':
            self.load_checkpoint_files('indels')


        self.report_stats()

        self.write_feature_table()
        self.write_unconserved_table()
        self.write_unpaired_table()
        self.write_sequences_table()
        self.write_trimmed_table()
        self.write_normalized_table()
        self.write_modified_table()

        # Write supplementary text files.
        self.write_nontrna_supplement()
        self.write_Tf_ends_supplement()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Total time elapsed (min)",
                                          time.time() - total_time_start,
                                          is_time_value=True))
            # Write an empty line to separate this run from any subsequent run starting from a
            # checkpoint writing to the same summary file.
            f.write("\n")


    def sanity_check(self):
        """Check `anvi-trnaseq` arguments."""
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

        self.parsed_3prime_termini = self.threeprime_termini_sanity_check()
        trnaidentifier.TRNAFeatureParameterizer.set_threeprime_termini(self.parsed_3prime_termini)
        # The following variable is only used as part of a heuristic in `find_indels`.
        self.max_length_3prime_terminus = max([len(t) for t in self.parsed_3prime_termini])

        self.run.info("Input FASTA file", self.input_fasta_path, nl_after=1)

        if not self.skip_fasta_check and not self.load_checkpoint:
            self.progress.new("Checking input FASTA defline format")
            self.progress.update("...")

            check_fasta_id_formatting(self.input_fasta_path)

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
        checkpoint_dir = self.checkpoint_dir
        load_checkpoint = self.load_checkpoint
        intermed_file_path_dict = self.intermed_file_path_dict

        checkpoint_subdir = os.path.join(checkpoint_dir, load_checkpoint.upper())
        if not os.path.exists(checkpoint_subdir):
            raise ConfigError("Intermediate files needed for running `anvi-trnaseq` with `--load-checkpoint` "
                              f"should be located in {checkpoint_subdir}, but this directory path does not exist. "
                              "You should probably run `anvi-trnaseq` from the beginning without `--load-checkpoint`. "
                              "To generate necessary intermediate files for future use of `--load-checkpoint`, use the flag `--write-checkpoints`.")

        missing_intermed_files = []
        for intermed_file_path in intermed_file_path_dict[load_checkpoint].values():
            if not os.path.exists(intermed_file_path):
                missing_intermed_files.append(intermed_file_path)
        if missing_intermed_files:
            raise ConfigError(f"Intermediate files needed for running `anvi-trnaseq` "
                              f"with `--load-checkpoint {load_checkpoint}` are missing: {', '.join(missing_intermed_files)}. "
                              "You should probably run `anvi-trnaseq` from the beginning without `--load-checkpoint`. "
                              "To generate necessary intermediate files for future use of `--load-checkpoint`, use the flag `--write-checkpoints`.")


    def threeprime_termini_sanity_check(self):
        """Check validity of provided tRNA 3' termini, returning a list of terminus strings."""
        valid_3prime_termini = []
        invalid_3prime_termini = []
        for terminus_3prime in self.param_3prime_termini.split(','):
            if terminus_3prime == '_':
                valid_3prime_termini.append('')
                continue

            for nt in terminus_3prime:
                if nt not in ALL_NTS:
                    invalid_3prime_termini.append(terminus_3prime)
                    break
            valid_3prime_termini.append(terminus_3prime)

        if invalid_3prime_termini:
            raise ConfigError(f"3' termini can consist of A, C, G, T, and N (any nucleotide) "
                              "or the discriminator nucleotide with no extension, symbolized by a single underscore, \"_\". "
                              f"The following invalid 3' sequence parameterizations were provided: {', '.join(invalid_3prime_termini)}")

        return valid_3prime_termini


    def create_trnaseq_database(self):
        """Create an empty tRNA-seq database."""
        meta_values = {'sample_id': self.sample_id,
                       'treatment': self.treatment,
                       'description': self.descrip if self.descrip else '_No description is provided_',
                       'INDELs_profiled': not self.skip_indel_profiling}
        dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True).create(meta_values)
        self.run.info("New tRNA-seq db", self.trnaseq_db_path, nl_after=1)


    def report_profiling_parameters(self):
        """Add profiling parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        set_meta_value = trnaseq_db.db.set_meta_value
        parameterizer = trnaidentifier.TRNAFeatureParameterizer()
        for param_tuple in parameterizer.list_accessible_param_tuples():
            set_meta_value(param_tuple[0], param_tuple[1])
        set_meta_value('min_length_long_5prime_extension', MIN_LENGTH_LONG_5PRIME_EXTENSION)
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
            # Use `param_3prime_termini` rather than `parsed_3prime_termini` here to write "_"
            # rather than "" for an absent 3' terminus.
            f.write(get_summary_line("Allowed 3' termini", ",".join(self.param_3prime_termini)))
            f.write(get_summary_line("Min length of \"long\" 5' extension", MIN_LENGTH_LONG_5PRIME_EXTENSION))


    def report_fragment_mapping_parameters(self):
        """Add fragment mapping parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db.set_meta_value('min_map_trna_fragment_size', self.min_trna_frag_size)
        trnaseq_db.disconnect()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Min length of mapped tRNA fragment", self.min_trna_frag_size))


    def report_substitution_analysis_parameters(self):
        """Add modification-induced substitution analysis parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        trnaseq_db.db.set_meta_value('agglomeration_max_mismatch_freq', self.agglom_max_mismatch_freq)
        trnaseq_db.disconnect()

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("Agglomeration max mismatch frequency", self.agglom_max_mismatch_freq))


    def report_indel_analysis_parameters(self):
        """Add modification-induced indel analysis parameters to the database."""
        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        set_meta_value = trnaseq_db.db.set_meta_value
        set_meta_value('max_indel_freq', self.max_indel_freq)
        set_meta_value('left_indel_buffer', self.left_indel_buffer)
        set_meta_value('right_indel_buffer', self.right_indel_buffer)
        trnaseq_db.disconnect()

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("INDELs profiled", not self.skip_indel_profiling))
            f.write(get_summary_line("Max indel frequency", self.max_indel_freq))
            f.write(get_summary_line("Left indel buffer", self.left_indel_buffer))
            f.write(get_summary_line("Right indel buffer", self.right_indel_buffer))


    def get_summary_line(self, label, value, is_time_value=False, padding=68):
        """Return a string formatted to be written to the summary statistics file."""
        # Report elapsed time in seconds in minutes.
        if is_time_value:
            value = "%.2f" % round(value / 60, 2)
        return '%s%s\t%s\n' % (label, ' ' + '.' * (padding - len(label)), value)


    def profile_trna(self):
        """Profile tRNA features in reads, finding Uf, Uc, and Un."""
        uniq_read_infos = self.unique_reads()

        start_time = time.time()

        pid = "Profiling tRNA features in unique reads"
        self.progress.new(pid)
        self.progress.update("...")

        # Count the number of reads and unique reads that have been added to the multiprocessing
        # input queue.
        total_read_count = 0
        total_uniq_count = len(uniq_read_infos)

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()
        profiler = trnaidentifier.Profiler()
        processes = [multiprocessing.Process(target=profile_worker, args=(input_queue, output_queue, profiler))
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
        dict_Uf = self.dict_Uf
        dict_Uc_nontrna = self.dict_Uc_nontrna
        dict_Un = self.dict_Un
        pp_total_uniq_count = pp(len(uniq_read_infos))
        while fetched_profile_count < total_uniq_count:
            self.progress.update_pid(pid)
            self.progress.update(f"{pp(input_count + 1)}-{pp(interval_stop)}/{pp_total_uniq_count}")

            while input_count < interval_stop:
                for uniq_read_info in uniq_read_infos[interval_start: interval_stop]:
                    input_queue.put(uniq_read_info)
                    total_read_count += uniq_read_info[2]
                    input_count += 1

            while fetched_profile_count < interval_stop:
                profile, read_count = output_queue.get()
                fetched_profile_count += 1

                name = profile.name
                if profile.is_predicted_trna:
                    dict_Uf[name] = UniqueFullProfileSequence(profile.input_seq, name, read_count, profile)
                else:
                    if profile.trunc_profile_index:
                        dict_Uc_nontrna[name] = UniqueTruncatedProfileSequence(profile.input_seq, name, read_count, profile)
                    else:
                        dict_Un[name] = UniqueSequence(profile.input_seq, name, read_count)

            interval_start = interval_stop
            interval_stop += profiling_chunk_size if interval_stop + profiling_chunk_size < total_uniq_count else total_uniq_count - interval_stop

        for p in processes:
            p.terminate()
            p.join()

        # Profiled seqs were added to the output queue as they were processed, so sort by name.
        self.dict_Uf = {name: seq for name, seq in sorted(dict_Uf.items())}
        self.dict_Uc_nontrna = {name: seq for name, seq in sorted(dict_Uc_nontrna.items())}
        self.dict_Un = {name: seq for name, seq in sorted(dict_Un.items())}

        get_summary_line = self.get_summary_line
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("Time elapsed profiling tRNA (min)", time.time() - start_time, is_time_value=True))
            f.write(get_summary_line("Reads processed", total_read_count))
            f.write(get_summary_line("Unique seqs processed", total_uniq_count))

        self.progress.end()

        self.run.info("Reads processed", total_read_count, mc='green')
        self.run.info("Unique seqs processed", total_uniq_count, mc='green')

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        set_meta_value = trnaseq_db.db.set_meta_value
        set_meta_value('input_reads', total_read_count)
        set_meta_value('input_U', total_uniq_count)
        trnaseq_db.disconnect()


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
        uniq_read_infos = []
        for cluster in Dereplicator(names, seqs).full_length_dereplicate():
            uniq_read_infos.append((cluster.member_seqs[0], cluster.member_names[0], len(cluster.member_names)))

        self.progress.end()
        return uniq_read_infos


    def report_profile_stats(self):
        """Report to terminal stats on Uf, Uc, and Un immediately after profiling."""
        seq_count_Uf = len(self.dict_Uf)
        read_count_Uf = 0
        seq_anticodon_count_Uf = 0
        read_anticodon_count_Uf = 0
        seq_complete_count_Uf = 0
        read_complete_count_Uf = 0
        max_reads_Uf = 0
        mean_read_3prime_length_Uf = 0
        seq_short_5prime_count_Uf = 0
        read_short_5prime_count_Uf = 0
        seq_long_5prime_count_Uf = 0
        read_long_5prime_count_Uf = 0
        mean_read_5prime_length_Uf = 0
        mean_seq_profiled_freq_Uf = 0
        mean_read_profiled_freq_Uf = 0
        mean_seq_unconserved_freq_Uf = 0
        mean_read_unconserved_freq_Uf = 0
        mean_seq_unpaired_freq_Uf = 0
        mean_read_unpaired_freq_Uf = 0
        mean_seq_extrap_freq_Uf = 0
        mean_read_extrap_freq_Uf = 0
        for seq_Uf in self.dict_Uf.values():
            read_count = seq_Uf.read_count
            read_count_Uf += read_count
            if seq_Uf.anticodon_string:
                seq_anticodon_count_Uf += 1
                read_anticodon_count_Uf += read_count
                if seq_Uf.has_complete_feature_set:
                    seq_complete_count_Uf += 1
                    read_complete_count_Uf += read_count
            if read_count > max_reads_Uf:
                max_reads_Uf = read_count

            mean_read_3prime_length_Uf += read_count * seq_Uf.length_3prime_terminus
            if seq_Uf.xtra_5prime_length:
                if seq_Uf.xtra_5prime_length < MIN_LENGTH_LONG_5PRIME_EXTENSION:
                    seq_short_5prime_count_Uf += 1
                    read_short_5prime_count_Uf += read_count
                else:
                    seq_long_5prime_count_Uf += 1
                    read_long_5prime_count_Uf += read_count
                    mean_read_5prime_length_Uf += read_count * seq_Uf.xtra_5prime_length

            profiled_freq = seq_Uf.profiled_seq_length / len(seq_Uf.string)
            mean_seq_profiled_freq_Uf += profiled_freq
            mean_read_profiled_freq_Uf += read_count * profiled_freq
            unconserved_freq = seq_Uf.num_unconserved / seq_Uf.profiled_seq_length
            mean_seq_unconserved_freq_Uf += unconserved_freq
            mean_read_unconserved_freq_Uf += read_count * unconserved_freq
            unpaired_freq = seq_Uf.num_unpaired / (seq_Uf.num_paired + seq_Uf.num_unpaired)
            mean_seq_unpaired_freq_Uf += unpaired_freq
            mean_read_unpaired_freq_Uf += read_count * unpaired_freq
            mean_seq_extrap_freq_Uf += seq_Uf.num_extrap_5prime_nts
            mean_read_extrap_freq_Uf += read_count * seq_Uf.num_extrap_5prime_nts
        mean_reads_Uf = read_count_Uf / seq_count_Uf
        mean_read_3prime_length_Uf /= read_count_Uf
        mean_read_5prime_length_Uf /= read_long_5prime_count_Uf
        mean_seq_profiled_freq_Uf /= seq_count_Uf
        mean_read_profiled_freq_Uf /= read_count_Uf
        mean_seq_unconserved_freq_Uf /= seq_count_Uf
        mean_read_unconserved_freq_Uf /= read_count_Uf
        mean_seq_unpaired_freq_Uf /= seq_count_Uf
        mean_read_unpaired_freq_Uf /= read_count_Uf
        mean_seq_extrap_freq_Uf /= seq_count_Uf
        mean_read_extrap_freq_Uf /= read_count_Uf

        seq_count_Uc = len(self.dict_Uc_nontrna)
        read_count_Uc = 0
        seq_anticodon_count_Uc = 0
        read_anticodon_count_Uc = 0
        max_reads_Uc = 0
        mean_read_3prime_length_Uc = 0
        mean_seq_profiled_freq_Uc = 0
        mean_read_profiled_freq_Uc = 0
        mean_seq_unconserved_freq_Uc = 0
        mean_read_unconserved_freq_Uc = 0
        mean_seq_unpaired_freq_Uc = 0
        mean_read_unpaired_freq_Uc = 0
        for seq_Uc in self.dict_Uc_nontrna.values():
            read_count = seq_Uc.read_count
            read_count_Uc += read_count
            if seq_Uc.anticodon_string:
                seq_anticodon_count_Uc += 1
                read_anticodon_count_Uc += read_count
            if read_count > max_reads_Uc:
                max_reads_Uc = read_count

            mean_read_3prime_length_Uc += read_count * seq_Uc.length_3prime_terminus

            profiled_freq = seq_Uc.profiled_seq_length / len(seq_Uc.string)
            mean_seq_profiled_freq_Uc += profiled_freq
            mean_read_profiled_freq_Uc += read_count * profiled_freq
            unconserved_freq = seq_Uc.num_unconserved / seq_Uc.profiled_seq_length
            mean_seq_unconserved_freq_Uc += unconserved_freq
            mean_read_unconserved_freq_Uc += read_count * unconserved_freq
            unpaired_freq = seq_Uc.num_unpaired / (seq_Uc.num_paired + seq_Uc.num_unpaired)
            mean_seq_unpaired_freq_Uc += unpaired_freq
            mean_read_unpaired_freq_Uc += read_count * unpaired_freq
        mean_reads_Uc = read_count_Uc / seq_count_Uc
        mean_read_3prime_length_Uc /= read_count_Uc
        mean_seq_profiled_freq_Uc /= seq_count_Uc
        mean_read_profiled_freq_Uc /= read_count_Uc
        mean_seq_unconserved_freq_Uc /= seq_count_Uc
        mean_read_unconserved_freq_Uc /= read_count_Uc
        mean_seq_unpaired_freq_Uc /= seq_count_Uc
        mean_read_unpaired_freq_Uc /= read_count_Uc

        seq_count_Un = len(self.dict_Un)
        read_count_Un = 0
        max_reads_Un = 0
        for seq_Un in self.dict_Un.values():
            read_count = seq_Un.read_count
            read_count_Un += read_count
            if read_count > max_reads_Un:
                max_reads_Un = read_count
        mean_reads_Un = read_count_Un / seq_count_Un

        warning = self.run.warning
        info_single = self.run.info_single
        info = self.run.info

        warning(None, "PROFILING RESULTS", lc='green', nl_before=1)
        info_single("subject to change -- see summary output file for final results")

        warning(None, "Unique seq counts", lc='cyan')
        info("tRNA profile", seq_count_Uf)
        info("Truncated tRNA profile", seq_count_Uc)
        info("No tRNA profile", seq_count_Un)

        warning(None, "Read counts", lc='cyan')
        info("tRNA profile", read_count_Uf)
        info("Truncated tRNA profile", read_count_Uc)
        info("No tRNA profile", read_count_Un)

        warning(None, "Unique seqs with tRNA profile", lc='cyan')
        info("Count with anticodon", seq_anticodon_count_Uf)
        info("Count with complete feature set", seq_complete_count_Uf)
        info("Mean reads per seq", round(mean_reads_Uf, 1))
        info("Max reads per seq", max_reads_Uf)
        info(f"Count with 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION - 1} extra 5' nts", seq_short_5prime_count_Uf)
        info(f"Count with ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", seq_long_5prime_count_Uf)
        info("Mean profiled nt freq", round(mean_seq_profiled_freq_Uf, 3))
        info("Mean freq of unconserved in profiled nts", round(mean_seq_unconserved_freq_Uf, 4))
        info("Mean freq of unpaired in stem nts", round(mean_seq_unpaired_freq_Uf, 4))
        info("Mean extrapolated 5' nt freq", round(mean_seq_extrap_freq_Uf, 3))

        warning(None, "Unique seqs with truncated tRNA profile", lc='cyan')
        info("Count with anticodon", seq_anticodon_count_Uc)
        info("Mean reads per seq", round(mean_reads_Uc, 1))
        info("Max reads per seq", max_reads_Uc)
        info("Mean profiled nt freq", round(mean_seq_profiled_freq_Uc, 3))
        info("Mean freq of unconserved in profiled nts", round(mean_seq_unconserved_freq_Uc, 4))
        info("Mean freq of unpaired in stem nts", round(mean_seq_unpaired_freq_Uc, 4))

        warning(None, "Unique seqs with no tRNA profile", lc='cyan')
        info("Mean reads per seq", round(mean_reads_Un, 1))
        info("Max reads per seq", max_reads_Un)

        warning(None, "Reads with tRNA profile", lc='cyan')
        info("Count with anticodon", read_anticodon_count_Uf)
        info("Count with complete feature set", read_complete_count_Uf)
        info("Mean length 3' terminus", round(mean_read_3prime_length_Uf, 1))
        info(f"Count with 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION - 1} extra 5' nts", read_short_5prime_count_Uf)
        info(f"Count with ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", read_long_5prime_count_Uf)
        info(f"Mean length ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} nt extension", round(mean_read_5prime_length_Uf, 1))
        info("Mean profiled nt freq", round(mean_read_profiled_freq_Uf, 3))
        info("Mean freq of unconserved in profiled nts", round(mean_read_unconserved_freq_Uf, 4))
        info("Mean freq of unpaired in stem nts", round(mean_read_unpaired_freq_Uf, 4))
        info("Mean extrapolated 5' nt freq", round(mean_read_extrap_freq_Uf, 3))

        warning(None, "Reads with truncated tRNA profile", lc='cyan')
        info("Spans anticodon", read_anticodon_count_Uc)
        info("Mean length 3' terminus", round(mean_read_3prime_length_Uc, 1))
        info("Mean profiled nt freq", round(mean_read_profiled_freq_Uc, 3))
        info("Mean freq of unconserved in profiled nts", round(mean_read_unconserved_freq_Uc, 4))
        info("Mean freq of unpaired in stem nts", round(mean_read_unpaired_freq_Uc, 4), nl_after=2 if self.write_checkpoints else 1)


    def trim_trna_ends(self):
        """Trim any nts 5' of the acceptor stem and 3' of the discriminator from Uf, forming Tf."""
        start_time = time.time()
        self.progress.new("Trimming the 3' and 5' ends of profiled tRNA")
        self.progress.update("...")

        seqs_Tf = self.get_trimmed_seqs([seq_Uf for seq_Uf in self.dict_Uf.values()], TrimmedFullProfileSequence)
        dict_Tf = self.dict_Tf
        for seq_Tf in sorted(seqs_Tf, key=lambda seq_Tf: seq_Tf.name):
            dict_Tf[seq_Tf.name] = seq_Tf

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed trimming profiled tRNA (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def get_trimmed_seqs(self, seqs_U, class_T):
        """Find Tf or Tc from Uf or Uc, respectively."""
        names = [seq_U.name for seq_U in seqs_U]
        if class_T == TrimmedFullProfileSequence:
            strings_T = [seq_U.string[seq_U.xtra_5prime_length: len(seq_U.string) - seq_U.length_3prime_terminus] for seq_U in seqs_U]
        elif class_T == TrimmedTruncatedProfileSequence:
            strings_T = [seq_U.string[: len(seq_U.string) - seq_U.length_3prime_terminus] for seq_U in seqs_U]

        clusters = Dereplicator(names, strings_T, extras=seqs_U).full_length_dereplicate()

        seqs_T = [class_T(cluster.member_seqs[0], cluster.member_extras) for cluster in clusters]
        return seqs_T


    def trim_truncated_profile_ends(self):
        """Trim any nts 3' of the discriminator from Uc, forming Tc."""
        start_time = time.time()
        self.progress.new("Trimming the 3' ends of seqs with truncated tRNA profiles")
        self.progress.update("...")

        seqs_Tc = self.get_trimmed_seqs([seq_Uc for seq_Uc in self.dict_Uc_nontrna.values()], TrimmedTruncatedProfileSequence)
        dict_Tc_nontrna = self.dict_Tc_nontrna
        for seq_Tc in sorted(seqs_Tc, key=lambda seq_Tc: seq_Tc.name):
            dict_Tc_nontrna[seq_Tc.name] = seq_Tc

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed trimming seqs with truncated feature profile (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def report_trim_stats(self):
        """Report to terminal stats to the terminal on Tf and Tc immediately after trimming steps."""
        count_Tf = len(self.dict_Tf)
        anticodon_count_Tf = 0
        complete_count_Tf = 0
        mean_uniq_seqs_Tf = len(self.dict_Uf) / count_Tf
        single_count_Tf = 0
        mean_reads_Tf = 0
        max_reads_Tf = 0
        long_5prime_count_Tf = 0
        for seq_Tf in self.dict_Tf.values():
            if seq_Tf.contains_anticodon:
                anticodon_count_Tf += 1
                if seq_Tf.has_complete_feature_set:
                    complete_count_Tf += 1
            if len(seq_Tf.names_U) == 1:
                single_count_Tf += 1
            read_count = seq_Tf.read_count
            mean_reads_Tf += read_count
            if read_count > max_reads_Tf:
                max_reads_Tf = read_count
            if seq_Tf.long_5prime_extension_dict:
                long_5prime_count_Tf += 1
        mean_reads_Tf /= count_Tf

        count_Tc = len(self.dict_Tc_nontrna)
        anticodon_count_Tc = 0
        mean_uniq_seqs_Tc = len(self.dict_Uc_nontrna) / count_Tc
        single_uniq_seq_count_Tc = 0
        mean_reads_Tc = 0
        max_reads_Tc = 0
        for seq_Tc in self.dict_Tc_nontrna.values():
            if seq_Tc.contains_anticodon:
                anticodon_count_Tc += 1
            if len(seq_Tc.names_U) == 1:
                single_uniq_seq_count_Tc += 1
            read_count = seq_Tc.read_count
            mean_reads_Tc += read_count
            if read_count > max_reads_Tc:
                max_reads_Tc = read_count
        mean_reads_Tc /= count_Tc

        warning = self.run.warning
        info_single = self.run.info_single
        info = self.run.info

        warning(None, "TRIMMING RESULTS", lc='green', nl_before=1)
        info_single("subject to change -- see summary output file for final results")

        warning(None, "Trimmed seq counts", lc='cyan')
        info("tRNA profile", count_Tf)
        info("Truncated tRNA profile", count_Tc)

        warning(None, "Trimmed seqs with tRNA profile", lc='cyan')
        info("Count with anticodon", anticodon_count_Tf)
        info("Count with complete feature set", complete_count_Tf)
        info("Mean unique seqs per seq", round(mean_uniq_seqs_Tf, 1))
        info("Count with single unique seq", single_count_Tf)
        info("Mean reads per seq", round(mean_reads_Tf, 1))
        info("Max reads per seq", max_reads_Tf)
        info(f"Count with ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", long_5prime_count_Tf)

        warning(None, "Trimmed seqs with truncated tRNA profile", lc='cyan')
        info("Count with anticodon", anticodon_count_Tc)
        info("Mean unique seqs per seq", round(mean_uniq_seqs_Tc, 1))
        info("Count with single unique seq", single_uniq_seq_count_Tc)
        info("Mean reads per seq", round(mean_reads_Tc, 1))
        info("Max reads per seq", max_reads_Tc, nl_after=1)


    def threeprime_dereplicate_profiled_trna(self):
        """Dereplicate Tf from the 3' end of longer Tf.

        EXAMPLE:
        Nf (Tf 1): TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Tf 2     :                       AATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Tf 3     :                                     GCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        """
        start_time = time.time()
        pid = "Dereplicating trimmed tRNA seqs from the 3' end"
        self.progress.new(pid)
        self.progress.update("...")

        # Prefix dereplicate Tf from the 3' end.
        names = []
        reverse_strings = []
        seqs_Tf = []
        for name, seq_Tf in self.dict_Tf.items():
            names.append(name)
            reverse_strings.append(seq_Tf.string[::-1])
            seqs_Tf.append(seq_Tf)
        clusters = Dereplicator(names, reverse_strings, extras=seqs_Tf).prefix_dereplicate()

        # Profiling may have found multiple Tf that would here be 3'-dereplicated as having
        # complete, but different, feature profiles. Assume that the shortest "completely profiled"
        # Tf in the cluster has the correct profile (with an exception explained with an asterisk
        # below). Reclassify discrepant 5' nts from Uf in longer Tf as extra 5' nts. Transfer the
        # profile information from the representative Uf of the shortest Tf to Uf of longer Tf,
        # replacing these longer Uf with Us objects. Produce a new Tf object from these Uf and Us.

        # Similarly, the longest Tf in the cluster may have an erroneous "incomplete profile." If
        # there is a shorter Tf in the cluster with a complete profile, then any longer Tf with an
        # incomplete profile can be reevaluated in the same manner.

        # This step, which likely affects a tiny number of Tf, helps reduce the number of extra,
        # wrong nts at the 5' end of downstream seed sequences. In debug mode, consolidated Tf are
        # written to a file.

        # Nota bene: A complete profile may be extrapolated at the 5' end of the 5' strand of the
        # acceptor stem. By default, only 1 nt may be extrapolated in the stem when all other nts
        # form base pairs. So it is hard to see how an extrapolated complete profile is more likely
        # to be inaccurate than the removed profiles of longer sequences in the cluster.

        # Do not check for feature-by-feature agreement among clustered profiles here, as 3' tRNA
        # fragments occasionally have some incorrect feature positions due to a paucity of sequence
        # information available in fragment profiling. It is conceivable that the seed sequence of
        # the cluster is assigned a wrong complete or incomplete profile while all of the shorter
        # sequences in the cluster are assigned correct incomplete profiles -- presumably a rare
        # inaccuracy that goes unchecked.
        self.progress.update_pid(pid)
        self.progress.update("Inspecting normalized seq clusters")

        dict_Nf = self.dict_Nf
        dict_Uf = self.dict_Uf

        # It is possible that trimmed sequences from multiple clusters can consolidate.
        dict_consol_Tf = {}

        # This dict is for an edge case explained below.
        dict_Tf_His = {}

        unrepresent_complete_profile_Tf_names = []
        for cluster in clusters:
            # Skip initialization of Nf objects, as additional Tf members are later added to the
            # objects after dereplicating Tc and mapping Tm.
            seqs_Tf = cluster.member_extras
            represent_name = seqs_Tf[0].name
            if len(seqs_Tf) == 1:
                dict_Nf[represent_name] = NormalizedFullProfileSequence(seqs_Tf)
                continue

            # Check that there are no shorter Tf in the cluster with a "complete profile".
            complete_profile_indices = []
            for index_Tf, seq_Tf in enumerate(seqs_Tf):
                if seq_Tf.has_complete_feature_set:
                    complete_profile_indices.append(index_Tf)

            if not complete_profile_indices:
                dict_Nf[represent_name] = NormalizedFullProfileSequence(seqs_Tf)
                continue

            if len(complete_profile_indices) == 1:
                # * Most frequently, the longest sequence has a complete profile. An edge case can
                # occur when the longest sequence does not have a complete profile while a shorter
                # sequence does. The shorter sequence may appear as a member of another cluster and
                # become the representative sequence due to consolidation. Prevent that
                # consolidation from happening, forming a normalized sequence as per usual from the
                # trimmed sequences, with the longest being representative. *
                dict_Nf[represent_name] = NormalizedFullProfileSequence(seqs_Tf)
                if complete_profile_indices != [0]:
                    unrepresent_complete_profile_Tf_names.append(seqs_Tf[complete_profile_indices[0]].name)
                continue

            # Reaching this point means that there are multiple Tf with "complete profiles" in the
            # cluster.

            # If the two shortest Tf with complete feature profiles differ by the
            # post-transcriptionally added 5'-G of tRNA-His, then they should both be maintained as
            # separate Tf.
            if seqs_Tf[complete_profile_indices[-2]].has_His_G:
                if seqs_Tf[complete_profile_indices[-1]].string == seqs_Tf[complete_profile_indices[-2]].string[1: ]:

                    if len(complete_profile_indices) == 2:
                        assert complete_profile_indices[-2] == 0
                        dict_Tf_His[seqs_Tf[1].name] = seqs_Tf
                        continue

                    # Perhaps more than two Tf in the cluster have "complete" profiles, though this
                    # has not been checked. In this case, consolidate the Tf, retaining the profile
                    # of the shortest.

            short_seq_Tf = seqs_Tf[complete_profile_indices[-1]]
            if short_seq_Tf.name in dict_consol_Tf:
                # Tf from multiple clusters consolidate with the same shorter Tf with a complete
                # profile.

                # Some of the longer Tf with rejected "complete" profiles occur in multiple
                # clusters. Clearly, the longest seed Tf, which is also being consolidated, must
                # differ between the clusters. Tf in the clusters that are shorter than the Tf with
                # the selected complete profile must be the same in each cluster.
                long_seqs_Tf = dict_consol_Tf[short_seq_Tf.name]['long_seqs_Tf']
                encountered_long_Tf_names = [seq_Tf.name for seq_Tf in long_seqs_Tf]
                for long_seq_Tf in seqs_Tf[: complete_profile_indices[-1]]:
                    if long_seq_Tf.name not in encountered_long_Tf_names:
                        long_seqs_Tf.append(long_seq_Tf)
            else:
                # This is the first time the shortest Tf with a complete profile has been processed
                # from a cluster.
                dict_consol_Tf[short_seq_Tf.name] = {'short_seq_Tf': short_seq_Tf,
                                                     'long_seqs_Tf': seqs_Tf[: complete_profile_indices[-1]],
                                                     'Nf_members': seqs_Tf[complete_profile_indices[-1] + 1: ]}

        # Consider the following edge case. One cluster had two Tf with complete profiles, Tf1 and
        # Tf2, so the two were consolidated, forming Nf1.
        # Tf1:        GGTGGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # Tf2 == Nf1:    GGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # Another cluster had two Tf, Tf3 and the same Tf2, only differing by a supposed tRNA-His
        # 5'-G.
        # Tf3:          GGGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # Tf2:           GGGAGAATTCCCGAGTGGCCAAGGGGGGCAGACTGTGTATCTGTTGCGTTTCGCTTCGATGGTTCGAATCCATCTTCTCCCA
        # Rather than creating another Nf, Nf2, seeded by Tf3, consolidate Tf from the two clusters.
        # This avoids producing two Nf, one of which is a 3' subsequence of the other.
        for name, seqs_Tf in dict_Tf_His.items():
            if name in dict_consol_Tf:
                dict_consol_Tf[name]['long_seqs_Tf'].append(seqs_Tf[0])
            else:
                dict_Nf[seqs_Tf[0].name] = NormalizedFullProfileSequence(seqs_Tf)

        consol_seqs_with_inconsis_profiles_file = open(self.consol_seqs_with_inconsis_profiles_path, 'w')
        consol_seqs_with_inconsis_profiles_file.write("Index\tTrimmed (0) or Unique (1)\tSequence\n")

        count_consol_Tf = 0
        dict_Us = self.dict_Us
        for subdict_consol_Tf in dict_consol_Tf.values():
            if subdict_consol_Tf['short_seq_Tf'].name in unrepresent_complete_profile_Tf_names:
                long_seqs_Tf = subdict_consol_Tf['long_seqs_Tf']
                dict_Nf[long_seqs_Tf[0].name] = NormalizedFullProfileSequence(long_seqs_Tf + [subdict_consol_Tf['short_seq_Tf']])
                continue

            count_consol_Tf += 1
            consol_seq_Tf = self.consolidate_trimmed_sequences(subdict_consol_Tf['short_seq_Tf'], subdict_consol_Tf['long_seqs_Tf'])
            dict_Nf[consol_seq_Tf.name] = NormalizedFullProfileSequence([consol_seq_Tf] + subdict_consol_Tf['Nf_members'])

            # Report consolidated Tf with different complete feature profiles.
            first_field = str(count_consol_Tf) + "\t"
            for seq_Tf in subdict_consol_Tf['long_seqs_Tf']:
                consol_seqs_with_inconsis_profiles_file.write(first_field)
                consol_seqs_with_inconsis_profiles_file.write("0\t")
                consol_seqs_with_inconsis_profiles_file.write(seq_Tf.string + "\n")
                for name_Us in seq_Tf.names_U:
                    consol_seqs_with_inconsis_profiles_file.write(first_field)
                    consol_seqs_with_inconsis_profiles_file.write("1\t")
                    consol_seqs_with_inconsis_profiles_file.write(dict_Us[name_Us].string + "\n")
            seq_Tf = subdict_consol_Tf['short_seq_Tf']
            consol_seqs_with_inconsis_profiles_file.write(first_field)
            consol_seqs_with_inconsis_profiles_file.write("0\t")
            consol_seqs_with_inconsis_profiles_file.write(seq_Tf.string + "\n")
            for name_Uf in seq_Tf.names_U:
                consol_seqs_with_inconsis_profiles_file.write(first_field)
                consol_seqs_with_inconsis_profiles_file.write("1\t")
                consol_seqs_with_inconsis_profiles_file.write(dict_Uf[name_Uf].string + "\n")

        self.count_consol_Tf = count_consol_Tf
        consol_seqs_with_inconsis_profiles_file.close()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed 3'-dereplicating trimmed profiled seqs", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def consolidate_trimmed_sequences(self, short_seq_Tf, long_seqs_Tf):
        """Consolidate longer Tf with a shorter Tf by pooling all of their Uf -- changing their
        profile information -- and generating a new Tf object."""
        short_seq_Tf_string = short_seq_Tf.string
        dict_Uf = self.dict_Uf
        short_seq_Uf = dict_Uf[short_seq_Tf.names_U[0]]

        # To transfer profile information from the representative short Uf to longer Uf, determine
        # where features are relative to the 3' end of the short Uf, which is found in the other Uf.
        replacement_info_dict = {'string_T': short_seq_Tf.string}
        feature_index_adjustment = -short_seq_Uf.xtra_5prime_length - len(short_seq_Tf_string)
        feature_starts_from_T_3prime = []
        for feature_start in short_seq_Uf.feature_starts:
            if isinstance(feature_start, int):
                feature_starts_from_T_3prime.append(feature_start + feature_index_adjustment)
            else:
                feature_starts_from_T_3prime.append(tuple([strand_start + feature_index_adjustment for strand_start in feature_start]))
        replacement_info_dict['feature_starts_from_T_3prime'] = feature_starts_from_T_3prime
        feature_stops_from_T_3prime = []
        for feature_stop in short_seq_Uf.feature_stops:
            if isinstance(feature_stop, int):
                feature_stops_from_T_3prime.append(feature_stop + feature_index_adjustment)
            else:
                feature_stops_from_T_3prime.append(tuple([strand_stop + feature_index_adjustment for strand_stop in feature_stop]))
        replacement_info_dict['feature_stops_from_T_3prime'] = feature_stops_from_T_3prime
        replacement_info_dict['has_His_G'] = short_seq_Uf.has_His_G
        replacement_info_dict['alpha_start_from_T_3prime'] = None if short_seq_Uf.alpha_start is None else short_seq_Uf.alpha_start + feature_index_adjustment
        replacement_info_dict['alpha_stop_from_T_3prime'] = None if short_seq_Uf.alpha_stop is None else short_seq_Uf.alpha_stop + feature_index_adjustment
        replacement_info_dict['beta_start_from_T_3prime'] = None if short_seq_Uf.beta_start is None else short_seq_Uf.beta_start + feature_index_adjustment
        replacement_info_dict['beta_stop_from_T_3prime'] = None if short_seq_Uf.beta_stop is None else short_seq_Uf.beta_stop + feature_index_adjustment
        replacement_info_dict['anticodon_string'] = short_seq_Uf.anticodon_string
        replacement_info_dict['anticodon_aa'] = short_seq_Uf.anticodon_aa
        replacement_info_dict['contains_anticodon'] = short_seq_Uf.contains_anticodon
        replacement_info_dict['num_conserved'] = short_seq_Uf.num_conserved
        replacement_info_dict['num_unconserved'] = short_seq_Uf.num_unconserved
        replacement_info_dict['num_paired'] = short_seq_Uf.num_paired
        replacement_info_dict['num_unpaired'] = short_seq_Uf.num_unpaired
        unconserved_info_from_T_3prime = []
        for unconserved_tuple in short_seq_Uf.unconserved_info:
            unconserved_info_from_T_3prime.append((unconserved_tuple[0] + feature_index_adjustment,
                                                   unconserved_tuple[1],
                                                   unconserved_tuple[2]))
        replacement_info_dict['unconserved_info_from_T_3prime'] = unconserved_info_from_T_3prime
        unpaired_info_from_T_3prime = []
        for unpaired_tuple in short_seq_Uf.unpaired_info:
            unpaired_info_from_T_3prime.append((unpaired_tuple[0] + feature_index_adjustment,
                                                unpaired_tuple[1] + feature_index_adjustment,
                                                unpaired_tuple[2],
                                                unpaired_tuple[3]))
        replacement_info_dict['unpaired_info_from_T_3prime'] = unpaired_info_from_T_3prime
        replacement_info_dict['profiled_seq_without_terminus_length'] = short_seq_Uf.profiled_seq_length - short_seq_Uf.length_3prime_terminus

        dict_Tf = self.dict_Tf
        dict_Tf.pop(short_seq_Tf.name)
        seqs_U = [dict_Uf[name_Uf] for name_Uf in short_seq_Tf.names_U]
        for seq_Uf in seqs_U:
            seq_Uf.name_T = None
        dict_Us = self.dict_Us
        for long_seq_Tf in long_seqs_Tf:
            dict_Tf.pop(long_seq_Tf.name)

            for seq_Uf in [dict_Uf[name_Uf] for name_Uf in long_seq_Tf.names_U]:
                seq_Uf.name_T = None
                dict_Uf.pop(seq_Uf.name)
                seq_Us = UniqueTransferredProfileSequence(seq_Uf, replacement_info_dict)
                seqs_U.append(seq_Us)
                dict_Us[seq_Us.name] = seq_Us

        consol_seqs_Tf = self.get_trimmed_seqs(seqs_U, TrimmedFullProfileSequence)

        if len(consol_seqs_Tf) > 1:
            raise ConfigError(f"Consolidation should have produced only 1 trimmed profiled tRNA sequence, not {len(consol_seqs_Tf)}.")

        consol_seq_Tf = consol_seqs_Tf[0]
        dict_Tf[consol_seq_Tf.name] = consol_seq_Tf

        return consol_seq_Tf


    def threeprime_dereplicate_truncated_sequences(self):
        """Recover Tc that are found to be 3' subseqs of Nf and thus legitimate 3' tRNA fragments.
        These Tc are folded into Nf, while unrecovered Tc are themselves 3'-dereplicated, forming
        another pool of Nc."""
        start_time = time.time()
        self.progress.new("Dereplicating trimmed seqs with a truncated feature profile")
        self.progress.update("...")

        # Prefix dereplicate both Tc and Nf from the 3' end.
        names_Tc = []
        reverse_Tc_strings = []
        seqs_Tc = []
        for name, seq_Tc in self.dict_Tc_nontrna.items():
            names_Tc.append(name)
            reverse_Tc_strings.append(seq_Tc.string[::-1])
            seqs_Tc.append(seq_Tc)
        names_Nf = []
        reverse_Nf_strings = []
        seqs_Nf = []
        for name, seq_Nf in self.dict_Nf.items():
            names_Nf.append(name)
            reverse_Nf_strings.append(seq_Nf.string[::-1])
            seqs_Nf.append(seq_Nf)
        clusters = Dereplicator(names_Tc + names_Nf,
                                reverse_Tc_strings + reverse_Nf_strings,
                                extras=seqs_Tc + seqs_Nf).prefix_dereplicate()

        # Associate each Tc with any Nf that contain it as a 3' subseq. Since a Tc can be a 3'
        # subseq of multiple Nf, do not reconstruct a feature profile for the Tc. (Similarly, the
        # seed Tf in an Nf need not have the same profile as other Tf in the Nf.) Uc in the Tc are
        # therefore not included in the features table of the tRNA-seq database.

        # Clusters cannot contain > 1 Nf, as these have already been 3'-dereplicated (by
        # definition). Nf can seed clusters and also be members of clusters seeded by a Tc. In the
        # latter case, only Tc that are shorter than the Nf in the cluster (3' subseqs of the Nf)
        # are incorporated as members of the Nf.

        # There are three types of cluster: 1. clusters consisting of a single Nf (ignore), 2.
        # clusters containing an Nf as seed or member with shorter Tc members, and 3. clusters
        # seeded by Tc. If a Tc is found in group 2 (part of one or more longer Nf) then ignore it
        # in group 3 (do not include it in Nc formed from group 3 clusters). The alternatives do not
        # make sense -- including the Tc in Nc but not Nf; or withholding the Tc from both Nc and
        # Nf, perhaps as a new category of sequence.

        # This dict relates Tc to Nf containing them.
        dict_Tc_Nf = defaultdict(list)
        # This dict relates Tc to other Tc found to be subseqs of the former.
        dict_Tc_Tc = {}
        for cluster in clusters:
            if len(cluster.member_names) == 1:
                if isinstance(cluster.member_extras[0], NormalizedSequence):
                    continue

            seq_Nf = None
            seed_seq_Tc = cluster.member_extras[0] if isinstance(cluster.member_extras[0], TrimmedTruncatedProfileSequence) else None
            if seed_seq_Tc:
                members_seed_seq_Tc = []
                dict_Tc_Tc[seed_seq_Tc.name] = (seed_seq_Tc, members_seed_seq_Tc)
            for seq in cluster.member_extras:
                # Members of each cluster are pre-sorted in descending order of seq length. There
                # cannot be an Nf and Tc of the same length.
                if seq_Nf:
                    if isinstance(seq, TrimmedTruncatedProfileSequence):
                        dict_Tc_Nf[seq.name].append(seq_Nf)
                        continue
                    else:
                        raise ConfigError("It appears that a cluster in the 3' dereplication "
                                          "of trimmed sequences with truncated profiles and normalized sequences with full profiles "
                                          "contains >1 normalized sequence, when it should only contain 0 or 1.")

                if isinstance(seq, NormalizedSequence):
                    seq_Nf = seq
                    continue

                # The cluster is seeded by a Tc.
                members_seed_seq_Tc.append(seq)

        # Add Tc to matching Nf.
        dict_Tc_nontrna = self.dict_Tc_nontrna
        dict_Tc_trna = self.dict_Tc_trna
        dict_Uc_nontrna = self.dict_Uc_nontrna
        dict_Uc_trna = self.dict_Uc_trna
        # It is important to determine whether Tc contain an anticodon to enable the later
        # measurement of isoacceptor abundances. A truncated profile may stop 3' of the anticodon,
        # but the anticodon may still be in the sequence. The presence of the anticodon in Tc is
        # therefore inferred from the longest Tf in the matching Nf.
        dict_Nf_anticodon = {} # This dict saves time finding the position of the anticodon relative to the 3' terminus of Nf
        RELATIVE_ANTICODON_LOOP_INDEX = self.RELATIVE_ANTICODON_LOOP_INDEX
        dict_Tf = self.dict_Tf
        dict_Uf = self.dict_Uf
        for name_Tc, seqs_Nf in dict_Tc_Nf.items():
            seq_Tc = dict_Tc_nontrna.pop(name_Tc)
            # Tc has been confirmed as tRNA, so transfer the object between dicts.
            seq_Tc.category = 'trna'
            dict_Tc_trna[seq_Tc.name] = seq_Tc
            length_Tc = len(seq_Tc.string)

            seqs_Uc = []
            for name_Uc in seq_Tc.names_U:
                # Uc has been confirmed as tRNA.
                seq_Uc = dict_Uc_nontrna.pop(name_Uc)
                dict_Uc_trna[name_Uc] = seq_Uc
                seqs_Uc.append(seq_Uc)

            for seq_Nf in seqs_Nf:
                seq_Nf.names_T.append(seq_Tc.name)
                seq_Nf.categories_T.append('Tc_trna')
                length_Nf = len(seq_Nf.string)
                seq_Nf.starts_T_in_N.append(length_Nf - length_Tc)
                seq_Nf.stops_T_in_N.append(length_Nf)
                seq_Tc.names_N.append(seq_Nf.name)

                if not seq_Tc.contains_anticodon:
                    # Determine from the first Nf in which Tc is found whether Tc contains the
                    # anticodon.
                    seq_Tf = dict_Tf[seq_Nf.names_T[0]]
                    try:
                        # The position of the anticodon in Nf has already been found.
                        anticodon_start_relative_to_3prime_terminus = dict_Nf_anticodon[seq_Nf.name]
                    except KeyError:
                        try:
                            anticodon_loop_start = seq_Tf.feature_starts[RELATIVE_ANTICODON_LOOP_INDEX]
                        except IndexError:
                            # The anticodon loop was not reached in the profile.
                            anticodon_loop_start = -1
                        if anticodon_loop_start > -1:
                            # The anticodon loop was profiled.
                            anticodon_start = anticodon_loop_start + 2
                            # The position of the anticodon relative to the 3' terminus is a negative number.
                            anticodon_start_relative_to_3prime_terminus = anticodon_start - dict_Uf[seq_Tf.names_U[0]].feature_starts[-1]
                        else:
                            # The anticodon loop was not profiled, indicated by a positive number.
                            anticodon_start_relative_to_3prime_terminus = 1
                        dict_Nf_anticodon[seq_Nf.name] = anticodon_start_relative_to_3prime_terminus
                    if anticodon_start_relative_to_3prime_terminus == 1:
                        continue
                    if length_Tc + anticodon_start_relative_to_3prime_terminus >= 0:
                        seq_Tc.contains_anticodon = True
                        for seq_Uc in seqs_Uc:
                            seq_Uc.contains_anticodon = True

        # Tc that don't match Nf are grouped into Nc.
        dict_Nc = self.dict_Nc
        for seed_Tc_name, entry in dict_Tc_Tc.items():
            seed_seq_Tc, member_seqs_Tc = entry

            # If a Tc is a 3' subseq of an Nf, then it and all shorter Tc should be excluded from
            # the new Nc.
            for index_Tc, seq_Tc in enumerate(member_seqs_Tc):
                if seq_Tc.name in dict_Tc_Nf:
                    member_seqs_Tc = member_seqs_Tc[: index_Tc]
                    break

            dict_Nc[seed_Tc_name] = NormalizedTruncatedProfileSequence(member_seqs_Tc)

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed recovering tRNA with truncated feature profile (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()


    def report_threeprime_dereplication_statistics(self):
        """Report to terminal stats regarding 3'-dereplication immediately after these steps."""
        count_Nf = len(self.dict_Nf)
        anticodon_count_Nf = 0
        complete_count_Nf = 0
        mean_spec_T_Nf = 0
        mean_nonspec_T_Nf = 0
        mean_spec_U_Nf = 0
        mean_nonspec_U_Nf = 0
        spec_reads_Nf = 0
        nonspec_reads_Nf = 0
        max_spec_reads_Nf = 0
        max_nonspec_reads_Nf = 0
        max_total_reads_Nf = 0
        dict_Uc_trna = self.dict_Uc_trna
        for seq_Nf in self.dict_Nf.values():
            if seq_Nf.contains_anticodon:
                anticodon_count_Nf += 1
                if seq_Nf.has_complete_feature_set:
                    complete_count_Nf += 1
            spec_reads = 0
            nonspec_reads = 0
            for category_T, name_T in zip(seq_Nf.categories_T, seq_Nf.names_T):
                seq_T = getattr(self, 'dict_' + category_T)[name_T]
                if category_T == 'Tf':
                    read_counts = []
                    for name_U, category_U in zip(seq_T.names_U, seq_T.categories_U):
                        read_counts.append(getattr(self, 'dict_' + category_U)[name_U].read_count)
                elif category_T == 'Tc_trna':
                    read_counts = [dict_Uc_trna[name_U].read_count for name_U in seq_T.names_U]
                else:
                    raise Exception(f"The trimmed seq ({seq_T.name}) in the normalized seq ({seq_Nf.name}) has the unexpected class, `{type(seq_T)}`.")
                if len(seq_T.names_N) == 1:
                    mean_spec_T_Nf += 1
                    for read_count in read_counts:
                        mean_spec_U_Nf += 1
                        spec_reads += read_count
                        spec_reads_Nf += read_count
                else:
                    mean_nonspec_T_Nf += 1
                    for read_count in read_counts:
                        mean_nonspec_U_Nf += 1
                        nonspec_reads += read_count
                        nonspec_reads_Nf += read_count
            if spec_reads > max_spec_reads_Nf:
                max_spec_reads_Nf = spec_reads
            if nonspec_reads > max_nonspec_reads_Nf:
                max_nonspec_reads_Nf = nonspec_reads
            if spec_reads + nonspec_reads > max_total_reads_Nf:
                max_total_reads_Nf = spec_reads + nonspec_reads
        mean_spec_T_Nf /= count_Nf
        mean_nonspec_T_Nf /= count_Nf
        mean_spec_U_Nf /= count_Nf
        mean_nonspec_U_Nf /= count_Nf
        spec_reads_Nf /= count_Nf
        nonspec_reads_Nf /= count_Nf
        count_Tc_trna_Nf = len(self.dict_Tc_trna)
        count_Uc_trna_Nf = len(self.dict_Uc_trna)
        read_count_Uc_trna_Nf = 0
        for seq_Uc_trna in self.dict_Uc_trna.values():
            read_count_Uc_trna_Nf += seq_Uc_trna.read_count
        mean_Tc_trna_Nf = count_Tc_trna_Nf / count_Nf
        mean_Uc_trna_Nf = count_Uc_trna_Nf / count_Nf
        mean_Uc_reads_Nf = read_count_Uc_trna_Nf / count_Nf

        count_Nc = len(self.dict_Nc)
        anticodon_count_Nc = 0
        mean_spec_T_Nc = 0
        mean_nonspec_T_Nc = 0
        mean_spec_U_Nc = 0
        mean_nonspec_U_Nc = 0
        mean_spec_reads_Nc = 0
        mean_nonspec_reads_Nc = 0
        dict_Tc_nontrna = self.dict_Tc_nontrna
        dict_Uc_nontrna = self.dict_Uc_nontrna
        for seq_Nc in self.dict_Nc.values():
            if seq_Nc.contains_anticodon:
                anticodon_count_Nc += 1
            for name_T in seq_Nc.names_T:
                seq_T = dict_Tc_nontrna[name_T]
                if len(seq_T.names_N) == 1:
                    mean_spec_T_Nc += 1
                    for name_U in seq_T.names_U:
                        mean_spec_U_Nc += 1
                        mean_spec_reads_Nc += dict_Uc_nontrna[name_U].read_count
                else:
                    mean_nonspec_T_Nc += 1
                    for name_U in seq_T.names_U:
                        mean_nonspec_U_Nc += 1
                        mean_nonspec_reads_Nc += dict_Uc_nontrna[name_U].read_count
        mean_spec_T_Nc /= count_Nc
        mean_nonspec_T_Nc /= count_Nc
        mean_spec_U_Nc /= count_Nc
        mean_nonspec_U_Nc /= count_Nc
        mean_spec_reads_Nc /= count_Nc
        mean_nonspec_reads_Nc /= count_Nc

        warning = self.run.warning
        info_single = self.run.info_single
        info = self.run.info

        warning(None, "3' DEREPLICATION RESULTS", lc='green')
        info_single("subject to change -- see summary output file for final results")

        warning(None, "Normalized seq counts", lc='cyan')
        info("tRNA profile", count_Nf)
        info("Truncated tRNA profile", count_Nc)

        warning(None, "Normalized seqs with tRNA profile", lc='cyan')
        info("Containing anticodon", anticodon_count_Nf)
        info("Containing complete feature set", complete_count_Nf)
        info("Mean specific trimmed seqs per seq", round(mean_spec_T_Nf, 1))
        info("Mean nonspecific trimmed seqs per seq", round(mean_nonspec_T_Nf, 1))
        info("Mean specific unique seqs per seq", round(mean_spec_U_Nf, 1))
        info("Mean nonspecific unique seqs per seq", round(mean_nonspec_U_Nf, 1))
        info("Mean specific reads per seq", round(spec_reads_Nf, 1))
        info("Mean nonspecific reads per seq", round(nonspec_reads_Nf, 1))
        info("Max specific reads per seq", max_spec_reads_Nf)
        info("Max nonspecific reads per seq", max_nonspec_reads_Nf)
        info("Max total reads per seq", max_total_reads_Nf)
        info("Recovered trimmed seqs with truncated profile", count_Tc_trna_Nf)
        info("Recovered unique seqs with truncated profile", count_Uc_trna_Nf)
        info("Recovered reads with truncated profile", read_count_Uc_trna_Nf)
        info("Mean recovered trunc trimmed seqs per seq", round(mean_Tc_trna_Nf, 2))
        info("Mean recovered trunc unique seqs per seq", round(mean_Uc_trna_Nf, 2))
        info("Mean recovered trunc reads per seq", round(mean_Uc_reads_Nf, 2))
        info("Consolidated trimmed tRNA seqs", self.count_consol_Tf, nl_after=1)

        warning(None, "Normalized seqs with truncated tRNA profile", lc='cyan')
        info("Containing anticodon", anticodon_count_Nc)
        info("Mean specific trimmed seqs per seq", round(mean_spec_T_Nc, 1))
        info("Mean nonspecific trimmed seqs per seq", round(mean_nonspec_T_Nc, 1))
        info("Mean specific unique seqs per seq", round(mean_spec_U_Nc, 1))
        info("Mean nonspecific unique seqs per seq", round(mean_nonspec_U_Nc, 1))
        info("Mean specific reads per seq", round(mean_spec_reads_Nc, 1))
        info("Mean nonspecific reads per seq", round(mean_nonspec_reads_Nc, 1), nl_after=2 if self.write_checkpoints else 1)


    def write_checkpoint_files(self, checkpoint_name):
        self.progress.new(f"Writing intermediate files for the \"{checkpoint_name}\" checkpoint")

        if not os.path.exists(self.checkpoint_dir):
            os.mkdir(self.checkpoint_dir)

        checkpoint_subdir_path = self.checkpoint_subdir_dict[checkpoint_name]
        if not os.path.exists(checkpoint_subdir_path):
            os.mkdir(checkpoint_subdir_path)

        for intermed_file_key, intermed_file_path in self.intermed_file_path_dict[checkpoint_name].items():
            self.progress.update(f"{os.path.basename(intermed_file_path)}")
            # The key, e.g., "dict_Uf", corresponds to the attribute to be saved to file.
            with open(intermed_file_path, 'wb') as intermed_file:
                pkl.dump(getattr(self, intermed_file_key), intermed_file, protocol=pkl.HIGHEST_PROTOCOL)

        self.progress.end()
        self.run.info_single(f"Wrote \"{checkpoint_name}\" checkpoint intermediate files to {checkpoint_subdir_path}", cut_after=200)


    def load_checkpoint_files(self, checkpoint_name):
        pid = f"Loading intermediate files at the \"{checkpoint_name}\" checkpoint"
        self.progress.new(pid)

        for intermed_file_key, intermed_file_path in self.intermed_file_path_dict[checkpoint_name].items():
            self.progress.update_pid(pid)
            self.progress.update(f"{os.path.basename(intermed_file_path)}")
            with open(intermed_file_path, 'rb') as f:
                setattr(self, intermed_file_key, pkl.load(f))

        with open(self.analysis_summary_path, 'a') as f:
            f.write(f"\nAnalysis restarted from the \"{checkpoint_name}\" checkpoint\n")

        self.progress.end()
        self.run.info_single(f"Loaded \"{checkpoint_name}\" checkpoint intermediate files from {self.checkpoint_subdir_dict[checkpoint_name]}")


    def threeprime_dereplicate_sequences_without_terminus(self):
        """Find tRNA sequences missing a 3' terminus. By default, U required a 3' terminus to have
        been profiled as tRNA. Un are searched against Nf. Un is recovered as tRNA when it is a 3'
        subseq of Nf or is longer than Nf with a complete profile and thus is shown to have a 5'
        extension. Recovered Un each generate a Um object."""
        self.progress.new("Dereplicating tRNA seqs ending in discriminator nt")
        self.progress.update("...")

        # 3'-dereplicate Nf and Un.
        names = []
        reverse_strings = []
        extras = []
        # Nf are added first to the dereplicator so that they appear first in the clusters. This
        # allows Un that are identical to the Nf (due to prior trimming of the 3' terminus from the
        # Nf) to always be recovered.
        for name_Nf, seq_Nf in self.dict_Nf.items():
            names.append(name_Nf)
            reverse_strings.append(seq_Nf.string[::-1])
            extras.append(seq_Nf)
        for name_Un, seq_Un in self.dict_Un.items():
            if len(seq_Un.string) >= self.min_trna_frag_size:
                names.append(name_Un)
                reverse_strings.append(seq_Un.string[::-1])
                extras.append(seq_Un)
        clusters = Dereplicator(names, reverse_strings, extras=extras).prefix_dereplicate()

        # Search clusters for Nf and qualifying Un. Some Un may be 3' subseqs of > 1 Nf. Un cannot
        # be longer than > 1 Nf, as Nf have already been dereplicated. The same Nf can be found in
        # multiple clusters as subseqs of different longer Un.
        dict_Un_Nf = defaultdict(list)
        dict_Un = self.dict_Un
        dict_Um = self.dict_Um
        dict_Tm = self.dict_Tm
        for cluster in clusters:
            if len(cluster.member_seqs) == 1:
                continue

            # Check that there is an Nf in the cluster -- there cannot be more than one.
            cluster_Nf_index = None
            seq_Nf = None
            length_Nf = None
            complete_feature_set_in_Nf = None
            for member_index, seq in enumerate(cluster.member_extras):
                if isinstance(seq, NormalizedFullProfileSequence):
                    cluster_Nf_index = member_index
                    seq_Nf = seq
                    length_Nf = len(seq_Nf.string)
                    if member_index > 0:
                        complete_feature_set_in_Nf = seq_Nf.has_complete_feature_set
                    break
            else:
                continue

            # To reach this point, an Nf must have been found in the cluster. Now process any longer
            # Un.
            for seq_Un in cluster.member_extras[: cluster_Nf_index]:
                # If the Nf has a complete feature profile (is a full-length tRNA), then the
                # overhanging 5' bases in the recovered Un can be trimmed as "extra" 5' bases.
                # Otherwise, it is possible the overhanging bases are part of an artifact joined to
                # the 5' end of a tRNA fragment, so conservatively ignore the Un.
                if not complete_feature_set_in_Nf:
                    break

                name_Un = seq_Un.name
                try:
                    dict_Un.pop(name_Un)
                except KeyError:
                    # The Un has already been recovered in another cluster, necessarily as a
                    # supersequence of the same Nf.
                    continue

                seq_Um = UniqueMappedSequence(seq_Un.string, name_Un, seq_Un.read_count, xtra_5prime_length=len(seq_Un.string) - length_Nf)
                dict_Um[name_Un] = seq_Um

                seq_Tm = TrimmedMappedSequence(seq_Um)
                seq_Nf.names_T.append(name_Un)
                seq_Nf.categories_T.append('Tm')
                seq_Tm.names_N.append(seq_Nf.name)
                seq_Nf.starts_T_in_N.append(0)
                seq_Nf.stops_T_in_N.append(length_Nf)
                dict_Tm[name_Un] = seq_Tm

            # Find all Un in the cluster ≤ Nf length.
            for seq_Un in cluster.member_extras[cluster_Nf_index + 1: ]:
                dict_Un_Nf[seq_Un.name].append(seq_Nf)

        for name_Un, seqs_Nf in dict_Un_Nf.items():
            seq_Un = dict_Un.pop(name_Un)

            seq_Um = UniqueMappedSequence(seq_Un.string, name_Un, seq_Un.read_count)
            dict_Um[name_Un] = seq_Um

            seq_Tm = TrimmedMappedSequence(seq_Um)
            length_Um = len(seq_Um.string)
            # The same Nf can be found in multiple clusters, so it can be represented multiple times
            # in `seqs_Nf`.
            for seq_Nf in set(seqs_Nf):
                seq_Nf.names_T.append(name_Un)
                seq_Nf.categories_T.append('Tm')
                seq_Tm.names_N.append(seq_Nf.name)
                length_Nf = len(seq_Nf.string)
                seq_Nf.starts_T_in_N.append(length_Nf - length_Um)
                seq_Nf.stops_T_in_N.append(length_Nf)
            dict_Tm[name_Un] = seq_Tm

        self.progress.end()


    def map_fragments(self):
        """Map unprofiled tRNA fragments to longer profiled tRNA sequences. Fragments only missing a
        3' terminus were already found with `threeprime_dereplicate_sequences_without_terminus` or
        by profiling if '' was an accepted 3' terminus.

        EXAMPLE:
        Nf:                   (GT)TCCGTGATAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGGTTCAATTCCCCGTCGCGGAG
        Um1 (extra 5' bases) :  T TCCGTGATAGTTTAATGGTCAGAATGG
        Um2 (interior)       :           TAGTTTAATGGTCAGAATGGGCGCTTGTCGCGTGCCAGATCGGGG
        """
        start_time = time.time()

        pid = "Set up search of unprofiled reads to profiled tRNA"
        self.progress.new(pid)

        self.progress.update("Getting queries from unprofiled reads")
        temp_dir_path = filesnpaths.get_temp_directory_path()
        query_fasta_path = os.path.join(temp_dir_path, 'query.fa')

        dict_Un = self.dict_Un
        min_trna_frag_size = self.min_trna_frag_size
        query_count = 0
        with open(query_fasta_path, 'w') as query_fasta:
            for seq_Un in [seq_Un for seq_Un in dict_Un.values() if len(seq_Un.string) >= min_trna_frag_size]:
                # Include Un length in the defline for the purposes of parsing vmatch output.
                query_fasta.write(f">{seq_Un.name}-{len(seq_Un.string)}\n{seq_Un.string}\n")
                query_count += 1


        self.progress.update_pid(pid)
        self.progress.update("Getting targets from profiled tRNAs")
        # Un are mapped to Nf with extra 5' bases added when present in underlying U. Multiple
        # targets for each Nf are therefore produced for different 5' extensions.
        target_fasta_path = os.path.join(temp_dir_path, 'target.fa')
        dict_Tf = self.dict_Tf
        dict_Uf = self.dict_Uf
        dict_Us = self.dict_Us
        with open(target_fasta_path, 'w') as target_fasta:
            for seq_Nf in self.dict_Nf.values():
                string_Nf = seq_Nf.string
                # The longest Tf (the first in the list) is by design the only one of the Tf forming
                # Nf that may have extra 5' bases.
                longest_Tf = dict_Tf[seq_Nf.names_T[0]]
                if longest_Tf.uniq_with_xtra_5prime_count > 0:
                    set_5prime_string = set()
                    for name_U in longest_Tf.names_U:
                        try:
                            seq_U = dict_Uf[name_U]
                        except KeyError:
                            seq_U = dict_Us[name_U]
                        if seq_U.xtra_5prime_length > 0:
                            set_5prime_string.add(seq_U.string[: seq_U.xtra_5prime_length])

                    # Avoid creating superfluous target seqs that are subseqs of other target seqs
                    # due to a 5' extension of an Nf being a subseq of a longer 5' extension of the
                    # same Nf.
                    strings_5prime = sorted(set_5prime_string, key=lambda string_5prime: -len(string_5prime))
                    string_5prime_additions = [strings_5prime[0]]
                    for string_5prime in strings_5prime[1: ]:
                        length_5prime = len(string_5prime)
                        for string_5prime_addition in string_5prime_additions:
                            if string_5prime == string_5prime_addition[-length_5prime: ]:
                                break
                        else:
                            string_5prime_additions.append(string_5prime)

                    for index_5prime, string_5prime in enumerate(strings_5prime):
                        # Use an index to distinguish otherwise equivalent targets with different 5'
                        # extensions of the same length.
                        target_fasta.write(f">{seq_Nf.name}-{len(string_5prime)}-{index_5prime}\n{string_5prime}{string_Nf}\n")
                else:
                    target_fasta.write(f">{seq_Nf.name}-0-0\n{string_Nf}\n") # no extra 5' bases
        self.progress.end()


        # Use a 10x bigger query chunk size than the Vmatch default, as the rather conservative
        # default is tailored to searches with mismatches/indels. It takes longer to process each
        # chunk in these searches, and these searches may generate more alignments per chunk.
        query_chunk_size_default = 10 * Vmatch.QUERY_CHUNK_SIZE_DEFAULT
        match_df = Vmatch(argparse.Namespace(match_mode='exact_query_substring',
                                             fasta_db_file=target_fasta_path,
                                             fasta_query_file=query_fasta_path,
                                             num_threads=self.num_threads,
                                             query_chunk_size=query_count // self.num_threads + 1 if query_count < query_chunk_size_default else query_chunk_size_default // self.num_threads,
                                             temp_dir=temp_dir_path)).search_queries()

        pid = "Filtering matches"
        self.progress.new(pid)
        self.progress.update("...")

        self.restructure_fragment_match_table(match_df)

        # Process each Un match. Each Un can match more than one Nf.
        match_gb = match_df.groupby('query_name')
        del match_df
        gc.collect()

        fragment_filter_progress_interval = 25000
        total_matched_queries = len(match_gb)
        pp_total_matched_queries = pp(total_matched_queries)
        num_filtered_queries = -1
        dict_Um = self.dict_Um
        dict_Tm = self.dict_Tm
        dict_Nf = self.dict_Nf
        for name_Un, query_match_df in match_gb:
            num_filtered_queries += 1
            if num_filtered_queries % fragment_filter_progress_interval == 0:
                pp_progress_interval_end = pp(total_matched_queries if num_filtered_queries + fragment_filter_progress_interval > total_matched_queries else num_filtered_queries + fragment_filter_progress_interval)
                self.progress.update_pid(pid)
                self.progress.update(f"Queries {pp(num_filtered_queries + 1)}-{pp_progress_interval_end}/{pp_total_matched_queries}")

            # Each Un with a validated match will yield a Um and Tm.
            seq_Um = None
            seq_Tm = None

            for name_Nf, length_target_5prime, query_start, length_Un in zip(query_match_df['target_name'],
                                                                             query_match_df['length_5prime'],
                                                                             query_match_df['query_start_in_target'],
                                                                             query_match_df['query_length']):
                query_stop = query_start + length_Un
                stop_Tm_in_Nf = query_stop - length_target_5prime

                if stop_Tm_in_Nf <= 0:
                    # Ignore queries that align entirely to extra 5' bases. Un mapping exclusively
                    # to the 5' extension that are long enough to fulfill the minimum length
                    # requirement may be mapping to an artifactual chimeric sequence.
                    continue

                seq_Nf = dict_Nf[name_Nf]

                if not seq_Um:
                    # Enter this block the first time the Un query validly matches an Nf.
                    seq_Un = dict_Un.pop(name_Un)

                    # Assume that 5' extensions are the same for the query regardless of the reference.
                    # This could be false when
                    # 1. tRNA profiling erroneously identified the end of the acceptor stem
                    # or 2. the query mapped to different places at the end of the acceptor stem in different tRNAs.
                    if length_target_5prime - query_start > 0:
                        length_query_5prime = length_target_5prime - query_start
                        start_Tm_in_Nf = 0
                    else:
                        length_query_5prime = 0
                        start_Tm_in_Nf = query_start - length_target_5prime

                    seq_Um = UniqueMappedSequence(seq_Un.string, name_Un, seq_Un.read_count, xtra_5prime_length=length_query_5prime)
                    dict_Um[name_Un] = seq_Um

                    seq_Tm = TrimmedMappedSequence(seq_Um)
                    seq_Nf.names_T.append(name_Un)
                    seq_Nf.categories_T.append('Tm')
                    seq_Tm.names_N.append(seq_Nf.name)
                    seq_Nf.starts_T_in_N.append(start_Tm_in_Nf)
                    seq_Nf.stops_T_in_N.append(stop_Tm_in_Nf)
                    dict_Tm[name_Un] = seq_Tm
                else:
                    for prev_name_T, prev_category_T in zip(seq_Nf.names_T[::-1], seq_Nf.categories_T[::-1]):
                        # Ensure that Tm maps to Nf only once. Multiple targets can be created from
                        # the same Nf for different 5' extensions. Tm are added after Tf and Tc to
                        # the list of T in Nf.
                        if seq_Tm.name == prev_name_T:
                            break

                        if prev_category_T != 'Tm':
                            seq_Nf.names_T.append(seq_Tm.name)
                            seq_Nf.categories_T.append('Tm')
                            seq_Tm.names_N.append(seq_Nf.name)
                            if length_target_5prime - query_start > 0:
                                start_Tm_in_Nf = 0
                            else:
                                start_Tm_in_Nf = query_start - length_target_5prime
                            seq_Nf.starts_T_in_N.append(start_Tm_in_Nf)
                            seq_Nf.stops_T_in_N.append(stop_Tm_in_Nf)
                            break
        self.progress.end()

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed mapping tRNA fragments (min)", time.time() - start_time, is_time_value=True))


    def restructure_fragment_match_table(self, match_df):
        """Helper method for `map_fragments`."""
        names_U = []
        for query_name in match_df['query_name']:
            # Sequence names in anvi'o cannot contain a hyphen.
            name_U, length_U = query_name.split('-')
            names_U.append(name_U)
        match_df.loc[:, 'query_name'] = names_U

        names_N = []
        lengths_5prime = []
        for target_name in match_df['target_name']:
            name_N, length_5prime, index_5prime = target_name.split('-')
            names_N.append(name_N)
            lengths_5prime.append(int(length_5prime))
        match_df.loc[:, 'target_name'] = names_N
        match_df['length_5prime'] = lengths_5prime


    def report_mapping_statistics(self):
        """Report to terminal stats on fragment mapping immediately after these steps."""
        count_spec_Nf = 0
        count_nonspec_Nf = 0
        count_any_Nf = 0
        mean_spec_Tm_Nf = 0
        mean_nonspec_Tm_Nf = 0
        spec_reads_Nf = 0
        nonspec_reads_Nf = 0
        absent_3prime_terminus_seqs_Tm = 0
        absent_3prime_terminus_reads_Tm = 0
        for seq_Nf in self.dict_Nf.values():
            if seq_Nf.spec_map_seq_count:
                count_spec_Nf += 1
                mean_spec_Tm_Nf += seq_Nf.spec_map_seq_count
                spec_reads_Nf += seq_Nf.spec_map_read_count
                if not PROFILE_ABSENT_3PRIME_TERMINUS:
                    absent_3prime_terminus_seqs_Tm += seq_Nf.absent_3prime_terminus_seq_count
                    absent_3prime_terminus_reads_Tm += seq_Nf.absent_3prime_terminus_read_count
            if seq_Nf.nonspec_map_seq_count:
                count_nonspec_Nf += 1
                mean_nonspec_Tm_Nf += seq_Nf.nonspec_map_seq_count
                nonspec_reads_Nf += seq_Nf.nonspec_map_read_count
            if seq_Nf.spec_map_seq_count or seq_Nf.nonspec_map_seq_count:
                count_any_Nf += 1
        count_Nf = len(self.dict_Nf)
        mean_spec_Tm_Nf /= count_Nf
        mean_nonspec_Tm_Nf /= count_Nf
        spec_reads_Nf /= count_Nf
        nonspec_reads_Nf /= count_Nf

        count_spec_Tm = 0
        reads_spec_Tm = 0
        spec_short_5prime_seq_Tm = 0
        spec_long_5prime_seq_Tm = 0
        spec_short_5prime_read_Tm = 0
        spec_long_5prime_read_Tm = 0
        spec_short_5prime_Nf_names = []
        spec_long_5prime_Nf_names = []
        nonspec_short_5prime_Nf_names = []
        nonspec_long_5prime_Nf_names = []
        for seq_Tm in self.dict_Tm.values():
            if len(seq_Tm.names_N) == 1:
                count_spec_Tm += 1
                reads_spec_Tm += seq_Tm.read_count
                if seq_Tm.read_with_xtra_5prime_count:
                    if seq_Tm.long_5prime_extension_dict:
                        spec_long_5prime_seq_Tm += 1
                        spec_long_5prime_read_Tm += seq_Tm.read_count
                        spec_long_5prime_Nf_names.append(seq_Tm.names_N[0])
                    else:
                        spec_short_5prime_seq_Tm += 1
                        spec_short_5prime_read_Tm += seq_Tm.read_count
                        spec_short_5prime_Nf_names.append(seq_Tm.names_N[0])
            else:
                if seq_Tm.read_with_xtra_5prime_count:
                    if seq_Tm.long_5prime_extension_dict:
                        nonspec_long_5prime_Nf_names.extend(seq_Tm.names_N)
                    else:
                        nonspec_short_5prime_Nf_names.extend(seq_Tm.names_N)
        spec_short_5prime_Nf = len(set(spec_short_5prime_Nf_names))
        spec_long_5prime_Nf = len(set(spec_long_5prime_Nf_names))
        nonspec_short_5prime_Nf = len(set(nonspec_short_5prime_Nf_names))
        nonspec_long_5prime_Nf = len(set(nonspec_long_5prime_Nf_names))

        warning = self.run.warning
        info_single = self.run.info_single
        info = self.run.info

        warning(None, "FRAGMENT MAPPING RESULTS", lc='green', nl_before=1 if self.write_checkpoints else 0)
        info_single("subject to change -- see summary output file for final results")

        warning(None, "Normalized seqs with tRNA profile", lc='cyan')
        info("With specific mapping", count_spec_Nf)
        info("With nonspecific mapping", count_nonspec_Nf)
        info("With any mapping", count_any_Nf)
        info("Mean specific mapped seqs per seq", round(mean_spec_Tm_Nf, 2))
        info("Mean nonspecific mapped seqs per seq", round(mean_nonspec_Tm_Nf, 2))
        info("Mean specific mapped reads per seq", round(spec_reads_Nf, 2))
        info("Mean nonspecific mapped reads per seq", round(nonspec_reads_Nf, 2))
        info(f"With specific mapping to 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_short_5prime_Nf)
        info(f"With specific mapping to ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_long_5prime_Nf)
        info(f"With nonspecific mapping to 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", nonspec_short_5prime_Nf)
        info(f"With nonspecific mapping to ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", nonspec_long_5prime_Nf)

        warning(None, "Mapped seq counts", lc='cyan')
        info("Specific seqs", count_spec_Tm)
        info("Specific reads", reads_spec_Tm)
        info(f"Specific seqs with 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_short_5prime_seq_Tm)
        info(f"Specific seqs with ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_long_5prime_seq_Tm)
        info(f"Specific reads with 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_short_5prime_read_Tm)
        info(f"Specific reads with ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_long_5prime_read_Tm)
        if not PROFILE_ABSENT_3PRIME_TERMINUS:
            info("Seqs only missing a 3' terminus", absent_3prime_terminus_seqs_Tm)
            info("Reads only missing a 3' terminus", absent_3prime_terminus_reads_Tm)
            info_single("Consider including an absent 3' terminus (by using '_') "
                        "in the `anvi-trnaseq` parameterization of allowed 3' termini "
                        "if the number of mapped seqs identical to a normalized seq but missing a 3' terminus seems high.",
                        mc='red')


    def report_initialized_normalized_sequence_coverage_statistics(self):
        """Report to terminal stats on N coverages immediately after N initialization."""
        spec_read_Nf = 0
        mean_spec_cov_Nf = 0
        mean_nonspec_cov_Nf = 0
        total_length_Nf = 0
        max_spec_cov_Nf = 0
        max_nonspec_cov_Nf = 0
        max_total_cov_Nf = 0
        for seq_Nf in self.dict_Nf.values():
            spec_read_Nf += seq_Nf.spec_read_count
            length_Nf = len(seq_Nf.string)
            mean_spec_cov = seq_Nf.mean_spec_cov
            mean_nonspec_cov = seq_Nf.mean_nonspec_cov
            mean_spec_cov_Nf += mean_spec_cov * length_Nf
            mean_nonspec_cov_Nf += mean_nonspec_cov * length_Nf
            total_length_Nf += length_Nf
            if mean_spec_cov > max_spec_cov_Nf:
                max_spec_cov_Nf = mean_spec_cov
            if mean_nonspec_cov > max_nonspec_cov_Nf:
                max_nonspec_cov_Nf = mean_nonspec_cov
            if mean_spec_cov + mean_nonspec_cov > max_total_cov_Nf:
                max_total_cov_Nf = mean_spec_cov + mean_nonspec_cov
        mean_spec_cov_Nf /= total_length_Nf
        mean_nonspec_cov_Nf /= total_length_Nf

        spec_read_Nc = 0
        mean_spec_cov_Nc = 0
        mean_nonspec_cov_Nc = 0
        total_length_Nc = 0
        max_spec_cov_Nc = 0
        max_nonspec_cov_Nc = 0
        max_total_cov_Nc = 0
        for seq_Nc in self.dict_Nc.values():
            spec_read_Nc += seq_Nc.spec_read_count
            length_Nc = len(seq_Nc.string)
            mean_spec_cov = seq_Nc.mean_spec_cov
            mean_nonspec_cov = seq_Nc.mean_nonspec_cov
            mean_spec_cov_Nc += mean_spec_cov * length_Nc
            mean_nonspec_cov_Nc += mean_nonspec_cov * length_Nc
            total_length_Nc += length_Nc
            if mean_spec_cov > max_spec_cov_Nc:
                max_spec_cov_Nc = mean_spec_cov
            if mean_nonspec_cov > max_nonspec_cov_Nc:
                max_nonspec_cov_Nc = mean_nonspec_cov
            if mean_spec_cov + mean_nonspec_cov > max_total_cov_Nc:
                max_total_cov_Nc = mean_spec_cov + mean_nonspec_cov
        mean_spec_cov_Nc /= total_length_Nc
        mean_nonspec_cov_Nc /= total_length_Nc

        warning = self.run.warning
        info_single = self.run.info_single
        info = self.run.info

        warning(None, "NORMALIZATION RESULTS", lc='green', nl_before=1)
        info_single("subject to change -- see summary output file for final results")

        warning(None, "Normalized seqs with tRNA profile", lc='cyan')
        info("Specific reads", spec_read_Nf)
        info("Mean specific coverage", round(mean_spec_cov_Nf, 2))
        info("Mean nonspecific coverage", round(mean_nonspec_cov_Nf, 2))
        info("Max specific coverage", round(max_spec_cov_Nf, 2))
        info("Max nonspecific coverage", round(max_nonspec_cov_Nf, 2))
        info("Max total coverage", round(max_total_cov_Nf, 2))

        warning(None, "Normalized seqs with truncated tRNA profile", lc='cyan')
        info("Specific reads", spec_read_Nc)
        info("Mean specific coverage", round(mean_spec_cov_Nc, 2))
        info("Mean nonspecific coverage", round(mean_nonspec_cov_Nc, 2))
        info("Max specific coverage", round(max_spec_cov_Nc, 2))
        info("Max nonspecific coverage", round(max_nonspec_cov_Nc, 2))
        info("Max total coverage", round(max_total_cov_Nc, 2), nl_after=2)


    def find_substitutions(self):
        """Find sites of potential modification-induced substitutions."""
        start_time = time.time()
        pid = "Finding modification-induced substitutions"
        self.progress.new(pid)
        self.progress.update("...")

        # Cluster Nf. Clusters agglomerate Nf that differ from at least one other Nf in the cluster
        # by no more than 3 nts in 100 (by default) in a gapless end-to-end alignment with no
        # clipping.
        dict_Nf = self.dict_Nf
        names_Nf = []
        strings_Nf = []
        dict_Nf_feature_completeness = {}
        for name_Nf, seq_Nf in dict_Nf.items():
            names_Nf.append(name_Nf)
            strings_Nf.append(seq_Nf.string)
            dict_Nf_feature_completeness[name_Nf] = seq_Nf.has_complete_feature_set
        self.progress.end()

        agglomerator = Agglomerator(names_Nf, strings_Nf, num_threads=self.num_threads)
        # Provide a priority function for seeding clusters that favors, in order:
        # 1. Nf with a complete set of tRNA features,
        # 2. longer Nf,
        # 3. Nf with more alignments in the all-against-all search,
        # 4. alphanumeric order of the Nf name.
        agglomerator.agglomerate(max_mismatch_freq=self.agglom_max_mismatch_freq,
                                 priority_function=lambda aligned_ref: (-dict_Nf_feature_completeness[aligned_ref.name],
                                                                        -len(aligned_ref.seq_string),
                                                                        -len(aligned_ref.alignments),
                                                                        aligned_ref.name))
        agglom_aligned_ref_dict = agglomerator.agglom_aligned_ref_dict

        pid = "Decomposing clusters"
        self.progress.new(pid)
        self.progress.update("...")

        excluded_Nf_names = [] # used to exclude Nf from being considered as aligned queries in clusters (see below)
        represent_Nb_names = [] # used to prevent the same M from being created twice
        dict_name_Nb_M = defaultdict(list)
        dict_M = self.dict_M
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

            # A mod requires at least 3 different nts to be detected, and each Nf differs by at
            # least 1 nt (mismatch or gap), so for a cluster to form an M, it must contain at least
            # 3 Nf.
            if len(aligned_ref.alignments) < 2:
                continue

            aligned_ref_length = len(aligned_ref.seq_string)
            valid_aligned_queries = []
            for alignment in aligned_ref.alignments:
                # Nf should only align at the 3' end. Alignments to the interior of Nf can
                # theoretically occur when the reference is a tRNA-tRNA chimera.
                if aligned_ref_length != alignment.target_start + alignment.alignment_length:
                    continue

                query_name = alignment.aligned_query.name
                # The Nf query may have formed a M already. If the Nf had a complete feature
                # profile, or if it was the same length as such a sequence, then it should not be
                # able to form a longer M that would have 5' nts beyond the end of a complete
                # feature profile. This is only relevant when nonspecific Nf membership of M is
                # allowed, which is not currently the case.
                if query_name in excluded_Nf_names:
                    continue

                valid_aligned_queries.append(dict_Nf[query_name])

            # Confirm that 2 or more queries passed the filters, so at least 3 Nf are still in the
            # cluster.
            if len(valid_aligned_queries) < 2:
                continue

            valid_aligned_queries.sort(key=lambda seq_Nf: (-len(seq_Nf.string), -seq_Nf.has_complete_feature_set, seq_Nf.name))
            seqs_Nf = np.array([dict_Nf[ref_name]] + valid_aligned_queries)

            self.decompose_substitution_cluster(seqs_Nf, represent_Nb_names, excluded_Nf_names, dict_name_Nb_M, dict_M)

        # Remove nonspecific Nf from M. Treat the remaining Nf in M as a cluster, from which a new M
        # may be generated.
        dict_M_with_nonspec_Nb = {}
        count_nonspec_Nb = 0
        while dict_name_Nb_M:
            name_Nb, seqs_M = dict_name_Nb_M.popitem()
            if len(seqs_M) > 1:
                count_nonspec_Nb += 1
                for seq_M in seqs_M:
                    seq_M.names_Nb.remove(name_Nb)
                    dict_M_with_nonspec_Nb[seq_M.name] = seq_M

        count_M_with_nonspec_Nb = 0
        count_rejected_M = 0
        for seq_M in dict_M_with_nonspec_Nb.values():
            count_M_with_nonspec_Nb += 1
            seqs_Nf = []
            for name_Nb in seq_M.names_Nb:
                seqs_Nf.append(dict_Nf[name_Nb])

            if not seqs_Nf:
                # All Nb were nonspecific.
                count_rejected_M += 1
                dict_M.pop(seq_M.name)
                continue

            try:
                # Nf is no longer a representative Nf of a M ... for the time being -- M might be
                # vindicated.
                represent_Nb_names.remove(seqs_Nf[0].name)
            except ValueError:
                pass

            dict_M.pop(seq_M.name)
            found_M = self.decompose_substitution_cluster(np.array(seqs_Nf), represent_Nb_names, excluded_Nf_names, dict_name_Nb_M, dict_M)
            if not found_M:
                count_rejected_M += 1

        with open(self.analysis_summary_path, 'a') as f:
            f.write(self.get_summary_line("Time elapsed finding modification-induced substitutions (min)", time.time() - start_time, is_time_value=True))

        self.progress.end()

        if count_nonspec_Nb:
            self.run.info_single(f"{pp(count_nonspec_Nb)} nonspecific norm seqs "
                                 f"were found in {pp(count_M_with_nonspec_Nb)} candidate mod seqs, "
                                 f"resulting in {pp(count_rejected_M)} rejected candidates",
                                 cut_after=100, nl_before=2 if self.write_checkpoints else 0)


    def decompose_substitution_cluster(self, input_seqs_Nf, represent_Nb_names, excluded_Nf_names, dict_name_Nb_M, dict_M):
        max_length_Nf = len(input_seqs_Nf[0].string)
        seq_array = np.zeros((len(input_seqs_Nf), max_length_Nf), dtype=int)
        # Rather than using the ASCII representation of each character, which saves some time in
        # converting the sequence string to a numpy array, constrain the integer representation to
        # the smallest possible range of integers to speed up the bincount method used to determine
        # the number of unique nts at an alignment position.
        for seq_index, seq_Nf in enumerate(input_seqs_Nf):
            seq_array[seq_index, max_length_Nf - len(seq_Nf.string): ] += [NT_INT_DICT[nt] for nt in seq_Nf.string]

        # Find positions in the alignment with nt variability.
        alignment_pos_uniq_nt_counts = (
            np.bincount(
                (seq_array + np.arange(max_length_Nf, dtype=int) * NUM_NT_BINS).ravel(),
                minlength=max_length_Nf * NUM_NT_BINS
            ).reshape(-1, NUM_NT_BINS)[:, 1:] != 0
        ).sum(axis=1)
        alignment_positions_3_4_nts = (alignment_pos_uniq_nt_counts > 2).nonzero()[0]

        # Modification-induced substitutions must have ≥ 3 nts.
        if not alignment_positions_3_4_nts.size:
            return False

        alignment_positions_2_nts = (alignment_pos_uniq_nt_counts == 2).nonzero()[0]
        clusters = deque(((seq_array, input_seqs_Nf, alignment_positions_3_4_nts), ))
        for alignment_pos in alignment_positions_2_nts:
            next_clusters = deque() # Make a new object with each iteration rather than clearing the same one

            while clusters:
                seq_array, seqs_Nf, alignment_positions_3_4_nts = clusters.pop()

                # A sub requires ≥ 3 different nts to be detected, and each Nf differs by ≥ 1 nt, so
                # for a cluster to form an M it must contain ≥ 3 Nf.
                if seqs_Nf.size < 3:
                    continue

                aligned_nts = seq_array[:, alignment_pos]
                nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]

                if (nt_counts != 0).sum() < 2:
                    # There are now < 2 nts at the alignment position in the cluster under
                    # consideration. 2 different nts are needed to distinguish SNVs.
                    next_clusters.appendleft((seq_array, seqs_Nf, alignment_positions_3_4_nts))
                    continue

                # Add a new cluster for each SNV to the stack of clusters to process if: 1. the new
                # cluster contains ≥ 3 Nf and 2. the longest Nf (with a complete feature profile, if
                # applicable) in the new cluster has not yet formed an M.
                represented_nts = nt_counts.nonzero()[0] + 1
                for nt in represented_nts:
                    split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]

                    if split_cluster_seq_indices.size > 2:
                        split_cluster_Nf_seqs = seqs_Nf[split_cluster_seq_indices]

                        if split_cluster_Nf_seqs[0].name in represent_Nb_names:
                            # Nf already seeded an M, which would be the same M, as the same
                            # subcluster of Nf can be found in different agglomerations.
                            continue

                        next_clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                  split_cluster_Nf_seqs,
                                                  alignment_positions_3_4_nts.copy()))

            if next_clusters:
                clusters = next_clusters
            else:
                return False

        # Check alignment positions previously found to have 3-4 nts. Further split clusters when
        # positions now have 2 nts.
        next_clusters = deque()
        while clusters:
            seq_array, seqs_Nf, alignment_positions_3_4_nts = clusters.pop()
            candidates_to_remove = []

            for i, alignment_pos in enumerate(alignment_positions_3_4_nts):
                aligned_nts = seq_array[:, alignment_pos]
                nt_counts = np.bincount(aligned_nts, minlength=NUM_NT_BINS)[1: ]
                # At least 3 different nts are needed at a position to predict a mod.
                represented_nts = nt_counts.nonzero()[0] + 1
                if represented_nts.size < 2:
                    candidates_to_remove.append(i)
                elif represented_nts.size == 2:
                    candidates_to_remove.append(i)
                    for nt in represented_nts:
                        split_cluster_seq_indices = (aligned_nts == nt).nonzero()[0]

                        # At least 3 Nf are needed, and the split cluster cannot have already formed
                        # an M.
                        if split_cluster_seq_indices.size > 2:
                            split_cluster_Nf_seqs = seqs_Nf[split_cluster_seq_indices]

                            if split_cluster_Nf_seqs[0].name in represent_Nb_names:
                                continue

                            clusters.appendleft((seq_array[split_cluster_seq_indices, :],
                                                 split_cluster_Nf_seqs,
                                                 np.delete(alignment_positions_3_4_nts, candidates_to_remove)))
                    # Reevaluate previous alignment positions in the split clusters.
                    break
            else:
                # At least 1 position was discounted as no longer having 3-4 different nts, but
                # these positions had < 2 nts, and so did not cause the cluster to be split into new
                # clusters. Therefore, do not cycle through the remaining positions again to find
                # any more with < 3 nts.
                if candidates_to_remove:
                    next_clusters.appendleft((seqs_Nf, np.delete(alignment_positions_3_4_nts, candidates_to_remove)))
                else:
                    next_clusters.appendleft((seqs_Nf, alignment_positions_3_4_nts))

        if not next_clusters:
            return False
        clusters = next_clusters

        while clusters:
            seqs_Nf, sub_positions = clusters.pop() # Nf should have retained their order
            seqs_Nf = list(seqs_Nf) # Turn the array into a list
            represent_Nf_seq = seqs_Nf[0]

            represent_Nf_length = len(represent_Nf_seq.string)
            represent_Nf_start_in_array = max_length_Nf - represent_Nf_length
            sub_positions -= represent_Nf_start_in_array
            seq_M = ModifiedSequence(seqs_Nf, tuple(sub_positions))

            represent_Nb_names.append(represent_Nf_seq.name)
            if represent_Nf_seq.has_complete_feature_set:
                for seq_Nf in seqs_Nf:
                    if len(seq_Nf.string) < represent_Nf_length:
                        break
                    excluded_Nf_names.append(seq_Nf.name)

            for seq_Nf in seqs_Nf:
                dict_name_Nb_M[seq_Nf.name].append(seq_M)

            dict_M[seq_M.name] = seq_M
            return True


    def report_sub_stats(self):
        """Report to terminal stats on potential modification-induced substitutions."""
        count_M = len(self.dict_M)
        dict_Nf = self.dict_Nf
        total_sub_count = 0
        total_length_M = 0
        for seq_M in self.dict_M.values():
            length_M = len(dict_Nf[seq_M.name].string)
            total_sub_count += len(seq_M.sub_positions)
            total_length_M += length_M
        mean_sub_per_seq = total_sub_count / count_M
        mean_sub_per_nt = total_sub_count / total_length_M

        warning = self.run.warning
        info = self.run.info

        warning(None, "SUBSTITUTION SEARCH RESULTS", lc='green', nl_before=1)
        info("Modified seqs", count_M)
        info("Mean (*potential*) subs per modified seq", round(mean_sub_per_seq, 1))
        info("Mean subs per nt in modified seq", round(mean_sub_per_nt, 3), nl_after=2)


    def find_indels(self):
        """Find mod-induced indels among Nq normalized seqs not known to be modified, notated *Nq*.
        These seqs form *Ni* objects, which are incorporated into corresponding mod seq objects,
        notated *M*.

        Nq are aligned to N with potential mod-induced subs comprising M, notated *Nb*. Vmatch
        alignments are conducted with Nq as queries and Nb as targets and vice versa, with Nb as
        queries and Nq as targets. Query seqs must be found fully in the target seq. M may be
        shorter than Nq due to unknown, untrimmed 5' and 3' nts in Nq.

        Two pools of Nq are searched for indels:
        1. normalized seqs with a full feature profile and not assigned to M, notated *Nqf*, and
        2. normalized "non-tRNA" seqs with truncated tRNA profiles, notated *Nc*.

        Why these two pools?
        1. Why are Nqf considered at all, when they have been successfully profiled, and therefore,
           presumably, do not have indels interrupting the profile? Indels can be erroneously
           accommodated by flexibility in feature lengths. For example, an indel associated with a
           mod in the D loop can cause the variable-length α or β sections of the D loop to be
           assigned one fewer or one more nt than is correct in order to optimize the profile.
        2. Why not all "non-tRNAs," why just those with a truncated feature profile (Nc)? Indels can
           cause truncation of the profile. "Non-tRNAs" without even a truncated profile (seqs that
           were not profiled past the min length threshold of the T arm) also have fewer
           opportunities for mod-induced mutations.

        Nq seqs rather than constituent T or U are searched for the sake of speed and simplicity.
        Ideally, Ni would be further processed, finding which of their constituent T and U actually
        contain the indels. However, *nonspecific* T and U in Nq are, by definition, in other Nq,
        theoretically permitting the ambiguity that a T would be marked as having an indel in one
        but not another Nq. This would not necessarily be an error, as identical underlying reads
        could theoretically originate from different cDNA seqs, with some containing an indel, and
        others, representing a different tRNA, not containing it."""
        pid = "Finding seqs with mod-induced indels"
        self.progress.new(pid)

        # Write FASTA files of queries and targets to a temp dir used in running Vmatch. Do not
        # allow the Vmatch driver to automatically remove the dir, as the FASTA file of Nq is used
        # in multiple searches. The file of parsed output generated by the Vmatch driver must be
        # removed between searches, as it would otherwise be appended by the next search.
        temp_dir_path = filesnpaths.get_temp_directory_path()
        parsed_output_path = os.path.join(temp_dir_path, 'parsed_output.tsv')
        fasta_path_Nb = os.path.join(temp_dir_path, 'Nb.fa')
        fasta_path_Nqf = os.path.join(temp_dir_path, 'Nqf.fa')

        # Write a FASTA file of Nb.
        self.progress.update_pid(pid)
        self.progress.update("Writing FASTA of norm tRNA seqs with mod-induced subs")
        count_Nb, max_length_M = self.write_fasta_Nb(fasta_path_Nb)

        # Write a FASTA file of Nqf.
        self.progress.update_pid(pid)
        self.progress.update("Writing FASTA of norm tRNA seqs without mod-induced subs")
        count_Nqf, max_length_Nqf = self.write_fasta_Nqf(fasta_path_Nqf)

        # Search Nqf against Nb.
        self.progress.update_pid(pid)
        self.progress.update("Searching for norm tRNA seqs within mod tRNA seqs")
        match_df = Vmatch(argparse.Namespace(match_mode='query_substring_with_indels',
                                             fasta_db_file=fasta_path_Nb,
                                             fasta_query_file=fasta_path_Nqf,
                                             num_threads=self.num_threads,
                                             query_chunk_size=count_Nqf // self.num_threads + 1 if count_Nqf < Vmatch.QUERY_CHUNK_SIZE_DEFAULT else 0,
                                             max_edit_dist=math.ceil(max_length_Nqf * self.max_indel_freq),
                                             min_ident=int(100 - 100 * self.max_indel_freq),
                                             align_output_length=10, # This value is chosen to speed up alignment parsing in Vmatch output.
                                             temp_dir=temp_dir_path,
                                             keep_temp_dir=True,
                                             edit_left_buffer=self.left_indel_buffer,
                                             edit_right_buffer=self.right_indel_buffer)).search_queries()
        self.organize_vmatch_driver_output(match_df, False)
        os.remove(parsed_output_path)

        results_dict = {}
        # The following method updates `results_dict`.
        count_Nqf_with_indels = self.process_Nq_with_indels(match_df, self.dict_Nf, results_dict, False)
        self.progress.end()
        self.run.info_single("Completed indel search Stage 1/4: norm tRNA seqs within mod tRNA seqs", nl_before=2 if self.write_checkpoints else 0)


        # Search Nb against Nqf.
        self.progress.new(pid)
        if count_Nqf_with_indels:
            # Indels were found in some Nqf, so rewrite the FASTA file of Nqf to exclude these.
            self.progress.update("Writing FASTA of norm tRNA seqs without known mod-induced mutations")
            count_Nqf, max_length_Nqf = self.write_fasta_Nqf(fasta_path_Nqf)

        self.progress.update_pid(pid)
        self.progress.update("Searching for mod tRNA seqs within norm tRNA seqs")
        match_df = Vmatch(argparse.Namespace(match_mode='query_substring_with_indels',
                                             fasta_db_file=fasta_path_Nqf,
                                             fasta_query_file=fasta_path_Nb,
                                             num_threads=self.num_threads,
                                             query_chunk_size=count_Nb // self.num_threads + 1 if count_Nb < Vmatch.QUERY_CHUNK_SIZE_DEFAULT else 0,
                                             max_edit_dist=math.ceil(max_length_M * self.max_indel_freq),
                                             min_ident=int(100 - 100 * self.max_indel_freq),
                                             align_output_length=10,
                                             temp_dir=temp_dir_path,
                                             keep_temp_dir=True,
                                             edit_left_buffer=self.left_indel_buffer,
                                             edit_right_buffer=self.right_indel_buffer)).search_queries()
        self.organize_vmatch_driver_output(match_df, True)
        os.remove(parsed_output_path)

        count_Nqf_with_indels += self.process_Nq_with_indels(match_df, self.dict_Nf, results_dict, True)
        self.progress.end()
        self.run.info_single("Completed indel search Stage 2/4: mod tRNA seqs within norm tRNA seqs")


        # Write a FASTA file of Nc.
        self.progress.new(pid)
        self.progress.update("Writing FASTA of norm trunc seqs")
        fasta_path_Nc = os.path.join(temp_dir_path, 'Nc.fa')
        count_Nc, max_length_Nc = self.write_fasta_Nc(fasta_path_Nc)

        # Search Nc against Nb.
        self.progress.update_pid(pid)
        self.progress.update("Searching for trunc tRNA seqs within mod tRNA seqs")
        match_df = Vmatch(argparse.Namespace(match_mode='query_substring_with_indels',
                                             fasta_db_file=fasta_path_Nb,
                                             fasta_query_file=fasta_path_Nc,
                                             num_threads=self.num_threads,
                                             query_chunk_size=count_Nc // self.num_threads + 1 if count_Nc < Vmatch.QUERY_CHUNK_SIZE_DEFAULT else 0,
                                             max_edit_dist=math.ceil(max_length_Nc * self.max_indel_freq),
                                             min_ident=int(100 - 100 * self.max_indel_freq),
                                             align_output_length=10,
                                             temp_dir=temp_dir_path,
                                             keep_temp_dir=True,
                                             edit_left_buffer=self.left_indel_buffer,
                                             edit_right_buffer=self.right_indel_buffer)).search_queries()
        self.organize_vmatch_driver_output(match_df, False)
        os.remove(parsed_output_path)

        count_Nc_with_indels = self.process_Nq_with_indels(match_df, self.dict_Nc, results_dict, False)
        self.progress.end()
        self.run.info_single("Completed indel search Stage 3/4: trunc tRNA seqs within mod tRNA seqs")


        # Search Nb against Nc.
        self.progress.new(pid)
        if count_Nc_with_indels:
            # Indels were found in some Nc, so rewrite the FASTA file of Nc to exclude these.
            self.progress.update("Writing FASTA of norm trunc seqs without known mod-induced mutations")
            count_Nc, max_length_Nc = self.write_fasta_Nc(fasta_path_Nc)

        self.progress.update_pid(pid)
        self.progress.update("Searching for mod tRNA seqs within trunc tRNA seqs")
        match_df = Vmatch(argparse.Namespace(match_mode='query_substring_with_indels',
                                             fasta_db_file=fasta_path_Nc,
                                             fasta_query_file=fasta_path_Nb,
                                             num_threads=self.num_threads,
                                             query_chunk_size=count_Nb // self.num_threads + 1 if count_Nb < Vmatch.QUERY_CHUNK_SIZE_DEFAULT else 0,
                                             max_edit_dist=math.ceil(max_length_M * self.max_indel_freq),
                                             min_ident=int(100 - 100 * self.max_indel_freq),
                                             align_output_length=10,
                                             temp_dir=temp_dir_path,
                                             keep_temp_dir=True,
                                             edit_left_buffer=self.left_indel_buffer,
                                             edit_right_buffer=self.right_indel_buffer)).search_queries()
        self.organize_vmatch_driver_output(match_df, True)
        os.remove(parsed_output_path)

        count_Nc_with_indels += self.process_Nq_with_indels(match_df, self.dict_Nc, results_dict, True)
        self.progress.end()
        self.run.info_single("Completed indel search Stage 4/4: mod tRNA seqs within trunc tRNA seqs")


        # Consolidate Nq differing by 5' and 3' extensions into a new Ni object.
        self.progress.new(pid)
        self.progress.update("Finalizing norm seqs with indels")
        self.add_Ni_to_M(results_dict)
        self.progress.end()

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        set_meta_value = trnaseq_db.db.set_meta_value
        set_meta_value('count_Nqf_with_indels', count_Nqf_with_indels)
        set_meta_value('count_Nc_with_indels', count_Nc_with_indels)
        trnaseq_db.disconnect()


    def write_fasta_Nb(self, fasta_path):
        """This helper method for `find_indels` writes a FASTA file of Nb constituting M."""
        count_Nb = 0
        max_length_M = 0
        dict_Nf = self.dict_Nf
        with open(fasta_path, 'w') as fasta:
            for name_M, seq_M in self.dict_M.items():
                length_M = seq_M.length
                for num_Nb, name_Nb in enumerate(seq_M.names_Nb):
                    string_Nb = dict_Nf[name_Nb].string
                    length_Nb = len(string_Nb)
                    fasta.write(f">{name_M}_{num_Nb}_{length_M - length_Nb}\n{string_Nb}\n")
                    count_Nb += 1
                    if length_Nb > max_length_M:
                        max_length_M = length_Nb
        return count_Nb, max_length_M


    def write_fasta_Nqf(self, fasta_path):
        """This helper method for `find_indels` writes a FASTA file of Nqf."""
        count_Nqf = 0
        max_length_Nqf = 0
        with open(fasta_path, 'w') as fasta:
            for name_Nf, seq_Nf in self.dict_Nf.items():
                if not seq_Nf.names_M:
                    string_Nf = seq_Nf.string
                    length_Nf = len(string_Nf)
                    fasta.write(f">{name_Nf}_{length_Nf}\n{string_Nf}\n")
                    count_Nqf += 1
                    if length_Nf > max_length_Nqf:
                        max_length_Nqf = length_Nf
        return count_Nqf, max_length_Nqf


    def write_fasta_Nc(self, fasta_path):
        """This helper method for `find_indels` writes a FASTA file of Nc."""
        max_length_Nc = 0
        with open(fasta_path, 'w') as fasta:
            for name_Nc, seq_Nc in self.dict_Nc.items():
                string_Nc = seq_Nc.string
                length_Nc = len(string_Nc)
                fasta.write(f">{name_Nc}_{length_Nc}\n{string_Nc}\n")
                if length_Nc > max_length_Nc:
                    max_length_Nc = length_Nc
        return len(self.dict_Nc), max_length_Nc


    def organize_vmatch_driver_output(self, match_df, queries_are_Nb):
        """This helper method for `find_indels` organizes a table of alignment data for further
        analysis."""
        if queries_are_Nb:
            col_name_Nb = 'query_name'
            col_name_N = 'target_name'
        else:
            col_name_Nb = 'target_name'
            col_name_N = 'query_name'

        names_M = []
        starts_Nb_in_M = []
        for defline in match_df[col_name_Nb]:
            split_name = defline.split('_')
            names_M.append('_'.join(split_name[: -2]))
            starts_Nb_in_M.append(int(split_name[-1]))
        match_df['M_name'] = names_M
        match_df['Nb_start_in_M'] = starts_Nb_in_M
        match_df.drop(col_name_Nb, axis=1, inplace=True)

        names_N = []
        lengths_N = []
        for defline in match_df[col_name_N]:
            split_name = defline.split('_')
            names_N.append('_'.join(split_name[: -1]))
            lengths_N.append(int(defline.split('_')[-1]))
        match_df['Nq_name'] = names_N
        match_df['Nq_length'] = lengths_N
        match_df.drop(col_name_N, axis=1, inplace=True)

        if queries_are_Nb:
            match_df.rename({'query_start_in_target': 'Nb_start_in_Nq',
                             'del_lengths': 'insert_lengths',
                             'target_align_del_starts': 'Nq_align_insert_starts',
                             'query_align_del_starts': 'Nb_align_insert_starts',
                             'insert_lengths': 'del_lengths',
                             'target_align_insert_starts': 'Nq_align_del_starts',
                             'query_align_insert_starts': 'Nb_align_del_starts'}, axis=1, inplace=True)
        else:
            match_df.rename({'query_start_in_target': 'Nq_start_in_Nb',
                             'query_align_insert_starts': 'Nq_align_insert_starts',
                             'target_align_insert_starts': 'Nb_align_insert_starts',
                             'query_align_del_starts': 'Nq_align_del_starts',
                             'target_align_del_starts': 'Nb_align_del_starts'}, axis=1, inplace=True)


    def process_Nq_with_indels(self, match_df, dict_N, results_dict, queries_are_Nb):
        """This helper method for `find_indels` finds indels in Nq from parsed Vmatch output. This
        method is generalized to handle Nq (Nqf or Nc) as queries and Nb as targets (Nq length ≤ Nb
        length in alignments) and vice versa (Nb as queries, Nq as targets, Nb length ≤ Nq length in
        alignments).

        A mod-induced insertion in Nq is a gap in Nb (or M), whereas a del in Nq corresponds to nts
        in Nb (or M). The position of an insertion in Nq or a del in M is marked as the first nt of
        the insertion or del. The position of an insertion in M or a del in Nq is marked by the
        position of the adjacent 5' nt."""
        count_of_Nq_with_indels = 0
        dict_M = self.dict_M
        max_length_3prime_terminus = self.max_length_3prime_terminus
        for name_Nq, match_df_Nq in match_df.groupby('Nq_name'):
            names_M = match_df_Nq['M_name']
            if len(set(names_M)) > 1:
                # Ignore Nq with indels that can arise from multiple M for simplicity's sake. This
                # can exclude real molecules with indels. For instance, a mod occurring next to a
                # SNV that distinguishes two M could generate a del that removes the SNV in both M,
                # resulting in a single Nq derived from two M.
                continue

            # Ignore Nq that align to different places in M.
            if queries_are_Nb:
                # When searching Nb against Nq, the start position of Nq in M is ≤ 0.
                starts_Nq_in_M = match_df_Nq['Nb_start_in_M'] - match_df_Nq['Nb_start_in_Nq']
                if len(set(starts_Nq_in_M)) > 1:
                    continue
            else:
                # When searching Nq against Nb, the start position of Nq in M is ≥ 0.
                starts_Nq_in_M = match_df_Nq['Nq_start_in_Nb'] + match_df_Nq['Nb_start_in_M']
                if len(set(starts_Nq_in_M)) > 1:
                    continue

            start_Nq_in_M = starts_Nq_in_M.iat[0]

            insert_length_configs = []
            sum_insert_lengths = []
            for insert_length_config in match_df_Nq['insert_lengths']:
                if insert_length_config:
                    parsed_insert_length_config = tuple(map(int, insert_length_config.split(',')))
                    insert_length_configs.append(parsed_insert_length_config)
                    sum_insert_lengths.append(sum(parsed_insert_length_config))
                else:
                    insert_length_configs.append(tuple())
                    sum_insert_lengths.append(0)
            sum_insert_lengths = np.array(sum_insert_lengths)
            del_length_configs = []
            sum_del_lengths = []
            for del_length_config in match_df_Nq['del_lengths']:
                if del_length_config:
                    parsed_del_length_config = tuple(map(int, del_length_config.split(',')))
                    del_length_configs.append(parsed_del_length_config)
                    sum_del_lengths.append(sum(parsed_del_length_config))
                else:
                    del_length_configs.append(tuple())
                    sum_del_lengths.append(0)
            sum_del_lengths = np.array(sum_del_lengths)

            # If Nq length ≤ M length, Nq must align to the 3' end of M. If Nq length > M length, Nq
            # must have the same number of extra 3' nts in all alignments.
            stops_Nq_in_M = starts_Nq_in_M + match_df_Nq['Nq_length'] - sum_insert_lengths + sum_del_lengths

            if len(set(stops_Nq_in_M)) > 1:
                continue
            stop_Nq_in_M = stops_Nq_in_M.iat[0]

            name_M = names_M.iat[0]
            try:
                results_M_dict = results_dict[name_M]
                seq_M = results_M_dict['M_name']['M_seq']
            except KeyError:
                results_M_dict = None
                seq_M = dict_M[name_M]

            if queries_are_Nb:
                # Ensure that the 3' overhang of Nq in M does not exceed the maximum 3' terminus
                # length allowed in profiling. This heuristic mainly addresses an inconsistency
                # caused by chimeric Nc targets with a 5' part that is tRNA. The Nb query may align
                # to the 5' part of the target, resulting in numerous extra 3' nts in Nc. Were these
                # sequences allowed, this would be the only case in the workflow where the 5' part
                # of a chimera is counted as tRNA.
                if stop_Nq_in_M - seq_M.length > max_length_3prime_terminus:
                    continue
            else:
                if stop_Nq_in_M != seq_M.length:
                    continue

            # Nq is now known to have indels. Their locations are to be determined.
            count_of_Nq_with_indels += 1

            seq_Nq = dict_N.pop(name_Nq)
            if not results_M_dict:
                results_dict[name_M] = results_M_dict = {'M_seq': seq_M}
                results_M_dict['Nq_seqs'] = []
                results_M_dict['Nq_starts_in_M'] = []
                results_M_dict['Nq_stops_in_M'] = []
                results_M_dict['Nq_insert_starts'] = []
                results_M_dict['M_insert_starts'] = []
                results_M_dict['insert_lengths'] = []
                results_M_dict['Nq_del_starts'] = []
                results_M_dict['M_del_starts'] = []
                results_M_dict['del_lengths'] = []

            results_M_dict['Nq_seqs'].append(seq_Nq)
            results_M_dict['Nq_starts_in_M'].append(start_Nq_in_M)
            results_M_dict['Nq_stops_in_M'].append(stop_Nq_in_M)

            indel_configs = []
            if queries_are_Nb:
                align_starts_in_Nq = match_df_Nq['Nb_start_in_Nq']
                align_starts_in_Nb = [0] * len(match_df_Nq)
            else:
                align_starts_in_Nq = [0] * len(match_df_Nq)
                align_starts_in_Nb = match_df_Nq['Nq_start_in_Nb']
            for (align_start_in_Nq,
                 align_start_in_Nb,
                 start_Nb_in_M,
                 align_insert_starts_Nq,
                 align_insert_starts_Nb,
                 insert_lengths,
                 align_del_starts_Nq,
                 align_del_starts_Nb,
                 del_lengths) in zip(align_starts_in_Nq,
                                     align_starts_in_Nb,
                                     match_df_Nq['Nb_start_in_M'],
                                     match_df_Nq['Nq_align_insert_starts'],
                                     match_df_Nq['Nb_align_insert_starts'],
                                     insert_length_configs,
                                     match_df_Nq['Nq_align_del_starts'],
                                     match_df_Nq['Nb_align_del_starts'],
                                     del_length_configs):
                indel_config = []
                if align_insert_starts_Nq:
                    indel_config.append(tuple(
                        [align_insert_start_Nq + align_start_in_Nq for align_insert_start_Nq
                         in map(int, align_insert_starts_Nq.split(','))]))
                    indel_config.append(tuple(
                        [align_insert_start_Nb + align_start_in_Nb + start_Nb_in_M for align_insert_start_Nb
                         in map(int, align_insert_starts_Nb.split(','))]))
                    indel_config.append(insert_lengths)
                else:
                    indel_config.append(tuple())
                    indel_config.append(tuple())
                    indel_config.append(tuple())
                if align_del_starts_Nq:
                    indel_config.append(tuple(
                        [align_del_start_Nq + align_start_in_Nq for align_del_start_Nq
                         in map(int, align_del_starts_Nq.split(','))]))
                    indel_config.append(tuple(
                        [align_del_start_Nb + align_start_in_Nb for align_del_start_Nb
                         in map(int, align_del_starts_Nb.split(','))]))
                    indel_config.append(del_lengths)
                else:
                    indel_config.append(tuple())
                    indel_config.append(tuple())
                    indel_config.append(tuple())
                indel_configs.append(tuple(indel_config))
            indel_configs = list(set(tuple(indel_configs)))

            if len(indel_configs) > 1:
                # Determine the configuration of indels in Nq from the different possibilities
                # suggested by the alignments.
                indel_config = self.select_indel_config(indel_configs, seq_M.sub_positions)
            else:
                # Only one configuration of indels in N was found.
                indel_config = indel_configs.pop()

            results_M_dict['Nq_insert_starts'].append(indel_config[0])
            results_M_dict['M_insert_starts'].append(indel_config[1])
            results_M_dict['insert_lengths'].append(indel_config[2])
            results_M_dict['Nq_del_starts'].append(indel_config[3])
            results_M_dict['M_del_starts'].append(indel_config[4])
            results_M_dict['del_lengths'].append(indel_config[5])

        return count_of_Nq_with_indels


    def select_indel_config(self, indel_configs, sub_positions):
        """This helper method for `process_Nq_with_indels` is called when multiple possible indel
        configurations are found for a given N aligned with M. The configuration that best fits the
        substitution configuration of M is selected. Fit is determined by the distance of indels to
        the closest substitution sites and the positions of indels relative to those substitutions
        -- when choosing between two indel candidates the same distance from a substitution, choose
        the more 5' indel."""
        sum_sub_distances = []
        counts_5prime_indels = []
        # Count the number of indels that occur to the 5' rather than 3' side of the closest sub.
        count_5prime_indels = 0
        for indel_config in indel_configs:
            # Analyze insertions.
            sum_sub_dist = 0
            for insert_start_pos, insert_length in zip(indel_config[2], indel_config[3]):
                insert_midpoint_pos = insert_start_pos + 0.5
                sub_index = bisect_left(sub_positions, insert_midpoint_pos)
                if sub_index == 0:
                    sum_sub_dist += sub_positions[sub_index] - insert_midpoint_pos
                    count_5prime_indels += 1
                elif sub_index == len(sub_positions):
                    sum_sub_dist += insert_midpoint_pos - sub_positions[sub_index - 1]
                else:
                    sub_pos_5prime = sub_positions[sub_index - 1]
                    sub_pos_3prime = sub_positions[sub_index]
                    sub_dist_5prime = insert_midpoint_pos - sub_pos_5prime
                    sub_dist_3prime = sub_pos_3prime - insert_midpoint_pos
                    if sub_dist_5prime >= sub_dist_3prime:
                        sum_sub_dist += sub_dist_3prime
                        count_5prime_indels += 1
                    else:
                        sum_sub_dist += sub_dist_5prime

            # Analyze deletions.
            for del_start_pos, del_length in zip(indel_config[4], indel_config[5]):
                del_midpoint_pos = del_start_pos + (del_length - 1) / 2
                sub_index = bisect_left(sub_positions, del_midpoint_pos)
                if sub_index == 0:
                    sum_sub_dist += sub_positions[sub_index] - del_midpoint_pos
                    count_5prime_indels += 1
                elif sub_index == len(sub_positions):
                    sum_sub_dist += del_midpoint_pos - sub_positions[sub_index - 1]
                else:
                    sub_pos_5prime = sub_positions[sub_index - 1]
                    sub_pos_3prime = sub_positions[sub_index]
                    sub_dist_5prime = del_midpoint_pos - sub_pos_5prime
                    sub_dist_3prime = sub_pos_3prime - del_midpoint_pos
                    if sub_dist_5prime >= sub_dist_3prime:
                        sum_sub_dist += sub_dist_3prime
                        count_5prime_indels += 1
                    else:
                        sum_sub_dist += sub_dist_5prime
            sum_sub_distances.append(sum_sub_dist)
            counts_5prime_indels.append(count_5prime_indels)

        min_sum_sub_dist = min(sum_sub_distances)
        max_5prime_indel_count = -1
        selected_config_index = -1
        config_index = 0
        for sub_sum_dist, count_5prime_indels in zip(sum_sub_distances, counts_5prime_indels):
            if min_sum_sub_dist == sub_sum_dist:
                if count_5prime_indels > max_5prime_indel_count:
                    # If multiple indel configurations somehow have both the same distance of indels
                    # from subs AND the same count of indels to the 5' side of the closest
                    # subs, then the first configuration is selected.
                    selected_config_index = config_index
            config_index += 1

        return indel_configs[selected_config_index]


    def add_Ni_to_M(self, results_dict):
        """This helper method for `find_indels` generates Ni from Nq and adds them to M."""
        dict_Ni_string = {}
        for name_M, results_M_dict in results_dict.items():
            names_Ni, starts_Ni_in_M = self.make_Ni(results_M_dict, dict_Ni_string)
            seq_M = results_M_dict['M_seq']
            seq_M.names_Ni = tuple(names_Ni)
            seq_M.starts_Ni_in_M = tuple(starts_Ni_in_M)


    def make_Ni(self, results_M_dict, dict_Ni_string):
        """Given an M with mapped Nq containing indels, make Ni from Nq.

        Nq may be Nqf or Nc. Each Nqf must contain ≥ 1 Tf, each of which is formed from ≥ 1 Uf.
        Each Nc must contain ≥ 1 Tc, each of which is formed from 1 Uc. Some Nq may contain Tm,
        each of which can only contain 1 Um. Uf/Uc are adjusted for extra 5'/3' nts, forming new Uip
        objects. Uip are pooled and trimmed, forming new Tip. 3' dereplication of Tip strings
        produces clusters of Tip in each new Ni. Ni strings are *already known* from mapping Nq to
        M (Nb). Ni objects can then be instantiated from their component Tip objects. Tim objects,
        derived from Tm adjusted for extra 5'/3' nts, are then added to Ni. Tim membership in Ni
        was known from the relation of Nq to Ni: Tm are part of Nq, and each Tm yields one Tim."""
        # Indels have not been added to M, which only has subs.
        seq_M = results_M_dict['M_seq']
        # From the search against M, Nq may have been found to have extra 5'/3' nts. The following
        # lists store information on Ni created from Nq.
        strings_Ni = []
        xtra_5primes_Nq = []
        xtra_3primes_Nq = []
        starts_Ni_in_M = []
        length_M = seq_M.length

        # Relate Tip strings to lists of Uip objects.
        dict_Uip_in_Tip = defaultdict(list)
        # Relate Ni strings to Tim objects.
        dict_Uim_Tim_in_Ni = defaultdict(list)
        # Relate Ti strings to their start positions in Ni.
        dict_Ti_start_in_Ni = {}
        dict_Ti_strings_from_T_name = defaultdict(list)
        excluded_Ti_strings = []
        for seq_Nq, start_Nq_in_M, stop_Nq_in_M in zip(results_M_dict['Nq_seqs'],
                                                       results_M_dict['Nq_starts_in_M'],
                                                       results_M_dict['Nq_stops_in_M']):
            string_Nq = seq_Nq.string
            length_Nq = len(string_Nq)
            # If Nq has 5' nts overhanging M, the start position of Nq in M is negative.
            xtra_5prime_Nq = -start_Nq_in_M if start_Nq_in_M < 0 else 0
            # If Nq has 3' nts overhanging M, indicating an adjustment is needed in the 3'
            # terminus, the stop position of Nq in M lies beyond the length of M.
            xtra_3prime_Nq = stop_Nq_in_M - length_M if stop_Nq_in_M > length_M else 0
            string_Ni = string_Nq[xtra_5prime_Nq: length_Nq - xtra_3prime_Nq]

            strings_Ni.append(string_Ni)
            xtra_5primes_Nq.append(xtra_5prime_Nq)
            xtra_3primes_Nq.append(xtra_3prime_Nq)
            starts_Ni_in_M.append(start_Nq_in_M + xtra_5prime_Nq)

            self.process_T_in_Nq(seq_Nq,
                                 string_Ni,
                                 xtra_5prime_Nq,
                                 xtra_3prime_Nq,
                                 dict_Ti_start_in_Ni,
                                 dict_Uip_in_Tip,
                                 dict_Uim_Tim_in_Ni,
                                 dict_Ti_strings_from_T_name,
                                 excluded_Ti_strings)

        # 3' dereplicate Tip to form Ni clusters and then Ni objects.
        seqs_Tip = []
        dict_Ti = self.dict_Ti
        dict_Ui = self.dict_Ui
        for string_Tip, seqs_Uip in dict_Uip_in_Tip.items():
            if string_Tip in excluded_Ti_strings:
                # Different Tip strings were found to originate from the same T.
                pass
            else:
                seq_Tip = TrimmedIndelSequence(string_Tip, seqs_Uip)
                seqs_Tip.append(seq_Tip)
                dict_Ti[seq_Tip.name] = seq_Tip
                for seq_Uip in seqs_Uip:
                    dict_Ui[seq_Uip.name] = seq_Uip

        names_Tip = [seq_Tip.name for seq_Tip in seqs_Tip]
        reverse_Tip_strings = [seq_Tip.string[::-1] for seq_Tip in seqs_Tip]
        extras_Tip = [(seq_Tip, dict_Ti_start_in_Ni[seq_Tip.string]) for seq_Tip in seqs_Tip]
        clusters = Dereplicator(names_Tip, reverse_Tip_strings, extras=extras_Tip).prefix_dereplicate()
        names_Ni = []
        dict_Ni = self.dict_Ni
        for cluster in clusters:
            seqs_Tip = []
            starts_Tip_in_Ni = []
            for extra in cluster.member_extras:
                seqs_Tip.append(extra[0])
                starts_Tip_in_Ni.append(extra[1])
            string_Ni = seqs_Tip[0].string
            index_Ni = strings_Ni.index(string_Ni)
            xtra_5prime_Nq = xtra_5primes_Nq[index_Ni]
            # Find the positions of insertions in Ni adjusting for extra 5' nts.
            insert_starts_Ni = tuple([insert_start_Nq - xtra_5prime_Nq for insert_start_Nq in results_M_dict['Nq_insert_starts'][index_Ni]])
            # Find the positions of deletions in Ni adjusting for extra 5' nts.
            del_starts_Ni = tuple([del_start_Nq - xtra_5prime_Nq for del_start_Nq in results_M_dict['Nq_del_starts'][index_Ni]])
            contains_anticodon = self.check_normalized_indel_sequence_for_anticodon(seq_M, starts_Ni_in_M[index_Ni])
            seq_Ni = NormalizedIndelSequence(string_Ni,
                                             seqs_Tip,
                                             starts_Tip_in_Ni,
                                             seq_M.name,
                                             insert_starts_Ni,
                                             tuple(results_M_dict['M_insert_starts'][index_Ni]),
                                             tuple(results_M_dict['insert_lengths'][index_Ni]),
                                             del_starts_Ni,
                                             tuple(results_M_dict['M_del_starts'][index_Ni]),
                                             tuple(results_M_dict['del_lengths'][index_Ni]),
                                             contains_anticodon)
            # Add Tim to Ni.
            seqs_Tim = []
            try:
                for seq_Uim, seq_Tim, start_Tim_in_Ni in dict_Uim_Tim_in_Ni[string_Ni]:
                    if seq_Tim.string not in excluded_Ti_strings:
                        name_Tim = seq_Tim.name
                        seq_Ni.names_T.append(name_Tim)
                        seq_Ni.starts_T_in_N.append(start_Tim_in_Ni)
                        dict_Ui[name_Tim] = seq_Uim
                        dict_Ti[name_Tim] = seq_Tim
            except KeyError:
                # Ni does not contain any Tim.
                pass
            seq_Ni.init(seqs_Tip + seqs_Tim, self.dict_Ui)
            dict_Ni[seq_Ni.name] = seq_Ni
            names_Ni.append(seq_Ni.name)

        return names_Ni, starts_Ni_in_M


    def process_T_in_Nq(self,
                        seq_Nq,
                        string_Ni,
                        xtra_5prime_Nq,
                        xtra_3prime_Nq,
                        dict_Ti_start_in_Ni,
                        dict_Uip_in_Tip,
                        dict_Uim_Tim_in_Ni,
                        dict_Ti_strings_from_T_name,
                        excluded_Ti_strings):
        """Process T in an Nq mapped to M with indels. Generate Ui and Tim objects. Gather
        information to generate Tip objects, which can come from multiple Ni."""
        length_Nq = len(seq_Nq.string)
        dict_Um = self.dict_Um
        try:
            # Nq is Nf.
            categories_T = seq_Nq.categories_T
        except AttributeError:
            # Nq is Nc.
            categories_T = ['Tc_nontrna'] * len(seq_Nq.names_T)
        for name_T, category_T, start_T_in_Nq in zip(seq_Nq.names_T, categories_T, seq_Nq.starts_T_in_N):
            seq_T = getattr(self, 'dict_' + category_T)[name_T]
            string_T = seq_T.string
            # T has 5' nts overhanging M when T is within `xtra_5prime_Nq` nts of the start of Nq.
            # These overhanging nts are trimmed in Ti.
            xtra_5prime_T = xtra_5prime_Nq - start_T_in_Nq if start_T_in_Nq < xtra_5prime_Nq else 0
            type_T = type(seq_T)
            # Extra 3' nts are handled differently in Tm and Tp.
            if type_T is TrimmedMappedSequence:
                # Find the distance from the end of Tm to the end of N.
                delta_3prime = length_Nq - start_T_in_Nq - len(string_T)
                # If Tm has 3' nts overhanging M, the stop position of Tm in M lies beyond the
                # length of M.
                xtra_3prime_T = xtra_3prime_Nq - delta_3prime if delta_3prime < xtra_3prime_Nq else 0
            else:
                # Tp ends at the same point as Nq, so the number of extra 3' nts is the same in
                # each.
                xtra_3prime_T = xtra_3prime_Nq
            string_Ti = string_T[xtra_5prime_T: len(string_T) - xtra_3prime_T]

            # It is required that the same T generate the same Ti in alignments to different M.
            # Conflicting Ti strings from the same T are recorded and their Ti objects later purged.
            strings_Ti_from_T_name = dict_Ti_strings_from_T_name[seq_T.name]
            if len(strings_Ti_from_T_name) == 0:
                strings_Ti_from_T_name.append(string_Ti)
            elif len(strings_Ti_from_T_name) == 1:
                if string_Ti != strings_Ti_from_T_name[0]:
                    excluded_Ti_strings.append(strings_Ti_from_T_name[0])
                    excluded_Ti_strings.append(string_Ti)
                    strings_Ti_from_T_name.append(string_Ti)
            else:
                if string_Ti not in strings_Ti_from_T_name:
                    excluded_Ti_strings.append(string_Ti)
                    strings_Ti_from_T_name.append(string_Ti)

            start_Ti_in_Ni = 0 if start_T_in_Nq < xtra_5prime_Nq else start_T_in_Nq - xtra_5prime_Nq
            dict_Ti_start_in_Ni[string_Ti] = start_Ti_in_Ni

            # Make Ui objects.
            for index_U, name_U in enumerate(seq_T.names_U):
                if type_T == TrimmedMappedSequence:
                    seq_Ui = UniqueIndelSequence(dict_Um[name_U], xtra_3prime_T, xtra_5prime_T)
                    # There is only 1 Um per Tm, so the loop can be terminated after the first
                    # iteration.
                    break

                # Only Up are considered after this point in the loop. Due to 5'/3' adjustments, Up
                # from different Tp and Np can end up in the same Tip. Therefore, initialize Tip
                # after processing all Nq.
                if type_T == TrimmedTruncatedProfileSequence:
                    seq_U = getattr(self, 'dict_Uc_' + seq_T.category)[name_U]
                    # Uc are flush with Tc at the 5' end, so they undergo the same 5' adjustment
                    # when converted into Ui and Ti.
                    length_5prime_U = xtra_5prime_T
                else:
                    seq_U = getattr(self, 'dict_' + seq_T.categories_U[index_U])[name_U]
                    # Adjust the nts beyond the 5' tRNA terminus as needed in Up.
                    length_5prime_U = xtra_5prime_T + seq_U.xtra_5prime_length

                length_3prime_U = xtra_3prime_T + seq_U.length_3prime_terminus
                seq_Uip = UniqueIndelSequence(seq_U, length_3prime_U, length_5prime_U)
                dict_Uip_in_Tip[string_Ti].append(seq_Uip)
            else:
                # This point is reached for Tp.
                continue

            # This point is reached for Tm. Since this is a mapped sequence, Tim is a mirror of Uim
            # and the string is not trimmed.
            seq_Tim = TrimmedIndelSequence(seq_Ui.string, [seq_Ui])
            dict_Uim_Tim_in_Ni[string_Ni].append((seq_Ui, seq_Tim, start_Ti_in_Ni))


    def check_normalized_indel_sequence_for_anticodon(self, seq_M, start_Ni_in_M):
        # Determine whether Ni contains the anticodon.
        seq_Tf = self.dict_Tf[self.dict_Nf[seq_M.name].names_T[0]]
        try:
            anticodon_loop_start = seq_Tf.feature_starts[self.RELATIVE_ANTICODON_LOOP_INDEX]
        except IndexError:
            # The anticodon loop was not reached in the profile.
            return False
        if anticodon_loop_start >= start_Ni_in_M:
            return True
        else:
            return False


    def report_M_stats(self):
        """Report to terminal stats on M. Stats on subs were already reported."""
        spec_read_count_M = 0
        spec_short_5prime_read_count_M = 0
        spec_long_5prime_read_count_M = 0
        mean_spec_cov_M = 0
        mean_nonspec_cov_M = 0
        total_length_M = 0
        max_spec_cov_M = 0
        max_nonspec_cov_M = 0
        max_total_cov_M = 0
        for seq_M in self.dict_M.values():
            spec_read_count_M += seq_M.spec_read_count
            if seq_M.spec_read_xtra_5prime_count:
                if seq_M.spec_long_5prime_extension_dict:
                    spec_long_5prime_read_count = sum(seq_M.spec_long_5prime_extension_dict.values())
                    spec_long_5prime_read_count_M += spec_long_5prime_read_count
                    spec_short_5prime_read_count_M += seq_M.spec_read_xtra_5prime_count - spec_long_5prime_read_count
                else:
                    spec_short_5prime_read_count_M += seq_M.spec_read_xtra_5prime_count
            if seq_M.nonspec_read_xtra_5prime_count:
                if seq_M.nonspec_long_5prime_extension_dict:
                    nonspec_long_5prime_read_count = sum(seq_M.nonspec_long_5prime_extension_dict.values())
            length_M = len(seq_M.consensus_string)
            mean_spec_cov = seq_M.mean_spec_cov
            mean_nonspec_cov = seq_M.mean_nonspec_cov
            mean_spec_cov_M += mean_spec_cov * length_M
            mean_nonspec_cov_M += mean_nonspec_cov * length_M
            total_length_M += length_M
            if mean_spec_cov > max_spec_cov_M:
                max_spec_cov_M = mean_spec_cov
            if mean_nonspec_cov > max_nonspec_cov_M:
                max_nonspec_cov_M = mean_nonspec_cov
            if mean_spec_cov + mean_nonspec_cov > max_total_cov_M:
                max_total_cov_M = mean_spec_cov + mean_nonspec_cov
        mean_spec_cov_M /= total_length_M
        mean_nonspec_cov_M /= total_length_M

        if not self.skip_indel_profiling:
            count_indel_M = 0
            count_insert_M = 0
            count_del_M = 0
            for seq_M in self.dict_M.values():
                if seq_M.insert_starts:
                    count_insert_M += 1
                if seq_M.del_starts:
                    count_del_M += 1
                if seq_M.insert_starts or seq_M.del_starts:
                    count_indel_M += 1

        warning = self.run.warning
        info = self.run.info

        warning(None, "MODIFICATION ANALYSIS RESULTS", lc='green', nl_before=1)

        warning(None, "Modified seqs", lc='cyan')
        info("Specific reads", spec_read_count_M)
        info(f"Specific reads with 1-{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_short_5prime_read_count_M)
        info(f"Specific reads with ≥{MIN_LENGTH_LONG_5PRIME_EXTENSION} extra 5' nts", spec_long_5prime_read_count_M)
        info("Mean specific coverage", round(mean_spec_cov_M, 2))
        info("Mean nonspecific coverage", round(mean_nonspec_cov_M, 2))
        info("Max specific coverage", round(max_spec_cov_M, 2))
        info("Max nonspecific coverage", round(max_nonspec_cov_M, 2))
        info("Max total coverage", round(max_total_cov_M, 2), nl_after=2 if self.skip_indel_profiling else 0)

        if not self.skip_indel_profiling:
            warning(None, "Results of indel search", lc='cyan')
            info("Modified seqs with indels", count_indel_M)
            info("Modified seqs with insertions", count_insert_M)
            info("Modified seqs with deletions", count_del_M)

            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            get_meta_value = trnaseq_db.db.get_meta_value
            count_Nqf_with_indels = get_meta_value('count_Nqf_with_indels')
            count_Nc_with_indels = get_meta_value('count_Nc_with_indels')
            trnaseq_db.disconnect()
            info("Normalized tRNA seqs found to have indels", count_Nqf_with_indels)
            info("Normalized trunc seqs found to have indels", count_Nc_with_indels, nl_after=2)


    def report_stats(self):
        """Add final stats to the db and write them to a summary file."""
        # Define Ntrna to encompass Nf and Ni.
        anticodon_profiled_trna_reads = 0
        complete_profiled_trna_reads = 0
        spec_reads_Ntrna = 0
        nonspec_reads_Ntrna = 0

        for seq_Tf in self.dict_Tf.values():
            if seq_Tf.contains_anticodon:
                anticodon_profiled_trna_reads += seq_Tf.read_count
                if seq_Tf.has_complete_feature_set:
                    complete_profiled_trna_reads += seq_Tf.read_count
            if len(seq_Tf.names_N) == 1:
                spec_reads_Ntrna += seq_Tf.read_count
            else:
                nonspec_reads_Ntrna += seq_Tf.read_count

        for seq_Tc_trna in self.dict_Tc_trna.values():
            if seq_Tc_trna.contains_anticodon:
                anticodon_profiled_trna_reads += seq_Tc_trna.read_count
            if len(seq_Tc_trna.names_N) == 1:
                spec_reads_Ntrna += seq_Tc_trna.read_count
            else:
                nonspec_reads_Ntrna += seq_Tc_trna.read_count
        profiled_trna_reads = spec_reads_Ntrna + nonspec_reads_Ntrna

        for seq_Tm in self.dict_Tm.values():
            if len(seq_Tm.names_N) == 1:
                spec_reads_Ntrna += seq_Tm.read_count
            else:
                nonspec_reads_Ntrna += seq_Tm.read_count

        for seq_Ti in self.dict_Ti.values():
            if len(seq_Ti.names_N) == 1:
                spec_reads_Ntrna += seq_Ti.read_count
            else:
                nonspec_reads_Ntrna += seq_Ti.read_count
        trna_reads = spec_reads_Ntrna + nonspec_reads_Ntrna

        # Define Nfu to be Nf not in M.
        count_M = len(self.dict_M)
        count_Nfu = 0
        for seq_Nf in self.dict_Nf.values():
            if not seq_Nf.names_M:
                count_Nfu += 1

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        set_meta_value = trnaseq_db.db.set_meta_value
        set_meta_value('trna_reads', trna_reads)
        set_meta_value('profiled_trna_reads', profiled_trna_reads)
        set_meta_value('anticodon_profiled_trna_reads', anticodon_profiled_trna_reads)
        set_meta_value('complete_profiled_trna_reads', complete_profiled_trna_reads)
        set_meta_value('spec_reads_Ntrna', spec_reads_Ntrna)
        set_meta_value('nonspec_reads_Ntrna', nonspec_reads_Ntrna)
        set_meta_value('count_M', count_M)
        set_meta_value('count_Nfu', count_Nfu)

        get_summary_line = self.get_summary_line
        get_meta_value = trnaseq_db.db.get_meta_value
        with open(self.analysis_summary_path, 'a') as f:
            f.write(get_summary_line("Input reads", get_meta_value('input_reads')))
            f.write(get_summary_line("Input uniq seqs", get_meta_value('input_U')))
            f.write(get_summary_line("tRNA reads", trna_reads))
            f.write(get_summary_line("Profiled tRNA reads", profiled_trna_reads))
            f.write(get_summary_line("Profiled tRNA reads with anticodon", anticodon_profiled_trna_reads))
            f.write(get_summary_line("Profiled reads with complete features", complete_profiled_trna_reads))
            f.write(get_summary_line("Reads specific to normalized tRNA seqs", spec_reads_Ntrna))
            f.write(get_summary_line("Reads nonspecific to normalized tRNA seqs", nonspec_reads_Ntrna))
            f.write(get_summary_line("Potentially modified seqs", count_M))
            f.write(get_summary_line("Normalized tRNA seqs without detected modifications", count_Nfu))

        trnaseq_db.disconnect()

        self.run.info("Summary", self.analysis_summary_path, nl_before=2)


    def write_feature_table(self):
        self.progress.new("Writing tRNA-seq db table of profiled tRNA features")
        self.progress.update("...")

        rows = []
        for dict_name, dict_U in zip(('dict_Uf', 'dict_Us', 'dict_Uc_trna'),
                                     (self.dict_Uf, self.dict_Us, self.dict_Uc_trna)):
            for seq_U in dict_U.values():
                if dict_name == 'dict_Uc_trna':
                    has_complete_feature_set = False
                    num_extrap_5prime_nts = 0
                    xtra_5prime_length = 0
                else:
                    has_complete_feature_set = seq_U.has_complete_feature_set
                    num_extrap_5prime_nts = seq_U.num_extrap_5prime_nts
                    xtra_5prime_length = seq_U.xtra_5prime_length
                length_U = len(seq_U.string)

                row = (
                    (seq_U.name,
                     has_complete_feature_set,
                     seq_U.anticodon_string,
                     seq_U.anticodon_aa,
                     length_U,
                     length_U - seq_U.profiled_seq_length,
                     seq_U.num_conserved,
                     seq_U.num_unconserved,
                     seq_U.num_paired,
                     seq_U.num_unpaired,
                     num_extrap_5prime_nts,
                     xtra_5prime_length,
                     seq_U.length_3prime_terminus)
                    # When tRNA features are not found at the 5' end of the read, the start and stop
                    # positions of these features also are not found.
                    + tuple([None for _ in range((len(TRNA_FEATURE_NAMES) - len(seq_U.feature_starts)))]) * 2
                    + tuple(chain(*zip(
                        [str(start) if isinstance(start, int)
                         else ','.join(map(str, start))
                         for start in seq_U.feature_starts],
                        # Convert pythonic stop position to real stop position of feature.
                        [str(stop - 1) if isinstance(stop, int)
                         else ','.join(str(strand_stop - 1) for strand_stop in stop)
                         for stop in seq_U.feature_stops])))
                    # The α and β sections of the D loop are "subfeatures," not "features," so add them
                    # to the row after the features.
                    + (seq_U.alpha_start,
                       seq_U.alpha_stop - 1 if seq_U.alpha_stop else None,
                       seq_U.beta_start,
                       seq_U.beta_stop - 1 if seq_U.beta_stop else None)
                )
                rows.append(row)

        if rows:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many(
                f'''INSERT INTO feature VALUES ({','.join('?' * len(tables.trnaseq_feature_table_structure))})''',
                rows)
            trnaseq_db.disconnect()

        self.progress.end()


    def write_unconserved_table(self):
        self.progress.new("Writing tRNA-seq db table of unconserved nts in fully profiled tRNA")
        self.progress.update("...")

        rows = []
        for seq_Uf in self.dict_Uf.values():
            for unconserved_tuple in seq_Uf.unconserved_info:
                rows.append((seq_Uf.name, ) + unconserved_tuple)

        if rows:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many(
                f'''INSERT INTO unconserved VALUES ({','.join('?' * len(tables.trnaseq_unconserved_table_structure))})''',
                rows)
            trnaseq_db.disconnect()

        self.progress.end()


    def write_unpaired_table(self):
        self.progress.new("Writing tRNA-seq db table of unpaired nts in fully profiled tRNA")
        self.progress.update("...")

        rows = []
        for seq_Uf in self.dict_Uf.values():
            for unpaired_tuple in seq_Uf.unpaired_info:
                rows.append((seq_Uf.name, ) + unpaired_tuple)

        if rows:
            trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
            trnaseq_db.db._exec_many(
                f'''INSERT INTO unpaired VALUES ({','.join('?' * len(tables.trnaseq_unpaired_table_structure))})''',
                rows)
            trnaseq_db.disconnect()

        self.progress.end()


    def write_sequences_table(self):
        self.progress.new("Writing tRNA-seq db table of unique tRNA seqs")
        self.progress.update("...")

        rows = []
        for info_U, dict_U in zip(('full_profile', 'transferred_profile', 'truncated_profile', 'mapped', 'indel_aligned'),
                                  (self.dict_Uf, self.dict_Us, self.dict_Uc_trna, self.dict_Um, self.dict_Ui)):
            for seq_U in dict_U.values():
                rows.append((seq_U.name, seq_U.read_count, info_U, seq_U.string))

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('sequences')
            trnaseq_db.db.create_table('sequences',
                                       tables.trnaseq_sequences_table_structure,
                                       tables.trnaseq_sequences_table_types)
        if rows:
            trnaseq_db.db._exec_many(
                f'''INSERT INTO sequences VALUES ({','.join('?' * len(tables.trnaseq_sequences_table_structure))})''',
                rows)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_trimmed_table(self):
        self.progress.new("Writing tRNA-seq db table of trimmed tRNA seqs")
        self.progress.update("...")

        rows = []
        for info_T, dict_T in zip(('full_profile', 'truncated_profile', 'mapped', 'indel_aligned'),
                                  (self.dict_Tf, self.dict_Tc_trna, self.dict_Tm, self.dict_Ti)):
            for seq_T in dict_T.values():
                termini_3prime = ''
                read_counts_3prime_termini = ''
                if info_T != 'mapped':
                    for string_3prime_terminus, read_count in seq_T.read_3prime_terminus_count_dict.items():
                        termini_3prime += string_3prime_terminus + ','
                        read_counts_3prime_termini += str(read_count) + ','

                if info_T != 'truncated_profile':
                    uniq_with_xtra_5prime_count = seq_T.uniq_with_xtra_5prime_count
                    read_with_xtra_5prime_count = seq_T.read_with_xtra_5prime_count
                else:
                    uniq_with_xtra_5prime_count = 0
                    read_with_xtra_5prime_count = 0

                rows.append(
                    (seq_T.name,
                     len(seq_T.names_U),
                     seq_T.read_count,
                     info_T,
                     seq_T.string,
                     len(seq_T.names_N),
                     uniq_with_xtra_5prime_count,
                     read_with_xtra_5prime_count,
                     termini_3prime,
                     read_counts_3prime_termini)
                )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('trimmed')
            trnaseq_db.db.create_table('trimmed',
                                       tables.trnaseq_trimmed_table_structure,
                                       tables.trnaseq_trimmed_table_types)
        if rows:
            trnaseq_db.db._exec_many(
                f'''INSERT INTO trimmed VALUES ({','.join('?' * len(tables.trnaseq_trimmed_table_structure))})''',
                rows)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_normalized_table(self):
        self.progress.new("Writing tRNA-seq db table of fragment-dereplicated (\"normalized\") tRNA seqs")
        self.progress.update("...")

        rows = []
        for info_N, dict_N in zip(('full_profile', 'indel_aligned'),
                                  (self.dict_Nf, self.dict_Ni)):
            for seq_N in dict_N.values():
                spec_long_5prime_extensions = ''
                spec_long_5prime_read_counts = ''
                for string_5prime, read_count in sorted(seq_N.spec_long_5prime_extension_dict.items(),
                                                        key=lambda item: -len(item[0])):
                    spec_long_5prime_extensions += string_5prime + ','
                    spec_long_5prime_read_counts += str(read_count) + ','

                nonspec_long_5prime_extensions = ''
                nonspec_long_5prime_read_counts = ''
                for string_5prime, read_count in sorted(seq_N.nonspec_long_5prime_extension_dict.items(),
                                                        key=lambda item: -len(item[0])):
                    nonspec_long_5prime_extensions += string_5prime + ','
                    nonspec_long_5prime_read_counts += str(read_count) + ','

                spec_3prime_termini = ''
                spec_3prime_terminus_read_counts = ''
                for string_3prime_terminus, read_count in seq_N.spec_read_3prime_terminus_count_dict.items():
                    spec_3prime_termini += string_3prime_terminus + ','
                    spec_3prime_terminus_read_counts += str(read_count) + ','

                nonspec_3prime_termini = ''
                nonspec_3prime_terminus_read_counts = ''
                for string_3prime_terminus, read_count in seq_N.nonspec_read_3prime_terminus_count_dict.items():
                    nonspec_3prime_termini += string_3prime_terminus + ','
                    nonspec_3prime_terminus_read_counts += str(read_count) + ','

                rows.append(
                    (seq_N.name,
                     len(seq_N.names_T),
                     info_N,
                     seq_N.mean_spec_cov,
                     seq_N.mean_nonspec_cov,
                     ','.join(map(str, seq_N.spec_covs)) + ',',
                     ','.join(map(str, seq_N.nonspec_covs)) + ',',
                     seq_N.spec_read_count,
                     seq_N.nonspec_read_count,
                     seq_N.spec_read_xtra_5prime_count,
                     seq_N.nonspec_read_xtra_5prime_count,
                     seq_N.spec_map_read_count,
                     seq_N.nonspec_map_read_count,
                     spec_long_5prime_extensions,
                     spec_long_5prime_read_counts,
                     nonspec_long_5prime_extensions,
                     nonspec_long_5prime_read_counts,
                     spec_3prime_termini,
                     spec_3prime_terminus_read_counts,
                     nonspec_3prime_termini,
                     nonspec_3prime_terminus_read_counts)
                )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('normalized')
            trnaseq_db.db.create_table('normalized',
                                       tables.trnaseq_normalized_table_structure,
                                       tables.trnaseq_normalized_table_types)
        if rows:
            trnaseq_db.db._exec_many(
                f'''INSERT INTO normalized VALUES ({','.join('?' * len(tables.trnaseq_normalized_table_structure))})''',
                rows)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_modified_table(self):
        self.progress.new("Writing tRNA-seq db table of potentially modified tRNA seqs")
        self.progress.update("...")

        rows = []
        for seq_M in self.dict_M.values():
            spec_long_5prime_extensions = ''
            spec_long_5prime_read_counts = ''
            for string_5prime, read_count in sorted(seq_M.spec_long_5prime_extension_dict.items(),
                                                    key=lambda item: -len(item[0])):
                spec_long_5prime_extensions += string_5prime + ','
                spec_long_5prime_read_counts += str(read_count) + ','

            nonspec_long_5prime_extensions = ''
            nonspec_long_5prime_read_counts = ''
            for string_5prime, read_count in sorted(seq_M.nonspec_long_5prime_extension_dict.items(),
                                                    key=lambda item: -len(item[0])):
                nonspec_long_5prime_extensions += string_5prime + ','
                nonspec_long_5prime_read_counts += str(read_count) + ','

            spec_3prime_termini = ''
            spec_3prime_terminus_read_counts = ''
            for string_3prime_terminus, read_count in seq_M.spec_read_3prime_terminus_count_dict.items():
                spec_3prime_termini += string_3prime_terminus + ','
                spec_3prime_terminus_read_counts += str(read_count) + ','

            nonspec_3prime_termini = ''
            nonspec_3prime_terminus_read_counts = ''
            for string_3prime_terminus, read_count in seq_M.nonspec_read_3prime_terminus_count_dict.items():
                nonspec_3prime_termini += string_3prime_terminus + ','
                nonspec_3prime_terminus_read_counts += str(read_count) + ','

            rows.append(
                (seq_M.name,
                 seq_M.mean_spec_cov,
                 seq_M.mean_nonspec_cov,
                 ','.join(map(str, seq_M.spec_covs)) + ',',
                 ','.join(map(str, seq_M.nonspec_covs)) + ',',
                 ','.join([str(sub_pos) for sub_pos in seq_M.sub_positions]) + ',')
                + tuple([','.join(map(str, seq_M.spec_sub_covs[:, nt_index - 1])) + ',' for nt_index in INT_NT_DICT])
                + tuple([','.join(map(str, seq_M.nonspec_sub_covs[:, nt_index - 1])) + ',' for nt_index in INT_NT_DICT])
                + (','.join([str(insert_start) for insert_start in seq_M.insert_starts]) + ',',
                   ','.join([insert_string for insert_string in seq_M.insert_strings]) + ',',
                   ','.join([str(spec_insert_cov) for spec_insert_cov in seq_M.spec_insert_covs]) + ',',
                   ','.join([str(nonspec_insert_cov) for nonspec_insert_cov in seq_M.nonspec_insert_covs]) + ',',
                   ','.join([str(del_start) for del_start in seq_M.del_starts]) + ',',
                   ','.join([str(del_length) for del_length in seq_M.del_lengths]) + ',',
                   ','.join([str(spec_del_cov) for spec_del_cov in seq_M.spec_del_covs]) + ',',
                   ','.join([str(nonspec_del_cov) for nonspec_del_cov in seq_M.nonspec_del_covs]) + ',',
                   seq_M.consensus_string,
                   len(seq_M.names_Nb),
                   ','.join(seq_M.names_Nb),
                   len(seq_M.names_Ni),
                   ','.join(seq_M.names_Ni),
                   seq_M.spec_read_count,
                   seq_M.nonspec_read_count,
                   seq_M.spec_read_xtra_5prime_count,
                   seq_M.nonspec_read_xtra_5prime_count,
                   seq_M.spec_map_read_count,
                   seq_M.nonspec_map_read_count,
                   spec_long_5prime_extensions,
                   spec_long_5prime_read_counts,
                   nonspec_long_5prime_extensions,
                   nonspec_long_5prime_read_counts,
                   spec_3prime_termini,
                   spec_3prime_terminus_read_counts,
                   nonspec_3prime_termini,
                   nonspec_3prime_terminus_read_counts)
            )

        trnaseq_db = dbops.TRNASeqDatabase(self.trnaseq_db_path, quiet=True)
        # Overwrite the existing table if starting from a checkpoint.
        if self.load_checkpoint:
            trnaseq_db.db.drop_table('modified')
            trnaseq_db.db.create_table('modified',
                                       tables.trnaseq_modified_table_structure,
                                       tables.trnaseq_modified_table_types)
        if rows:
            trnaseq_db.db._exec_many(
                f'''INSERT INTO modified VALUES ({','.join('?' * len(tables.trnaseq_modified_table_structure))})''',
                rows)
        trnaseq_db.disconnect()

        self.progress.end()


    def write_nontrna_supplement(self):
        """Write a supplementary file on Un and Uc_nontrna."""
        self.progress.new("Writing file of unique seqs not identified as tRNA")
        self.progress.update("...")

        with open(self.path_Un_supplement, 'w') as file_Un_supplement:
            file_Un_supplement.write("\t".join(self.UNIQ_NONTRNA_HEADER) + "\n")
            for seq_Un in self.dict_Un.values():
                file_Un_supplement.write(seq_Un.name + "\t"
                                         + str(seq_Un.read_count) + "\t"
                                         + "\t"
                                         + seq_Un.string + "\n")

            for seq_Uc_nontrna in self.dict_Uc_nontrna.values():
                file_Un_supplement.write(seq_Uc_nontrna.name + "\t"
                                         + str(seq_Uc_nontrna.read_count) + "\t"
                                         + str(seq_Uc_nontrna.trunc_profile_index) + "\t"
                                         + seq_Uc_nontrna.string + "\n")

        self.progress.end()

        self.run.info("Unique non-tRNA supplement", self.path_Un_supplement)


    def write_Tf_ends_supplement(self):
        """Write a supplementary file showing the spectrum of 5'/3' extensions of Tf."""
        self.progress.new("Writing file showing 5'/3' ends of trimmed, fully profiled tRNA seqs")
        self.progress.update("...")

        with open(self.path_Tf_ends, 'w') as file_Tf_ends:
            file_Tf_ends.write("\t".join(self.TRIMMED_ENDS_HEADER) + "\n")
            for seq_Tf in sorted(self.dict_Tf.values(), key=lambda seq_Tf: -seq_Tf.read_count):
                name_Tf = seq_Tf.name
                seqs_U = [getattr(self, 'dict_' + category_U)[name_U] for category_U, name_U in zip(seq_Tf.categories_U, seq_Tf.names_U)]
                for seq_U in sorted(seqs_U, key=lambda seq_U: (-seq_U.xtra_5prime_length, -seq_U.length_3prime_terminus)):
                    file_Tf_ends.write(name_Tf + "\t"
                                       + seq_U.name + "\t"
                                       + seq_U.string[: seq_U.xtra_5prime_length] + "\t"
                                       + seq_U.string[len(seq_U.string) - seq_U.length_3prime_terminus: ] + "\t"
                                       + str(seq_U.read_count) + "\n")

        self.progress.end()

        self.run.info("Trimmed tRNA supplement", self.path_Tf_ends)


def profile_worker(input_queue, output_queue, profiler):
    """This client for `trnaidentifier.Profiler.profile` is located outside the `TRNASeqDataset`
    class to allow multiprocessing."""
    while True:
        seq_string, represent_name, read_count = input_queue.get()
        output_queue.put((profiler.profile(seq_string, name=represent_name), read_count))


class NormalizedSequenceSummary(object):
    """Relevant data from a normalized sequence stored in a tRNA-seq database is reloaded into this
    object."""

    __slots__ = (
        'name',
        'sample_id',
        'string',
        'anticodon_string',
        'feature_dict',
        'feature_threshold_start',
        'mean_spec_cov',
        'spec_covs',
        'nonspec_covs',
        'spec_nt_covs_dict',
        'nonspec_nt_covs_dict',
        'summary_M'
    )

    def __init__(self):
        for attr_name in NormalizedSequenceSummary.__slots__:
            setattr(self, attr_name, None)


class NormalizedIndelSequenceSummary(NormalizedSequenceSummary):
    """Relevant data from a normalized sequence with indels stored in a tRNA-seq database is
    reloaded into this object."""

    __slots__ = (
        'insert_starts',
        'insert_strings',
        'spec_insert_covs',
        'nonspec_insert_covs',
        'del_starts',
        'del_lengths',
        'spec_del_covs',
        'nonspec_del_covs'
    )

    def __init__(self):
        super().__init__()
        for attr_name in NormalizedIndelSequenceSummary.__slots__:
            setattr(self, attr_name, None)


class ModifiedSequenceSummary(object):
    """Relevant data from a modified sequence stored in a tRNA-seq database is reloaded into this
    object."""

    __slots__ = (
        'name',
        'sample_id',
        'consensus_string',
        'sub_positions',
        'spec_nt_covs_dict',
        'nonspec_nt_covs_dict',
        'spec_covs',
        'nonspec_covs',
        'insert_starts',
        'insert_strings',
        'spec_insert_covs',
        'nonspec_insert_covs',
        'del_starts',
        'del_lengths',
        'spec_del_covs',
        'nonspec_del_covs',
        'summaries_Nb',
        'summaries_Ni'
    )

    def __init__(self):
        for attr_name in ModifiedSequenceSummary.__slots__:
            setattr(self, attr_name, None)


class SeedSequence(object):

    __slots__ = (
        'name',
        'string',
        'meets_feature_threshold',
        'summaries_Nu',
        'summaries_M',
        'anticodon_string',
        'feature_dict',
        'total_spec_covs',
        'total_nonspec_covs',
        'total_mean_spec_cov',
        'total_mean_nonspec_cov',
        'sample_spec_covs_dict',
        'sample_nonspec_covs_dict',
        'sample_summed_covs_dict',
        'sample_spec_nt_covs_dict',
        'sample_nonspec_nt_covs_dict',
        'sample_summed_nt_covs_dict',
        'sample_mean_spec_cov_dict',
        'sample_mean_nonspec_cov_dict',
        'sample_mean_summed_cov_dict',
        'sample_std_spec_cov_dict',
        'sample_std_nonspec_cov_dict',
        'sample_std_summed_cov_dict',
        'sample_spec_abund_dict',
        'sample_nonspec_abund_dict',
        'sample_summed_abund_dict',
        'sample_spec_rel_abund_dict',
        'sample_nonspec_rel_abund_dict',
        'sample_summed_rel_abund_dict',
        'sample_spec_detection_dict',
        'sample_nonspec_detection_dict',
        'sample_summed_detection_dict',
        'sample_mean_Q2Q3_spec_cov_dict',
        'sample_mean_Q2Q3_nonspec_cov_dict',
        'sample_mean_Q2Q3_summed_cov_dict',
        'sample_normalized_mean_Q2Q3_spec_cov_dict',
        'sample_normalized_mean_Q2Q3_nonspec_cov_dict',
        'sample_normalized_mean_Q2Q3_summed_cov_dict',
        'sample_spec_max_normalized_ratio_dict',
        'sample_nonspec_max_normalized_ratio_dict',
        'sample_summed_max_normalized_ratio_dict',
        'gc_fraction',
        'sample_sub_positions_dict',
        'total_sub_positions',
        'sample_variability_dict',
        'sample_insert_dict',
        'sample_del_dict'
    )

    def __init__(self):
        for attr_name in SeedSequence.__slots__:
            setattr(self, attr_name, None)


class DatabaseMerger(object):
    """Merges tRNA-seq database(s) into contigs, auxiliary, and profile databases. "Contigs" in
    this context are tRNA seed sequences representing tRNA identified in the samples."""

    # The following constants are columns needed from tables of a tRNA-seq database.
    # Load all feature positional indices but the 3' terminus start and stop indices.
    FEATURE_INDEX_COLS_OF_INTEREST = list(chain(*zip([f + '_start' for f in TRNA_FEATURE_NAMES[: -1]],
                                                     [f + '_stop' for f in TRNA_FEATURE_NAMES[: -1]]))) + ['alpha_start', 'alpha_stop', 'beta_start', 'beta_stop']
    FEATURE_TABLE_COLS_OF_INTEREST = [
        'name',
        'anticodon_sequence',
        'num_extra_fiveprime',
    ] + FEATURE_INDEX_COLS_OF_INTEREST
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
        'names_of_normalized_seqs_without_indels',
        'names_of_normalized_seqs_with_indels',
        'substitution_positions',
        'substitution_A_specific_coverage',
        'substitution_C_specific_coverage',
        'substitution_G_specific_coverage',
        'substitution_T_specific_coverage',
        'substitution_A_nonspecific_coverage',
        'substitution_C_nonspecific_coverage',
        'substitution_G_nonspecific_coverage',
        'substitution_T_nonspecific_coverage',
        'insertion_starts',
        'insertion_seqs',
        'insertion_specific_coverages',
        'insertion_nonspecific_coverages',
        'deletion_starts',
        'deletion_lengths',
        'deletion_specific_coverages',
        'deletion_nonspecific_coverages',
        'consensus_sequence'
    ]

    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        # Argument group A: MANDATORY
        self.trnaseq_db_paths = A('input')
        self.out_dir = A('output_dir')
        self.project_name = A('project_name')

        # Argument group B: EXTRAS
        self.num_threads = A('num_threads')
        self.seed_limit = A('max_reported_trna_seeds')
        self.overwrite_out_dest = A('overwrite_output_destinations')
        self.descrip_path = os.path.abspath(A('description')) if A('description') else None

        # Argument group C: ADVANCED
        self.feature_threshold = A('feature_threshold')
        self.preferred_treatment = A('preferred_treatment')
        self.nonspec_output = A('nonspecific_output')
        self.min_variation = A('min_variation')
        self.min_third_fourth_nt = A('min_third_fourth_nt')
        self.min_indel_fraction = A('min_indel_fraction')
        self.distance = A('distance') or constants.distance_metric_default
        self.linkage = A('linkage') or constants.linkage_method_default

        if not self.project_name:
            raise ConfigError("Please specify a name for the collection of input tRNA-seq dbs using --project-name or -n.")
        if not self.out_dir:
            raise ConfigError("Please provide an output directory using --output-dir or -o.")

        self.contigs_db_path = None
        self.contigs_db_hash = None

        self.spec_out_dir = None
        self.spec_profile_db_path = None
        self.spec_auxiliary_db_path = None

        self.nonspec_out_dir = None
        self.nonspec_profile_db_path = None
        self.nonspec_auxiliary_db_path = None
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
        self.dict_summaries_Nu = defaultdict(list)
        self.dict_summaries_M = defaultdict(list)

        self.sample_total_mean_spec_cov_dict = None
        self.sample_total_discriminator_spec_cov_dict = None

        self.seeds = None
        self.total_seed_length = None

        self.sample_total_spec_cov_dict = None
        self.sample_total_nonspec_cov_dict = None
        self.sample_total_summed_cov_dict = None
        self.sample_overall_mean_spec_cov_dict = None
        self.sample_mean_nonspec_cov_dict = None
        self.sample_mean_summed_cov_dict = None
        self.sample_normalization_multiplier_dict = None

        self.overall_mean_spec_cov = None
        self.overall_mean_nonspec_cov = None

        self.variable_nts_table_entries = None
        self.spec_indels_table_entries = None
        self.nonspec_indels_table_entries = None
        self.summed_indels_table_entries = None


    def process(self):
        """Orchestrate the steps needed to create contigs, profile and auxiliary databases."""
        self.sanity_check()
        filesnpaths.gen_output_directory(self.out_dir, delete_if_exists=self.overwrite_out_dest)

        self.load_trnaseq_dbs()

        self.form_seeds()

        filesnpaths.gen_output_directory(self.spec_out_dir)
        self.generate_contigs_database()
        self.generate_auxiliary_database('specific')

        self.set_sample_total_coverages()
        self.set_sample_overall_mean_coverages()
        self.set_sample_mean_coverages()
        self.set_sample_coverage_standard_deviations()
        self.set_sample_abundances()
        self.set_sample_normalization_multipliers()
        self.set_sample_normalized_mean_Q2Q3_coverages()
        self.set_sample_detections()
        self.set_sample_relative_abundances()
        self.set_sample_max_normalized_ratios()
        self.set_variable_nucleotides_table_entries()
        self.set_indels_table_entries()
        self.generate_profile_database('specific')

        if self.nonspec_out_dir:
            filesnpaths.gen_output_directory(self.nonspec_out_dir, delete_if_exists=self.overwrite_out_dest)
            self.generate_auxiliary_database('nonspecific')
            self.generate_profile_database('nonspecific')
        if self.combined_out_dir:
            filesnpaths.gen_output_directory(self.combined_out_dir, delete_if_exists=self.overwrite_out_dest)
            self.generate_auxiliary_database('combined')
            self.generate_profile_database('combined')
        if self.summed_out_dir:
            filesnpaths.gen_output_directory(self.summed_out_dir, delete_if_exists=self.overwrite_out_dest)
            self.generate_auxiliary_database('summed')
            self.generate_profile_database('summed')


    def sanity_check(self):
        """Check `anvi-merge-trnaseq` arguments."""
        for trnaseq_db_path in self.trnaseq_db_paths:
            is_trnaseq_db(trnaseq_db_path)
        self.populate_trnaseq_dbs_info_dict()
        self.trnaseq_db_sample_ids = [trnaseq_db_info_dict['sample_id'] for trnaseq_db_info_dict in self.trnaseq_dbs_info_dict.values()]
        if len(self.trnaseq_dbs_info_dict) != len(set(self.trnaseq_db_sample_ids)):
            raise ConfigError("Sample IDs in each input tRNA-seq db must be unique. "
                              "This is not the case with your input. "
                              "Here are the sample names so you can see which ones occur more than once: "
                              f"'{', '.join(self.trnaseq_db_sample_ids)}'")
        self.num_trnaseq_dbs = len(self.trnaseq_db_sample_ids)

        self.check_trnaseq_db_versions()

        self.out_dir = filesnpaths.check_output_directory(self.out_dir, ok_if_exists=self.overwrite_out_dest)
        self.out_dir = os.path.abspath(self.out_dir)

        self.contigs_db_path = os.path.join(self.out_dir, 'CONTIGS.db')
        self.contigs_db_hash = 'hash' + str('%08x' % random.randrange(16**8))

        self.spec_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'SPECIFIC_COVERAGE'), ok_if_exists=self.overwrite_out_dest)
        self.spec_profile_db_path = os.path.join(self.spec_out_dir, 'PROFILE.db')
        self.spec_auxiliary_db_path = os.path.join(self.spec_out_dir, 'AUXILIARY-DATA.db')

        if not 1 <= self.num_threads <= multiprocessing.cpu_count():
            raise ConfigError(f"The number of threads to use must be a positive integer less than or equal to {multiprocessing.cpu_count()}. Try again!")

        self.set_treatment_preference()

        self.set_nonspecific_database_info()

        check_sample_id(self.project_name)

        if self.descrip_path:
            filesnpaths.is_file_plain_text(self.descrip_path)
            self.descrip_path = os.path.abspath(self.descrip_path)
            self.descrip = open(self.descrip_path).read()

        if self.seed_limit == -1:
            self.seed_limit = MAXSIZE
        elif self.seed_limit < 1:
            raise ConfigError(f"{self.seed_limit} is an invalid value for `--max-reported-seed-seqs`. "
                              "To remove the limit on tRNA seeds reported to the contigs db, provide a value of -1. "
                              "Otherwise provide an integer greater than 0.")

        self.run.info("Input tRNA-seq dbs", ", ".join(self.trnaseq_db_paths))
        if self.preferred_treatment:
            self.run.info("Databases preferred for seed formation",
                          ", ".join([trnaseq_db_path for trnaseq_db_num, trnaseq_db_path in enumerate(self.trnaseq_db_paths)
                                     if trnaseq_db_num in self.preferred_trnaseq_db_nums]))
        self.run.info("Output directory", self.out_dir)


    def populate_trnaseq_dbs_info_dict(self):
        """Get the meta-data from the input tRNA-seq databases."""
        for trnaseq_db_path in self.trnaseq_db_paths:
            trnaseq_db = dbops.TRNASeqDatabase(trnaseq_db_path)
            self.trnaseq_dbs_info_dict[trnaseq_db_path] = trnaseq_db.meta


    def check_trnaseq_db_versions(self):
        if len(set([trnaseq_db_info_dict['version'] for trnaseq_db_info_dict in self.trnaseq_dbs_info_dict.values()])) > 1:
            trnaseq_db_version_report = "\n".join(
                [trnaseq_db_path + " : " + trnaseq_db_info_dict['version']
                 for trnaseq_db_path, trnaseq_db_info_dict in self.trnaseq_dbs_info_dict.items()])
            if anvio.FORCE:
                self.run.warning("Not all input tRNA-seq dbs have the same version number, but since you have used the `--force` flag, "
                                 "`anvi-merge-trnaseq` will proceed though this is dangerous and may lead to errors. "
                                 f"Here is the version number of each database:\n{trnaseq_db_version_report}")
            else:
                raise ConfigError("Not all input tRNA-seq dbs have the same version number. "
                                  f"Here is the version number of each db:\n{trnaseq_db_version_report}")


    def set_treatment_preference(self):
        if not self.preferred_treatment:
            return

        input_treatments = [trnaseq_db_info_dict['treatment'] for trnaseq_db_info_dict in self.trnaseq_dbs_info_dict.values()]
        self.preferred_trnaseq_db_sample_ids = []
        self.preferred_trnaseq_db_nums = []
        if self.preferred_treatment not in input_treatments:
            raise ConfigError(f"You provided a preferred treatment type, {self.preferred_treatment}, "
                              "but it was not found in any of the input dbs, "
                              f"which were found to have the following treatments: {', '.join(input_treatments)}.")
        for trnaseq_db_num, treatment in enumerate(input_treatments):
            if self.preferred_treatment == treatment:
                self.preferred_trnaseq_db_sample_ids.append(self.trnaseq_db_sample_ids[trnaseq_db_num])
                self.preferred_trnaseq_db_nums.append(trnaseq_db_num)


    def set_nonspecific_database_info(self):
        self.nonspec_db_types = self.nonspec_output.split(',')
        for nonspec_db_type in self.nonspec_db_types:
            if nonspec_db_type not in ['nonspecific_db', 'combined_db', 'summed_db']:
                raise ConfigError("The nonspecific profile db types provided by `--nonspecific-output` are not recognized. "
                                  "The db types must be comma separated without spaces, e.g., 'nonspecific_db,combined_db,summed_db'. "
                                  f"Your argument was: {self.nonspec_output}'")

        if 'nonspecific_db' in self.nonspec_db_types:
            self.nonspec_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'NONSPECIFIC_COVERAGE'), ok_if_exists=self.overwrite_out_dest)
            self.nonspec_profile_db_path = os.path.join(self.nonspec_out_dir, 'PROFILE.db')
            self.nonspec_auxiliary_db_path = os.path.join(self.nonspec_out_dir, 'AUXILIARY-DATA.db')

        if 'combined_db' in self.nonspec_db_types:
            self.combined_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'COMBINED_COVERAGE'), ok_if_exists=self.overwrite_out_dest)
            self.combined_profile_db_path = os.path.join(self.combined_out_dir, 'PROFILE.db')
            self.combined_auxiliary_db_path = os.path.join(self.combined_out_dir, 'AUXILIARY-DATA.db')

        if 'summed_db' in self.nonspec_db_types:
            self.summed_out_dir = filesnpaths.check_output_directory(os.path.join(self.out_dir, 'SUMMED_COVERAGE'), ok_if_exists=self.overwrite_out_dest)
            self.summed_profile_db_path = os.path.join(self.summed_out_dir, 'PROFILE.db')
            self.summed_auxiliary_db_path = os.path.join(self.summed_out_dir, 'AUXILIARY-DATA.db')


    def load_trnaseq_dbs(self):
        """Load information from input tRNA-seq databases."""
        loaded_db_count = 0
        trnaseq_db_paths = self.trnaseq_db_paths
        num_trnaseq_db_paths = len(trnaseq_db_paths)
        pid = "Loading seq info from tRNA-seq dbs"
        self.progress.new(pid)
        self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} dbs loaded")

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue_Nu_summaries = manager.Queue()
        output_queue_M_summaries = manager.Queue()
        processes = [multiprocessing.Process(target=trnaseq_db_loader,
                                args=(input_queue, output_queue_Nu_summaries, output_queue_M_summaries, self))
                     for _ in range(self.num_threads)]
        for p in processes:
            p.start()

        for trnaseq_db_path in trnaseq_db_paths:
            input_queue.put(trnaseq_db_path)

        poison_pill_count = 0
        dict_summaries_Nu = self.dict_summaries_Nu
        dict_summaries_M = self.dict_summaries_M
        empty = queue.Empty
        db_completion_dict = {}
        while poison_pill_count < len(trnaseq_db_paths) * 2:
            # For each input db, one poison pill is put at the end of the Nu queue and another at
            # the end of the M queue.
            try:
                summary_item = output_queue_Nu_summaries.get_nowait()
                try:
                    trnaseq_db_path, summary_Nu = summary_item
                    dict_summaries_Nu[trnaseq_db_path].append(summary_Nu)
                except ValueError:
                    # The poison pill indicates there are no more Nu to be retrieved.
                    trnaseq_db_path = summary_item
                    poison_pill_count += 1
                    try:
                        # The poison pill was already returned in the db's M queue, so all Nu and M
                        # have been retrieved, and the db has been fully processed.
                        db_completion_dict[trnaseq_db_path] += 1
                        loaded_db_count += 1
                        self.progress.update_pid(pid)
                        self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} dbs loaded")
                    except KeyError:
                        # The poison pill has not been returned in the db's M queue.
                        db_completion_dict[trnaseq_db_path] = 1
            except empty:
                pass # waiting on items in the Nu queue

            try:
                summary_item = output_queue_M_summaries.get_nowait()
                try:
                    trnaseq_db_path, summary_M = summary_item
                    dict_summaries_M[trnaseq_db_path].append(summary_M)
                except ValueError:
                    trnaseq_db_path = summary_item
                    poison_pill_count += 1
                    try:
                        db_completion_dict[trnaseq_db_path] += 1
                        loaded_db_count += 1
                        self.progress.update_pid(pid)
                        self.progress.update(f"{loaded_db_count}/{num_trnaseq_db_paths} dbs loaded")
                    except KeyError:
                        db_completion_dict[trnaseq_db_path] = 1
            except empty:
                pass # waiting on items in the M queue

        for p in processes:
            p.terminate()
            p.join()

        self.progress.end()


    def load_trnaseq_database_sequence_summaries(self, trnaseq_db_path):
        """Load necessary tRNA sequence data from the input tRNA-seq database.

        Unmodified normalized sequences and "modified" sequences, comprising clustered normalized
        sequences, are stored in distinct data structures.
        """
        trnaseq_db_num = list(self.trnaseq_dbs_info_dict.keys()).index(trnaseq_db_path)
        sample_id = self.trnaseq_db_sample_ids[trnaseq_db_num]

        trnaseq_db = dbops.TRNASeqDatabase(trnaseq_db_path)

        dict_N_summary = {} # Used to link N to M summaries

        # Load columns from the "normalized", "features", and "trimmed" tables.
        df_N = pd.merge(
            pd.DataFrame(
                trnaseq_db.db.get_some_columns_from_table('normalized', ', '.join(self.NORM_TABLE_COLS_OF_INTEREST)),
                columns=self.NORM_TABLE_COLS_OF_INTEREST
            ).set_index('name'),
            pd.merge(pd.DataFrame(
                        trnaseq_db.db.get_some_columns_from_table('feature', ', '.join(self.FEATURE_TABLE_COLS_OF_INTEREST)),
                        columns=self.FEATURE_TABLE_COLS_OF_INTEREST
                     ).set_index('name'),
                     pd.DataFrame(
                        trnaseq_db.db.get_some_columns_from_table('trimmed', ', '.join(self.TRIMMED_TABLE_COLS_OF_INTEREST)),
                        columns=self.TRIMMED_TABLE_COLS_OF_INTEREST
                     ).set_index('name'),
                     left_index=True,
                     right_index=True),
            left_index=True,
            right_index=True)

        FEATURE_INDEX_COLS_OF_INTEREST = self.FEATURE_INDEX_COLS_OF_INTEREST
        threshold_feature = self.feature_threshold + '_start'
        # The starts of both strands of the stem are recorded, so pick the start of the 5' strand.
        is_threshold_feature_stem = True if 'stem' in threshold_feature else False
        for info_N in df_N.itertuples():
            if info_N.id_info == 'indel_aligned':
                # Ignore Ni. The coverage of indels themselves is recorded in the parent M, but the
                # contribution of Ni to nt coverage is ignored. Inclusion of Ni would produce
                # numerous complications (e.g., they don't have feature profiles).
                continue

            summary_N = NormalizedSequenceSummary()
            summary_N.name = info_N.Index
            summary_N.sample_id = sample_id
            summary_N.mean_spec_cov = info_N.mean_specific_coverage
            # There is a trailing comma in the coverage strings.
            summary_N.spec_covs = np.fromiter(map(int, info_N.specific_coverages.split(',')[: -1]), int)
            summary_N.nonspec_covs = np.fromiter(map(int, info_N.nonspecific_coverages.split(',')[: -1]), int)
            summary_N.string = info_N.sequence
            summary_N.anticodon_string = info_N.anticodon_sequence if info_N.anticodon_sequence else ''
            length_N = len(summary_N.string)
            summary_N.spec_nt_covs_dict = spec_nt_covs_dict = {nt: np.zeros(length_N, dtype=int) for nt in UNAMBIG_NTS}
            summary_N.nonspec_nt_covs_dict = nonspec_nt_covs_dict = {nt: np.zeros(length_N, dtype=int) for nt in UNAMBIG_NTS}
            for nt_pos, nt, spec_cov, nonspec_cov in zip(range(length_N), summary_N.string, summary_N.spec_covs, summary_N.nonspec_covs):
                spec_nt_covs_dict[nt][nt_pos] += spec_cov
                nonspec_nt_covs_dict[nt][nt_pos] += nonspec_cov

            length_5prime = info_N.num_extra_fiveprime
            summary_N.feature_dict = feature_dict = OrderedDict()
            for feature_index_col in FEATURE_INDEX_COLS_OF_INTEREST: # `trna_his_position_0_start`, etc.
                db_value = getattr(info_N, feature_index_col)
                if isinstance(db_value, str):
                    # The string contains the start indices of stem strands.
                    if length_5prime:
                        split_db_value = tuple(map(int, db_value.split(',')))
                        feature_dict[feature_index_col] = (split_db_value[0] - length_5prime, split_db_value[1] - length_5prime)
                    else:
                        feature_dict[feature_index_col] = tuple(map(int, db_value.split(',')))
                else:
                    try:
                        # N/A in numerical columns where the feature is not a stem load as `np.nan`.
                        # Represent missing feature indices as an impossibly negative number.
                        if np.isnan(db_value):
                            feature_dict[feature_index_col] = -1000
                        else:
                            feature_dict[feature_index_col] = int(db_value) - length_5prime
                    except TypeError:
                        # N/A in string columns where the feature is a stem load as `None`.
                        feature_dict[feature_index_col] = -1000

            db_value = getattr(info_N, threshold_feature)
            if is_threshold_feature_stem:
                summary_N.feature_threshold_start = int(db_value.split(',')[0]) - length_5prime if isinstance(db_value, str) else -1000
            else:
                summary_N.feature_threshold_start = -1000 if np.isnan(db_value) else db_value - length_5prime

            dict_N_summary[summary_N.name] = summary_N

        summaries_M = []
        for info_M in pd.DataFrame(
            trnaseq_db.db.get_some_columns_from_table('modified', ', '.join(self.MOD_TABLE_COLS_OF_INTEREST)),
            columns=self.MOD_TABLE_COLS_OF_INTEREST).set_index('name').itertuples():
            summary_M = ModifiedSequenceSummary()
            summary_M.name = info_M.Index
            summary_M.sample_id = sample_id
            # There is a trailing comma in the substitution and indel coverage and position strings.
            summary_M.sub_positions = np.fromiter(map(int, info_M.substitution_positions.split(',')[: -1]), int)
            summary_M.consensus_string = consensus_string = info_M.consensus_sequence

            # Make nt variability arrays covering every position in the seq. Start with arrays of
            # overall specific/nonspecific coverage with nonzero values for the nts found in the
            # consensus seq, and then correct the variable positions.
            seq_length = len(consensus_string)
            summary_M.spec_nt_covs_dict = spec_nt_covs_dict = {nt: np.zeros(seq_length, int) for nt in UNAMBIG_NTS}
            summary_M.nonspec_nt_covs_dict = nonspec_nt_covs_dict = {nt: np.zeros(seq_length, int) for nt in UNAMBIG_NTS}
            pos = 0
            for nt, spec_cov, nonspec_cov in zip(consensus_string,
                                                 info_M.specific_coverages.split(',')[: -1],
                                                 info_M.nonspecific_coverages.split(',')[: -1]):
                spec_nt_covs_dict[nt][pos] = spec_cov
                nonspec_nt_covs_dict[nt][pos] = nonspec_cov
                pos += 1
            for nt in UNAMBIG_NTS:
                sub_positions = summary_M.sub_positions
                for sub_pos, spec_nt_cov, nonspec_nt_cov in zip(
                    sub_positions,
                    map(int, getattr(info_M, 'substitution_' + nt + '_specific_coverage').split(',')[: -1]),
                    map(int, getattr(info_M, 'substitution_' + nt + '_nonspecific_coverage').split(',')[: -1])):
                    spec_nt_covs_dict[nt][sub_pos] = spec_nt_cov
                    nonspec_nt_covs_dict[nt][sub_pos] = nonspec_nt_cov

            if info_M.insertion_starts == ',':
                summary_M.insert_starts = np.zeros(0, dtype=int)
                summary_M.insert_strings = []
                summary_M.spec_insert_covs = np.zeros(0, dtype=int)
                summary_M.nonspec_insert_covs = np.zeros(0, dtype=int)
            else:
                summary_M.insert_starts = np.fromiter(map(int, info_M.insertion_starts.split(',')[: -1]), int)
                summary_M.insert_strings = list(info_M.insertion_seqs.split(',')[: -1])
                summary_M.spec_insert_covs = np.fromiter(map(int, info_M.insertion_specific_coverages.split(',')[: -1]), int)
                summary_M.nonspec_insert_covs = np.fromiter(map(int, info_M.insertion_nonspecific_coverages.split(',')[: -1]), int)

            if info_M.deletion_starts == ',':
                summary_M.del_starts = np.zeros(0, dtype=int)
                summary_M.del_lengths = []
                summary_M.spec_del_covs = np.zeros(0, dtype=int)
                summary_M.nonspec_del_covs = np.zeros(0, dtype=int)
            else:
                summary_M.del_starts = np.fromiter(map(int, info_M.deletion_starts.split(',')[: -1]), int)
                summary_M.del_lengths = list(map(int, info_M.deletion_lengths.split(',')[: -1]))
                summary_M.spec_del_covs = np.fromiter(map(int, info_M.deletion_specific_coverages.split(',')[: -1]), int)
                summary_M.nonspec_del_covs = np.fromiter(map(int, info_M.deletion_nonspecific_coverages.split(',')[: -1]), int)

            summary_M.summaries_Nb = summaries_Nb = []
            for name_Nb in info_M.names_of_normalized_seqs_without_indels.split(','):
                summary_Nb = dict_N_summary[name_Nb]

                # Cross-reference the M summary with summary objects of constituent Nb.
                summary_Nb.summary_M = summary_M
                summaries_Nb.append(summary_Nb)

                # Ensure that all constituent Nb have coverage arrays flush with those of M.
                if summary_Nb.spec_covs.size < seq_length:
                    elongation_5prime = np.zeros(seq_length - summary_Nb.spec_covs.size, int)
                    self.elongate_normalized_sequence_fiveprime(summary_Nb, elongation_5prime)

            summaries_M.append(summary_M)

        summaries_Nu = [summary_N for summary_N in dict_N_summary.values() if not summary_N.summary_M]

        return summaries_Nu, summaries_M


    def elongate_normalized_sequence_fiveprime(self, summary_N, elongation_5prime):
        """Seed sequences can be longer than the normalized sequences from the individual samples,
        requiring addition of empty positions in the normalized sequence coverage arrays at the 5'
        end and reindexing of features, as normalized (and modified) sequences are aligned from the
        3' end."""
        summary_N.spec_covs = np.concatenate([elongation_5prime, summary_N.spec_covs])
        summary_N.nonspec_covs = np.concatenate([elongation_5prime, summary_N.nonspec_covs])
        spec_nt_covs_dict = summary_N.spec_nt_covs_dict
        for nt, spec_covs in spec_nt_covs_dict.items():
            spec_nt_covs_dict[nt] = np.concatenate([elongation_5prime, spec_covs])
        nonspec_nt_covs_dict = summary_N.nonspec_nt_covs_dict
        for nt, nonspec_covs in nonspec_nt_covs_dict.items():
            nonspec_nt_covs_dict[nt] = np.concatenate([elongation_5prime, nonspec_covs])

        elongation_length = elongation_5prime.size
        new_feature_dict = OrderedDict()
        for feature, feature_index in summary_N.feature_dict.items():
            try:
                new_feature_dict[feature] = (feature_index[0] + elongation_length, feature_index[1] + elongation_length)
            except TypeError:
                # Missing features with an index of np.nan result in np.nan after addition.
                try:
                    new_feature_dict[feature] = feature_index + elongation_length
                except TypeError:
                    # `feature_index` is None rather than np.nan for a missing stem start index. For
                    # convenience, replace it with np.nan.
                    new_feature_dict[feature] = np.nan
        summary_N.feature_dict = new_feature_dict
        summary_N.feature_threshold_start + elongation_length


    def form_seeds(self):
        """Form tRNA seeds through comparison of sequences from the input samples.

        Normalized sequences (N) are compared. These include N underlying modified sequences (M),
        called Nm, and N that are not part of M, called Nu.

        Modification-induced mutations (substitutions and indels) complicate seed formation. If N is
        shared identically (not as a subsequence) between samples, then N from the samples and any M
        that N are part of are combined into a single seed. The seed sequence is that of the longest
        N and need not be found in every sample.

        It is a heuristic to exactly match N, rather than to check whether one N is a 3' subsequence
        of the other, or to instead compare underlying trimmed sequences making up N. This heuristic
        should not distort sample merging for the more abundant tRNA species, in particular, as
        these are most likely to be represented by reads spanning the full length of the tRNA,
        producing the same N."""
        pid = "Forming seed seqs from input samples"
        self.progress.new(pid)

        string_N_seed_dict = {}
        for trnaseq_db_num, trnaseq_db_path in enumerate(self.trnaseq_db_paths):
            sample_id = self.trnaseq_db_sample_ids[trnaseq_db_num]
            self.progress.update_pid(pid)
            self.progress.update(f"Adding {sample_id}")

            summaries_Nu = self.dict_summaries_Nu[trnaseq_db_path]
            summaries_M = self.dict_summaries_M[trnaseq_db_path]

            for summary_Nu in summaries_Nu:
                # Process Nu.
                string_Nu = summary_Nu.string
                try:
                    # Nu has already been found in another dataset.
                    seed = string_N_seed_dict[string_Nu]
                except KeyError:
                    # Create a new seed based on Nu.
                    seed = SeedSequence()
                    seed.name = summary_Nu.name + '_' + sample_id
                    seed.string = string_Nu
                    seed.meets_feature_threshold = True if summary_Nu.feature_threshold_start >= 0 else False
                    seed.summaries_Nu = [summary_Nu]
                    seed.summaries_M = []
                    string_N_seed_dict[string_Nu] = seed
                    continue

                if len(string_Nu) < len(seed.string):
                    # Nu is shorter than the existing seed. This implies that the seed was formed
                    # from a longer M in another dataset. Extend Nu coverage arrays the needed
                    # amount at the 5' end. (Note: It is impossible here for Nu to be longer than
                    # the seed.)
                    elongation_5prime = np.zeros(len(seed.string) - len(string_Nu), int)
                    self.elongate_normalized_sequence_fiveprime(summary_Nu, elongation_5prime)
                seed.summaries_Nu.append(summary_Nu)

            for summary_M in summaries_M:
                # Find seeds from other datasets containing any of Nb forming the M under
                # consideration. If >1 seed is identified, they are merged.
                seed_dict = {}
                for summary_Nb in summary_M.summaries_Nb:
                    try:
                        # Nb is represented in another dataset.
                        seed = string_N_seed_dict[summary_Nb.string]
                    except KeyError:
                        continue
                    seed_dict[seed.name] = seed

                if not seed_dict:
                    # Create a new seed based on M.
                    seed = SeedSequence()
                    seed.name = summary_M.name + '_' + sample_id
                    seed.string = summary_M.consensus_string
                    seed.summaries_Nu = []
                    seed.summaries_M = []
                    seed.summaries_M.append(summary_M)
                    for summary_Nb in summary_M.summaries_Nb:
                        if summary_Nb.feature_threshold_start >= 0:
                            seed.meets_feature_threshold = True
                            break
                    else:
                        seed.meets_feature_threshold = False
                    for summary_Nb in summary_M.summaries_Nb:
                        string_N_seed_dict[summary_Nb.string] = seed
                    continue

                if len(seed_dict) == 1:
                    # M shares ≥1 Nb with 1 seed.
                    seed_name, seed = seed_dict.popitem()
                    if len(summary_M.consensus_string) < len(seed.string):
                        # M is shorter than the seed, so its coverage arrays must be extended with
                        # zeros at the 5' end.
                        elongation_5prime = np.zeros(len(seed.string) - len(summary_M.consensus_string), int)
                        self.elongate_modified_sequence_fiveprime(summary_M, elongation_5prime)
                    elif len(summary_M.consensus_string) > len(seed.string):
                        # M is longer than the seed, so the coverage arrays of N forming the seed
                        # must be extended with zeros at the 5' end.
                        elongation_5prime = np.zeros(len(summary_M.consensus_string) - len(seed.string), int)
                        for summary_Nu in seed.summaries_Nu:
                            self.elongate_normalized_sequence_fiveprime(summary_Nu, elongation_5prime)
                        for other_summary_M in seed.summaries_M:
                            self.elongate_modified_sequence_fiveprime(other_summary_M, elongation_5prime)
                        seed.name = summary_M.name + '_' + sample_id
                        seed.string = summary_M.consensus_string
                        for summary_Nb in summary_M.summaries_Nb:
                            if summary_Nb.feature_threshold_start >= 0:
                                seed.meets_feature_threshold = True
                                break
                        else:
                            seed.meets_feature_threshold = False
                    seed.summaries_M.append(summary_M)
                    for summary_Nb in summary_M.summaries_Nb:
                        string_N_seed_dict[summary_Nb.string] = seed
                    continue

                # To reach this point, M must map to >1 seed.
                sorted_seeds = sorted([seed for seed in seed_dict.values()], key=lambda seed: -len(seed.string))
                max_seed_length = len(sorted_seeds[0].string)
                length_M = len(summary_M.consensus_string)
                if length_M < max_seed_length:
                    # Extend coverage arrays of M.
                    elongation_5prime = np.zeros(max_seed_length - length_M, int)
                    self.elongate_modified_sequence_fiveprime(summary_M, elongation_5prime)

                    # Extend coverage arrays of shorter seeds now grouped with a longer seed.
                    for seed in seed_dict.values():
                        if len(seed.string) < max_seed_length:
                            elongation_5prime = np.zeros(max_seed_length - len(seed.string), int)
                            for summary_Nu in seed.summaries_Nu:
                                self.elongate_normalized_sequence_fiveprime(summary_Nu, elongation_5prime)
                            for other_summary_M in seed.summaries_M:
                                self.elongate_modified_sequence_fiveprime(other_summary_M, elongation_5prime)

                    new_seed = SeedSequence()
                    longest_seed = sorted_seeds[0]
                    new_seed.name = longest_seed.name
                    new_seed.string = longest_seed.string
                    new_seed.meets_feature_threshold = longest_seed.meets_feature_threshold
                elif length_M > max_seed_length:
                    # Extend coverage arrays of seeds.
                    for seed in seed_dict.values():
                        elongation_5prime = np.zeros(length_M - len(seed.string), int)
                        for summary_Nu in seed.summaries_Nu:
                            self.elongate_normalized_sequence_fiveprime(summary_Nu, elongation_5prime)
                        for other_summary_M in seed.summaries_M:
                            self.elongate_modified_sequence_fiveprime(other_summary_M, elongation_5prime)

                    new_seed = SeedSequence()
                    new_seed.name = summary_M.name + '_' + sample_id
                    new_seed.string = summary_M.consensus_string
                    for seed in seed_dict.values():
                        if seed.meets_feature_threshold:
                            new_seed.meets_feature_threshold = True
                            break
                    else:
                        for summary_Nb in summary_M.summaries_Nb:
                            if summary_Nb.feature_threshold_start >= 0:
                                new_seed.meets_feature_threshold = True
                                break
                        else:
                            new_seed.meets_feature_threshold = False
                else:
                    # M is the same length as the longest seed. Extend coverage arrays of shorter
                    # seeds now grouped with a longer seed.
                    for seed in seed_dict.values():
                        if len(seed.string) < max_seed_length:
                            elongation_5prime = np.zeros(max_seed_length - len(seed.string), int)
                            for summary_Nu in seed.summaries_Nu:
                                self.elongate_normalized_sequence_fiveprime(summary_Nu, elongation_5prime)
                            for other_summary_M in seed.summaries_M:
                                self.elongate_modified_sequence_fiveprime(other_summary_M, elongation_5prime)

                    new_seed = SeedSequence()
                    new_seed.name = summary_M.name + '_' + sample_id
                    new_seed.string = summary_M.consensus_string
                    if sorted_seeds[0].meets_feature_threshold:
                        new_seed.meets_feature_threshold = True
                    else:
                        for summary_Nb in summary_M.summaries_Nb:
                            if summary_Nb.feature_threshold_start >= 0:
                                new_seed.meets_feature_threshold = True
                                break
                        else:
                            new_seed.meets_feature_threshold = False

                # Now that all of the coverage arrays are reconciled in length, M and constituent N
                # of the matching seeds can be added to the new seed.
                new_seed.summaries_Nu = []
                new_seed.summaries_M = []
                for seed in seed_dict.values():
                    new_seed.summaries_Nu += seed.summaries_Nu
                    new_seed.summaries_M += seed.summaries_M
                new_seed.summaries_M.append(summary_M)

                for summary_Nu in new_seed.summaries_Nu:
                    string_N_seed_dict[summary_Nu.string] = new_seed
                for summary_M in new_seed.summaries_M:
                    for summary_Nb in summary_M.summaries_Nb:
                        string_N_seed_dict[summary_Nb.string] = new_seed

        self.progress.update_pid(pid)
        self.progress.update("...")
        # The seed references in the dict need to be dereplicated.
        seeds = list({seed.name: seed for seed in string_N_seed_dict.values()}.values())

        # Disregard seeds that do not reach the 5' feature threshold.
        seeds = [seed for seed in seeds if seed.meets_feature_threshold]

        self.progress.update_pid(pid)
        self.progress.update("Checking feature profiles")
        # Reads with extra 5' nts beyond the acceptor stem sometimes generate false positive tRNA
        # feature profiles. These erroneous profiles shoehorn the 5' extra nts into the profile
        # through nt accommodation in variable-length sections of the profile, such as the D loop.
        # Make sure there are no N shorter than the seed sequence that have a full-length profile.

        # The algorithm requires seq summaries sorted in descending order of seq length.
        for seed in seeds:
            seed.summaries_Nu = sorted([summary_Nu for summary_Nu in seed.summaries_Nu],
                                       key=lambda summary_Nu: -len(summary_Nu.string))
            seed.summaries_M = sorted([summary_M for summary_M in seed.summaries_M],
                                      key=lambda summary_M: -len(summary_M.consensus_string))
            for summary_M in seed.summaries_M:
                summary_M.summaries_Nb = sorted([summary_Nb for summary_Nb in summary_M.summaries_Nb],
                                                key=lambda summary_Nb: -len(summary_Nb.string))

        names_seeds_with_conflict = set()
        seed_indices_below_feature_threshold = []
        for seed_index, seed in enumerate(seeds):
            while True:
                # Remove the longest N (Nu or Nb) with inconsistent full-length profiles until only
                # N with consistent full-length profiles remain. Exit the while loop at that point.
                seed_length = len(seed.string)
                longest_Nu_indices = []
                longest_Nb_indices = []
                # Get the feature profile of the first N with the same seq as the seed plus a
                # full-length profile.
                seed_features = None
                # The following variable is used for a specific circumstance. Nu are searched before
                # M. An Nu the same length as the seed with a full-length profile may not be found
                # while a shorter Nu with a full-length profile is found, apparently a conflict.
                # However, an M the same length as the seed with a full-length profile and a
                # tRNA-His 5'-G can subsequently explain the full-length profile of the shorter Nu.
                check_His_G = False
                # Search for a conflict between the full-length profile of a shorter seq and the
                # seed profile.
                found_conflict = False

                for index_Nu, summary_Nu in enumerate(seed.summaries_Nu):
                    if len(summary_Nu.string) == seed_length:
                        longest_Nu_indices.append(index_Nu)
                        if seed_features:
                            # Only record the first qualifying profile.
                            pass
                        elif summary_Nu.feature_dict['fiveprime_acceptor_stem_sequence_start'] >= -10:
                            # Extrapolated acceptor stems, e.g., with starts of -1, can be
                            # discredited.
                            seed_features = tuple(summary_Nu.feature_dict.values())
                        else:
                            # The seed does not have a full-length profile, so stop looking.
                            break
                    elif summary_Nu.feature_dict['fiveprime_acceptor_stem_sequence_start'] >= 0:
                        # A shorter N than the seed was assigned a full-length profile. This
                        # indicates a conflict (with one exception involving tRNA-His 5'-G), and the
                        # longest N will be removed from the seed.
                        if seed_features:
                            if seed_features[0] == 0:
                                if (summary_Nu.feature_dict['fiveprime_acceptor_stem_sequence_start'] == 1
                                    and ANTICODON_AA_DICT[summary_Nu.anticodon_string] == 'His'):
                                    # The longest N was also assigned a full-length profile and only
                                    # differs from the shorter profile by a 5'-G at the end of
                                    # tRNA-His, so a conflict was not found.
                                    break
                        if (not seed_features
                            and summary_Nu.feature_dict['fiveprime_acceptor_stem_sequence_start'] == 1
                            and ANTICODON_AA_DICT[summary_Nu.anticodon_string] == 'His'):
                            check_His_G = True
                        found_conflict = True
                        break

                for index_M, summary_M in enumerate(seed.summaries_M):
                    if len(summary_M.consensus_string) == seed_length:
                        for index_Nb, summary_Nb in enumerate(summary_M.summaries_Nb):
                            # Each M consists of Nb which may be of varying lengths.
                            if len(summary_Nb.string) == seed_length:
                                longest_Nb_indices.append((index_M, index_Nb))
                                if seed_features:
                                    pass
                                if summary_Nb.feature_dict['fiveprime_acceptor_stem_sequence_start'] >= -10:
                                    seed_features = tuple(summary_Nb.feature_dict.values())
                                    if check_His_G:
                                        if seed_features[0] == 0:
                                            found_conflict = False
                                else:
                                    # The seed does not have a full-length profile, so stop
                                    # searching M.
                                    break
                        else:
                            continue
                        # To reach this point, it was found that the seed does not have a
                        # full-length profile, so stop searching.
                        break
                    else:
                        for summary_Nb in summary_M.summaries_Nb:
                            if summary_Nb.feature_dict['fiveprime_acceptor_stem_sequence_start'] >= 0:
                                if seed_features:
                                    if seed_features[0] == 0:
                                        if (summary_Nb.feature_dict['fiveprime_acceptor_stem_sequence_start'] == 1
                                            and ANTICODON_AA_DICT[summary_Nb.anticodon_string] == 'His'):
                                            break
                                found_conflict = True
                                break

                if not found_conflict:
                    try:
                        max_length_Nu = len(seed.summaries_Nu[0].string)
                    except IndexError:
                        max_length_Nu = 0
                    try:
                        max_length_M = len(seed.summaries_M[0].consensus_string)
                    except IndexError:
                        max_length_M = 0
                    if not max_length_Nu:
                        selected_summary = seed.summaries_M[0].summaries_Nb[0]
                    elif not max_length_M:
                        selected_summary = seed.summaries_Nu[0]
                    else:
                        if max_length_Nu >= max_length_M:
                            selected_summary = seed.summaries_Nu[0]
                        else:
                            selected_summary = seed.summaries_M[0].summaries_Nb[0]
                    seed.feature_dict = selected_summary.feature_dict

                    if selected_summary.feature_threshold_start < 0:
                        seed.meets_feature_threshold = False
                        seed_indices_below_feature_threshold.append(seed_index)
                    break
                names_seeds_with_conflict.add(seed.name)

                # Remove the longest Nu from the seed.
                for longest_Nu_index in longest_Nu_indices[::-1]:
                    seed.summaries_Nu.pop(longest_Nu_index)

                altered_M_indices = set()
                for longest_Nb_index in longest_Nb_indices[::-1]:
                    seed.summaries_M[longest_Nb_index[0]].summaries_Nb.pop(longest_Nb_index[1])
                altered_M_indices = set([longest_Nb_index[0] for longest_Nb_index in longest_Nb_indices])

                # If all Nb are removed from an M, then remove M from the seed.
                altered_M_summaries = []
                if longest_Nb_indices:
                    empty_M_indices = []
                    for index_M, summary_M in enumerate(seed.summaries_M):
                        if not summary_M.summaries_Nb:
                            empty_M_indices.append(index_M)
                        elif index_M in altered_M_indices:
                            altered_M_summaries.append(summary_M)
                    for empty_M_index in empty_M_indices[::-1]:
                        seed.summaries_M.pop(empty_M_index)

                # Find the new length of the seed.
                max_length_Nu = 0
                if seed.summaries_Nu:
                    max_length_Nu = len(seed.summaries_Nu[0].string)
                max_length_M = 0
                for summary_M in seed.summaries_M:
                    length_M = len(summary_M.summaries_Nb[0].string)
                    if length_M > max_length_M:
                        max_length_M = length_M
                new_seed_length = max(max_length_Nu, max_length_M)
                seed_length_diff = seed_length - new_seed_length

                # Adjust the following attributes of Nu, M, and Nb summary objects that depend upon
                # seed length: coverage array lengths, feature indices, threshold feature start,
                # indel information (some indels may now be removed from M).
                for summary_Nu in seed.summaries_Nu:
                    self.shorten_normalized_sequence_fiveprime(summary_Nu, seed_length_diff)

                for summary_M in seed.summaries_M:
                    self.shorten_modified_sequence_fiveprime(summary_M, seed_length_diff)

                if max_length_Nu >= max_length_M:
                    seed.string = seed.summaries_Nu[0].string
                    selected_summary = seed.summaries_Nu[0]
                else:
                    seed.string = seed.summaries_M[0].consensus_string
                    selected_summary = seed.summaries_M[0].summaries_Nb[0]
                seed.feature_dict = selected_summary.feature_dict
        for seed_index in seed_indices_below_feature_threshold[::-1]:
            seeds.pop(seed_index)


        self.progress.update_pid(pid)
        self.progress.update("Assigning anticodons")
        # Assign the anticodon by comparing the mean specific coverage of N comprising the seed, a
        # simpler approximation of extracting the anticodon coverage from each N.
        for seed in seeds:
            if not seed.summaries_M:
                # The seed is comprised entirely of Nu without unlike nts.
                seed.anticodon_string = seed.summaries_Nu[0].anticodon_string
                continue

            anticodon_cov_dict = defaultdict(int)
            for summary_Nu in seed.summaries_Nu:
                # Mean specific coverage does not include any 5' padding of zero coverage from the
                # formation of the seed.
                anticodon_cov_dict[summary_Nu.anticodon_string] += summary_Nu.mean_spec_cov
            for summary_M in seed.summaries_M:
                for summary_Nb in summary_M.summaries_Nb:
                    anticodon_cov_dict[summary_Nb.anticodon_string] += summary_Nb.mean_spec_cov

            seed.anticodon_string = sorted(anticodon_cov_dict.items(), key=lambda anticodon_item: -anticodon_item[1])[0][0]
        seeds = [seed for seed in seeds if seed.anticodon_string]


        self.progress.update_pid(pid)
        self.progress.update("Calculating coverages")
        self.sample_total_mean_spec_cov_dict = sample_total_mean_spec_cov_dict = defaultdict(int)
        self.sample_total_discriminator_spec_cov_dict = sample_total_discriminator_spec_cov_dict = defaultdict(int)
        for seed in seeds:
            # Find specific coverages of M by summing specific coverages of nts.
            for summary_M in seed.summaries_M:
                summary_M.spec_covs = np.array([covs for covs in summary_M.spec_nt_covs_dict.values()]).sum(axis=0)

            # Sum specific coverages from each seq comprising a seed. Calculate the mean specific
            # coverage of the seed.
            seed.total_spec_covs = np.zeros(len(seed.string), int)
            for summary_Nu in seed.summaries_Nu:
                seed.total_spec_covs += summary_Nu.spec_covs
                sample_total_mean_spec_cov_dict[summary_Nu.sample_id] += summary_Nu.mean_spec_cov
                sample_total_discriminator_spec_cov_dict[summary_Nu.sample_id] += summary_Nu.spec_covs[-1]
            for summary_M in seed.summaries_M:
                seed.total_spec_covs += summary_M.spec_covs
                sample_total_mean_spec_cov_dict[summary_M.sample_id] += summary_M.spec_covs.mean()
                sample_total_discriminator_spec_cov_dict[summary_M.sample_id] += summary_M.spec_covs[-1]

            seed.total_mean_spec_cov = seed.total_spec_covs.mean()
        # Select the top seeds by specific coverage.
        seeds = sorted(seeds, key=lambda seed: -seed.total_mean_spec_cov)[: self.seed_limit]


        for seed in seeds:
            # Find nonspecific coverages of M by summing nonspecific coverage of nts.
            for summary_M in seed.summaries_M:
                summary_M.nonspec_covs = np.array([covs for covs in summary_M.nonspec_nt_covs_dict.values()]).sum(axis=0)

            # Sum nonspecific coverages from each seq comprising a seed. Calculate the mean
            # nonspecific coverage of the seed.
            seed.total_nonspec_covs = np.zeros(len(seed.string), int)
            for summary_Nu in seed.summaries_Nu:
                seed.total_nonspec_covs += summary_Nu.nonspec_covs
            for summary_M in seed.summaries_M:
                nonspec_nt_covs = np.array([covs for covs in summary_M.nonspec_nt_covs_dict.values()])
                seed.total_nonspec_covs += nonspec_nt_covs.sum(axis=0)

            seed.total_mean_nonspec_cov = seed.total_nonspec_covs.mean()


        # The consensus sequence of a seed consists of the nts with the maximum specific coverage
        # summed across constituent sequences. When certain tRNA-seq treatments are preferred (e.g.,
        # demethylase), the chosen nts at substitution sites are on the basis of seqs from preferred
        # samples.
        for seed in seeds:
            if not seed.summaries_M:
                # The seed is comprised entirely of Nu, which must be subsequences of one another,
                # so there is no variation in nt composition at any position. The sequence string
                # was already assigned and does not need to be altered.
                continue

            total_nt_cov_dict = {nt: np.zeros(len(seed.string), int) for nt in UNAMBIG_NTS}
            for summary_Nu in seed.summaries_Nu:
                for nt in UNAMBIG_NTS:
                    total_nt_cov_dict[nt] += summary_Nu.spec_nt_covs_dict[nt]
            for summary_M in seed.summaries_M:
                for nt in UNAMBIG_NTS:
                    total_nt_cov_dict[nt] += summary_M.spec_nt_covs_dict[nt]

            seed.string = ''.join([INT_NT_DICT[i + 1] for i in np.argmax(np.array([total_nt_cov_dict[nt] for nt in UNAMBIG_NTS]), axis=0)])

            # Set GC fraction of the seed.
            seed.gc_fraction = sum([1 for nt in seed.string if nt == 'C' or nt == 'G']) / len(seed.string)
        self.seeds = seeds
        self.total_seed_length = sum([len(seed.string) for seed in seeds])

        self.set_sample_covs()
        self.progress.update_pid(pid)
        self.progress.update("Setting mod-induced substitutions")
        self.set_substitutions()
        self.progress.update_pid(pid)
        self.progress.update("Setting mod-induced indels")
        self.set_sample_indels()

        self.progress.end()
        self.run.info("Candidate seeds with feature conflicts", len(names_seeds_with_conflict), nl_before=1, nl_after=1)


    def elongate_modified_sequence_fiveprime(self, summary_M, elongation_5prime):
        """Seed sequences can be longer than the normalized sequences from the individual samples,
        requiring addition of empty positions in the normalized sequence coverage arrays at the 5'
        end, as normalized (and modified) sequences are aligned from the 3' end."""
        for summary_Nb in summary_M.summaries_Nb:
            self.elongate_normalized_sequence_fiveprime(summary_Nb, elongation_5prime)

        elongation_length = elongation_5prime.size
        # The positions of substitutions are recorded in the seed sequence index.
        summary_M.sub_positions += elongation_length
        for nt in UNAMBIG_NTS:
            summary_M.spec_nt_covs_dict[nt] = np.concatenate([elongation_5prime, summary_M.spec_nt_covs_dict[nt]])
            summary_M.nonspec_nt_covs_dict[nt] = np.concatenate([elongation_5prime, summary_M.nonspec_nt_covs_dict[nt]])
        summary_M.insert_starts += elongation_length
        summary_M.del_starts += elongation_length


    def shorten_normalized_sequence_fiveprime(self, summary_N, reduction_5prime):
        """Remove some 5' nucleotides from the normalized sequence and recompute certain
        attributes."""
        new_feature_dict = {}
        for feature, feature_index in summary_N.feature_dict.items():
            try:
                new_feature_dict[feature] = feature_index - reduction_5prime
            except TypeError:
                new_feature_dict[feature] = (feature_index[0] - reduction_5prime, feature_index[1] - reduction_5prime)
        summary_N.feature_dict = new_feature_dict
        summary_N.feature_threshold_start -= reduction_5prime

        summary_N.spec_covs = summary_N.spec_covs[reduction_5prime: ]
        summary_N.nonspec_covs = summary_N.nonspec_covs[reduction_5prime: ]

        spec_nt_covs_dict = summary_N.spec_nt_covs_dict
        for nt, spec_covs in summary_N.spec_nt_covs_dict.items():
            spec_nt_covs_dict[nt] = spec_covs[reduction_5prime: ]
        nonspec_nt_covs_dict = summary_N.nonspec_nt_covs_dict
        for nt, spec_covs in summary_N.nonspec_nt_covs_dict.items():
            nonspec_nt_covs_dict[nt] = spec_covs[reduction_5prime: ]


    def shorten_modified_sequence_fiveprime(self, summary_M, reduction_5prime):
        """Remove some 5' nucleotides from the modified sequence and recompute certain
        attributes."""
        for summary_Nb in summary_M.summaries_Nb:
            self.shorten_normalized_sequence_fiveprime(summary_Nb, reduction_5prime)

        sub_positions = summary_M.sub_positions - reduction_5prime
        summary_M.sub_positions = sub_positions[np.where(sub_positions >= 0)[0][0]: ]

        # The following indel adjustments are very crude because Ni are not currently tracked.
        if summary_M.insert_starts.size:
            insert_starts = summary_M.insert_starts - reduction_5prime
            new_insert_starts = np.where(insert_starts >= 0)[0]
            summary_M.insert_starts = new_insert_starts
            if new_insert_starts.size:
                first_retained_insert_index = new_insert_starts[0]
                summary_M.insert_strings = summary_M.insert_strings[first_retained_insert_index: ]
                summary_M.spec_insert_covs = summary_M.spec_insert_covs[first_retained_insert_index: ]
                summary_M.nonspec_insert_covs = summary_M.nonspec_insert_covs[first_retained_insert_index: ]
            else:
                summary_M.insert_strings = []
                summary_M.spec_insert_covs = np.zeros(0, dtype=int)
                summary_M.nonspec_insert_covs = np.zeros(0, dtype=int)

        if summary_M.del_starts.size:
            del_starts = summary_M.del_starts - reduction_5prime
            new_del_starts = np.where(del_starts >= 0)[0]
            summary_M.del_starts = new_del_starts
            if new_del_starts.size:
                first_retained_del_index = new_del_starts[0]
                summary_M.del_lengths = summary_M.del_lengths[first_retained_del_index: ]
                summary_M.spec_del_covs = summary_M.spec_del_covs[first_retained_del_index: ]
                summary_M.nonspec_del_covs = summary_M.nonspec_del_covs[first_retained_del_index: ]
            else:
                summary_M.del_lengths = []
                summary_M.spec_del_covs = np.zeros(0, dtype=int)
                summary_M.nonspec_del_covs = np.zeros(0, dtype=int)

        # Find the nt coverages of M from the Nb remaining in M.
        seq_length = summary_M.spec_nt_covs_dict['A'].size - reduction_5prime
        summary_M.spec_nt_covs_dict = spec_nt_covs_dict = {nt: np.zeros(seq_length, int) for nt in UNAMBIG_NTS}
        summary_M.nonspec_nt_covs_dict = nonspec_nt_covs_dict = {nt: np.zeros(seq_length, int) for nt in UNAMBIG_NTS}
        for summary_Nb in summary_M.summaries_Nb:
            for nt in UNAMBIG_NTS:
                spec_nt_covs_dict[nt] += summary_Nb.spec_nt_covs_dict[nt]
                nonspec_nt_covs_dict[nt] += summary_Nb.nonspec_nt_covs_dict[nt]
        summary_M.spec_covs = sum(spec_nt_covs_dict.values())
        summary_M.nonspec_covs = sum(nonspec_nt_covs_dict.values())

        # Set a consensus seq using the nts with the highest specific cov at each sub position.
        summary_M.consensus_string = ''.join(
            [INT_NT_DICT[nt_int + 1] for nt_int in
             np.argmax(np.stack(tuple(spec_nt_covs_dict.values()), axis=0), axis=0)])


    def set_sample_covs(self):
        """Determine coverages of seeds in each sample. Specific, nonspecific and summed coverages
        are found for each nucleotide, as well as overall and for Q2-Q3."""
        for seed in self.seeds:
            sample_spec_covs_dict = {}
            sample_nonspec_covs_dict = {}
            sample_summed_covs_dict = {}
            sample_spec_nt_covs_dict = {}
            sample_nonspec_nt_covs_dict = {}
            sample_summed_nt_covs_dict = {}
            seed_length = len(seed.string)

            for sample_id in self.trnaseq_db_sample_ids:
                sample_spec_covs_dict[sample_id] = np.zeros(seed_length, int)
                sample_nonspec_covs_dict[sample_id] = np.zeros(seed_length, int)

                sample_spec_nt_covs_dict[sample_id] = [np.zeros(seed_length, int) for _ in UNAMBIG_NTS]
                sample_nonspec_nt_covs_dict[sample_id] = [np.zeros(seed_length, int) for _ in UNAMBIG_NTS]

            for summary_Nu in seed.summaries_Nu:
                sample_id = summary_Nu.sample_id

                sample_spec_covs_dict[sample_id] += summary_Nu.spec_covs
                sample_nonspec_covs_dict[sample_id] += summary_Nu.nonspec_covs

                for nt_num, nt in enumerate(UNAMBIG_NTS):
                    sample_spec_nt_covs_dict[sample_id][nt_num] += summary_Nu.spec_nt_covs_dict[nt]
                    sample_nonspec_nt_covs_dict[sample_id][nt_num] += summary_Nu.nonspec_nt_covs_dict[nt]

            for summary_M in seed.summaries_M:
                sample_id = summary_M.sample_id

                sample_spec_covs_dict[sample_id] += summary_M.spec_covs
                sample_nonspec_covs_dict[sample_id] += summary_M.nonspec_covs

                for nt_num, nt in enumerate(UNAMBIG_NTS):
                    sample_spec_nt_covs_dict[sample_id][nt_num] += summary_M.spec_nt_covs_dict[nt]
                    sample_nonspec_nt_covs_dict[sample_id][nt_num] += summary_M.nonspec_nt_covs_dict[nt]

            for sample_id in self.trnaseq_db_sample_ids:
                sample_summed_covs_dict[sample_id] = sample_spec_covs_dict[sample_id] + sample_nonspec_covs_dict[sample_id]

                sample_spec_nt_covs = sample_spec_nt_covs_dict[sample_id]
                sample_nonspec_nt_covs = sample_nonspec_nt_covs_dict[sample_id]
                sample_summed_nt_covs_dict[sample_id] = [spec_covs + nonspec_covs for spec_covs, nonspec_covs
                                                         in zip(sample_spec_nt_covs, sample_nonspec_nt_covs)]

            seed.sample_spec_covs_dict = sample_spec_covs_dict
            seed.sample_nonspec_covs_dict = sample_nonspec_covs_dict
            seed.sample_summed_covs_dict = sample_summed_covs_dict
            seed.sample_spec_nt_covs_dict = sample_spec_nt_covs_dict
            seed.sample_nonspec_nt_covs_dict = sample_nonspec_nt_covs_dict
            seed.sample_summed_nt_covs_dict = sample_summed_nt_covs_dict

            quartile = int(seed_length * 0.25)
            seed.sample_mean_Q2Q3_spec_cov_dict = {}
            seed.sample_mean_Q2Q3_nonspec_cov_dict = {}
            seed.sample_mean_Q2Q3_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                seed.sample_mean_Q2Q3_spec_cov_dict[sample_id] = np.mean(
                    sorted(sample_spec_covs_dict[sample_id])[quartile: -quartile])
                seed.sample_mean_Q2Q3_nonspec_cov_dict[sample_id] = np.mean(
                    sorted(sample_nonspec_covs_dict[sample_id])[quartile: -quartile])
                seed.sample_mean_Q2Q3_summed_cov_dict[sample_id] = np.mean(
                    sorted(sample_summed_covs_dict[sample_id])[quartile: -quartile])


    def set_substitutions(self):
        """Predict positions with modification-induced substitutions in the tRNA seed.

        A sub requires a certain level of third- and/or fourth-most abundant nts at the position in
        ≥1 sample. A sub in any particular sample additionally requires a certain level of second-
        through fourth-most abundant nts at the position.

        There is currently an idiosyncracy in how subs are set that results in the retention, but
        potential masking, of SNVs. If the position of a potential modification does not meet the
        coverage threshold for third- and fourth-most abundant nucleotides, the seed is not split
        into separate seeds around those SNVs, as occurs in `anvi-trnaseq`. Instead, the SNVs are
        simply not reported. This is a downside to imposing the aforementioned coverage
        threshold."""
        # Division by zero issues a numpy warning, but we handle it immediately by converting the
        # nan result to zero so that a warning is not produced. Unfortunately, this is the only way
        # to suppress the warning that is known by the humble developer.
        np.seterr(invalid='ignore')
        min_variation = self.min_variation
        min_third_fourth_nt = self.min_third_fourth_nt
        for seed in self.seeds:
            sample_sub_positions_dict = {}
            sample_variability_dict = {}
            sample_spec_nt_covs_dict = seed.sample_spec_nt_covs_dict
            seed_length = len(seed.string)
            sample_variations = []
            third_fourth_variations = np.zeros(seed_length)
            for sample_id, spec_nt_covs in sample_spec_nt_covs_dict.items():
                spec_nt_covs_array = np.array(spec_nt_covs)
                spec_nt_covs_array.sort(axis=0)
                first_covs = spec_nt_covs_array[-1, :]
                second_covs = spec_nt_covs_array[-2, :]
                summed_covs = spec_nt_covs_array.sum(axis=0)
                sample_variations.append(np.nan_to_num(1 - first_covs / summed_covs))
                third_fourth_variations += np.nan_to_num(1 - (first_covs + second_covs) / summed_covs) >= min_third_fourth_nt
            sample_variations = np.array(sample_variations)
            third_fourth_variations = (third_fourth_variations > 0)
            total_sub_positions = np.nonzero((sample_variations >= min_variation).any(axis=0) & third_fourth_variations)[0]
            sub_sample_variations = sample_variations[:, total_sub_positions]
            for sample_num, sample_id in enumerate(sample_spec_nt_covs_dict.keys()):
                sample_sub_positions = total_sub_positions[np.nonzero(sub_sample_variations[sample_num, :] >= min_variation)[0]]
                sample_sub_positions_dict[sample_id] = sample_sub_positions.tolist()
                sample_variability_dict[sample_id] = sample_sub_positions.size * 1000 / seed_length
            seed.total_sub_positions = total_sub_positions.tolist()
            seed.sample_sub_positions_dict = sample_sub_positions_dict
            seed.sample_variability_dict = sample_variability_dict
        np.seterr(invalid='warn')

        if self.preferred_treatment:
            self.set_consensus_substitution_nucleotides()


    def set_consensus_substitution_nucleotides(self):
        """Change predicted nucleotides in seed consensus sequences to those supported by the
        samples with the preferred treatment (e.g., demethylase splits) with the goal of increasing
        the accuracy of the underlying base call."""
        for seed in self.seeds:
            seed_string = seed.string
            preferred_nt_cov_dict = {nt: np.zeros(len(seed_string), int) for nt in UNAMBIG_NTS}
            for summary_Nu in seed.summaries_Nu:
                if summary_Nu.sample_id in self.preferred_trnaseq_db_sample_ids:
                    for nt in UNAMBIG_NTS:
                        preferred_nt_cov_dict[nt] += summary_Nu.spec_nt_covs_dict[nt]
            for summary_M in seed.summaries_M:
                if summary_M.sample_id in self.preferred_trnaseq_db_sample_ids:
                    for nt in UNAMBIG_NTS:
                        preferred_nt_cov_dict[nt] += summary_M.spec_nt_covs_dict[nt]

            preferred_nt_cov_array = np.array([preferred_nt_cov_dict[nt] for nt in UNAMBIG_NTS])
            for sub_pos in seed.total_sub_positions:
                sub_covs = preferred_nt_cov_array[:, sub_pos]
                if sub_covs.sum() == 0:
                    # The preferred treatments do not have specific coverage of the sub site.
                    continue
                seed_string = seed_string[: sub_pos] + INT_NT_DICT[np.argmax(sub_covs) + 1] + seed_string[sub_pos + 1: ]
            seed.string = seed_string


    def set_sample_indels(self):
        for seed in self.seeds:
            sample_insert_dict = {}
            sample_del_dict = {}

            for sample_id in self.trnaseq_db_sample_ids:
                sample_insert_dict[sample_id] = []
                sample_del_dict[sample_id] = []

            for summary_M in seed.summaries_M:
                sample_id = summary_M.sample_id

                sample_insert_info = sample_insert_dict[sample_id]
                for insert_start, insert_string, spec_insert_cov, nonspec_insert_cov in zip(summary_M.insert_starts,
                                                                                            summary_M.insert_strings,
                                                                                            summary_M.spec_insert_covs,
                                                                                            summary_M.nonspec_insert_covs):
                    sample_insert_info.append((insert_start, insert_string, spec_insert_cov, nonspec_insert_cov))

                sample_del_info = sample_del_dict[sample_id]
                for del_start, del_length, spec_del_cov, nonspec_del_cov in zip(summary_M.del_starts,
                                                                                summary_M.del_lengths,
                                                                                summary_M.spec_del_covs,
                                                                                summary_M.nonspec_del_covs):
                    sample_del_info.append((del_start, del_length, spec_del_cov, nonspec_del_cov))

            seed.sample_insert_dict = sample_insert_dict
            seed.sample_del_dict = sample_del_dict


    def generate_contigs_database(self):
        """Generate a contigs database of tRNA seeds. The create method of `dbops.ContigsDatabase`
        is not used because it tries to call genes, count kmers, and do other things that are
        irrelevant to tRNA-seq reads. There are no tRNA splits, but to satisfy the structure of the
        database, call every contig a split, and maintain tables for both contigs and splits."""
        self.progress.new("Generating a contigs db of tRNA seeds")
        self.progress.update("...")

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path)
        contigs_db.touch('trnaseq')

        set_meta_value = contigs_db.db.set_meta_value
        # Meta-values are set like in `dbops.ContigsDatabase.create`.
        set_meta_value('db_type', 'contigs')
        set_meta_value('db_variant', 'trnaseq')
        set_meta_value('project_name', self.project_name)
        set_meta_value('description', self.descrip if self.descrip else '_No description is provided_')
        set_meta_value('contigs_db_hash', self.contigs_db_hash)
        set_meta_value('split_length', 10000) # sys.maxsize
        set_meta_value('num_contigs', len(self.seeds))
        set_meta_value('num_splits', len(self.seeds))
        set_meta_value('total_length', self.total_seed_length)
        set_meta_value('kmer_size', 0)
        set_meta_value('gene_level_taxonomy_source', None)
        set_meta_value('gene_function_sources', 'Transfer_RNAs')
        set_meta_value('genes_are_called', True)
        set_meta_value('external_gene_calls', True)
        set_meta_value('external_gene_amino_acid_seqs', False)
        set_meta_value('skip_predict_frame', True)
        set_meta_value('splits_consider_gene_calls', False)
        set_meta_value('scg_taxonomy_was_run', False)
        set_meta_value('scg_taxonomy_database_version', None)
        set_meta_value('trna_taxonomy_was_run', False)
        set_meta_value('trna_taxonomy_database_version', None)
        set_meta_value('creation_date', time.time())

        insert_many = contigs_db.db.insert_many
        insert_many('contig_sequences', [(seed.name, seed.string) for seed in self.seeds])
        insert_many('contigs_basic_info', self.get_contigs_basic_info_table_entries())
        insert_many('splits_basic_info', self.get_splits_basic_info_table_entries())
        insert_many('hmm_hits', self.get_hmm_hits_table_entries())
        insert_many('hmm_hits_in_splits', self.get_hmm_hits_in_splits_table_entries())
        # tRNA predictions are treated like HMM or tRNAScan-SE hits. The blank columns of the HMM
        # hits info table are 'ref', 'search_type', 'domain' and 'genes'.
        contigs_db.db.insert('hmm_hits_info', ('Transfer_RNAs', '', 'Transfer_RNAs', None, ''))
        insert_many('genes_in_contigs', self.get_genes_in_contigs_table_entries())
        insert_many('gene_amino_acid_sequences', [(seed_num, '') for seed_num in range(len(self.seeds))])
        insert_many('genes_in_splits', self.get_genes_in_splits_table_entries())
        insert_many('gene_functions', self.get_gene_functions_table_entries())
        insert_many('trna_feature', self.get_trna_feature_table_entries())

        contigs_db.disconnect()

        self.progress.end()


    def get_contigs_basic_info_table_entries(self):
        entries = []
        for seed in self.seeds:
            seed_string = seed.string
            entries.append(
                (seed.name,
                 len(seed_string),
                 seed.gc_fraction,
                 1)
            )
        return entries


    def get_splits_basic_info_table_entries(self):
        entries = []
        for seed in self.seeds:
            entries.append(
                (seed.name + '_split_00001',
                 0, # Order of split in parent contig
                 0, # Start in contig
                 len(seed.string), # Stop in contig
                 len(seed.string), # Split length
                 seed.gc_fraction, # GC content of split
                 seed.gc_fraction, # GC content of parent contig
                 seed.name)
            )
        return entries


    def get_hmm_hits_table_entries(self):
        """tRNA seeds are analogous to tRNA gene predictions from a metagenomic contigs database."""
        entries = []
        for seed_num, seed in enumerate(self.seeds):
            entries.append(
                (seed_num, # Entry ID
                 'Transfer_RNAs', # Source, à la tRNA gene prediction via tRNAScan-SE
                 sha1(seed.string.encode('utf-8')).hexdigest(), # "Gene unique identifier"
                 seed_num, # "Gene callers ID"
                 ANTICODON_AA_DICT[seed.anticodon_string] + '_' + seed.anticodon_string, # "Gene name", à la tRNA gene prediction via tRNAScan-SE
                 '-', # "Gene HMM ID"
                 0.0) # "HMM E-value"
            )
        return entries


    def get_hmm_hits_in_splits_table_entries(self):
        entries = []
        for seed_num, seed in enumerate(self.seeds):
            entries.append(
                (seed_num, # Entry ID
                 seed.name + '_split_00001', # Split name
                 100, # Percentage of "HMM hit" in split
                 'Transfer_RNAs')
            )
        return entries


    def get_genes_in_contigs_table_entries(self):
        entries = []
        for seed_num, seed in enumerate(self.seeds):
            entries.append(
                (seed_num, # Gene callers ID
                 seed.name, # Contig name
                 0, # Gene start in contig
                 len(seed.string), # Gene stop in contig
                 'f', # Direction of gene call on contig
                 0, # Is partial gene call: for now, say all seeds are "full tRNAs"
                 2, # Call type: 1 = coding, 2 = noncoding, 3 = unknown
                 'anvi-trnaseq', # Gene caller
                 tables.trnaseq_db_version) # Version of caller
            )
        return entries


    def get_genes_in_splits_table_entries(self):
        entries = []
        for seed_num, seed in enumerate(self.seeds):
            entries.append(
                (seed.name + '_split_00001',
                 seed_num,
                 0,
                 len(seed.string),
                 100)
            )
        return entries


    def get_gene_functions_table_entries(self):
        entries = []
        for seed_num, seed in enumerate(self.seeds):
            entries.append(
                (seed_num,
                 'Transfer_RNAs',
                 '%s_%s_%d' % (ANTICODON_AA_DICT[seed.anticodon_string], seed.anticodon_string, seed_num),
                 'tRNA transcript',
                 0.0)
            )
        return entries


    def get_trna_feature_table_entries(self):
        entries = []
        for seed_num, seed in enumerate(self.seeds):
            entry = [seed_num]
            for feature, feature_index in seed.feature_dict.items():
                try:
                    if feature_index >= -10:
                        entry.append(feature_index)
                    else:
                        entry.append(None)
                except TypeError:
                    entry.append(','.join(map(str, feature_index)))
            entries.append(tuple(entry))
        return entries


    def generate_auxiliary_database(self, db_cov_type):
        if db_cov_type == 'specific':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.spec_auxiliary_db_path, self.contigs_db_hash, db_variant='trnaseq', create_new=True)
            for seed in self.seeds:
                split_name = seed.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name, sample_id, seed.sample_spec_covs_dict[sample_id])
        elif db_cov_type == 'nonspecific':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.nonspec_auxiliary_db_path, self.contigs_db_hash, db_variant='trnaseq', create_new=True)
            for seed in self.seeds:
                split_name = seed.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name, sample_id, seed.sample_nonspec_covs_dict[sample_id])
        elif db_cov_type == 'combined':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.combined_auxiliary_db_path, self.contigs_db_hash, db_variant='trnaseq', create_new=True)
            for seed in self.seeds:
                split_name = seed.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name, sample_id + '_specific', seed.sample_spec_covs_dict[sample_id])
                    auxiliary_db.append(split_name, sample_id + '_nonspecific', seed.sample_nonspec_covs_dict[sample_id])
        elif db_cov_type == 'summed':
            auxiliary_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.summed_auxiliary_db_path, self.contigs_db_hash, db_variant='trnaseq', create_new=True)
            for seed in self.seeds:
                split_name = seed.name + '_split_00001'
                for sample_id in self.trnaseq_db_sample_ids:
                    auxiliary_db.append(split_name, sample_id, (seed.sample_spec_covs_dict[sample_id] + seed.sample_nonspec_covs_dict[sample_id]))
        else:
            raise ConfigError(f"The type of profile database provided, {db_cov_type}, "
                              "is not among those that are recognized: 'specific', 'nonspecific', 'combined', and 'summed'.")

        auxiliary_db.store()
        auxiliary_db.close()


    def set_sample_total_coverages(self):
        """For each input sample, find the total specific, nonspecific and summed coverage of the
        seeds across all positions (single integers)."""
        sample_total_spec_cov_dict = {sample_id: 0 for sample_id in self.trnaseq_db_sample_ids}
        sample_total_nonspec_cov_dict = {sample_id: 0 for sample_id in self.trnaseq_db_sample_ids}
        for seed in self.seeds:
            for summary_Nu in seed.summaries_Nu:
                sample_id = summary_Nu.sample_id
                spec_cov = summary_Nu.spec_covs.sum()
                nonspec_cov = summary_Nu.nonspec_covs.sum()
                sample_total_spec_cov_dict[sample_id] += spec_cov
                sample_total_nonspec_cov_dict[sample_id] += nonspec_cov
            for summary_M in seed.summaries_M:
                sample_id = summary_M.sample_id
                spec_cov = summary_M.spec_covs.sum()
                nonspec_cov = summary_M.nonspec_covs.sum()
                sample_total_spec_cov_dict[sample_id] += spec_cov
                sample_total_nonspec_cov_dict[sample_id] += nonspec_cov

        sample_total_summed_cov_dict = {}
        for sample_id in self.trnaseq_db_sample_ids:
            sample_total_summed_cov_dict[sample_id] = sample_total_spec_cov_dict[sample_id] + sample_total_nonspec_cov_dict[sample_id]

        self.sample_total_spec_cov_dict = sample_total_spec_cov_dict
        self.sample_total_nonspec_cov_dict = sample_total_nonspec_cov_dict
        self.sample_total_summed_cov_dict = sample_total_summed_cov_dict


    def set_sample_overall_mean_coverages(self):
        """For each input sample, find the mean specific, nonspecific and summed coverage of all
        seeds across all positions (single numbers)."""
        sample_overall_mean_spec_cov_dict = {}
        sample_overall_mean_nonspec_cov_dict = {}
        sample_overall_mean_summed_cov_dict = {}
        for sample_id, total_spec_cov in self.sample_total_spec_cov_dict.items():
            sample_overall_mean_spec_cov_dict[sample_id] = total_spec_cov / self.total_seed_length
        for sample_id, total_nonspec_cov in self.sample_total_nonspec_cov_dict.items():
            sample_overall_mean_nonspec_cov_dict[sample_id] = total_nonspec_cov / self.total_seed_length
        for sample_id, total_summed_cov in self.sample_total_summed_cov_dict.items():
            sample_overall_mean_summed_cov_dict[sample_id] = total_summed_cov / self.total_seed_length

        self.sample_overall_mean_spec_cov_dict = sample_overall_mean_spec_cov_dict
        self.sample_overall_mean_nonspec_cov_dict = sample_overall_mean_nonspec_cov_dict
        self.sample_overall_mean_summed_cov_dict = sample_overall_mean_summed_cov_dict


    def set_sample_mean_coverages(self):
        """Set the mean coverage of each seed in a sample."""
        for seed in self.seeds:
            sample_mean_spec_cov_dict = {}
            sample_mean_nonspec_cov_dict = {}
            sample_mean_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_mean_spec_cov_dict[sample_id] = seed.sample_spec_covs_dict[sample_id].mean()
                sample_mean_nonspec_cov_dict[sample_id] = seed.sample_nonspec_covs_dict[sample_id].mean()
                sample_mean_summed_cov_dict[sample_id] = seed.sample_summed_covs_dict[sample_id].mean()
            seed.sample_mean_spec_cov_dict = sample_mean_spec_cov_dict
            seed.sample_mean_nonspec_cov_dict = sample_mean_nonspec_cov_dict
            seed.sample_mean_summed_cov_dict = sample_mean_summed_cov_dict


    def set_sample_coverage_standard_deviations(self):
        """Set the standard deviation of the coverage of each seed in a sample."""
        for seed in self.seeds:
            sample_std_spec_cov_dict = {}
            sample_std_nonspec_cov_dict = {}
            sample_std_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_std_spec_cov_dict[sample_id] = seed.sample_spec_covs_dict[sample_id].std()
                sample_std_nonspec_cov_dict[sample_id] = seed.sample_nonspec_covs_dict[sample_id].std()
                sample_std_summed_cov_dict[sample_id] = seed.sample_summed_covs_dict[sample_id].std()
            seed.sample_std_spec_cov_dict = sample_std_spec_cov_dict
            seed.sample_std_nonspec_cov_dict = sample_std_nonspec_cov_dict
            seed.sample_std_summed_cov_dict = sample_std_summed_cov_dict


    def set_sample_abundances(self):
        """For each sample, and for specific and nonspecific coverages, abundance is defined as the
        mean coverage of the seed divided by the mean total coverage of the sample across all
        seeds."""
        for seed in self.seeds:
            sample_spec_abund_dict = {}
            sample_nonspec_abund_dict = {}
            sample_summed_abund_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_spec_abund_dict[sample_id] = seed.sample_mean_spec_cov_dict[sample_id] / self.sample_overall_mean_spec_cov_dict[sample_id]
                sample_nonspec_abund_dict[sample_id] = seed.sample_mean_nonspec_cov_dict[sample_id] / self.sample_overall_mean_nonspec_cov_dict[sample_id]
                sample_summed_abund_dict[sample_id] = seed.sample_mean_summed_cov_dict[sample_id] / self.sample_overall_mean_summed_cov_dict[sample_id]
            seed.sample_spec_abund_dict = sample_spec_abund_dict
            seed.sample_nonspec_abund_dict = sample_nonspec_abund_dict
            seed.sample_summed_abund_dict = sample_summed_abund_dict


    def set_sample_normalization_multipliers(self):
        """Set a normalization constant for each sample to scale their coverages, allowing the
        relative abundance of seeds in a sample to be compared between samples. Normalization is
        based on the total specific coverage of each sample -- one can imagine other ways of doing
        this, including use of summed specific and nonspecific coverage, but this would require
        deconvoluting the multiple representation of nonspecific reads."""
        sample_normalization_multiplier_dict = {}
        min_total_spec_cov = min([v for v in self.sample_total_spec_cov_dict.values()])
        for sample_id, total_spec_cov in self.sample_total_spec_cov_dict.items():
            sample_normalization_multiplier_dict[sample_id] = min_total_spec_cov / total_spec_cov
        self.sample_normalization_multiplier_dict = sample_normalization_multiplier_dict


    def set_sample_normalized_mean_Q2Q3_coverages(self):
        """Scale mean coverages for comparison across samples."""
        for seed in self.seeds:
            sample_normalized_mean_Q2Q3_spec_cov_dict = {}
            sample_normalized_mean_Q2Q3_nonspec_cov_dict = {}
            sample_normalized_mean_Q2Q3_summed_cov_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_normalized_mean_Q2Q3_spec_cov_dict[sample_id] = seed.sample_mean_Q2Q3_spec_cov_dict[sample_id] * self.sample_normalization_multiplier_dict[sample_id]
                sample_normalized_mean_Q2Q3_nonspec_cov_dict[sample_id] = seed.sample_mean_Q2Q3_nonspec_cov_dict[sample_id] * self.sample_normalization_multiplier_dict[sample_id]
                sample_normalized_mean_Q2Q3_summed_cov_dict[sample_id] = seed.sample_mean_Q2Q3_summed_cov_dict[sample_id] * self.sample_normalization_multiplier_dict[sample_id]
            seed.sample_normalized_mean_Q2Q3_spec_cov_dict = sample_normalized_mean_Q2Q3_spec_cov_dict
            seed.sample_normalized_mean_Q2Q3_nonspec_cov_dict = sample_normalized_mean_Q2Q3_nonspec_cov_dict
            seed.sample_normalized_mean_Q2Q3_summed_cov_dict = sample_normalized_mean_Q2Q3_summed_cov_dict


    def set_sample_detections(self):
        """Find the proportion of each seed sequence covered by reads in a sample."""
        for seed in self.seeds:
            seed_length = len(seed.string)
            sample_spec_detection_dict = {}
            sample_nonspec_detection_dict = {}
            sample_summed_detection_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_spec_detection_dict[sample_id] = seed.sample_spec_covs_dict[sample_id].nonzero()[0].size / seed_length
                sample_nonspec_detection_dict[sample_id] = seed.sample_nonspec_covs_dict[sample_id].nonzero()[0].size / seed_length
                sample_summed_detection_dict[sample_id] = seed.sample_summed_covs_dict[sample_id].nonzero()[0].size / seed_length
            seed.sample_spec_detection_dict = sample_spec_detection_dict
            seed.sample_nonspec_detection_dict = sample_nonspec_detection_dict
            seed.sample_summed_detection_dict = sample_summed_detection_dict


    def set_sample_relative_abundances(self):
        """Relative abundance represents the coverage of the seed in one sample relative to the
        total coverage of the seed across samples -- relative abundances sum to one across
        samples."""
        np.seterr(invalid='ignore')
        for seed in self.seeds:
            sample_spec_rel_abund_dict = {}
            sample_nonspec_rel_abund_dict = {}
            sample_summed_rel_abund_dict = {}
            pansample_normalized_mean_Q2Q3_spec_cov = sum(seed.sample_normalized_mean_Q2Q3_spec_cov_dict.values())
            pansample_normalized_mean_Q2Q3_nonspec_cov = sum(seed.sample_normalized_mean_Q2Q3_nonspec_cov_dict.values())
            pansample_normalized_mean_Q2Q3_summed_cov = sum(seed.sample_normalized_mean_Q2Q3_summed_cov_dict.values())
            for sample_id in self.trnaseq_db_sample_ids:
                sample_spec_rel_abund_dict[sample_id] = seed.sample_normalized_mean_Q2Q3_spec_cov_dict[sample_id] / pansample_normalized_mean_Q2Q3_spec_cov
                sample_nonspec_rel_abund_dict[sample_id] = seed.sample_normalized_mean_Q2Q3_nonspec_cov_dict[sample_id] / pansample_normalized_mean_Q2Q3_nonspec_cov
                sample_summed_rel_abund_dict[sample_id] = seed.sample_normalized_mean_Q2Q3_summed_cov_dict[sample_id] / pansample_normalized_mean_Q2Q3_summed_cov
            seed.sample_spec_rel_abund_dict = sample_spec_rel_abund_dict
            seed.sample_nonspec_rel_abund_dict = sample_nonspec_rel_abund_dict
            seed.sample_summed_rel_abund_dict = sample_summed_rel_abund_dict
        np.seterr(invalid='warn')


    def set_sample_max_normalized_ratios(self):
        """The max normalized coverage ratio represents the coverage of the seed in one sample
        relative to the max coverage amongst the samples -- one sample will always have a value
        equal to 1."""
        for seed in self.seeds:
            sample_spec_max_normalized_ratio_dict = {}
            sample_nonspec_max_normalized_ratio_dict = {}
            sample_summed_max_normalized_ratio_dict = {}
            max_normalized_mean_Q2Q3_spec_cov = max(seed.sample_normalized_mean_Q2Q3_spec_cov_dict.values())
            max_normalized_mean_Q2Q3_nonspec_cov = max(seed.sample_normalized_mean_Q2Q3_nonspec_cov_dict.values())
            max_normalized_mean_Q2Q3_summed_cov = max(seed.sample_normalized_mean_Q2Q3_summed_cov_dict.values())
            for sample_id in self.trnaseq_db_sample_ids:
                sample_spec_max_normalized_ratio_dict[sample_id] = seed.sample_normalized_mean_Q2Q3_spec_cov_dict[sample_id] / max_normalized_mean_Q2Q3_spec_cov if max_normalized_mean_Q2Q3_spec_cov else 0
                sample_nonspec_max_normalized_ratio_dict[sample_id] = seed.sample_normalized_mean_Q2Q3_nonspec_cov_dict[sample_id] / max_normalized_mean_Q2Q3_nonspec_cov if max_normalized_mean_Q2Q3_nonspec_cov else 0
                sample_summed_max_normalized_ratio_dict[sample_id] = seed.sample_normalized_mean_Q2Q3_summed_cov_dict[sample_id] / max_normalized_mean_Q2Q3_summed_cov if max_normalized_mean_Q2Q3_summed_cov else 0
            seed.sample_spec_max_normalized_ratio_dict = sample_spec_max_normalized_ratio_dict
            seed.sample_nonspec_max_normalized_ratio_dict = sample_nonspec_max_normalized_ratio_dict
            seed.sample_summed_max_normalized_ratio_dict = sample_summed_max_normalized_ratio_dict


    def set_variable_nucleotides_table_entries(self):
        """Variable nucleotides in the profile databases are those with predicted
        modification-induced substitutions, not single nucleotide variants. Subs are determined from
        specific coverage and are currently only reported in the specific coverage profile database.
        (Therefore, `anvi-interactive` does not display subs with "combined" or "summed" coverage
        profile databases.)"""
        entries = []
        for sample_id in self.trnaseq_db_sample_ids:
            for seed_num, seed in enumerate(self.seeds):
                spec_covs = seed.sample_spec_covs_dict[sample_id]
                spec_nt_cov_arrays = seed.sample_spec_nt_covs_dict[sample_id]
                for pos in seed.sample_sub_positions_dict[sample_id]:
                    total_cov = spec_covs[pos]
                    spec_nt_covs = [arr[pos] for arr in spec_nt_cov_arrays]
                    max_nt_cov = max(spec_nt_covs)
                    sorted_nt_covs = sorted(zip(UNAMBIG_NTS, spec_nt_covs), key=lambda nt_item: -nt_item[1])
                    ref_nt = sorted_nt_covs[0][0]
                    secondary_nt = sorted_nt_covs[1][0]
                    entries.append((sample_id,
                                    seed.name + '_split_00001',
                                    pos, # Position in split
                                    pos, # Position in contig
                                    seed_num, # Corresponding gene call
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
                                    spec_nt_covs[0], # A coverage
                                    spec_nt_covs[1], # C coverage
                                    spec_nt_covs[2], # G coverage
                                    spec_nt_covs[3], # T coverage
                                    0))
        self.variable_nts_table_entries = entries


    def set_indels_table_entries(self):
        """Indels are determined separately from specific and nonspecific coverages."""
        spec_entries = []
        nonspec_entries = []
        summed_entries = []
        min_indel_fraction = self.min_indel_fraction
        for sample_id in self.trnaseq_db_sample_ids:
            for seed_num, seed in enumerate(self.seeds):
                insert_info = seed.sample_insert_dict[sample_id]
                del_info = seed.sample_del_dict[sample_id]
                spec_covs = seed.sample_spec_covs_dict[sample_id]
                nonspec_covs = seed.sample_nonspec_covs_dict[sample_id]

                for insert_start, insert_string, insert_spec_cov, insert_nonspec_cov in insert_info:
                    spec_cov = (spec_covs[insert_start] + spec_covs[insert_start + 1]) / 2

                    # A frequency of 1 only occurs when there is specific coverage of the insertion
                    # but not the reference nucleotide.
                    insert_freq = 1 if spec_cov == 0 else insert_spec_cov / spec_cov
                    if insert_freq >= min_indel_fraction:
                        spec_entries.append((sample_id,
                                             seed.name + '_split_00001',
                                             insert_start, # Position in split
                                             insert_start, # Position in contig
                                             seed_num, # Corresponding gene call
                                             1, # In noncoding gene call
                                             0, # In coding gene call
                                             0, # Base position in codon (0 for noncoding gene call)
                                             -1, # Codon order in gene (-1 for noncoding gene call)
                                             0, # Coverage outlier in split (0 or 1)
                                             0, # Coverage outlier in contig (0 or 1)
                                             seed.string[insert_start], # Reference nt
                                             'INS', # Type of indel
                                             insert_string, # Indel sequence ('' for deletion)
                                             len(insert_string), # Indel length
                                             insert_spec_cov, # Deletion count (coverage)
                                             spec_cov)) # Reference sequence coverage

                    nonspec_cov = (nonspec_covs[insert_start] + nonspec_covs[insert_start + 1]) / 2
                    insert_freq = 1 if nonspec_cov == 0 else insert_nonspec_cov / nonspec_cov
                    if insert_freq >= min_indel_fraction:
                        nonspec_entries.append((sample_id,
                                                seed.name + '_split_00001',
                                                insert_start,
                                                insert_start,
                                                seed_num,
                                                1,
                                                0,
                                                0,
                                                -1,
                                                0,
                                                0,
                                                seed.string[insert_start],
                                                'INS',
                                                insert_string,
                                                len(insert_string),
                                                insert_nonspec_cov,
                                                nonspec_cov))

                    sum_cov = spec_cov + nonspec_cov
                    insert_sum_cov = insert_spec_cov + insert_nonspec_cov
                    insert_freq = 1 if sum_cov == 0 else insert_sum_cov / sum_cov
                    if insert_freq >= min_indel_fraction:
                        summed_entries.append((sample_id,
                                               seed.name + '_split_00001',
                                               insert_start,
                                               insert_start,
                                               seed_num,
                                               1,
                                               0,
                                               0,
                                               -1,
                                               0,
                                               0,
                                               seed.string[insert_start],
                                               'INS',
                                               insert_string,
                                               len(insert_string),
                                               insert_sum_cov,
                                               sum_cov))

                for del_start, del_length, del_spec_cov, del_nonspec_cov in del_info:
                    spec_cov = spec_covs[del_start: del_start + del_length].mean()
                    del_freq = 1 if spec_cov == 0 else del_spec_cov / spec_cov
                    if del_freq >= min_indel_fraction:
                        spec_entries.append((sample_id,
                                             seed.name + '_split_00001',
                                             del_start, # Position in split
                                             del_start, # Position in contig
                                             seed_num, # Corresponding gene call
                                             1, # In noncoding gene call
                                             0, # In coding gene call
                                             0, # Base position in codon (0 for noncoding gene call)
                                             -1, # Codon order in gene (-1 for noncoding gene call)
                                             0, # Coverage outlier in split (0 or 1)
                                             0, # Coverage outlier in contig (0 or 1)
                                             seed.string[del_start], # Reference nt
                                             'DEL', # Type of indel
                                             '', # Indel sequence ('' for deletion)
                                             del_length, # Indel length
                                             del_spec_cov, # Deletion count (coverage)
                                             spec_cov)) # Reference sequence coverage

                    nonspec_cov = nonspec_covs[del_start: del_start + del_length].mean()
                    del_freq = 1 if nonspec_cov == 0 else del_nonspec_cov / nonspec_cov
                    if del_freq >= min_indel_fraction:
                        nonspec_entries.append((sample_id,
                                                seed.name + '_split_00001',
                                                del_start,
                                                del_start,
                                                seed_num,
                                                1,
                                                0,
                                                0,
                                                -1,
                                                0,
                                                0,
                                                seed.string[del_start],
                                                'DEL',
                                                '',
                                                del_length,
                                                del_nonspec_cov,
                                                nonspec_cov))

                    sum_cov = spec_cov + nonspec_cov
                    del_sum_cov = del_spec_cov + del_nonspec_cov
                    del_freq = 1 if sum_cov == 0 else del_sum_cov / sum_cov
                    if del_freq >= min_indel_fraction:
                        summed_entries.append((sample_id,
                                               seed.name + '_split_00001',
                                               del_start,
                                               del_start,
                                               seed_num,
                                               1,
                                               0,
                                               0,
                                               -1,
                                               0,
                                               0,
                                               seed.string[del_start],
                                               'DEL',
                                               '',
                                               del_length,
                                               del_sum_cov,
                                               sum_cov))
        self.spec_indels_table_entries = spec_entries
        self.nonspec_indels_table_entries = nonspec_entries
        self.summed_indels_table_entries = summed_entries


    def generate_profile_database(self, db_cov_type):
        self.progress.new(f"Generating {db_cov_type} profile db")
        self.progress.update("...")

        if db_cov_type == 'specific':
            profile_db_path = self.spec_profile_db_path
        elif db_cov_type == 'nonspecific':
            profile_db_path = self.nonspec_profile_db_path
        elif db_cov_type == 'combined':
            profile_db_path = self.combined_profile_db_path
        elif db_cov_type == 'summed':
            profile_db_path = self.summed_profile_db_path
        else:
            raise ConfigError(f"The tRNA-seq coverage type, {db_cov_type}, is not recognized. "
                              "The only valid options are specific, nonspecific, summed, and combined.")

        profile_db = dbops.ProfileDatabase(profile_db_path)
        profile_db.touch()

        set_meta_value = profile_db.db.set_meta_value
        # Profile database meta-values are set in a parallel fashion to `merger.MultipleRuns.merge`.
        set_meta_value('creation_date', time.time())
        set_meta_value('db_type', 'profile')
        set_meta_value('db_variant', 'trnaseq')

        contigs_db = dbops.ContigsDatabase(self.contigs_db_path, quiet=True)
        set_meta_value('contigs_db_hash', contigs_db.meta['contigs_db_hash'])
        set_meta_value('sample_id', contigs_db.meta['project_name'])
        contigs_db.disconnect()

        if db_cov_type == 'combined':
            set_meta_value('samples', ', '.join([sample_id + '_' + cov_type
                                                 for sample_id in self.trnaseq_db_sample_ids
                                                 for cov_type in ('specific', 'nonspecific')]))
        else:
            set_meta_value('samples', ', '.join([sample_id for sample_id in self.trnaseq_db_sample_ids]))
        if db_cov_type == 'specific':
            set_meta_value('sample_total_mean_specific_coverage', ', '.join(map(str, [round(self.sample_total_mean_spec_cov_dict[sample_id], 1)
                                                                                      for sample_id in self.trnaseq_db_sample_ids])))
            set_meta_value('sample_total_discriminator_specific_coverage', ', '.join(map(str, [self.sample_total_discriminator_spec_cov_dict[sample_id]
                                                                                               for sample_id in self.trnaseq_db_sample_ids])))
        # The total number of reads "mapped" is not calculated due to various complexities.
        # set_meta_value('total_reads_mapped', -1)
        set_meta_value('merged', True)
        set_meta_value('blank', False)
        set_meta_value('default_view', 'mean_coverage')
        set_meta_value('min_contig_length', 1)
        set_meta_value('max_contig_length', MAXSIZE)
        set_meta_value('SNVs_profiled', False)
        set_meta_value('SCVs_profiled', False)
        set_meta_value('INDELs_profiled', False)
        set_meta_value('num_contigs', len(self.seeds))
        set_meta_value('num_splits', len(self.seeds))
        set_meta_value('total_length', self.total_seed_length)
        set_meta_value('min_coverage_for_variability', 1)
        set_meta_value('min_indel_fraction', 0)
        set_meta_value('report_variability_full', False)
        set_meta_value('description', self.descrip if self.descrip else '_No description is provided_')
        # set_meta_value('min_percent_identity', -1)

        # Whereas variability in metagenomics refers to SNVs, here it refers to inferred
        # modification-induced substitutions. Subs are only identified from specific coverage.
        if db_cov_type == 'specific':
            profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('variable_nucleotides', ','.join('?' * len(tables.variable_nts_table_structure))),
                                     self.variable_nts_table_entries)
            profile_db.db.commit()

        if db_cov_type == 'specific' or db_cov_type == 'nonspecific' or db_cov_type == 'summed':
            tables_to_create = [('sample_mean_' + db_cov_type.replace('specific', 'spec') + '_cov_dict', 'mean_coverage'),
                                ('sample_std_' + db_cov_type.replace('specific', 'spec') + '_cov_dict', 'std_coverage'),
                                ('sample_' + db_cov_type.replace('specific', 'spec') + '_abund_dict', 'abundance'),
                                ('sample_' + db_cov_type.replace('specific', 'spec') + '_detection_dict', 'detection'),
                                ('sample_mean_Q2Q3_' + db_cov_type.replace('specific', 'spec') + '_cov_dict', 'mean_coverage_Q2Q3')]

            for attr, table_basename in tables_to_create:
                data_dict = self.get_specific_nonspecific_or_summed_data_dict(attr)
                self.create_contigs_and_splits_tables(profile_db_path, table_basename, data_dict)
            # Variability is the measure of the frequency of inferred modification-induced
            # substitutions in seeds. Subs are only calculated from specific coverage -- nonspecific
            # coverage is ignored.
            variability_data_dict = self.get_specific_nonspecific_or_summed_data_dict('sample_variability_dict')
            self.create_contigs_and_splits_tables(profile_db_path, 'variability', data_dict)

            profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('indels', ','.join('?' * len(tables.indels_table_structure))),
                                     getattr(self, db_cov_type.replace('specific', 'spec') + '_indels_table_entries'))
        elif db_cov_type == 'combined':
            tables_to_create = [('sample_mean_spec_cov_dict', 'sample_mean_nonspec_cov_dict', 'mean_coverage'),
                                ('sample_std_spec_cov_dict', 'sample_std_nonspec_cov_dict', 'std_coverage'),
                                ('sample_spec_abund_dict', 'sample_nonspec_abund_dict', 'abundance'),
                                ('sample_spec_detection_dict', 'sample_nonspec_detection_dict', 'detection'),
                                ('sample_mean_Q2Q3_spec_cov_dict', 'sample_mean_Q2Q3_nonspec_cov_dict', 'mean_coverage_Q2Q3')]

            for spec_attr, nonspec_attr, table_basename in tables_to_create:
                data_dict = self.get_combined_data_dict(spec_attr, nonspec_attr)
                self.create_contigs_and_splits_tables(profile_db_path, table_basename, data_dict)
            # Variability is the measure of the frequency of inferred modification-induced
            # substitutions in seeds. Subs are only calculated from specific coverage -- nonspecific
            # coverage is ignored.
            variability_data_dict = self.get_combined_data_dict('sample_variability_dict', 'sample_variability_dict')
            self.create_contigs_and_splits_tables(profile_db_path, 'variability', data_dict)

            combined_indels_table_entries = []
            for entry in self.spec_indels_table_entries:
                combined_indels_table_entries.append((entry[0] + '_specific', ) + entry[1: ])
            for entry in self.nonspec_indels_table_entries:
                combined_indels_table_entries.append((entry[0] + '_nonspecific', ) + entry[1: ])
            profile_db.db._exec_many('''INSERT INTO %s VALUES (%s)'''
                                     % ('indels', ','.join('?' * len(tables.indels_table_structure))),
                                     combined_indels_table_entries)

        profile_db.db.commit()
        self.progress.end()

        # Add layers for anticodon and corresponding amino acid.
        items_additional_data_table = miscdata.MiscDataTableFactory(argparse.Namespace(profile_db=profile_db_path, target_data_table='items'))
        data_dict = {}
        for seed in self.seeds:
            data_dict[seed.name + '_split_00001'] = {'anticodon': seed.anticodon_string,
                                                     'amino_acid': ANTICODON_AA_DICT[seed.anticodon_string]}
        items_additional_data_table.add(data_dict, ['anticodon', 'amino_acid'])

        # Cluster tRNA seeds to form the central dendrogram in anvi-interactive.
        dbops.do_hierarchical_clustering_of_items(profile_db_path,
                                                  constants.clustering_configs['trnaseq'],
                                                  [seed.name + '_split_00001' for seed in self.seeds],
                                                  {'CONTIGS.db': self.contigs_db_path, 'PROFILE.db': profile_db_path},
                                                  input_directory=os.path.dirname(profile_db_path),
                                                  default_clustering_config=constants.trnaseq_default,
                                                  distance=self.distance,
                                                  linkage=self.linkage,
                                                  run=self.run,
                                                  progress=self.progress)
        set_meta_value('items_ordered', True)
        profile_db.db.disconnect()

        # Cluster samples by "view" data to find possible sample layer orderings.
        profile_db_super = dbops.ProfileSuperclass(argparse.Namespace(profile_db=profile_db_path))
        profile_db_super.load_views(omit_parent_column=True)
        layer_orders_data_dict = {}
        failed_attempts = []
        for essential_field in constants.essential_data_fields_for_anvio_profiles:
            try:
                data_value = clustering.get_newick_tree_data_for_dict(profile_db_super.views[essential_field]['dict'],
                                                                      distance=self.distance,
                                                                      linkage=self.linkage,
                                                                      transpose=True)
                layer_orders_data_dict[essential_field] = {'data_value': data_value, 'data_type': 'newick'}
            except:
                failed_attempts.append(essential_field)
        if not len(layer_orders_data_dict):
            self.run.warning("This may or may not be important: "
                             "anvi'o attempted to generate orders for your samples based on the view data. It failed :/")
            return
        if len(failed_attempts):
            self.run.warning(f"While anvi'o was trying to generate clusterings of samples based on view data available in the {db_cov_type} profile, "
                             "clustering of some of the essential data failed. "
                             "It is likely not a very big deal, but you shall be the judge of it. "
                             "Anvi'o now proceeds to store layers order information for those view items the clustering in fact worked. "
                             f"Here is the list of stuff that failed: '{', '.join(failed_attempts)}'")
        # Add the layer orders quietly.
        TableForLayerOrders(argparse.Namespace(profile_db=profile_db_path), r=terminal.Run(verbose=False)).add(layer_orders_data_dict)


    def get_specific_nonspecific_or_summed_data_dict(self, seed_attr):
        """Get data from seeds to generate a table in a specific, nonspecific, or summed profile
        database."""
        data_dict = {}
        for seed in self.seeds:
            data_dict[seed.name + '_split_00001'] = sample_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_dict[sample_id] = getattr(seed, seed_attr)[sample_id]
        return data_dict


    def get_combined_data_dict(self, spec_seed_attr, nonspec_seed_attr):
        """Get data from seeds to generate a table in a combined profile database."""
        data_dict = {}
        for seed in self.seeds:
            data_dict[seed.name + '_split_00001'] = sample_dict = {}
            for sample_id in self.trnaseq_db_sample_ids:
                sample_dict[sample_id + '_specific'] = getattr(seed, spec_seed_attr)[sample_id]
                sample_dict[sample_id + '_nonspecific'] = getattr(seed, nonspec_seed_attr)[sample_id]
        return data_dict


    def create_contigs_and_splits_tables(self, profile_db_path, table_basename, data_dict):
        """Create a pair of tables in a profile database. Contigs and splits tables contain the same
        information since tRNA, unlike a metagenomic contig, is not long enough to be split."""
        TablesForViews(profile_db_path).create_new_view(
            view_data=data_dict,
            table_name=table_basename + '_contigs',
            view_name=None,
            from_matrix_form=True)
        TablesForViews(profile_db_path).create_new_view(
            view_data=data_dict,
            table_name=table_basename + '_splits',
            view_name=table_basename,
            from_matrix_form=True)


def trnaseq_db_loader(input_queue, output_queue_Nu_summaries, output_queue_M_summaries, db_merger):
    """This client for `DatabaseMerger.load_trnaseq_database_sequence_summaries` is located outside the
    `DatabaseMerger` class to allow multiprocessing."""
    while True:
        trnaseq_db_path = input_queue.get()
        summaries_Nu, summaries_M = db_merger.load_trnaseq_database_sequence_summaries(trnaseq_db_path)
        for summary_Nu in summaries_Nu:
            output_queue_Nu_summaries.put((trnaseq_db_path, summary_Nu))
        for summary_M in summaries_M:
            output_queue_M_summaries.put((trnaseq_db_path, summary_M))
        output_queue_Nu_summaries.put(trnaseq_db_path)
        output_queue_M_summaries.put(trnaseq_db_path)


class ResultTabulator(object):

    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # Group 1: MANDATORY
        self.contigs_db_path = A('contigs_db')
        self.spec_profile_db_path = A('specific_profile_db')
        # Group 2: OPTIONAL
        self.nonspec_profile_db_path = A('nonspecific_profile_db')
        self.out_dir = os.path.abspath(A('output_dir')) if A('output_dir') else os.path.dirname(self.contigs_db_path)
        self.overwrite_out_dest = A('overwrite_output_destinations')

        if not self.contigs_db_path:
            raise ConfigError("Please provide the path to a `trnaseq`-variant contigs database using `--contigs-db` or `-c`.")
        if not self.spec_profile_db_path:
            raise ConfigError("Please provide the path to a profile database containing specific coverage information on tRNA seeds using `--specific-profile-db` or `-s`.")


    def process(self):
        self.sanity_check()
        self.contigs_db = dbops.ContigsDatabase(self.contigs_db_path, quiet=True)
        self.contig_names = self.contigs_db.db.get_single_column_from_table(tables.genes_in_contigs_table_name, 'contig')

        self.set_ordinal_features()
        self.generate_seed_output()
        self.generate_modification_output()

        self.contigs_db.disconnect()


    def sanity_check(self):
        """Check `anvi-tabulate-trnaseq` arguments."""
        self.contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')
        if self.contigs_db_info.variant != 'trnaseq':
            raise ConfigError("The input contigs database is not of the `trnaseq` variant...")

        self.spec_profile_db_info = DBInfo(self.spec_profile_db_path, expecting='profile')
        self.nonspec_profile_db_info = DBInfo(self.nonspec_profile_db_path, expecting='profile') if self.nonspec_profile_db_path else None
        self.spec_aux_db_path = os.path.join(os.path.dirname(self.spec_profile_db_path), 'AUXILIARY-DATA.db')
        self.spec_aux_db_info = DBInfo(self.spec_aux_db_path, expecting='auxiliary data for coverages')
        self.nonspec_aux_db_path = os.path.join(os.path.dirname(self.nonspec_profile_db_path), 'AUXILIARY-DATA.db') if self.nonspec_profile_db_path else None
        if self.nonspec_aux_db_path:
            nonspec_aux_db_info = DBInfo(self.nonspec_aux_db_path, expecting='auxiliary data for coverages')

        filesnpaths.is_output_dir_writable(self.out_dir)

        self.spec_seed_out_path = os.path.join(self.out_dir, 'SEEDS_SPECIFIC.txt')
        filesnpaths.is_output_file_writable(self.spec_seed_out_path, ok_if_exists=self.overwrite_out_dest)
        if os.path.exists(self.spec_seed_out_path) and self.overwrite_out_dest:
            os.remove(self.spec_seed_out_path)

        self.nonspec_seed_out_path = os.path.join(self.out_dir, 'SEEDS_NONSPECIFIC.txt') if self.nonspec_profile_db_path else None
        if self.nonspec_seed_out_path:
            filesnpaths.is_output_file_writable(self.nonspec_seed_out_path, ok_if_exists=self.overwrite_out_dest)
            if os.path.exists(self.nonspec_seed_out_path) and self.overwrite_out_dest:
                os.remove(self.nonspec_seed_out_path)

        self.mod_out_path = os.path.join(self.out_dir, 'MODIFICATIONS.txt')
        filesnpaths.is_output_file_writable(self.mod_out_path, ok_if_exists=self.overwrite_out_dest)
        if os.path.exists(self.mod_out_path) and self.overwrite_out_dest:
            os.remove(self.mod_out_path)


    def set_ordinal_features(self):
        # The term "feature" is used loosely in the following dict, as it separates out "subfeatures" of the
        # features just listed, e.g., `d_loop` is split into `d_loop_prealpha`, `d_loop_alpha`, etc.
        self.pretty_ordinal_extrema_dict = OrderedDict([('trna_his_position_0_start', 101),
                                                        ('trna_his_position_0_stop', 101),
                                                        ('fiveprime_acceptor_stem_sequence_start', 201),
                                                        ('fiveprime_acceptor_stem_sequence_stop', 207),
                                                        ('position_8_start', 301),
                                                        ('position_8_stop', 301),
                                                        ('position_9_start', 401),
                                                        ('position_9_stop', 401),
                                                        ('fiveprime_d_stem_sequence_start', 501),
                                                        ('fiveprime_d_stem_sequence_stop', 504),
                                                        ('d_loop_prealpha_start', 601),
                                                        ('d_loop_prealpha_stop', 602),
                                                        ('d_loop_alpha_start', 701),
                                                        ('d_loop_alpha_stop', 703), # variable: default α length range 1-3
                                                        ('d_loop_postalpha_start', 801),
                                                        ('d_loop_postalpha_stop', 802),
                                                        ('d_loop_beta_start', 901),
                                                        ('d_loop_beta_stop', 903), # variable: default β length range 1-3
                                                        ('d_loop_postbeta_start', 1001),
                                                        ('d_loop_postbeta_stop', 1001),
                                                        ('threeprime_d_stem_sequence_start', 1101),
                                                        ('threeprime_d_stem_sequence_stop', 1104),
                                                        ('position_26_start', 1201),
                                                        ('position_26_stop', 1201),
                                                        ('fiveprime_anticodon_stem_sequence_start', 1301),
                                                        ('fiveprime_anticodon_stem_sequence_stop', 1305),
                                                        ('anticodon_loop_start', 1401),
                                                        ('anticodon_loop_stop', 1407),
                                                        ('threeprime_anticodon_stem_sequence_start', 1501),
                                                        ('threeprime_anticodon_stem_sequence_stop', 1505),
                                                        ('v_loop_start', 1601),
                                                        ('v_loop_stop', 1623), # variable: default V loop length range 4-5; 9-23
                                                        ('fiveprime_t_stem_sequence_start', 1701),
                                                        ('fiveprime_t_stem_sequence_stop', 1705),
                                                        ('t_loop_start', 1801),
                                                        ('t_loop_stop', 1807),
                                                        ('threeprime_t_stem_sequence_start', 1901),
                                                        ('threeprime_t_stem_sequence_stop', 1905),
                                                        ('threeprime_acceptor_stem_sequence_start', 2001),
                                                        ('threeprime_acceptor_stem_sequence_stop', 2007),
                                                        ('discriminator_start', 2101),
                                                        ('discriminator_stop', 2101)])
        self.ordinal_extrema_dict = ordinal_extrema_dict = OrderedDict()
        self.ordinal_dict = ordinal_dict = OrderedDict()
        rank = 1
        for feature_item in zip(list(self.pretty_ordinal_extrema_dict.items())[::2],
                                list(self.pretty_ordinal_extrema_dict.items())[1::2]):
            ordinal_extrema_dict[feature_item[0][0]] = rank
            ordinal_extrema_dict[feature_item[1][0]] = rank + feature_item[1][1] - feature_item[0][1]
            feature_name = feature_item[0][0].split('_start')[0]
            for feature_index in range(1, feature_item[1][1] - feature_item[0][1] + 2):
                ordinal_dict[feature_name + '_' + str(feature_index)] = rank
                rank += 1
        self.reverse_ordinal_dict = OrderedDict([(rank, ordinal_name) for ordinal_name, rank in ordinal_dict.items()])

        d_loop_alpha_entries = [('d_loop_alpha_' + str(feature_index), str(feature_index + 15)) for feature_index
                                in range(1, ordinal_extrema_dict['d_loop_alpha_stop'] - ordinal_extrema_dict['d_loop_alpha_start'] + 2)]
        if len(d_loop_alpha_entries) > 2:
            d_loop_alpha_entries[2] = ('d_loop_alpha_3', '17a')
        if len(d_loop_alpha_entries) > 3:
            for entry_index, entry in d_loop_alpha_entries[3: ]:
                d_loop_alpha_entries[entry_index] = (entry[0], '')
        d_loop_beta_entries = [('d_loop_beta_' + str(feature_index), str(feature_index + 19)) for feature_index
                                in range(1, ordinal_extrema_dict['d_loop_beta_stop'] - ordinal_extrema_dict['d_loop_beta_start'] + 2)]
        if len(d_loop_beta_entries) > 1:
            d_loop_beta_entries[1] = ('d_loop_beta_2', '20a')
        if len(d_loop_beta_entries) > 2:
            d_loop_beta_entries[2] = ('d_loop_beta_3', '20b')
        if len(d_loop_beta_entries) > 3:
            for entry_index, entry in d_loop_beta_entries[3: ]:
                d_loop_beta_entries[entry_index] = (entry[0], '')
        self.canonical_dict = OrderedDict([('trna_his_position_0_1', '0')] +
                                          [('fiveprime_acceptor_stem_sequence_' + str(feature_index), str(feature_index)) for feature_index
                                           in range(1, ordinal_extrema_dict['fiveprime_acceptor_stem_sequence_stop'] - ordinal_extrema_dict['fiveprime_acceptor_stem_sequence_start'] + 2)] +
                                          [('position_8_1', '8'),
                                           ('position_9_1', '9')] +
                                          [('fiveprime_d_stem_sequence_' + str(feature_index), str(feature_index + 9)) for feature_index
                                           in range(1, ordinal_extrema_dict['fiveprime_d_stem_sequence_stop'] - ordinal_extrema_dict['fiveprime_d_stem_sequence_start'] + 2)] +
                                          [('d_loop_prealpha_' + str(feature_index), str(feature_index + 13)) for feature_index
                                           in range(1, ordinal_extrema_dict['d_loop_prealpha_stop'] - ordinal_extrema_dict['d_loop_prealpha_start'] + 2)] +
                                          d_loop_alpha_entries +
                                          [('d_loop_postalpha_' + str(feature_index), str(feature_index + 17)) for feature_index
                                           in range(1, ordinal_extrema_dict['d_loop_postalpha_stop'] - ordinal_extrema_dict['d_loop_postalpha_start'] + 2)] +
                                          d_loop_beta_entries +
                                          [('d_loop_postbeta_' + str(feature_index), str(feature_index + 20)) for feature_index
                                           in range(1, ordinal_extrema_dict['d_loop_postbeta_stop'] - ordinal_extrema_dict['d_loop_postbeta_start'] + 2)] +
                                          [('threeprime_d_stem_sequence_' + str(feature_index), str(feature_index + 21)) for feature_index
                                           in range(1, ordinal_extrema_dict['threeprime_d_stem_sequence_stop'] - ordinal_extrema_dict['threeprime_d_stem_sequence_start'] + 2)] +
                                          [('position_26_1', '26')] +
                                          [('fiveprime_anticodon_stem_sequence_' + str(feature_index), str(feature_index + 26)) for feature_index
                                           in range(1, ordinal_extrema_dict['fiveprime_anticodon_stem_sequence_stop'] - ordinal_extrema_dict['fiveprime_anticodon_stem_sequence_start'] + 2)] +
                                          [('anticodon_loop_' + str(feature_index), str(feature_index + 31)) for feature_index
                                           in range(1, ordinal_extrema_dict['anticodon_loop_stop'] - ordinal_extrema_dict['anticodon_loop_start'] + 2)] +
                                          [('threeprime_anticodon_stem_sequence_' + str(feature_index), str(feature_index + 38)) for feature_index
                                           in range(1, ordinal_extrema_dict['threeprime_anticodon_stem_sequence_stop'] - ordinal_extrema_dict['threeprime_anticodon_stem_sequence_start'] + 2)] +
                                          [('v_loop_' + str(feature_index), '') for feature_index in range(1, ordinal_extrema_dict['v_loop_stop'] - ordinal_extrema_dict['v_loop_start'] + 2)] +
                                          [('fiveprime_t_stem_sequence_' + str(feature_index), str(feature_index + 48)) for feature_index
                                           in range(1, ordinal_extrema_dict['fiveprime_t_stem_sequence_stop'] - ordinal_extrema_dict['fiveprime_t_stem_sequence_start'] + 2)] +
                                          [('t_loop_' + str(feature_index), str(feature_index + 53)) for feature_index
                                           in range(1, ordinal_extrema_dict['t_loop_stop'] - ordinal_extrema_dict['t_loop_start'] + 2)] +
                                          [('threeprime_t_stem_sequence_' + str(feature_index), str(feature_index + 60)) for feature_index
                                           in range(1, ordinal_extrema_dict['threeprime_t_stem_sequence_stop'] - ordinal_extrema_dict['threeprime_t_stem_sequence_start'] + 2)] +
                                          [('threeprime_acceptor_stem_sequence_' + str(feature_index), str(feature_index + 65)) for feature_index
                                           in range(1, ordinal_extrema_dict['threeprime_acceptor_stem_sequence_stop'] - ordinal_extrema_dict['threeprime_acceptor_stem_sequence_start'] + 2)] +
                                          [('discriminator_1', '73')])


    def generate_seed_output(self):
        contig_name_iter = iter(self.contig_names)
        feature_table_dict = self.contigs_db.db.get_table_as_dict(tables.trna_seed_feature_table_name)
        self.taxonomy_table_dict = taxonomy_table_dict = {}
        for taxonomy_item in self.contigs_db.db.get_table_as_dataframe(tables.trna_taxonomy_table_name,
                                                                       error_if_no_data=False,
                                                                       columns_of_interest=['gene_callers_id',
                                                                                            't_domain',
                                                                                            't_phylum',
                                                                                            't_class',
                                                                                            't_order',
                                                                                            't_family',
                                                                                            't_genus',
                                                                                            't_species',
                                                                                            'percent_identity']).itertuples():
            taxonomy_table_dict[taxonomy_item.gene_callers_id] = (taxonomy_item.t_domain if taxonomy_item.t_domain else '',
                                                                  taxonomy_item.t_phylum if taxonomy_item.t_phylum else '',
                                                                  taxonomy_item.t_class if taxonomy_item.t_class else '',
                                                                  taxonomy_item.t_order if taxonomy_item.t_order else '',
                                                                  taxonomy_item.t_family if taxonomy_item.t_family else '',
                                                                  taxonomy_item.t_genus if taxonomy_item.t_genus else '',
                                                                  taxonomy_item.t_species if taxonomy_item.t_species else '',
                                                                  str(taxonomy_item.percent_identity) if taxonomy_item.percent_identity else '')

        spec_profile_db = dbops.ProfileDatabase(self.spec_profile_db_path, quiet=True)
        get_meta_value = spec_profile_db.db.get_meta_value
        self.sample_names = sample_names = get_meta_value('samples').split(', ')
        sample_total_mean_spec_covs = tuple(map(float, get_meta_value('sample_total_mean_specific_coverage').split(', ')))
        sample_total_discriminator_spec_covs = tuple(map(int, get_meta_value('sample_total_discriminator_specific_coverage').split(', ')))
        mean_spec_cov_df = spec_profile_db.db.get_table_as_dataframe('mean_coverage_contigs')
        spec_profile_db.disconnect()
        mean_spec_cov_df['contig_name'] = mean_spec_cov_df['item'].apply(lambda s: s.split('_split_00001')[0])
        mean_spec_cov_df = mean_spec_cov_df.drop(['item', 'layer'], axis=1)
        mean_spec_cov_df = mean_spec_cov_df.rename({'value': 'mean_spec_cov'}, axis=1)
        mean_spec_cov_dict = {}
        for contig_name, contig_df in mean_spec_cov_df.groupby('contig_name'):
            mean_spec_cov_dict[contig_name] = tuple(contig_df['mean_spec_cov'])

        do_nonspec = True if self.nonspec_profile_db_path else False
        if do_nonspec:
            nonspec_profile_db = dbops.ProfileDatabase(self.nonspec_profile_db_path, quiet=True)
            get_meta_value = nonspec_profile_db.db.get_meta_value
            self.sample_names = sample_names = get_meta_value('samples').split(', ')
            mean_nonspec_cov_df = nonspec_profile_db.db.get_table_as_dataframe('mean_coverage_contigs')
            nonspec_profile_db.disconnect()
            mean_nonspec_cov_df['contig_name'] = mean_nonspec_cov_df['item'].apply(lambda s: s.split('_split_00001')[0])
            mean_nonspec_cov_df = mean_nonspec_cov_df.drop(['item', 'layer'], axis=1)
            mean_nonspec_cov_df = mean_nonspec_cov_df.rename({'value': 'mean_nonspec_cov'}, axis=1)
            mean_nonspec_cov_dict = {}
            for contig_name, contig_df in mean_nonspec_cov_df.groupby('contig_name'):
                mean_nonspec_cov_dict[contig_name] = tuple(contig_df['mean_nonspec_cov'])

        anticodon_aa_items = [(anticodon, aa) for aa, anticodon in
                              [anticodon_aa_item.split('_') for anticodon_aa_item in
                               self.contigs_db.db.get_single_column_from_table(tables.hmm_hits_table_name, 'gene_name')]]
        self.gene_callers_id_anticodon_aa_dict = gene_callers_id_anticodon_aa_dict = {}
        for gene_callers_id, anticodon_aa_item in enumerate(anticodon_aa_items):
            gene_callers_id_anticodon_aa_dict[gene_callers_id] = anticodon_aa_item
        anticodon_aa_iter = iter(anticodon_aa_items)

        # Coverages are stored as bytes or as a single integer when coverage is zero for the sample.
        # Convert coverage entries back to numpy arrays accordingly.
        convert_binary_blob_to_numpy_array = partial(convert_binary_blob_to_numpy_array, dtype=auxiliarydataops.TRNASEQ_COVERAGE_DTYPE)
        convert_int_to_numpy_array = lambda x: np.zeros(x, dtype=auxiliarydataops.TRNASEQ_COVERAGE_DTYPE)
        def convert_coverage_entry_to_numpy_array(cov_entry):
            try:
                cov_array = convert_binary_blob_to_numpy_array(cov_entry)
            except TypeError:
                cov_array = convert_int_to_numpy_array(cov_entry)
            return cov_array

        spec_aux_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.spec_aux_db_path, self.contigs_db_info.hash, db_variant='trnaseq')
        spec_aux_df = spec_aux_db.db.get_table_as_dataframe(tables.split_coverages_table_name)
        spec_aux_db.close()
        spec_aux_df['contig_name'] = spec_aux_df['split_name'].apply(lambda s: s.split('_split_00001')[0])
        spec_aux_df = spec_aux_df.drop('split_name', axis=1)
        spec_aux_df['coverages'] = spec_aux_df['coverages'].apply(convert_coverage_entry_to_numpy_array)

        if do_nonspec:
            nonspec_aux_db = auxiliarydataops.AuxiliaryDataForSplitCoverages(self.nonspec_aux_db_path, self.contigs_db_info.hash, db_variant='trnaseq')
            nonspec_aux_df = nonspec_aux_db.db.get_table_as_dataframe(tables.split_coverages_table_name)
            nonspec_aux_db.close()
            nonspec_aux_df['contig_name'] = spec_aux_df['contig_name'].values
            nonspec_aux_df = nonspec_aux_df.drop('split_name', axis=1)
            nonspec_aux_df['coverages'] = nonspec_aux_df['coverages'].apply(convert_coverage_entry_to_numpy_array)

        spec_covs_dict = {}
        sample_rel_discriminator_spec_cov_dict = {}
        for contig_name, contig_df in spec_aux_df.groupby('contig_name'):
            spec_covs_dict[contig_name] = dict(zip(contig_df['sample_name'], contig_df['coverages']))
            sample_rel_discriminator_spec_cov_dict[contig_name] = tuple([round(covs[-1] / total_discriminator_cov, 8) for covs, total_discriminator_cov
                                                                         in zip(contig_df['coverages'], sample_total_discriminator_spec_covs)])

        if do_nonspec:
            nonspec_covs_dict = {}
            for contig_name, contig_df in nonspec_aux_df.groupby('contig_name'):
                nonspec_covs_dict[contig_name] = dict(zip(contig_df['sample_name'], contig_df['coverages']))

        canonical_dict = self.canonical_dict
        spec_top_header = ("gene_callers_id",
                           "contig_name",
                           "anticodon",
                           "aa",
                           "domain",
                           "phylum",
                           "class",
                           "order",
                           "family",
                           "genus",
                           "species",
                           "taxon_percent_id",
                           "sample_name",
                           "mean_coverage",
                           "relative_mean_coverage",
                           "relative_discriminator_coverage") + tuple(self.ordinal_dict.keys())
        spec_middle_header = "\t".join(("", ) * (len(spec_top_header) - len(self.ordinal_dict))
                                       + tuple(map(str, range(1, len(self.ordinal_dict) + 1)))) + "\n"
        spec_bottom_header = "\t".join(("", ) * (len(spec_top_header) - len(self.canonical_dict))
                                       + tuple([canonical_dict[ordinal_name] for ordinal_name in self.ordinal_dict])) + "\n"
        spec_top_header = "\t".join(spec_top_header) + "\n"

        if do_nonspec:
            nonspec_top_header = ("gene_callers_id",
                                  "contig_name",
                                  "anticodon",
                                  "aa",
                                  "domain",
                                  "phylum",
                                  "class",
                                  "order",
                                  "family",
                                  "genus",
                                  "species",
                                  "taxon_percent_id",
                                  "sample_name",
                                  "mean_coverage") + tuple(self.ordinal_dict.keys())
            nonspec_middle_header = "\t".join(("", ) * (len(spec_top_header) - len(self.ordinal_dict))
                                              + tuple(map(str, range(1, len(self.ordinal_dict) + 1)))) + "\n"
            nonspec_bottom_header = "\t".join(("", ) * (len(spec_top_header) - len(self.canonical_dict))
                                              + tuple([canonical_dict[ordinal_name] for ordinal_name in self.ordinal_dict])) + "\n"
            nonspec_top_header = "\t".join(nonspec_top_header) + "\n"


        spec_out_file = open(self.spec_seed_out_path, 'a')
        spec_out_file.write(spec_top_header + spec_middle_header + spec_bottom_header)
        if do_nonspec:
            nonspec_out_file = open(self.nonspec_seed_out_path, 'a')
            nonspec_out_file.write(nonspec_top_header + nonspec_middle_header + nonspec_bottom_header)

        reverse_ordinal_dict = self.reverse_ordinal_dict
        self.contig_name_gene_callers_id_dict = contig_name_gene_callers_id_dict = {}
        self.contig_feature_coord_dict = contig_feature_coord_dict = {}
        for gene_callers_id, feature_dict in feature_table_dict.items():
            contig_name = next(contig_name_iter)
            anticodon, aa = next(anticodon_aa_iter)
            contig_name_gene_callers_id_dict[contig_name] = gene_callers_id

            try:
                t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species, t_percent_id = taxonomy_table_dict[gene_callers_id]
            except KeyError:
                t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species, t_percent_id = ('', ) * 8

            sample_mean_spec_covs = mean_spec_cov_dict[contig_name]
            sample_rel_mean_spec_covs = [round(mean_spec_cov / total_mean_spec_cov, 8) for mean_spec_cov, total_mean_spec_cov
                                         in zip(sample_mean_spec_covs, sample_total_mean_spec_covs)]
            sample_rel_discriminator_spec_covs = sample_rel_discriminator_spec_cov_dict[contig_name]

            if do_nonspec:
                sample_mean_nonspec_covs = mean_nonspec_cov_dict[contig_name]

            contig_spec_covs_dict = {sample_name: iter(covs) for sample_name, covs in spec_covs_dict[contig_name].items()}
            if do_nonspec:
                contig_nonspec_covs_dict = {sample_name: iter(covs) for sample_name, covs in nonspec_covs_dict[contig_name].items()}
            sample_range = range(len(contig_spec_covs_dict))
            spec_contig_output_rows = []
            for sample_name, mean_spec_cov, rel_mean_cov, rel_discriminator_cov in zip(sample_names,
                                                                                       sample_mean_spec_covs,
                                                                                       sample_rel_mean_spec_covs,
                                                                                       sample_rel_discriminator_spec_covs):
                spec_contig_output_rows.append([str(gene_callers_id),
                                                contig_name,
                                                anticodon,
                                                aa,
                                                t_domain,
                                                t_phylum,
                                                t_class,
                                                t_order,
                                                t_family,
                                                t_genus,
                                                t_species,
                                                t_percent_id,
                                                sample_name,
                                                str(round(mean_spec_cov, 1)),
                                                str(rel_mean_cov),
                                                str(rel_discriminator_cov)])

            if do_nonspec:
                nonspec_contig_output_rows = []
                for spec_contig_output_row, mean_nonspec_cov in zip(spec_contig_output_rows, sample_mean_nonspec_covs):
                    nonspec_contig_output_rows.append(spec_contig_output_row[: -3] + [str(round(mean_nonspec_cov, 1))])

            feature_start_in_seq = None
            feature_stop_in_seq = None
            feature_ordinal_start = None

            tabulated_contig_spec_covs = [[] for sample_index in sample_range]
            if do_nonspec:
                tabulated_contig_nonspec_covs = [[] for sample_index in sample_range]
            feature_coord_map = {}
            ordinal_extremum_index = 0
            for ordinal_extremum_name, ordinal_extremum in self.ordinal_extrema_dict.items():
                try:
                    extremum_in_seq = feature_dict[ordinal_extremum_name]
                except KeyError:
                    if ordinal_extremum_name == 'd_loop_prealpha_start':
                        extremum_in_seq = feature_dict['d_loop_start']
                    elif ordinal_extremum_name == 'd_loop_prealpha_stop':
                        alpha_start_in_seq = feature_dict['alpha_start']
                        extremum_in_seq = alpha_start_in_seq - 1 if pd.notnull(alpha_start_in_seq) else np.nan
                    elif ordinal_extremum_name == 'd_loop_alpha_start':
                        extremum_in_seq = alpha_start_in_seq = feature_dict['alpha_start']
                    elif ordinal_extremum_name == 'd_loop_alpha_stop':
                        extremum_in_seq = alpha_stop_in_seq = feature_dict['alpha_stop']
                    elif ordinal_extremum_name == 'd_loop_postalpha_start':
                        extremum_in_seq = alpha_stop_in_seq + 1 if pd.notnull(alpha_stop_in_seq) else np.nan
                    elif ordinal_extremum_name == 'd_loop_postalpha_stop':
                        beta_start_in_seq = feature_dict['beta_start']
                        extremum_in_seq = beta_start_in_seq - 1 if pd.notnull(beta_start_in_seq) else np.nan
                    elif ordinal_extremum_name == 'd_loop_beta_start':
                        extremum_in_seq = beta_start_in_seq
                    elif ordinal_extremum_name == 'd_loop_beta_stop':
                        extremum_in_seq = beta_stop_in_seq = feature_dict['beta_stop']
                    elif ordinal_extremum_name == 'd_loop_postbeta_start':
                        extremum_in_seq = beta_stop_in_seq + 1 if pd.notnull(beta_stop_in_seq) else np.nan
                    elif ordinal_extremum_name == 'd_loop_postbeta_stop':
                        extremum_in_seq = feature_dict['d_loop_stop']

                if pd.notnull(feature_start_in_seq):
                    feature_stop_in_seq = extremum_in_seq
                    feature_ordinal_max_stop = ordinal_extremum
                    for ordinal_index, seq_index in enumerate(range(feature_start_in_seq, feature_stop_in_seq + 1)):
                        if seq_index >= 0:
                            for sample_num, spec_cov_iter in enumerate(contig_spec_covs_dict.values()):
                                tabulated_contig_spec_covs[sample_num].append(next(spec_cov_iter))
                            if do_nonspec:
                                for sample_num, nonspec_cov_iter in enumerate(contig_nonspec_covs_dict.values()):
                                    tabulated_contig_nonspec_covs[sample_num].append(next(nonspec_cov_iter))
                            rank = feature_ordinal_start + ordinal_index
                            feature_coord_map[seq_index] = (reverse_ordinal_dict[rank], rank)
                        else:
                            for sample_num in sample_range:
                                tabulated_contig_spec_covs[sample_num].append('')
                                if do_nonspec:
                                    tabulated_contig_nonspec_covs[sample_num].append('')
                    uncovered_ordinal_length = feature_ordinal_max_stop - feature_ordinal_start - ordinal_index
                    if uncovered_ordinal_length:
                        cov_extension = [''] * uncovered_ordinal_length
                        for sample_num in sample_range:
                            tabulated_contig_spec_covs[sample_num].extend(cov_extension)
                            if do_nonspec:
                                tabulated_contig_nonspec_covs[sample_num].extend(cov_extension)
                    feature_start_in_seq = None
                    feature_stop_in_seq = None
                    feature_ordinal_start = None
                    feature_ordinal_max_stop = None
                elif ordinal_extremum_index % 2 == 1:
                    cov_extension = [''] * (ordinal_extremum + 1 - feature_ordinal_start)
                    for sample_num in sample_range:
                        tabulated_contig_spec_covs[sample_num].extend(cov_extension)
                        if do_nonspec:
                            tabulated_contig_nonspec_covs[sample_num].extend(cov_extension)
                    feature_start_in_seq = None
                    feature_ordinal_start = None
                else:
                    feature_start_in_seq = extremum_in_seq
                    feature_ordinal_start = ordinal_extremum
                ordinal_extremum_index += 1

            contig_feature_coord_dict[contig_name] = feature_coord_map

            for spec_contig_output_row, contig_spec_covs_row in zip(spec_contig_output_rows, tabulated_contig_spec_covs):
                spec_out_file.write("\t".join(spec_contig_output_row + list(map(str, contig_spec_covs_row))) + "\n")
            if do_nonspec:
                for nonspec_contig_output_row, contig_nonspec_covs_row in zip(nonspec_contig_output_rows, tabulated_contig_nonspec_covs):
                    nonspec_out_file.write("\t".join(nonspec_contig_output_row + list(map(str, contig_nonspec_covs_row))) + "\n")


    def generate_modification_output(self):
        profile_db = dbops.ProfileDatabase(self.spec_profile_db_path)
        var_nts_df = profile_db.db.get_table_as_dataframe('variable_nucleotides')
        profile_db.disconnect()
        var_nts_df['contig_name'] = var_nts_df['split_name'].apply(lambda s: s.split('_split_00001')[0])
        contig_name_gene_callers_id_dict = self.contig_name_gene_callers_id_dict
        var_nts_df['gene_callers_id'] = var_nts_df['contig_name'].apply(lambda c: contig_name_gene_callers_id_dict[c])
        var_nts_df = var_nts_df.drop(['pos_in_contig',
                                      'corresponding_gene_call',
                                      'in_noncoding_gene_call',
                                      'in_coding_gene_call',
                                      'base_pos_in_codon',
                                      'codon_order_in_gene',
                                      'cov_outlier_in_split',
                                      'departure_from_reference',
                                      'competing_nts'],
                                     axis=1)
        var_nts_df = var_nts_df.sort_values('gene_callers_id')

        header = "\t".join(("gene_callers_id",
                            "contig_name",
                            "anticodon",
                            "aa",
                            "domain",
                            "phylum",
                            "class",
                            "order",
                            "family",
                            "genus",
                            "species",
                            "taxon_percent_id",
                            "seed_position",
                            "ordinal_name",
                            "ordinal_position",
                            "canonical_position",
                            "reference",
                            "sample_name")
                           + UNAMBIG_NTS) + "\n"
        out_file = open(self.mod_out_path, 'a')
        out_file.write(header)

        gene_callers_id_anticodon_aa_dict = self.gene_callers_id_anticodon_aa_dict
        taxonomy_table_dict = self.taxonomy_table_dict
        contig_feature_coord_dict = self.contig_feature_coord_dict
        canonical_dict = self.canonical_dict
        sample_names = self.sample_names
        for contig_index, contig_df in var_nts_df.groupby(['gene_callers_id', 'contig_name']):
            gene_callers_id = contig_index[0]
            contig_name = contig_index[1]
            anticodon, aa = gene_callers_id_anticodon_aa_dict[gene_callers_id]
            try:
                t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species, t_percent_id = taxonomy_table_dict[gene_callers_id]
            except KeyError:
                t_domain, t_phylum, t_class, t_order, t_family, t_genus, t_species, t_percent_id = ('', ) * 8
            gene_callers_id = str(gene_callers_id)
            contig_df = contig_df.sort_values('pos')
            feature_coord_dict = contig_feature_coord_dict[contig_name]
            for seed_pos, sub_df in contig_df.groupby('pos'):
                try:
                    ordinal_name, ordinal_pos = feature_coord_dict[seed_pos]
                except KeyError:
                    # The modification position may be 5' of the 5'-most feature defined in the
                    # seed.
                    ordinal_name = 'NA'
                    ordinal_pos = 'NA'
                str_seed_pos = str(seed_pos)
                ordinal_pos = str(ordinal_pos)
                try:
                    canonical_pos = str(canonical_dict[ordinal_name])
                except KeyError:
                    canonical_pos = 'NA'
                sample_row_dict = {}
                for sample_row in sub_df.itertuples():
                    out_row = [gene_callers_id,
                               contig_name,
                               anticodon,
                               aa,
                               t_domain,
                               t_phylum,
                               t_class,
                               t_order,
                               t_family,
                               t_genus,
                               t_species,
                               t_percent_id,
                               str_seed_pos,
                               ordinal_name,
                               ordinal_pos,
                               canonical_pos,
                               sample_row.reference,
                               sample_row.sample_id]
                    sample_row_dict[out_row[-1]] = out_row
                    out_row.extend([str(getattr(sample_row, nt)) for nt in UNAMBIG_NTS])
                for sample_name in sample_names:
                    try:
                        out_row = sample_row_dict[sample_name]
                    except KeyError:
                        out_row = [gene_callers_id,
                                   contig_name,
                                   anticodon,
                                   aa,
                                   t_domain,
                                   t_phylum,
                                   t_class,
                                   t_order,
                                   t_family,
                                   t_genus,
                                   t_species,
                                   t_percent_id,
                                   str_seed_pos,
                                   ordinal_name,
                                   ordinal_pos,
                                   canonical_pos,
                                   sample_row.reference,
                                   sample_name] + [""] * len(UNAMBIG_NTS)
                    out_file.write("\t".join(out_row) + "\n")


class ResultPlotter(object):

    NT_COLORS = ('blue', 'green', 'orange', 'red')
    RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    DEFAULT_FORMAT_PARAM_DICT = {
        'nt_key_y': 1.07,
        'nt_key_font_size': 3,
        'group_id_y': 1.07,
        'group_id_font_size': 6,
        'seed_count_y': 1.07,
        'seed_count_font_size': 6,
        'cov_tick_font_size': 3,
        'mod_frac_tick_font_size': 2,
        'sample_rel_abund_font_size': 3,
        'ordinal_pos_tick_font_size': 4
    }

    def __init__(self, args, run=terminal.Run(width=100)):
        self.run = run

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # Group 1: MANDATORY
        self.contigs_db_path = A('contigs_db')
        self.spec_txt_path = A('seeds_specific_txt')
        self.mod_txt_path = A('modifications_txt')
        self.out_dir = A('output_dir')

        # Group 2: OPTIONAL
        self.overwrite_out_dest = A('overwrite_output_destinations')
        # The following parameters are not currently available to the user.
        self.nonspec_txt_path = None
        self.cov_cutoff = 1

        if not self.contigs_db_path:
            raise ConfigError("Please provide the path to a `trnaseq`-variant contigs database using `--contigs-db` or `-c`.")
        if not self.spec_txt_path:
            raise ConfigError("Please provide the path to a tRNA seeds table of data including specific coverages using `--seeds-specific-txt` or `-s`.")
        if not self.mod_txt_path:
            raise ConfigError("Please provide the path to a tRNA modifications table using `--modifications-txt` or `-m`.")

        self.format_param_dict = self.DEFAULT_FORMAT_PARAM_DICT


    def go(self):
        """Load data and get standard input from user to produce plots."""
        self.sanity_check()
        self.load_data()
        self.print_help()

        RANKS = self.RANKS
        DEFAULT_FORMAT_PARAM_DICT = self.DEFAULT_FORMAT_PARAM_DICT
        warning = self.run.warning
        while True:
            entry = input("\033[36;1;3mInput -> \033[0m")

            if entry == 'q':
                return
            elif entry == 'h':
                self.print_help()
                continue
            elif entry == 'f':
                self.print_formatting_options()
                continue

            loop_ranks = []
            format_command = False
            taxon_rank_filter = None
            taxon_filter = None
            single_aa = None
            single_anticodon = None
            out_dir = self.out_dir
            for field_index, field in enumerate([field.strip() for field in entry.split(';')]):
                split_field = field.split(',')

                if len(split_field) == 1:
                    if split_field[0] == '':
                        warning("An empty semicolon-delimited field was given.", "INVALID FIELD")
                        break

                    single_field = split_field[0]
                    if single_field == 'f':
                        if field_index > 0:
                            warning("A format command must start with f; and not include plot commands.", "INVALID FIELD")
                            break
                        format_command = True
                    elif single_field in RANKS:
                        loop_ranks.append(single_field)
                    elif single_field in AMINO_ACIDS:
                        if single_aa:
                            warning("Multiple amino acids were given.", "INVALID FIELD")
                            break
                        single_aa = single_field
                    elif single_field in ANTICODONS:
                        if single_anticodon:
                            warning("Multiple anticodons were given.", "INVALID FIELD")
                            break
                        single_anticodon = single_field
                    elif format_command:
                        if single_field == 'reset':
                            self.format_param_dict = DEFAULT_FORMAT_PARAM_DICT
                        else:
                            warning("A format option must have a value of 'reset' or be similar to 'nt_key_y, 1.1'", "INVALID FIELD")
                    else:
                        warning(f"{single_field} is not recognized as a taxonomic rank, amino acid, or anticodon.", "INVALID FIELD")
                        break
                elif len(split_field) == 2:
                    if format_command:
                        format_param = split_field[0].strip()
                        if format_param not in DEFAULT_FORMAT_PARAM_DICT:
                            warning(f"{format_param} is not recognized. Type f for format options.", "INVALID FIELD")
                        self.format_param_dict[format_param] = float(split_field[1].strip())
                        continue

                    if split_field[0].strip() == 'outdir':
                        out_dir = split_field[1].strip()
                        if os.path.isdir(out_dir):
                            if not os.access(os.path.abspath(out_dir), os.W_OK):
                                warning(f"Permission denied to create files in {out_dir}", "INVALID DIRECTORY PATH")
                                break
                        else:
                            try:
                                os.mkdir(out_dir)
                            except FileNotFoundError:
                                warning(f"{out_dir}", "INVALID DIRECTORY PATH")
                                break
                        continue

                    if taxon_filter:
                        warning("Multiple fields with a comma, interpreted as a taxon, were given.", "INVALID FIELD")
                        break
                    taxon_rank_filter = split_field[0].strip()

                    if taxon_rank_filter not in RANKS:
                        warning(f"{taxon_rank_filter} is not recognized as a taxonomic rank.", "INVALID FIELD")
                        break
                    taxon_filter = split_field[1].strip()
                else:
                    warning(f"A field with {len(split_field)} commas was given: '{field}'.", "INVALID FIELD")
                    break
            else:
                if format_command:
                    continue

                if not taxon_rank_filter and not loop_ranks:
                    warning(f"A taxon or at least one rank must be given.", "PROGRAM REQUIREMENT")
                    continue

                loop_ranks = sorted(loop_ranks, key=lambda rank: RANKS.index(rank))
                if taxon_rank_filter and loop_ranks:
                    if RANKS.index(taxon_rank_filter) >= RANKS.index(loop_ranks[0]):
                        warning(f"A taxon ('{taxon_rank_filter} {taxon_filter}') lower than a rank ('{loop_ranks[0]}') was given.", "INVALID FIELD")
                        continue

                if single_aa and single_anticodon:
                    if ANTICODON_AA_DICT[single_anticodon] == single_aa:
                        warning(f"An amino acid is not needed with the anticodon.", "UNNECESSARY FIELD")
                        single_aa = None
                    else:
                        warning(f"The anticodon ('{single_anticodon}') does not decode the amino acid ('{single_aa}').", "INVALID FIELD")
                        continue

                spec_df = self.spec_df
                mod_df = self.mod_df
                nonspec_df = self.nonspec_df if self.nonspec_txt_path else None

                if single_anticodon:
                    spec_df = spec_df[(spec_df['anticodon'] == single_anticodon).squeeze()]
                    mod_df = mod_df[(mod_df['anticodon'] == single_anticodon)]
                    nonspec_df = nonspec_df[(nonspec_df['anticodon'] == single_anticodon).squeeze()] if self.nonspec_txt_path else None

                if taxon_rank_filter:
                    # Subset the data to the taxon of interest.
                    spec_df = spec_df[(spec_df[taxon_rank_filter] == taxon_filter).squeeze()]
                    mod_df = mod_df[(mod_df[taxon_rank_filter] == taxon_filter)]
                    nonspec_df = nonspec_df[(nonspec_df[taxon_rank_filter] == taxon_filter).squeeze()] if self.nonspec_txt_path else None

                if single_anticodon:
                    if taxon_filter:
                        if loop_ranks: # (potentially) multiple taxa at ranks below the taxon filter
                            for rank in loop_ranks:
                                self.plot_rank(rank,
                                               spec_df,
                                               mod_df,
                                               ANTICODON_AA_DICT[single_anticodon],
                                               single_anticodon,
                                               out_dir,
                                               nonspec_df=nonspec_df)
                        else: # single taxon
                            self.plot_items(spec_df,
                                            mod_df,
                                            ANTICODON_AA_DICT[single_anticodon],
                                            single_anticodon,
                                            taxon_rank_filter,
                                            taxon_filter,
                                            out_dir,
                                            nonspec_df=nonspec_df)
                    else: # (potentially) multiple taxa
                        for rank in loop_ranks:
                            self.plot_rank(rank,
                                           spec_df,
                                           mod_df,
                                           ANTICODON_AA_DICT[single_anticodon],
                                           single_anticodon,
                                           out_dir,
                                           nonspec_df=nonspec_df)
                else: # (potentially) multiple anticodons
                    if single_aa:
                        for anticodon in AA_ANTICODON_DICT[single_aa]:
                            anticodon_spec_df = spec_df[(spec_df['anticodon'] == anticodon).squeeze()]
                            anticodon_mod_df = mod_df[(mod_df['anticodon'] == anticodon)]
                            anticodon_nonspec_df = nonspec_df[(nonspec_df['anticodon'] == anticodon).squeeze()] if self.nonspec_txt_path else None
                            if taxon_filter:
                                if loop_ranks: # (potentially) multiple taxa at ranks below the taxon filter
                                    for rank in loop_ranks:
                                        self.plot_rank(rank,
                                                       anticodon_spec_df,
                                                       anticodon_mod_df,
                                                       single_aa,
                                                       anticodon,
                                                       out_dir,
                                                       nonspec_df=anticodon_nonspec_df)
                                else: # single taxon
                                    self.plot_items(anticodon_spec_df,
                                                    anticodon_mod_df,
                                                    single_aa,
                                                    anticodon,
                                                    taxon_rank_filter,
                                                    taxon_filter,
                                                    out_dir,
                                                    nonspec_df=anticodon_nonspec_df)
                            else: # (potentially) multiple taxa
                                for rank in loop_ranks:
                                    self.plot_rank(rank,
                                                   anticodon_spec_df,
                                                   anticodon_mod_df,
                                                   single_aa,
                                                   anticodon,
                                                   out_dir,
                                                   nonspec_df=anticodon_nonspec_df)
                    else: # all anticodons
                        for aa, anticodons in AA_ANTICODON_DICT.items():
                            for anticodon in anticodons:
                                anticodon_spec_df = spec_df[(spec_df['anticodon'] == anticodon).squeeze()]
                                anticodon_mod_df = mod_df[(mod_df['anticodon'] == anticodon)]
                                anticodon_nonspec_df = nonspec_df[(nonspec_df['anticodon'] == anticodon).squeeze()] if self.nonspec_txt_path else None
                                if taxon_filter:
                                    if loop_ranks: # (potentially) multiple taxa at ranks below the taxon filter
                                        for rank in loop_ranks:
                                            self.plot_rank(rank,
                                                           anticodon_spec_df,
                                                           anticodon_mod_df,
                                                           aa,
                                                           anticodon,
                                                           out_dir,
                                                           nonspec_df=anticodon_nonspec_df)
                                    else: # single taxon
                                        self.plot_items(anticodon_spec_df,
                                                        anticodon_mod_df,
                                                        aa,
                                                        anticodon,
                                                        taxon_rank_filter,
                                                        taxon_filter,
                                                        out_dir,
                                                        nonspec_df=anticodon_nonspec_df)
                                else: # (potentially) multiple taxa
                                    for rank in loop_ranks:
                                        self.plot_rank(rank,
                                                       anticodon_spec_df,
                                                       anticodon_mod_df,
                                                       aa,
                                                       anticodon,
                                                       out_dir,
                                                       nonspec_df=anticodon_nonspec_df)
                print("\n")


    def sanity_check(self):
        """Check `anvi-plot-trnaseq` arguments."""
        self.contigs_db_info = DBInfo(self.contigs_db_path, expecting='contigs')
        if self.contigs_db_info.variant != 'trnaseq':
            raise ConfigError("The input contigs database is not of the `trnaseq` variant...")

        filesnpaths.is_file_tab_delimited(self.spec_txt_path)
        if self.nonspec_txt_path:
            filesnpaths.is_file_tab_delimited(self.nonspec_txt_path)
        filesnpaths.is_file_tab_delimited(self.mod_txt_path)

        if self.out_dir:
            filesnpaths.is_output_dir_writable(self.out_dir)


    def load_data(self):
        # Load the tables of seed specific coverages and modifications.
        dtype_dict = {rank: str for rank in self.RANKS}
        self.spec_df = spec_df = pd.read_csv(self.spec_txt_path, sep='\t', header=[0, 1, 2], dtype=dtype_dict)

        # The column multiindex has 3 levels. Some of the labels in the levels should be empty
        # strings, but upon loading the table, these are replaced with default values that need to
        # be changed to empty strings. The column multiindex for the positional coverage columns
        # also needs to be isolated.
        new_header = []
        cov_header = []
        in_cov_cols = False
        self.cov_x_labels = cov_x_labels = []
        for col_names in spec_df.columns:
            new_col_names = []
            for col_name in col_names:
                if 'Unnamed: ' in col_name:
                    new_col_names.append('')
                else:
                    new_col_names.append(col_name)
                if col_name == 'trna_his_position_0_1':
                    in_cov_cols = True
            new_header.append(new_col_names)
            if in_cov_cols:
                cov_header.append(tuple(new_col_names))
                cov_x_labels.append(new_col_names[2])
        spec_df.columns = pd.MultiIndex.from_tuples(new_header)
        self.cov_index = pd.MultiIndex.from_tuples(cov_header)

        self.cov_length = len(self.cov_index.levels[0])
        self.cov_x_values = np.arange(self.cov_length)

        self.mod_y_ticks = np.arange(0, 1.25, 0.25)

        # Each seed in a seeds coverage table has rows for all samples in the same order. Get sample
        # names from the first seed in the table.
        self.num_samples = len(set(spec_df['sample_name'].squeeze()))
        self.sample_names = spec_df['sample_name'].squeeze().iloc[: self.num_samples].tolist()

        dtype_dict = {rank: str for rank in self.RANKS}
        dtype_dict['ordinal_position'] = pd.Int16Dtype()
        dtype_dict['canonical_position'] = str
        self.mod_df = pd.read_csv(self.mod_txt_path, sep='\t', header=0, dtype=dtype_dict)


    def print_help(self):
        self.run.warning("", header="INPUT HELP", lc='green')
        info_single = self.run.info_single
        info_single("This program generates plots of seed coverage and modification levels across samples.", cut_after=100)
        info_single("Use a semicolon to separate fields and a comma to separate rank and taxon within a field.", cut_after=100, nl_after=1)

        info_single("Write a single plot by specifying a taxon and anticodon.", mc='cyan')
        info_single("This input will plot Clostridia Arg-TCT seeds: class, Clostridia; TCT", mc='magenta', level=2)
        info_single("Redundant but accepted: class, Clostridia; Arg; TCT", mc='magenta', level=2, nl_after=1)

        info_single("Write plots at a single taxonomic rank.", mc='cyan')
        info_single("Every class + every anticodon: class", mc='magenta', level=2)
        info_single("Every class + Arg-TCT: class; TCT", mc='magenta', level=2, nl_after=1)

        info_single("Write plots for every anticodon decoding an amino acid.", mc='cyan')
        info_single("One class + every Arg anticodon: class, Clostridia; Arg", mc='magenta', level=2)
        info_single("Every class + every Arg anticodon: class; Arg", mc='magenta', level=2, nl_after=1)

        info_single("Write plots for multiple taxonomic ranks.", mc='cyan')
        info_single("Every phylum and class + every anticodon: phylum; class", mc='magenta', level=2)
        info_single("Every order in class Clostridia + every anticodon: class, Clostridia; order", mc='magenta', level=2)
        info_single("Every genus in class Clostridia + Arg-TCT: class, Clostridia; genus; TCT", mc='magenta', level=2)
        info_single("Every order and family in class Clostridia + every Arg anticodon: class, Clostridia; order; family; Arg", cut_after=120, mc='magenta', level=2)
        info_single("Invalid: phylum; class, Clostridia", mc='red', level=2, nl_after=1)

        info_single("The output directory can be changed to a new or existing directory.", mc='cyan', cut_after=100)
        info_single("Example: class, Clostridia; Arg; outdir, plots/class/Clostridia/aa/Arg", mc='magenta', level=2, nl_after=1)

        info_single("Only one taxon, amino acid, or anticodon can be given at a time.", mc='cyan')
        info_single("Invalid: phylum, Firmicutes; class, Clostridia", mc='red', level=2)
        info_single("Invalid: class, Clostridia; class, Bacteroidia", mc='red', level=2)
        info_single("Invalid: Arg; Asp", mc='red', level=2, nl_after=1)

        info_single("Alter plot appearance with format options.", mc='cyan')
        info_single("The available parameters, such as tick label size, may need to be tuned to make everything fit.", cut_after=120, mc='magenta', level=2)
        info_single("A format command affects all subsequent plots.", mc='magenta', level=2)
        info_single("A format command starts with f;", mc='magenta', level=2)
        info_single("Example: f; sample_rel_abund_size, 2; group_id_y, 1.1; seed_count_y, 1.1", mc='magenta', cut_after=120, level=2)
        info_single("A format command should not include a command to make a plot.", mc='magenta', level=2)
        info_single("Invalid: f; sample_rel_abund_size, 2; class", mc='red', level=2)
        info_single("Reset format options to their default values: f; reset", mc='magenta', level=2, nl_after=1)

        info_single("Type f for format options.")
        info_single("Type h for this help message.")
        info_single("Type q to quit.", nl_after=2)


    def print_formatting_options(self):
        self.run.warning("", header="FORMAT OPTIONS", lc='green')
        info_single = self.run.info_single
        DEFAULT_FORMAT_PARAM_DICT = self.DEFAULT_FORMAT_PARAM_DICT
        info_single("The name and default value of the option is given.", mc='magenta', level=1)
        for param, value in self.DEFAULT_FORMAT_PARAM_DICT.items():
            info_single(f"{param}, {value}", mc='yellow', level=2)
        print("\n")


    def plot_rank(self, rank, spec_df, mod_df, aa, anticodon, out_dir, nonspec_df=None):
        spec_gb = spec_df.groupby(rank)
        mod_gb = mod_df.groupby(rank)
        nonspec_gb = nonspec_df.groupby(rank) if self.nonspec_txt_path else None

        for taxon, taxon_spec_df in spec_gb:
            try:
                taxon_mod_df = mod_gb.get_group(taxon)
            except KeyError:
                taxon_mod_df = pd.DataFrame(columns=self.mod_df.columns)
            taxon_nonspec_df = nonspec_gb.get_group(taxon) if self.nonspec_txt_path else None

            self.plot_items(taxon_spec_df, taxon_mod_df, aa, anticodon, rank, taxon, out_dir, nonspec_df=taxon_nonspec_df)


    def plot_items(self, spec_df, mod_df, aa, anticodon, rank, taxon, out_dir, nonspec_df=None):
        if spec_df.empty:
            if self.nonspec_txt_path:
                if nonspec_df.empty:
                    self.run.info_single(f"No values for {rank} {taxon} {aa}-{anticodon}")
                    return
            else:
                self.run.info_single(f"No values for {rank} {taxon} {aa}-{anticodon}")
                return

        sample_spec_gb = spec_df.groupby(('sample_name', '', ''))
        sample_mod_gb = mod_df.groupby('sample_name')
        sample_nonspec_gb = nonspec_df.groupby(('sample_name', '', '')) if self.nonspec_txt_path else None

        fig = plt.figure(1)
        num_samples = self.num_samples
        GridSpec(num_samples, 1)

        bottom_cov_ax = plt.subplot2grid((num_samples, 1), (num_samples - 1, 0))
        cov_axs = []
        for grid_row in range(0, num_samples - 1):
            cov_ax = plt.subplot2grid((num_samples, 1), (grid_row, 0), sharex=bottom_cov_ax)
            cov_ax.tick_params(bottom=False, labelbottom=False)
            cov_axs.append(cov_ax)
        cov_axs.append(bottom_cov_ax)

        plt.subplots_adjust(hspace=0)

        format_param_dict = self.format_param_dict
        cov_index = self.cov_index
        cov_cutoff = self.cov_cutoff
        cov_length = self.cov_length
        cov_x_values = self.cov_x_values
        cov_x_labels = self.cov_x_labels
        NT_COLORS = self.NT_COLORS
        mod_y_ticks = self.mod_y_ticks
        ax_count = -1
        for sample_name, cov_ax in zip(self.sample_names, cov_axs):
            ax_count += 1
            try:
                sample_spec_df = sample_spec_gb.get_group(sample_name)
            except KeyError:
                raise ConfigError("The input seeds specific coverage table appears not to have been generated correctly -- "
                                  "not every seed has a row for every sample. "
                                  "If you (wisely) generated the table with `anvi-tabulate-trnaseq` please contact the developer.")

            sum_rel_mean_spec_cov = sample_spec_df['relative_mean_coverage'].sum(axis=0)
            if sum_rel_mean_spec_cov == 0:
                cov_ax.tick_params(left=False, right=False, bottom=False, labelleft=False, labelright=False)
                cov_ax.annotate(f"{sample_name}",
                                xy=(1.04, 1),
                                xycoords=('axes fraction', 'axes fraction'),
                                va='top',
                                fontsize=format_param_dict['sample_rel_abund_font_size'])
                continue

            # Display seed coverages from bottom to top in descending order of mean specific coverage.
            sample_spec_df = sample_spec_df.sort_values('mean_coverage', ascending=False)

            spec_cov_array = sample_spec_df[cov_index].values
            cov_y_values = np.zeros(self.cov_length)
            for spec_covs, mean_spec_cov in zip(spec_cov_array, sample_spec_df['mean_coverage']):
                cov_y_values += spec_covs
                cov_ax.plot(cov_x_values, cov_y_values, linewidth=0.15, color='gray')

            # Plot a line showing the total specific coverage of all seeds under consideration.
            cov_ax.plot(cov_x_values, spec_cov_array.sum(axis=0), linewidth=0.15, color='purple')

            cov_ax.set_ylim(0)
            cov_ax.yaxis.set_tick_params(labelsize=format_param_dict['cov_tick_font_size'], length=2, pad=0)
            if ax_count == 0 and ax_count == len(cov_axs) - 1:
                cov_ax.yaxis.set_major_locator(MaxNLocator(nbins=3, steps=[5], integer=True, min_n_ticks=2))
            elif ax_count == 0:
                cov_ax.yaxis.set_major_locator(MaxNLocator(nbins=3, steps=[5], integer=True, prune='lower', min_n_ticks=2))
            elif ax_count == len(cov_axs) - 1:
                cov_ax.yaxis.set_major_locator(MaxNLocator(nbins=3, steps=[5], integer=True, prune='upper', min_n_ticks=2))
            else:
                cov_ax.yaxis.set_major_locator(MaxNLocator(nbins=3, steps=[5], integer=True, prune='both', min_n_ticks=2))

            sum_rel_discriminator_spec_cov = sample_spec_df['relative_discriminator_coverage'].sum(axis=0)
            cov_ax.annotate(f"{sample_name}\n\nRel abunds\nmean: {sum_rel_mean_spec_cov:.1e}\n3': {sum_rel_discriminator_spec_cov:.1e}",
                            xy=(1.04, 0.98),
                            xycoords=('axes fraction', 'axes fraction'),
                            va='top',
                            fontsize=format_param_dict['sample_rel_abund_font_size'])

            mod_ax = cov_ax.twinx()
            try:
                sample_mod_df = sample_mod_gb.get_group(sample_name)
            except KeyError:
                # The seeds under consideration have no modifications.
                mod_ax.tick_params(right=False, labelright=False)
                continue

            for ordinal_pos, pos_df in sample_mod_df.groupby('ordinal_position'):
                try:
                    total_pos_cov = spec_cov_array[:, ordinal_pos - 1].sum()
                except IndexError:
                    # The ordinal position is NaN when the modification position is 5' of the
                    # 5'-most feature defined in the seed.
                    continue
                mod_nt_covs = pos_df[UNAMBIG_NTS_LIST].sum(axis=0)
                nt_order = np.argsort(mod_nt_covs)[::-1]
                nt_colors = [NT_COLORS[nt_int] for nt_int in nt_order]
                prev_mod_nt_freq = 0
                for mod_nt_cov, nt_color in zip(mod_nt_covs[nt_order], nt_colors):
                    mod_nt_freq = mod_nt_cov / total_pos_cov
                    mod_ax.bar([ordinal_pos - 1], prev_mod_nt_freq + mod_nt_freq, width=0.5, bottom=prev_mod_nt_freq, color=nt_color)
                    prev_mod_nt_freq += mod_nt_freq

            mod_ax.set_ylim(0, 1)
            if ax_count == 0 and ax_count == len(cov_axs) - 1:
                mod_ax.set_yticks(mod_y_ticks)
            elif ax_count == 0:
                mod_ax.set_yticks(mod_y_ticks[1: ])
            elif ax_count == len(cov_axs) - 1:
                mod_ax.set_yticks(mod_y_ticks[: -1])
            else:
                mod_ax.set_yticks(mod_y_ticks[1: -1])
            mod_ax.yaxis.set_tick_params(labelsize=format_param_dict['mod_frac_tick_font_size'], length=2, pad=0)

        bottom_cov_ax.set_xlim(-1, cov_length)
        bottom_cov_ax.set_xticks(cov_x_values)
        bottom_cov_ax.set_xticklabels(cov_x_labels, fontsize=format_param_dict['ordinal_pos_tick_font_size'], rotation='vertical')
        bottom_cov_ax.xaxis.set_tick_params(length=2, pad=0)

        top_cov_ax = cov_axs[0]
        # Place a color key in the upper left corner for modification nt coverages.
        NT_COLORS = self.NT_COLORS
        nt_num = 0
        for nt, nt_color in zip(UNAMBIG_NTS, NT_COLORS):
            top_cov_ax.annotate(f"{nt}",
                                xy=(0.01 + 0.007 * nt_num, format_param_dict['nt_key_y']),
                                xycoords=('axes fraction', 'axes fraction'),
                                fontsize=format_param_dict['nt_key_font_size'],
                                ha='center',
                                color=nt_color)
            nt_num += 1

        top_cov_ax.annotate(f"{aa}-{anticodon}\n{rank.capitalize()} {taxon}",
                            xy=(0.05, format_param_dict['group_id_y']),
                            xycoords=('axes fraction', 'axes fraction'),
                            fontsize=format_param_dict['group_id_font_size'])

        top_cov_ax.annotate(f"Seeds: {pp(len(sample_spec_df))}",
                            xy=(1, format_param_dict['seed_count_y']),
                            xycoords=('axes fraction', 'axes fraction'),
                            ha='right',
                            fontsize=format_param_dict['seed_count_font_size'])

        filename = f"{aa}_{anticodon}_{rank}_{taxon}.png"
        out_path = os.path.join(out_dir, filename)
        overwritten = True if os.path.exists(out_path) else False
        plt.savefig(out_path, dpi=600, bbox_inches='tight')
        self.run.info_single(f"{out_path}{' overwritten' if overwritten else ''}", mc='green')
