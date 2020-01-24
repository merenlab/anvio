# -*- coding: utf-8
# pylint: disable=line-too-long
"""tRNA identification from tRNA-seq sequences -- classes adapted from tRNA-seq-tools"""

import anvio.db as db
import anvio.dbops as dbops
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths
import anvio.tables as tables
import anvio.terminal as terminal
import anvio.utils as utils

import itertools
import os.path

from anvio.constants import WC_plus_wobble_base_pairs as WC_PLUS_WOBBLE_BASE_PAIRS
from anvio.errors import ConfigError
from collections import Counter, OrderedDict

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


# The following are comprehensive guidelines for describing tRNA
ACCEPTOR_ARM_GUIDELINES = {
    'tRNA-His position -1 base': 'G',
    'acceptor stem length': 7,
    'acceptor position 1 base': 'C',
    'acceptor position 2 base': 'C',
    'acceptor position 3 base': 'A',
    'max acceptor stem pair mismatches': 1}
D_ARM_GUIDELINES = {
    'D stem length range': (3, 4),
    'D loop length range': (7, 11),
    'position 14 base': 'A',
    'position 15 bases': ('A', 'G'),
    'alpha positions length range': (1, 3),
    'position 18 base': 'G',
    'position 19 base': 'G',
    'beta positions length range': (1, 3),
    'position 21 bases': ('A', 'G'),
    'max D stem pair mismatches': 1,
    'max D loop anomalies': 1,}
ANTICODON_ARM_GUIDELINES = {
    'anticodon stem length': 5,
    'anticodon loop length': 7,
    'position 33 base': 'T',
    'position 37 bases': ('A', 'G'),
    'max anticodon arm anomalies': 1}
V_ARM_GUIDELINES = {
    'Type I length range': (4, 5),
    'Type II length range': (12, 23)}
T_ARM_GUIDELINES = {
    'T stem length': 5,
    'T loop length': 7,
    'position 53 base': 'G',
    'position 54 base': 'T',
    'position 55 base': 'T',
    'position 56 base': 'C',
    'position 57 bases': ('A', 'G'),
    'position 58 base': 'A',
    'position 61 base': 'C',
    'max T stem pair mismatches': 1,
    'max T loop plus positions 53 and 61 anomalies': 1}
GENERAL_GUIDELINES = {
    'min length': 24}

GUIDELINE_DESCRIPTIONS = {
    'tRNA-His position -1 base': ("If the tRNA has a histidine anticodon, "
                                  "there is an extra G at the start of the tRNA."),
    'acceptor stem length': "There are 7 base pairs before the discriminator in the acceptor stem.",
    'acceptor position 1 base': "The first base of the acceptor is C.",
    'acceptor position 2 base': "The second base of the acceptor is C.",
    'acceptor position 3 base': "The third base of the acceptor is A.",
    'max acceptor stem pair mismatches': "At most 1 of 7 base pairs in the acceptor stem should be missing.",
    'D stem length range': "The D stem contains 3-4 base pairs.",
    'D loop length range': "The D loop contains 7-11 bases.",
    'position 14 base': "Position 14 is G.",
    'position 15 bases': "Position 15 is a purine.",
    'alpha positions length range': "The alpha site of the D loop contains 1-3 bases (positions 16, 17, and 17A).",
    'position 18 base': "Position 18 is G.",
    'position 19 base': "Position 19 is G.",
    'beta positions length range': "The beta site of the D loop contains 1-3 bases (positions 20, 20A, and 20B).",
    'position 21 bases': "Position 21 is a purine.",
    'max D stem pair mismatches': "At most 1 of 3-4 base pairs in the D stem should be missing.",
    'max D loop anomalies': "At most 1 of 5 conserved bases in the D loop should differ from expectation.",
    'anticodon stem length': "There are 5 base pairs in the anticodon stem.",
    'anticodon loop length': "There are 7 bases in the anticodon loop.",
    'position 33 base': "Position 33 is T.",
    'position 37 bases': "Position 37 is a purine.",
    'max anticodon arm anomalies': ("At most 1 of 5 base pairs in the anticodon stem should be missing "
                                    "or at most 1 of the 2 conserved bases in the anticodon loop should differ from expectation."),
    'Type I length range': "The V arm is 4-5 nucleotides in Type I tRNA.",
    'Type II length range': "The V arm is 12-23 nucleotides in Type II tRNA.",
    'T stem length': "There are 5 base pairs in the T stem.",
    'T loop length': "There are 7 bases in the T loop.",
    'position 53 base': "Position 53 is G.",
    'position 54 base': "Position 54 is T.",
    'position 55 base': "Position 55 is T.",
    'position 56 base': "Position 56 is C.",
    'position 57 bases': "Position 57 is a purine.",
    'position 58 base': "Position 58 is A.",
    'position 61 base': "Position 61 is C.",
    'max T stem pair mismatches': "At most 1 of 5 base pairs in the T stem should be missing.",
    'max T loop plus positions 53 and 61 anomalies': ("At most 1 base of 7 conserved bases in the T loop "
                                                      "and the flanking bases at positions 53 and 61 should differ from expectation."),
    'min length': ("To be identifiable as tRNA, a tRNA-seq sequence should be at least 24 bases long, "
                   "encompassing the conserved base at position 53.")}

# Distances between benchmark bases in tRNA
# All combinations of variable lengths within the D arm
D_ARM_VARIABLE_LENGTH_COMBOS = tuple(itertools.product(
    D_ARM_GUIDELINES['D stem length range'],
    D_ARM_GUIDELINES['alpha positions length range'],
    D_ARM_GUIDELINES['beta positions length range']))
# Range of possible distances from position 27 at the start of the anticodon arm
# to position 49 at the start of the T arm
POSITION_27_TO_POSITION_49_RANGE = (
    (ANTICODON_ARM_GUIDELINES['anticodon stem length']
     + ANTICODON_ARM_GUIDELINES['anticodon loop length']
     + ANTICODON_ARM_GUIDELINES['anticodon stem length']
     + V_ARM_GUIDELINES['Type I length range'][0]),
    (ANTICODON_ARM_GUIDELINES['anticodon stem length']
     + ANTICODON_ARM_GUIDELINES['anticodon loop length']
     + ANTICODON_ARM_GUIDELINES['anticodon stem length']
     + V_ARM_GUIDELINES['Type I length range'][1]))
# Distance from the first conserved position in the T arm to the end of the tRNA
POSITION_53_TO_THREEPRIME_DISTANCE = (
    T_ARM_GUIDELINES['T loop length']
    + T_ARM_GUIDELINES['T stem length']
    + ACCEPTOR_ARM_GUIDELINES['acceptor stem length']
    + 4)


class Sorter:
    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.sample_name = A('sample_name')
        self.input_fasta_path = A('input_fasta')
        self.output_db_path = A('output_db_path')

        self.run = terminal.Run()
        self.progress = terminal.Progress()

        self.stats_dict = Counter()

        self.extractor = Extractor()
        self.db = None
        self.seq_count_dict = {}

    def process(self):
        self.screen_user_args()

        # Empty tRNA db
        tRNA_db = dbops.tRNADatabase(self.output_db_path)
        tRNA_db.create(meta_values={'sample_name': self.sample_name})

        # list buffer for results
        results_buffer = []

        # Arbitrary max size to store and reset buffer
        memory_max = 2000000

        # Create a REJECTED_SEQUENCES tabular output.
        reject_txt_path = os.path.join(
            os.path.dirname(self.output_db_path),
            filesnpaths.get_name_from_file_path(self.output_db_path) + 'REJECTED_SEQUENCES.txt')

        tRNA_profiler = tRNAProfiler(reject_txt_path)

        input_fasta = fastalib.SequenceSource(self.input_fasta_path)

        self.run.info("Sample name", self.sample_name)
        self.run.info("Input FASTA", self.input_fasta_path)

        self.progress.new("Identifying tRNA features")
        self.progress.update("...")
        while next(input_fasta):

            self.stats_dict['total_seqs'] += 1
            seq = input_fasta.seq

            current_tRNA_features = tRNAFeatures()
            is_tRNA = tRNA_profiler.check_if_tRNA(seq, input_fasta.id, current_tRNA_features)
            if is_tRNA:
                tRNA_profiler.find_remaining_tRNA_features(seq, current_tRNA_features)
            else:
                self.stats_dict['total_rejected'] += 1

            seq_length = len(seq)
            current_tRNA_features = tRNAFeatures(tRNA_profiler.filters)
            current_tRNA_features.seq_length = seq_length


    def screen_user_args():
        if not self.sample_name:
            raise ConfigError("You must provide a sample name"
                              "(which should uniquely describe this dataset).")
        if '-' in self.sample_name:
            self.sample_name = self.sample_name.replace('-', '_')
            self.run.warning("I just replaced all '-' characters with '_' characters in your sample name."
                             "This program does not like '-' characters in sample names.")
        utils.check_sample_id(self.sample_name)

        filesnpaths.is_file_fasta_formatted(self.input_fasta_path)
        self.input_fasta_path = os.path.abspath(self.input_fasta_path)

        if not self.output_db_path:
            raise ConfigError("You must provide an output database file path"
                              "(which should end with extension '.db').")
        if not self.output_db_path.endswith('.db'):
            raise ConfigError("The output database file name must end with '.db'.")
        if filesnpaths.is_file_exists(self.output_db_path, dont_raise=True):
            raise ConfigError("The output file already exists."
                              "We don't like overwriting stuff here :/")
        filesnpaths.is_output_file_writable(self.output_db_path)


class Extractor:
    def __init__(self):
        self.extractor_stats_file = ""
        self.extractor_stats = ExtractorStats()


class ExtractorStats:
    def __init__(self):
        self.total_seqs = 0
        self.type_I_seqs = 0
        self.type_II_seqs = 0
        self.type_I_V_arm_lengths = dict(
            [(length, 0) for length in V_ARM_GUIDELINES['Type I length range']])
        self.type_II_V_arm_lengths = dict(
            [(length, 0) for length in V_ARM_GUIDELINES['Type II length range']])


def check_position(seq, index, bases, allow_neg_index=True):
    if not allow_neg_index:
        if index < 0:
            return False
    try:
        if seq[index] in bases:
            return True
    except IndexError:
        return False
    return False


def find_stem_pair_mismatches(seq, stem_start_index, stem_length, loop_length, base_pairings):
    stem_mismatch_indices = []
    for fiveprime_index, threeprime_index in zip(
        range(stem_start_index, stem_start_index + stem_length),
        range(stem_start_index + stem_length + loop_length,
              stem_start_index + stem_length + loop_length + stem_length)):
        if seq[fiveprime_index] not in base_pairings[seq[threeprime_index]]:
            stem_mismatch_indices.append((fiveprime_index, threeprime_index))
    return stem_mismatch_indices


def find_acceptor_anomalies(acceptor_seq, start_index):
    acceptor_seq = 'N' * start_index + acceptor_seq
    acceptor_anomalies = []
    if acceptor_seq[0] != ACCEPTOR_ARM_GUIDELINES['acceptor position 1 base']:
        acceptor_anomalies.append(0)
    if acceptor_seq[1] != ACCEPTOR_ARM_GUIDELINES['acceptor position 2 base']:
        acceptor_anomalies.append(1)
    if acceptor_seq[2] != ACCEPTOR_ARM_GUIDELINES['acceptor position 3 base']:
        acceptor_anomalies.append(2)
    return acceptor_anomalies


def find_T_loop_anomalies(T_loop_seq, start_index=0):
    T_loop_seq = start_index * 'N' + T_loop_seq
    T_loop_anomalies = []
    if not check_position(T_loop_seq, 0, T_ARM_GUIDELINES['position 54 base']):
        T_loop_anomalies.append(0)
    if not check_position(T_loop_seq, 1, T_ARM_GUIDELINES['position 55 base']):
        T_loop_anomalies.append(1)
    if not check_position(T_loop_seq, 2, T_ARM_GUIDELINES['position 56 base']):
        T_loop_anomalies.append(2)
    if not check_position(T_loop_seq, 3, T_ARM_GUIDELINES['position 57 bases']):
        T_loop_anomalies.append(3)
    if not check_position(T_loop_seq, 4, T_ARM_GUIDELINES['position 58 bases']):
        T_loop_anomalies.append(4)
    if not check_position(T_loop_seq, 7, T_ARM_GUIDELINES['position 61 bases']):
        T_loop_anomalies.append(7)
    return T_loop_anomalies


def find_anticodon_loop_anomalies(anticodon_loop_seq):
    anticodon_loop_anomalies = []
    if not check_position(anticodon_loop_seq, 1, ANTICODON_ARM_GUIDELINES['position 33 base']):
        anticodon_loop_anomalies.append(1)
    if not check_position(anticodon_loop_seq, 5, ANTICODON_ARM_GUIDELINES['position 37 base']):
        anticodon_loop_anomalies.append(5)
    return anticodon_loop_anomalies


def find_D_loop_anomalies(D_loop_seq, alpha_length, beta_length, start_index=0):
    D_loop_seq = start_index * 'N' + D_loop_seq
    D_loop_anomalies = []
    if not check_position(D_loop_seq, 0, D_ARM_GUIDELINES['position 14 base']):
        D_loop_anomalies.append(0)
    if not check_position(D_loop_seq, 1, D_ARM_GUIDELINES['position 15 bases']):
        D_loop_anomalies.append(1)
    if not check_position(D_loop_seq, 2 + alpha_length, D_ARM_GUIDELINES['position 18 base']):
        D_loop_anomalies.append(2 + alpha_length)
    if not check_position(D_loop_seq, 3 + alpha_length, D_ARM_GUIDELINES['position 19 base']):
        D_loop_anomalies.append(3 + alpha_length)
    if not check_position(D_loop_seq, 4 + alpha_length + beta_length, D_ARM_GUIDELINES['position 21 bases']):
        D_loop_anomalies.append(4 + alpha_length + beta_length)
    return D_loop_anomalies


class tRNAProfiler:
    def __init__(self, reject_txt_path):
        # Sequence filters to screen for tRNA
        self.FILTERS = {}
        self.FILTERS['min length'] = GENERAL_GUIDELINES['min length']
        self.FILTERS['acceptor stem length'] = ACCEPTOR_ARM_GUIDELINES['acceptor stem length']
        self.FILTERS['acceptor position 1 base'] = ACCEPTOR_ARM_GUIDELINES['acceptor position 1 base']
        self.FILTERS['acceptor position 2 base'] = ACCEPTOR_ARM_GUIDELINES['acceptor position 2 base']
        self.FILTERS['acceptor position 3 base'] = ACCEPTOR_ARM_GUIDELINES['acceptor position 3 base']
        self.FILTERS['T stem length'] = T_ARM_GUIDELINES['T stem length']
        self.FILTERS['max T loop plus positions 53 and 61 anomalies'] = T_ARM_GUIDELINES[
            'max T loop plus positions 53 and 61 anomalies']

        # Table of rejection information
        self.reject_txt_path = reject_txt_path


    def check_if_tRNA(self, seq, seq_name, current_tRNA_features):
        is_tRNA = True
        # Store the default value if the criterion is met.
        # Otherwise store the aberrant value.
        reject_statuses = {}


        # General filter
        if len(seq) >= self.FILTERS['min length']:
            reject_statuses['min length'] = self.FILTERS['min length']
        else:
            is_tRNA = False
            reject_statuses['min length'] = len(seq)


        # Acceptor filter
        acceptor_seq = seq[max([-len(seq), -3]):]
        acceptor_start = 3 - min(len(seq), 3)
        acceptor_anomalies = find_acceptor_anomalies(acceptor_seq, acceptor_start)
        if acceptor_anomalies:
            is_tRNA = False
            # First acceptor position
            if len(acceptor_seq) < 3:
                reject_statuses['acceptor position 1 base'] = 'NA'
            else:
                if acceptor_seq[0] == self.FILTERS['acceptor position 1 base']:
                    reject_statuses['acceptor position 1 base'] = self.FILTERS['acceptor position 1 base']
                else:
                    reject_statuses['acceptor position 1 base'] = acceptor_seq[0]
            # Second acceptor position
            if len(acceptor_seq) < 2:
                reject_statuses['acceptor position 2 base'] = 'NA'
            else:
                if acceptor_seq[1] == self.FILTERS['acceptor position 2 base']:
                    reject_statuses['acceptor position 2 base'] = self.FILTERS['acceptor position 2 base']
                else:
                    reject_statuses['acceptor position 2 base'] = acceptor_seq[1]
            # Third acceptor position
            if acceptor_seq[2] == self.FILTERS['acceptor position 3 base']:
                reject_statuses['acceptor position 3 base'] = self.FILTERS['acceptor position 3 base']
            else:
                reject_statuses['acceptor position 3 base'] = acceptor_seq[2]
        else:
            reject_statuses['acceptor position 1 base'] = self.FILTERS['acceptor position 1 base']
            reject_statuses['acceptor position 2 base'] = self.FILTERS['acceptor position 2 base']
            reject_statuses['acceptor position 3 base'] = self.FILTERS['acceptor position 3 base']


        # T loop filter
        # In the process of checking the T loop,
        # the positions of T loop bases in the sequence are recorded in current_tRNA_features.
        position_53_index = len(seq) - POSITION_53_TO_THREEPRIME_DISTANCE - 1
        has_position_53_anomaly = 1 - check_position(
            seq, position_53_index, T_ARM_GUIDELINES['position 53 base'], allow_neg_index=False)
        T_loop_seq = seq[max(position_53_index + 1, 0): max(position_53_index + 8, 0)]
        T_loop_start = abs(min(max(position_53_index + 1, -7), 0))
        T_loop_anomalies = find_T_loop_anomalies(T_loop_seq, T_loop_start)
        has_position_61_anomaly = 1 - check_position(
            seq, position_53_index + 8, T_ARM_GUIDELINES['position 61 base'], allow_neg_index=False)

        if (has_position_53_anomaly + len(T_loop_anomalies) + has_position_61_anomaly
            > self.FILTERS['max T loop plus positions 53 and 61 anomalies']):
            is_tRNA = False
            reject_statuses['max T loop plus positions 53 and 61 anomalies'] = (
                has_position_53_anomaly + len(T_loop_anomalies) + has_position_61_anomaly)
        else:
            # Congratulations, the sequence is a tRNA!
            current_tRNA_features.seq = seq
            current_tRNA_features.seq_length = len(seq)
            current_tRNA_features.T_loop_plus_positions_53_and_61_anomaly_indices = []
            if has_position_53_anomaly and position_53_index >= 0:
                current_tRNA_features.T_loop_plus_positions_53_and_61_anomaly_indices.append(position_53_index)
            if T_loop_anomalies:
                current_tRNA_features.T_loop_plus_positions_53_and_61_anomaly_indices.extend(
                    [position_53_index + 1 + T_loop_index for T_loop_index in T_loop_anomalies
                     if position_53_index + 1 + T_loop_index >= 0])
            if has_position_61_anomaly and position_61_index >= 0:
                current_tRNA_features.T_loop_plus_positions_53_and_61_anomaly_indices.append(position_53_index + 8)
            if position_53_index + 1 >= 0:
                current_tRNA_features.T_loop_bounds = (position_53_index + 1, position_53_index + 7)
            reject_statuses['max T loop plus positions 53 and 61 anomalies'] = FILTERS[
                'max T loop plus positions 53 and 61 anomalies']

        if position_53_index < 0:
            reject_statuses['position 53 base'] = 'NA'
        else:
            if has_position_53_anomaly:
                reject_statuses['position 53 base'] = seq[position_53_index]
            else:
                reject_statuses['position 53 base'] = FILTERS['position 53 base']
        if len(T_loop_seq) < 7:
            reject_statuses['position 54 base'] = 'NA'
        else:
            if 0 in T_loop_anomalies:
                reject_statuses['position 54 base'] = T_loop_seq[0]
            else:
                reject_statuses['position 54 base'] = FILTERS['position 54 base']
        if len(T_loop_seq) < 6:
            reject_statuses['position 55 base'] = 'NA'
        else:
            if 1 in T_loop_anomalies:
                reject_statuses['position 55 base'] = T_loop_seq[1]
            else:
                reject_statuses['position 55 base'] = FILTERS['position 55 base']
        if len(T_loop_seq) < 5:
            reject_statuses['position 56 base'] = 'NA'
        else:
            if 2 in T_loop_anomalies:
                reject_statuses['position 56 base'] = T_loop_seq[2]
            else:
                reject_statuses['position 56 base'] = FILTERS['position 56 base']
        if len(T_loop_seq) < 4:
            reject_statuses['position 57 base'] = 'NA'
        else:
            if 3 in T_loop_anomalies:
                reject_statuses['position 57 base'] = T_loop_seq[3]
            else:
                reject_statuses['position 57 base'] = FILTERS['position 57 bases']
        if len(T_loop_seq) < 3:
            reject_statuses['position 58 base'] = 'NA'
        else:
            if 4 in T_loop_anomalies:
                reject_statuses['position 58 base'] = T_loop_seq[4]
            else:
                reject_statuses['position 58 base'] = FILTERS['position 58 base']
        if position_53_index + 8 < 0:
            reject_statuses['position 61 base'] = 'NA'
        else:
            if has_position_61_anomaly:
                reject_statuses['position 61 base'] = seq[position_53_index + 8]
            else:
                reject_statuses['position 61 base'] = FILTERS['position 61 base']


        # Write rejection information.
        reject_info_row = seq_name + "\t" + seq + "\t"
        for reject_guideline, correct_value in self.FILTERS.items():
            if reject_statuses[reject_guideline] == correct_value:
                reject_info_row += "\t"
            else:
                reject_info_row += reject_statuses[reject_guideline] + "\t"
        reject_info_row += "\n"

        with open(self.reject_txt_path, 'a') as f:
            f.write(reject_info_row)

        return is_tRNA


    def find_remaining_tRNA_features(self, seq, current_tRNA_features):
        # The sequence string, sequence lenth, bounds of the T loop,
        # and anomalies in the T loop and positions 53 and 61
        # were added as features by check_if_tRNA.

        # T stem
        position_49_index = current_tRNA_features.T_loop_bounds[0] - 5
        # Ensure that the sequence is long enough to include the full stem.
        if position_49_index < 0:
            return
        T_stem_mismatch_indices = find_stem_pair_mismatches(
            seq, position_49_index, 5, 7, WC_PLUS_WOBBLE_BASE_PAIRS)
        if len(T_stem_mismatch_indices) <= T_ARM_GUIDELINES['max T stem pair mismatches']:
            # Record mismatched indices in order from 5' to 3' along the sequence.
            current_tRNA_features.T_stem_mismatch_indices.extend(
                [mismatch_pair[0] for mismatch_pair in stem_mismatch_indices])
            current_tRNA_features.T_stem_mismatch_indices.extend(
                [mismatch_pair[1] for mismatch_pair in reversed(stem_mismatch_indices)])
            # Bounds of 5' half of stem
            current_tRNA_features.T_stem_bounds.append((position_49_index, position_49_index + T_STEM_LENGTH - 1))
            # Bounds of 3' half of stem
            current_tRNA_features.T_stem_bounds.append(
                (current_tRNA_features.T_loop_bounds[1] + 1,
                 current_tRNA_features.T_loop_bounds[1] + T_STEM_LENGTH))
        # Continue searching for features before the T stem even if the T stem was not identified.


        # Anticodon arm
        # Since there are few landmark bases in the anticodon loop,
        # it can only be reliably identified along with the anticodon stem.
        if position_49_index - POSITION_27_TO_POSITION_49_RANGE[0] < 0:
            # The anticodon arm cannot be contained within this fragmentary tRNA sequence.
            return
        if position_49_index - POSITION_27_TO_POSITION_49_RANGE[1] < 0:
            position_27_min_index = 0
        else:
            position_27_min_index = position_49_index - POSITION_27_TO_POSITION_49_RANGE[1]
        position_27_max_index = position_49_index - POSITION_27_TO_POSITION_49_RANGE[0]

        # Find candidate anticodon arms.
        position_27_index_possibilities = []
        for position_27_index in range(position_27_min_index, position_27_max_index):
            anticodon_stem_mismatch_indices = find_stem_pair_mismatches(
                seq, position_27_index, 5, 7, WC_PLUS_WOBBLE_BASE_PAIRS)
            anticodon_loop_anomalies = find_anticodon_loop_anomalies(
                seq[position_27_index + 5: position_27_index + 12])
            if (len(anticodon_stem_mismatch_indices) + len(anticodon_loop_anomalies)
                <= ANTICODON_ARM_GUIDELINES['max anticodon arm anomalies']):
                position_27_index_possibilities.append(position_27_index)

        if len(position_27_index_possibilities) == 0:
            # No anticodon arm candidate was found.
            return
        elif len(position_27_index_possibilities) == 1:
            # One anticodon arm candidate was found -- keep it.
            # Bounds of 5' half of stem
            current_tRNA_features.anticodon_stem_bounds.append((position_27_index, position_27_index + 4))
            # Bounds of 3' half of stem
            current_tRNA_features.anticodon_stem_bounds.append((position_27_index + 12, position_27_index + 16))
            current_tRNA_features.anticodon_loop_bounds.append((position_27_index + 5, position_27_index + 11))
            current_tRNA_features.V_arm_bounds.append((position_27_index + 17, current_tRNA_features.T_stem_bounds[0][0] - 1))
            # Record mismatched indices in order from 5' to 3' along the sequence.
            current_tRNA_features.anticodon_stem_mismatch_indices.extend(
                [mismatch_pair[0] for mismatch_pair in anticodon_stem_mismatch_indices])
            current_tRNA_features.anticodon_stem_mismatch_indices.extend(
                [mismatch_pair[1] for mismatch_pair in anticodon_stem_mismatch_indices])
            current_tRNA_features.anticodon_anomaly_indices.extend(
                [position_27_index + 5 + anticodon_loop_index for anticodon_loop_index in anticodon_loop_anomalies])
        else:
            # Multiple anticodon arm candidates were found -- report them.
            progress.update("Multiple anticodon arm candidates "
                            "starting at indices %s were found in the following sequence: "
                            "%s" % (', '.join(map(str, position_27_index_possibilities)), seq))
            return


        # D arm
        # Check that a D loop can be accommodated by the sequence
        # by looking for the shortest possible distance from position 14 to position 27.
        position_27_index = current_tRNA_features.anticodon_stem_bounds[0][0]
        if position_27_index - 16 < 0:
            return

        # Find candidate D loops.
        D_loop_candidate_dict = {}
        for D_arm_variable_lengths in D_ARM_VARIABLE_LENGTH_COMBOS:
            position_14_index = position_27_index - sum(D_arm_variable_lengths) - 6
            alpha_length = D_arm_variable_lengths[1]
            beta_length = D_arm_variable_lengths[2]
            position_22_index = position_14_index + alpha_length + beta_length + 5
            D_loop_seq = seq[max(position_14_index, 0): max(position_22_index, 0)]
            D_loop_start = abs(max(min(position_14_index, 0), -(alpha_length + beta_length + 5)))
            D_loop_anomalies = find_D_loop_anomalies(
                D_loop_seq, alpha_length, beta_length, D_loop_start)
            if len(D_loop_anomalies) <= D_ARM_GUIDELINES['max D loop anomalies']:
                D_loop_candidate_dict[D_arm_variable_lengths] = D_loop_anomalies

        # Find candidate D stems.
        D_stem_candidate_dict = {}
        # Check that a D stem can be accommodated by the sequence.
        if position_27_index - 19 >= 0:
            for D_arm_variable_lengths in D_ARM_VARIABLE_LENGTH_COMBOS:
                D_stem_length = D_arm_variable_lengths[0]
                alpha_length = D_arm_variable_lengths[1]
                beta_length = D_arm_variable_lengths[2]
                position_10_index = position_27_index - sum(D_arm_variable_lengths) - D_stem_length
                D_stem_mismatch_indices = find_stem_pair_mismatches(
                    seq, position_10_index, D_stem_length, alpha_length + beta_length + 5, WC_PLUS_WOBBLE_BASE_PAIRS)
                if len(D_stem_mismatch_indices) <= D_ARM_GUIDELINES['max D stem pair mismatches']:
                    D_stem_candidate_dict[D_arm_variable_lengths] = D_stem_mismatch_indices

        # Try to select a D arm that contains an identifiable stem and loop and has the fewest anomalies.
        D_arm_candidates = []
        D_arm_min_anomalies = (
            D_ARM_GUIDELINES['max D loop anomalies'] + D_ARM_GUIDELINES['max D stem pair mismatches'] + 1)
        for D_arm_variable_lengths in D_ARM_VARIABLE_LENGTH_COMBOS:
            if D_arm_variable_length in D_loop_candidate_dict and D_arm_variable_length in D_stem_candidate_dict:
                if (D_loop_candidate_dict[D_arm_variable_length] + D_stem_candidate_dict[D_arm_variable_length]
                    <= D_arm_min_anomalies):
                    D_arm_candidates.append(D_arm_variable_lengths)
                    D_arm_min_anomalies = D_loop_candidate_dict[D_arm_variable_length] + D_stem_candidate_dict[D_arm_variable_length]
        if len(D_arm_candidates) == 0:
            # No D arm candidate.
        elif len(D_arm_candidates) == 1:
            # One clear D arm candidate.
            D_arm_variable_lengths = D_arm_candidates[0]
            current_tRNA_features.D_stem_bounds.append(
                (position_27_index - sum(D_arm_variable_lengths) - D_arm_variable_lengths[0] - 6,
                 position_27_index - sum(D_arm_variable_lengths[0]) - 6))
            current_tRNA_features.D_stem_bounds.append(
                (position_27_index - D_arm_variable_lengths[0] - 1, position_27_index - 2))
            current_tRNA_features.D_loop_bounds.append(
                (position_27_index - sum(D_arm_variable_lengths) - 6,
                 position_27_index - D_arm_variable_lengths[0] - 2))
            current_tRNA_features.D_stem_mismatch_indices.append(
                [mismatch_pair[0] for mismatch_pair in D_stem_mismatch_indices])
            current_tRNA_features.D_stem_mismatch_indices.append(
                [mismatch_pair[1] for mismatch_pair in D_stem_mismatch_indices])
            current_tRNA_features.D_loop_anomaly_indices.extend(
                [position_27_index - sum(D_arm_variable_lengths) - 6 + D_loop_index
                 for D_loop_index in D_loop_anomalies])
        else:
            # Multiple D arm candidates.
            position_10_index_candidates = [
                position_27_index - sum(D_arm_variable_lengths) - D_arm_variable_lengths[0] - 6
                for D_arm_variable_lengths in D_arm_candidates]
            progress.update(
                "Multiple D arm candidates starting at indices %s "
                "were found in the following sequence: "
                "%s" % (', '.join(map(str, position_10_index_candidates)), seq))


        D_STEM_LENGTH_RANGE = D_ARM_GUIDELINES['D stem length range']
        MAX_D_STEM_PAIR_MISMATCHES = D_ARM_GUIDELINES['max D stem pair mismatches']
        MAX_D_LOOP_ANOMALIES = D_ARM_GUIDELINES['max D loop anomalies']
        POSITION_14_BASE = D_ARM_GUIDELINES['position 14 base']
        POSITION_15_BASES = D_ARM_GUIDELINES['position 15 bases']
        ALPHA_POSITIONS_LENGTH_RANGE = D_ARM_GUIDELINES['alpha positions length range']
        POSITION_18_BASE = D_ARM_GUIDELINES['position 18 base']
        POSITION_19_BASE = D_ARM_GUIDELINES['position 19 base']
        BETA_POSITIONS_LENGTH_RANGE = D_ARM_GUIDELINES['beta positions length range']
        POSITION_21_BASES = D_ARM_GUIDELINES['position 21 bases']
        # Check all possible combinations for the D loop,
        # as there is a small probability of mistaken D loop identification.
        D_loop_anomaly_indices_dict = {}
        for D_arm_variable_lengths in self.D_ARM_VARIABLE_LENGTHS:
            # There are 5 "fixed length" bases in the D loop
            # and 1 base between the D stem and anticodon stem,
            # for 6 "fixed length" bases.
            position_14_index = position_27_index - sum(D_arm_variable_lengths) - 6
            if position_14_index < 0:
                # A D loop of the length being considered cannot be contained within this fragmentary tRNA sequence.
                continue

            D_loop_anomaly_indices = []
            stem_corresponds_to_loop = False
            # Check each conserved position in the D loop.
            if seq[position_14_index] != POSITION_14_BASE:
                D_loop_anomaly_indices.append(position_14_index)

            position_15_index = position_14_index + 1
            position_15_base = seq[position_15_index]
            if position_15_base != POSITION_15_BASES[0] and position_15_base != POSITION_15_BASES[1]:
                D_loop_anomaly_indices.append(position_15_index)

            position_18_index = position_15_index + D_arm_variable_lengths[1]
            if seq[position_18_index] != POSITION_18_BASE:
                D_loop_anomaly_indices.append(position_18_index)

            position_19_index = position_18_index + 1
            if seq[position_19_index] != POSITION_19_BASE:
                D_loop_anomaly_indices.append(position_19_index)

            position_21_index = position_19_index + D_arm_variable_lengths[2]
            position_21_base = seq[position_21_index]
            if position_21_base != POSITION_21_BASES[0] and position_21_base != POSITION_21_BASES[1]:
                D_loop_anomaly_indices.append(position_21_index)

            if len(D_loop_anomaly_indices) <= MAX_D_LOOP_ANOMALIES:
                D_loop_anomaly_indices_dict[D_arm_variable_lengths] = D_loop_anomaly_indices

        # For now, stop annotating tRNA if multiple possible D loops were found.
        if len(D_loop_anomaly_indices_dict) > 0:
            self.progress.update("Multiple possible D loops were found in %s" % seq)
            return
        D_arm_variable_lengths, D_loop_anomaly_indices = D_loop_anomaly_indices_dict.items()[0]
        current_tRNA_features.D_loop_bounds = (
            position_27_index - sum(D_arm_variable_lengths) - 6,
            position_27_index - 1 - D_arm_variable_lengths[0])
        current_tRNA_features.D_loop_anomaly_indices = D_loop_anomaly_indices

        # Look for the D stem.
        position_10_index = position_14_index - D_arm_variable_lengths[0]
        if position_10_index >= 0:
            D_stem_mismatch_indices = find_stem_pair_mismatches(
                seq,
                position_10_index,
                D_arm_variable_lengths[0],
                5 + D_arm_variable_lengths[1] + D_arm_variable_lengths[2])
        if len(D_stem_mismatch_indices) <= MAX_D_STEM_PAIR_MISMATCHES:
            current_tRNA_features.D_stem_mismatch_indices = D_stem_mismatch_indices
        else:
            return




class tRNAFeatures:
    def __init__(self):
        self.seq = ''
        self.seq_length = 0
        self.is_full_length = False

        # Feature indices start at the 3' end of the sequence
        # Acceptor stem features
        self.acceptor_stem_bounds = []
        self.fiveprime_leader_bounds = []
        self.acceptor_stem_mismatch_indices = []
        self.has_tRNAHis_fiveprime_G = False

        # T stem loop features
        self.T_stem_bounds = []
        self.T_loop_bounds = []
        self.T_stem_mismatch_indices = []
        self.T_loop_anomaly_indices = []

        # V loop features
        self.V_arm_bounds = []

        # Anticodon stem loop features
        self.anticodon_bounds = []
        self.anticodon_stem_bounds = []
        self.anticodon_loop_bounds = []
        self.anticodon_stem_mismatch_indices = []
        self.anticodon_loop_anomaly_indices = []

        # D stem loop features
        self.D_stem_bounds = []
        self.D_loop_bounds = []
        self.D_stem_mismatch_indices = []
        self.D_loop_anomaly_indices = []