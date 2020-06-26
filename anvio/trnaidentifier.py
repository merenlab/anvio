# -*- coding: utf-8
# pylint: disable=line-too-long
"""tRNA identification from tRNA-seq reads"""

import itertools

from anvio.constants import WC_plus_wobble_base_pairs as WC_PLUS_WOBBLE_BASE_PAIRS
from anvio.constants import anticodon_to_AA as ANTICODON_TO_AA
from anvio.errors import TransferRNAIdentifierError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


THREEPRIME_VARIANT_LIST = [
    'CCA',
    'C', 'CC',
    'CCC', 'CCG', 'CCT',
    'CAA', 'CGA', 'CTA',
    'ACA', 'GCA', 'TCA',
    'CCAA', 'CCAC', 'CCAG', 'CCAT',
    'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT',
    'CCACA', 'CCACC', 'CCACG', 'CCACT',
    'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT',
    'CCATA', 'CCATC', 'CCATG', 'CCATT']


class _TransferRNAFeature:
    conserved_nucleotides = ({}, )

    def __init__(
        self,
        string_components, # ex. ('CCA', )
        num_allowed_unconserved=-1, # ex. 0
        cautious=False):

        if cautious:
            if type(string_components) != tuple:
                raise TransferRNAIdentifierError(
                    "`string_components` must be in the form of a tuple, "
                    "e.g., ('ACTGG', 'CCAGT'). "
                    "Your `string_components` were %s"
                    % (string_components, ))
        self.string_components = string_components

        self.num_nucleotides = sum(map(len, string_components))

        # By default, base conservation is not enforced.
        if num_allowed_unconserved == -1:
            self.num_allowed_unconserved = sum(len(d) for d in self.conserved_nucleotides)
        else:
            self.num_allowed_unconserved = num_allowed_unconserved


    def check_conserved_nucleotides(self):
        num_conserved = 0
        num_unconserved = 0 # Can include N "padding" in extrapolated 5' feature
        conserved_status = []
        for substring, nuc_dict in zip(self.string_components, self.conserved_nucleotides):
            substring_statuses = []
            conserved_status.append(substring_statuses)
            for pos, expected_nucleotides in nuc_dict.items():
                try:
                    observed_nucleotide = substring[pos]
                except IndexError:
                    # This occurs when charging status is measured.
                    # When the acceptor is charged, the acceptor is changed to 'CC' in the read.
                    # The acceptor's conserved nucleotides include the 3' A, which must be ignored.
                    break
                if observed_nucleotide in expected_nucleotides:
                    num_conserved += 1
                    substring_statuses.append(
                        (pos, True, observed_nucleotide, expected_nucleotides))
                else:
                    num_unconserved += 1
                    substring_statuses.append(
                        (pos, False, observed_nucleotide, expected_nucleotides))
        if num_unconserved > self.num_allowed_unconserved:
            meets_conserved_thresh = False
        else:
            meets_conserved_thresh = True
        return (meets_conserved_thresh, num_conserved, num_unconserved, conserved_status)


    @staticmethod
    def list_all_tRNA_features():
        return [
            TransferRNAHisPositionZero,
            AcceptorStem,
            FiveprimeAcceptorStemSeq,
            PositionEight,
            PositionNine,
            DArm,
            DStem,
            FiveprimeDStemSeq,
            DLoop,
            ThreeprimeDStemSeq,
            PositionTwentySix,
            AnticodonArm,
            AnticodonStem,
            FiveprimeAnticodonStemSeq,
            AnticodonLoop,
            ThreeprimeAnticodonStemSeq,
            VLoop,
            TArm,
            TStem,
            FiveprimeTStemSeq,
            TLoop,
            ThreeprimeTStemSeq,
            ThreeprimeAcceptorStemSeq,
            Discriminator,
            Acceptor]


class _Nucleotide(_TransferRNAFeature):
    allowed_input_lengths = ((1, ), )
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))

    def __init__(
        self,
        string, # must be a string of length 1
        num_allowed_unconserved=True,
        start_index=None,
        stop_index=None,
        cautious=False):

        self.string = string
        self.start_index = start_index
        self.stop_index = stop_index

        super().__init__(
            (string, ),
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved,
         self.conserved_status) = self.check_conserved_nucleotides()


class _Sequence(_TransferRNAFeature):
    def __init__(
        self,
        substrings, # must be a string, tuple of strings, or tuple of _Nucleotide/_Sequence objects
        num_allowed_unconserved=-1,
        start_index=None,
        stop_index=None,
        cautious=False):

        if type(substrings) == str:
            string_components = (substrings, )
        elif all([type(s) == str for s in substrings]):
            string_components = substrings
        elif all([type(s) == _Nucleotide or type(s) == _Sequence for s in substrings]):
            string_components = tuple(s.string for s in substrings)
        self.string = ''.join(substrings)
        self.start_index = start_index
        self.stop_index = stop_index

        super().__init__(
            string_components,
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved, # Can include N "padding" in extrapolated 5' feature
         self.conserved_status) = self.check_conserved_nucleotides()


class _Loop(_Sequence):
    def __init__(
        self,
        substrings, # must be a string, tuple of strings, or tuple of _Nucleotide/_Sequence objects
        num_allowed_unconserved=-1,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class _Stem(_TransferRNAFeature):
    def __init__(
        self,
        fiveprime_seq,
        threeprime_seq,
        num_allowed_unpaired=0,
        num_allowed_unconserved=-1,
        cautious=False):

        if cautious:
            if type(fiveprime_seq) != _Sequence or type(threeprime_seq) != _Sequence:
                raise TransferRNAIdentifierError(
                    "You can only define a _Stem from _Sequence objects.")
        self.fiveprime_seq = fiveprime_seq
        self.threeprime_seq = threeprime_seq

        self.canonical_positions=(*self.fiveprime_seq.canonical_positions,
                                  *self.threeprime_seq.canonical_positions)
        self.conserved_nucleotides = (*self.fiveprime_seq.conserved_nucleotides,
                                      *self.threeprime_seq.conserved_nucleotides)

        if cautious:
            if (tuple(map(len, self.fiveprime_seq.string_components))
                != tuple(map(len, self.threeprime_seq.string_components[::-1]))):
                raise TransferRNAIdentifierError(
                    "The two _Sequence objects, %s and %s, "
                    "that were used to define your _Stem are not the same length."
                    % (self.fiveprime_seq.string_components, threeprime_seq.string_components))
            length = sum(map(len, self.fiveprime_seq.string_components))
            if num_allowed_unpaired > length:
                raise TransferRNAIdentifierError(
                    "You tried to leave at most %d base pairs unpaired, "
                    "but there are only %d base pairs in the stem."
                    % (num_allowed_unpaired, length))
        self.num_allowed_unpaired = num_allowed_unpaired

        self.start_indices = (self.fiveprime_seq.start_index, self.threeprime_seq.start_index)
        self.stop_indices = (self.fiveprime_seq.stop_index, self.threeprime_seq.stop_index)

        super().__init__(
            (*self.fiveprime_seq.string_components,
             *self.threeprime_seq.string_components),
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)

        (self.meets_pair_thresh,
         self.num_paired,
         self.num_unpaired, # Can include N "padding" in extrapolated 5' feature
         self.paired_status) = self.check_pairs()


    def check_pairs(self):
        num_paired = 0
        num_unpaired = 0 # Can include N "padding" in extrapolated 5' feature
        paired_status = []
        for fiveprime_nuc, threeprime_nuc in zip(
            self.fiveprime_seq.string, self.threeprime_seq.string[::-1]):
            if fiveprime_nuc in WC_PLUS_WOBBLE_BASE_PAIRS[threeprime_nuc]:
                num_paired += 1
                paired_status.append((True, fiveprime_nuc, threeprime_nuc))
            else:
                num_unpaired += 1
                paired_status.append((False, fiveprime_nuc, threeprime_nuc))
        if num_unpaired > self.num_allowed_unpaired:
            meets_pair_thresh = False
        else:
            meets_pair_thresh = True
        return (meets_pair_thresh, num_paired, num_unpaired, paired_status)


class _Arm(_TransferRNAFeature):
    def __init__(
        self,
        stem, # must be _Stem object
        loop, # must be _Loop object
        num_allowed_unconserved=-1,
        cautious=False):

        if cautious:
            if type(stem) != _Stem or type(loop) != _Loop:
                raise TransferRNAIdentifierError(
                    "A `_Stem` and a `_Loop` are required input to create an `_Arm`.")
        self.stem = stem
        self.loop = loop

        self.canonical_positions = (*stem.fiveprime_seq.canonical_positions,
                                    *loop.canonical_positions,
                                    *stem.threeprime_seq.canonical_positions)
        if cautious:
            if (tuple(pos for component_positions in canonical_positions
                      for pos in component_positions)
                != tuple(range(canonical_positions[0][0],
                               (canonical_positions[0][0]
                                + len(canonical_positions[0])
                                + len(canonical_positions[1])
                                + len(canonical_positions[2]))))):
                raise TransferRNAIdentifierError(
                    "The canonical positions in an `_Arm` must be contiguous. "
                    "These were yours: %s. "
                    "These came from the canonical positions in _Stem, %s, "
                    "and the canonical positions in _Loop, %s."
                    % (canonical_positions,
                       stem.canonical_positions,
                       loop.canonical_positions))

        self.conserved_nucleotides=(*stem.fiveprime_seq.conserved_nucleotides,
                                    *loop.conserved_nucleotides,
                                    *stem.threeprime_seq.conserved_nucleotides)

        self.start_index = self.stem.start_indices[0]
        self.stop_index = self.stem.stop_indices[1]

        super().__init__(
            (*stem.fiveprime_seq.string_components,
             *loop.string_components,
             *stem.threeprime_seq.string_components),
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)

        if (self.stem.fiveprime_seq.num_unconserved
            + self.loop.num_unconserved
            + self.stem.threeprime_seq.num_unconserved) > self.num_allowed_unconserved:
            self.meets_conserved_thresh = False
        else:
            self.meets_conserved_thresh = True


class TransferRNAHisPositionZero(_Nucleotide):
    name = 'tRNA-His Position 0'
    canonical_positions = ((-1, ), )
    conserved_nucleotides = ({0: 'G'}, )

    def __init__(
        self,
        string,
        num_allowed_unconserved=0,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            string,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class AcceptorStem(_Stem):
    # For our purposes, the acceptor stem only includes the base-paired nucleotides in the stem.
    name = 'Acceptor Stem'

    def __init__(
        self,
        fiveprime_seq,
        threeprime_seq,
        num_allowed_unpaired=1,
        cautious=False):

        super().__init__(
            fiveprime_seq,
            threeprime_seq,
            num_allowed_unpaired=num_allowed_unpaired,
            cautious=cautious)


class FiveprimeAcceptorStemSeq(_Sequence):
    name = '5\' Acceptor Stem Sequence'
    canonical_positions = ((1, 2, 3, 4, 5, 6, 7), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({}, )
    stem_class = AcceptorStem

    def __init__(
        self,
        substrings,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)

        if len(self.string) != 7:
            raise TransferRNAIdentifierError(
                "Your `FiveprimeAcceptorSeq` was not the required 7 bases long: %s" % self.string)


class PositionEight(_Nucleotide):
    name = 'Position 8'
    canonical_positions = ((8, ), )
    conserved_nucleotides = ({0: 'T'}, )

    def __init__(
        self,
        string,
        num_allowed_unconserved=1,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            string,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class PositionNine(_Nucleotide):
    name = 'Position 9'
    canonical_positions = ((9, ), )

    def __init__(
        self,
        string,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            string,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class DArm(_Arm):
    name = 'D Arm'

    def __init__(
        self,
        stem,
        loop,
        num_allowed_unconserved=2,
        cautious=False):

        super().__init__(
            stem,
            loop,
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)


class DStem(_Stem):
    name = 'D Stem'
    arm_class = DArm

    def __init__(
        self,
        fiveprime_seq,
        threeprime_seq,
        num_allowed_unpaired=1,
        cautious=False):

        super().__init__(
            fiveprime_seq,
            threeprime_seq,
            num_allowed_unpaired=num_allowed_unpaired,
            cautious=cautious)


class FiveprimeDStemSeq(_Sequence):
    name = '5\' D Stem Sequence'
    canonical_positions = ((10, 11, 12), (13, ))
    allowed_input_lengths = tuple(itertools.product((3, ), (0, 1)))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucleotides = ({}, {})
    arm_class = DArm
    stem_class = DStem

    def __init__(
        self,
        positions_10_to_12_string,
        position_13_string='',
        start_index=None,
        stop_index=None,
        cautious=False):

        if cautious:
            if len(positions_10_to_12_string) != 3:
                raise TransferRNAIdentifierError(
                    "Your `positions_10_to_12_string` was not the required 3 bases long: %s"
                    % positions_10_to_12_string)
            if not 0 <= len(position_13_string) <= 1:
                raise TransferRNAIdentifierError(
                    "Your `position_13_string` was not the required 0 or 1 bases long: %s"
                    % position_13_string)

        super().__init__(
            (positions_10_to_12_string,
             position_13_string),
            start_index=start_index,
            stop_index=stop_index,
            cautious=False)


class DLoop(_Loop):
    name = 'D Loop'
    canonical_positions = ((14, 15), (16, 17), (18, 19), (20, ), (21, ))
    allowed_input_lengths = tuple(itertools.product((2, ), (1, 2, 3), (2, ), (1, 2, 3), (1, )))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucleotides = (
        {0: 'A', 1: ('A', 'G')}, {}, {0: 'G', 1: 'G'}, {}, {0: ('A', 'G')})
    arm_class = DArm
    stem_class = DStem

    def __init__(
        self,
        positions_14_to_15_string,
        alpha_positions_string,
        positions_18_to_19_string,
        beta_positions_string,
        position_21_string,
        num_allowed_unconserved=2,
        start_index=None,
        stop_index=None,
        cautious=False):

        if cautious:
            if len(positions_14_to_15_string) != 1:
                raise TransferRNAIdentifierError(
                    "Your `positions_14_to_15_string` was not the required 1 base long: %s"
                    % positions_14_to_15_string)
            if not 1 <= len(alpha_positions_string) <= 3:
                raise TransferRNAIdentifierError(
                    "Your `alpha_positions_string` was not the required 1 to 3 bases long: %s"
                    % alpha_positions_string)
            if len(positions_18_to_19_string) != 2:
                raise TransferRNAIdentifierError(
                    "Your `positions_18_to_19_string` was not the required 2 bases long: %s"
                    % positions_18_to_19_string)
            if not 1 <= len(beta_positions_string) <= 3:
                raise TransferRNAIdentifierError(
                    "Your `beta_positions_string` was not the required 1 to 3 bases long: %s"
                    % beta_positions_string)
            if len(position_21_string) != 1:
                raise TransferRNAIdentifierError(
                    "Your `position_21_string` was not the required 1 base long: %s"
                    % position_21_string)
        alpha_start_index = 1
        alpha_stop_index = alpha_start_index + len(alpha_positions_string)
        self.alpha_seq = _Sequence(
            alpha_positions_string, start_index=alpha_start_index, stop_index=alpha_stop_index)
        beta_start_index = alpha_stop_index + 2
        beta_stop_index = beta_start_index + len(beta_positions_string)
        self.beta_seq = _Sequence(
            beta_positions_string, start_index=beta_start_index, stop_index=beta_stop_index)

        super().__init__(
            (positions_14_to_15_string,
             alpha_positions_string,
             positions_18_to_19_string,
             beta_positions_string,
             position_21_string),
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class ThreeprimeDStemSeq(_Sequence):
    name = '3\' D Stem Sequence'
    canonical_positions = ((22, ), (23, 24, 25))
    allowed_input_lengths = tuple(itertools.product((0, 1), (3, )))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucleotides = ({}, {})
    arm_class = DArm
    stem_class = DStem

    def __init__(
        self,
        position_22_string,
        positions_23_to_25_string='',
        start_index=None,
        stop_index=None,
        cautious=False):

        if cautious:
            if not 0 <= len(position_22_string) <= 1:
                raise TransferRNAIdentifierError(
                    "Your `position_22_string` was not the required 1 base long: %s"
                    % position_22_string)
            if len(positions_23_to_25_string) != 3:
                raise TransferRNAIdentifierError(
                    "Your `positions_23_to_25_string` was not the required 1 to 3 bases long: %s"
                    % positions_23_to_25_string)

        super().__init__(
            (position_22_string, positions_23_to_25_string),
            start_index=start_index,
            stop_index=stop_index,
            cautious=False)


class PositionTwentySix(_Nucleotide):
    name = 'Position 26'
    canonical_positions = ((26, ), )

    def __init__(
        self,
        string,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            string,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class AnticodonArm(_Arm):
    name = 'Anticodon Arm'

    def __init__(
        self,
        stem,
        loop,
        num_allowed_unconserved=1,
        cautious=False):

        self.anticodon = loop.anticodon

        super().__init__(
            stem,
            loop,
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)


class AnticodonStem(_Stem):
    name = 'Anticodon Stem'
    arm_class = AnticodonArm

    def __init__(
        self,
        fiveprime_seq,
        threeprime_seq,
        num_allowed_unpaired=1,
        cautious=False):

        super().__init__(
            fiveprime_seq,
            threeprime_seq,
            num_allowed_unpaired=num_allowed_unpaired,
            cautious=cautious)


class FiveprimeAnticodonStemSeq(_Sequence):
    name = '5\' Anticodon Stem Sequence'
    canonical_positions = ((27, 28, 29, 30, 31), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({}, )
    stem_class = AnticodonStem
    arm_class = AnticodonArm

    def __init__(
        self,
        substrings,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class Anticodon(_Sequence):
    name = 'Anticodon'
    canonical_positions = ((34, 35, 36))
    allowed_input_lengths = ((3, ), )
    summed_input_lengths = (3, )
    conserved_nucleotides = ({}, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(
        self,
        substrings,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)

        try:
            self.aa_string = ANTICODON_TO_AA[self.string]
        except KeyError:
            self.aa_string = 'NA'


class AnticodonLoop(_Loop):
    name = 'Anticodon Loop'
    canonical_positions = ((32, 33, 34, 35, 36, 37, 38), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({1: 'T', 5: ('A', 'G')}, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(
        self,
        substrings,
        num_allowed_unconserved=1,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)

        if len(self.string) != 7:
            raise TransferRNAIdentifierError(
                "Your `AnticodonLoop` was not the required 7 bases long: %s" % self.string)

        self.anticodon = Anticodon(self.string[2: 5])


class ThreeprimeAnticodonStemSeq(_Sequence):
    name = '3\' Anticodon Stem Sequence'
    canonical_positions = ((39, 40, 41, 42, 43), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({}, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(
        self,
        substrings,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class VLoop(_Loop):
    name = 'V Loop'
    canonical_positions = ((44, 45, 46, 47, 48), )
    allowed_input_lengths = tuple((l, ) for l in range(4, 24))
    summed_input_lengths = tuple(range(4, 24))
    conserved_nucleotides = ({}, )

    def __init__(
        self,
        substrings,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)

        if 4 <= len(self.string) <= 5:
            self.type = 'I'
        elif 12 <= len(self.string) <= 23:
            self.type = 'II'
        else:
            self.type = 'NA'


class TArm(_Arm):
    name = 'T Arm'
    def __init__(
        self,
        stem,
        loop,
        num_allowed_unconserved=2,
        cautious=False):

        super().__init__(
            stem,
            loop,
            num_allowed_unconserved=num_allowed_unconserved,
            cautious=cautious)


class TStem(_Stem):
    name = 'T Stem'
    arm_class = TArm

    def __init__(
        self,
        fiveprime_seq,
        threeprime_seq,
        num_allowed_unpaired=1,
        cautious=False):

        super().__init__(
            fiveprime_seq,
            threeprime_seq,
            num_allowed_unpaired=num_allowed_unpaired,
            cautious=cautious)


class FiveprimeTStemSeq(_Sequence):
    name = '5\' T Stem Sequence'
    canonical_positions = ((49, 50, 51, 52, 53), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({4: 'G'}, )
    arm_class = TArm
    stem_class = TStem

    def __init__(
        self,
        substrings,
        num_allowed_unconserved=1,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class TLoop(_Loop):
    name = 'T Loop'
    canonical_positions = ((54, 55, 56, 57, 58, 59, 60), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({0: 'T', 1: 'T', 2: 'C', 3: ('A', 'G'), 4: 'A'}, )
    arm_class = TArm

    def __init__(
        self,
        substrings,
        num_allowed_unconserved=2,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class ThreeprimeTStemSeq(_Sequence):
    name = '3\' T Stem Sequence'
    canonical_positions = ((61, 62, 63, 64, 65), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({0: 'C'}, )
    arm_class = TArm
    stem_class = TStem

    def __init__(
        self,
        substrings,
        num_allowed_unconserved=1,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class ThreeprimeAcceptorStemSeq(_Sequence):
    name = '3\' Acceptor Stem Sequence'
    canonical_positions=((66, 67, 68, 69, 70, 71, 72), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({}, )
    stem_class = AcceptorStem

    def __init__(
        self,
        substrings,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)

        if len(self.string) != 7:
            raise TransferRNAIdentifierError(
                "Your `ThreeprimeAcceptorSeq` was not the required 7 bases long: %s" % self.string)


class Discriminator(_Nucleotide):
    name = 'Discriminator'
    canonical_positions = ((73, ), )

    def __init__(
        self,
        string,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            string,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)


class Acceptor(_Sequence):
    name = 'Acceptor'
    canonical_positions = ((74, 75, 76), )
    allowed_input_lengths = ((3, ), (2, ), (1, )) # a minority of reads end CC or C instead of CCA
    summed_input_lengths = (3, 2, 1)
    conserved_nucleotides = ({0: 'C', 1: 'C', 2: 'A'}, )
    charging_recorded = False

    def __init__(
        self,
        substrings,
        num_allowed_unconserved=0,
        num_extra_threeprime=0,
        start_index=None,
        stop_index=None,
        cautious=False):

        super().__init__(
            substrings,
            num_allowed_unconserved=num_allowed_unconserved,
            start_index=start_index,
            stop_index=stop_index,
            cautious=cautious)

        # A relatively common read "error" is an erroneous base in the acceptor.
        # Errors seem to occur most frequently at the 5'-C or 3'-A of the CCA sequence.
        # Possible acceptor errors are not factored into Acceptor.conserved_nucleotides,
        # so errors are recorded as "unconserved" nucleotides.
        # To allow for an error while num_allowed_unconserved is 0,
        # Acceptor.meets_conserved_thresh is set to True.
        # The offending nucleotide and its position can be found in Acceptor.conserved_status.
        # An unconserved base is only allowed in a full-length acceptor at the 3' end of the read
        # (as in, without extra 3' nucleotides beyond the acceptor)
        # to avoid too much leeway resulting in false positive profiles of shorter sequences.
        if self.num_unconserved == 1:
            if self.num_nucleotides == 3:
                if num_extra_threeprime == 0:
                    self.meets_conserved_thresh = True

        if self.charging_recorded:
            if self.string == 'CC' or self.string == 'C':
                self.charged = False
            else:
                self.charged = True # the loaded amino acid "protects" the A
        else:
            self.charged = None


def set_up_charging_analysis():
    Acceptor.charging_recorded = True


def _get_max_fiveprime_lengths(threeprime_to_fiveprime_feature_classes):
    fiveprime_max_lengths = [threeprime_to_fiveprime_feature_classes[-1].summed_input_lengths[-1]]
    for feature_class in threeprime_to_fiveprime_feature_classes[::-1][1:]:
        if hasattr(feature_class, 'summed_input_lengths'):
            fiveprime_max_lengths.insert(
                0, feature_class.summed_input_lengths[-1] + fiveprime_max_lengths[0])
        else:
            fiveprime_max_lengths.insert(0, fiveprime_max_lengths[0])
    return fiveprime_max_lengths


def profile_wrapper(input_queue, output_queue):
    ''' Iterates a Queue of FASTA records and builds a Queue of tRNA profiles '''
    while True:
        # A record is a tuple of read name and sequence.
        name, read = input_queue.get(True)
        output_queue.put(Profile(read, name))


class Profile:
    fiveprime_to_threeprime_feature_classes = _TransferRNAFeature.list_all_tRNA_features()
    threeprime_to_fiveprime_feature_classes = fiveprime_to_threeprime_feature_classes[::-1]
    stem_formation_triggers = [
        FiveprimeTStemSeq, FiveprimeAnticodonStemSeq, FiveprimeDStemSeq, FiveprimeAcceptorStemSeq]
    arm_formation_triggers = [TStem, AnticodonStem, DStem]
    threeprime_stem_seq_indices = {
        TStem: threeprime_to_fiveprime_feature_classes.index(ThreeprimeTStemSeq),
        AnticodonStem: threeprime_to_fiveprime_feature_classes.index(ThreeprimeAnticodonStemSeq),
        DStem: threeprime_to_fiveprime_feature_classes.index(ThreeprimeDStemSeq),
        AcceptorStem: threeprime_to_fiveprime_feature_classes.index(ThreeprimeAcceptorStemSeq)}
    D_loop_index = threeprime_to_fiveprime_feature_classes.index(DLoop)
    loop_indices = {
        TArm: threeprime_to_fiveprime_feature_classes.index(TLoop),
        AnticodonArm: threeprime_to_fiveprime_feature_classes.index(AnticodonLoop),
        DArm: D_loop_index}
    anticodon_loop_index = threeprime_to_fiveprime_feature_classes.index(AnticodonLoop)
    extrapolation_ineligible_features = [
        ThreeprimeAcceptorStemSeq,
        ThreeprimeTStemSeq,
        VLoop,
        ThreeprimeAnticodonStemSeq,
        ThreeprimeDStemSeq]
    fiveprime_max_lengths = _get_max_fiveprime_lengths(threeprime_to_fiveprime_feature_classes)
    long_read_unprofiled_thresh = 5


    def __init__(self, read, name=''):
        self.read = read
        self.name = name

        (self.profiled_seq,
        self.features,
        self.num_conserved,
        self.num_unconserved,
        self.num_paired,
        self.num_unpaired,
        self.num_extra_threeprime,
        self.num_in_extrapolated_fiveprime_feature,
        self.is_mature,
        self.is_long_read) = self.get_profile(self.read[::-1], '', [])

        self.feature_names = [f.name for f in self.features]

        if self.features:
            # The variable lengths of the alpha and beta subsequences of the D loop
            # can create sequence alignment problems downstream.
            # Therefore, their bounds are tracked explicitly.
            if self.D_loop_index < len(self.features):
                D_loop = self.features[-self.D_loop_index - 1]
                alpha_seq = D_loop.alpha_seq
                beta_seq = D_loop.beta_seq
                self.alpha_start = D_loop.start_index + alpha_seq.start_index
                self.alpha_stop = D_loop.start_index + alpha_seq.stop_index
                self.beta_start = D_loop.start_index + beta_seq.start_index
                self.beta_stop = D_loop.start_index + beta_seq.stop_index
            else:
                self.alpha_start = None
                self.alpha_stop = None
                self.beta_start = None
                self.beta_stop = None

            if self.anticodon_loop_index < len(self.features):
                anticodon = self.features[-self.anticodon_loop_index - 1].anticodon
                self.anticodon_seq = anticodon.string
                self.anticodon_aa = anticodon.aa_string
            else:
                self.anticodon_seq = ''
                self.anticodon_aa = ''

            self.charged = self.features[-1].charged # charge state of the acceptor
        else:
            self.alpha_start = None
            self.alpha_stop = None
            self.beta_start = None
            self.beta_stop = None

            self.anticodon_seq = ''
            self.anticodon_aa = ''

            self.charged = None

        self.is_fully_profiled = (self.read == self.profiled_seq)


    def get_profile(
        self,
        unprofiled_read, # string in 3' to 5' direction
        profiled_read, # string in 5' to 3' direction
        features, # listed in the 5' to 3' direction
        num_conserved=0,
        num_unconserved=0,
        num_paired=0,
        num_unpaired=0,
        num_extra_threeprime=0,
        feature_index=0,
        is_mature=False):

        profile_candidates = []

        if feature_index == len(self.threeprime_to_fiveprime_feature_classes):
            # All tRNA features were profiled, including the tRNA-His 5' G.
            if len(unprofiled_read) > self.long_read_unprofiled_thresh:
                is_long_read = True # read is longer than full-length tRNA
            else:
                is_long_read = False
            return (
                profiled_read, # string in 5' to 3' direction
                features, # listed in the 5' to 3' direction
                num_conserved,
                num_unconserved,
                num_paired,
                num_unpaired,
                num_extra_threeprime,
                0, # number of nucleotides in extrapolated 5' feature
                is_mature,
                is_long_read)
        if not unprofiled_read:
            # The full length of the read was profiled.
            return (
                profiled_read,
                features, # listed in the 5' to 3' direction
                num_conserved,
                num_unconserved,
                num_paired,
                num_unpaired,
                num_extra_threeprime,
                0, # number of nucleotides in extrapolated 5' feature
                is_mature,
                False) # read is not longer than full-length tRNA


        feature_class = self.threeprime_to_fiveprime_feature_classes[feature_index]
        if feature_class in self.stem_formation_triggers: # 3' stem seq triggers stem formation
            # Prepare to form a stem as well as the 3' sequence.
            make_stem = True
            stem_class = feature_class.stem_class
            threeprime_stem_seq = features[-self.threeprime_stem_seq_indices[stem_class] - 1]
            if stem_class in self.arm_formation_triggers:
                # Prepare to form an arm (stem + loop) as well as the stem.
                make_arm = True
                arm_class = stem_class.arm_class
                loop = features[-self.loop_indices[arm_class] - 1]
            else:
                make_arm = False
        else:
            make_stem = False
            make_arm = False
            if feature_class == TransferRNAHisPositionZero:
                # Check for tRNA-His based on the anticodon sequence.
                # tRNA-His uniquely has an extra nucleotide (G) at the 5' end.
                anticodon_string = features[-self.anticodon_loop_index - 1].anticodon.string
                try:
                    aa_string = ANTICODON_TO_AA[anticodon_string]
                except KeyError:
                    aa_string = 'NA'
                if aa_string != 'His':
                    if len(unprofiled_read) > self.long_read_unprofiled_thresh:
                        is_long_read = True # read is longer than full-length tRNA
                    else:
                        is_long_read = False
                    return (
                        profiled_read,
                        features, # listed in the 5' to 3' direction
                        num_conserved,
                        num_unconserved,
                        num_paired,
                        num_unpaired,
                        num_extra_threeprime,
                        0, # number of nucleotides in extrapolated 5' feature
                        is_mature,
                        is_long_read)


        # The list stores the result of a recursive function call finding subsequent 5' features.
        incremental_profile_candidates = []

        # Each primary sequence feature takes (sub)sequence inputs, which can be of varying length.
        # Consider each possible combination of input lengths for the feature,
        # e.g., the D loop contains alpha and beta subsequences of variable length.
        for input_lengths, summed_input_length in zip(
            feature_class.allowed_input_lengths, feature_class.summed_input_lengths):

            # Strands of unequal length cannot form a stem.
            if make_stem:
                # Compare the lengths of the 5' and 3' sequence components.
                if input_lengths != tuple(map(len, threeprime_stem_seq.string_components[::-1])):
                    continue

            # Determine whether there is enough information in the remaining 5' end of the read
            # to assign it to a feature despite the incompleteness of the feature sequence.
            if len(unprofiled_read) < summed_input_length:

                # Features that lack conserved positions
                # or that do not form base pairs with a previously profiled stem sequence
                # do not contain any information for extrapolation of an incomplete sequence.
                if feature_class in self.extrapolation_ineligible_features:
                    continue
                # The unprofiled sequence must be at least 6 nucleotides long
                # to span the 2 conserved positions in the anticodon loop (positions 33 and 37).
                # Conservation of these positions is here considered to be the minimum information
                # needed to identify the anticodon loop.
                elif feature_class == AnticodonLoop:
                    if len(unprofiled_read) < 6:
                        continue

                string_components = [] # feature input substrings
                num_processed_bases = 0
                for input_length in input_lengths[::-1]: # create substrings from 3' to 5'
                    threeprime_to_fiveprime_string = ''
                    for _ in range(input_length):
                        if num_processed_bases >= len(unprofiled_read):
                            threeprime_to_fiveprime_string += 'N' # pad the string at the 5' end
                        else:
                            threeprime_to_fiveprime_string += unprofiled_read[num_processed_bases]
                        num_processed_bases += 1
                    string_components.insert(0, threeprime_to_fiveprime_string[::-1]) # 5' to 3'
                # WITH THE N PADDING, THE START INDEX OF THE FEATURE IN THE READ IS NEGATIVE.
                feature = feature_class(
                    *string_components,
                    start_index=len(self.read) - len(profiled_read) - num_processed_bases,
                    stop_index=len(self.read) - len(profiled_read))
                # The sequence is valid if it doesn't have too many unconserved bases.
                if feature.meets_conserved_thresh:
                    if make_stem:
                        stem = stem_class(feature, threeprime_stem_seq)
                        # The stem is valid if it doesn't have too many unpaired bases.
                        if stem.meets_pair_thresh:
                            if make_arm:
                                arm = arm_class(stem, loop)
                                # The arm is valid if it doesn't have too many unconserved bases.
                                if arm.meets_conserved_thresh:
                                    incremental_profile_candidates.append((
                                        unprofiled_read[::-1], # flip orientation to 5' to 3'
                                        [arm, stem, feature],
                                        feature.num_conserved,
                                        feature.num_unconserved,
                                        stem.num_paired,
                                        stem.num_unpaired,
                                        # number of nucleotides in extrapolated 5' feature
                                        summed_input_length - len(unprofiled_read)))
                                    continue
                            else:
                                incremental_profile_candidates.append((
                                    unprofiled_read[::-1], # flip orientation to 5' to 3'
                                    [stem, feature],
                                    feature.num_conserved,
                                    feature.num_unconserved,
                                    stem.num_paired,
                                    stem.num_unpaired,
                                    # number of nucleotides in extrapolated 5' feature
                                    summed_input_length - len(unprofiled_read)))
                                continue
                    else:
                        incremental_profile_candidates.append((
                            unprofiled_read[::-1], # flip orientation to 5' to 3'
                            [feature],
                            feature.num_conserved,
                            feature.num_unconserved,
                            0, # number of paired nucleotides (not considering a stem)
                            0, # number of unpaired nucleotides (not considering a stem)
                            # number of nucleotides in extrapolated 5' feature
                            summed_input_length - len(unprofiled_read)))
                        continue

            # The procedure for assigning full-length features
            # is similar to the prior precedure for partial-length features
            # with a few efficiencies and special consideration of the acceptor.
            else:
                string_components = [] # feature input substrings
                num_processed_bases = 0
                for input_length in input_lengths[::-1]: # create substrings from 3' to 5'
                    string_components.insert(
                        0,
                        unprofiled_read[num_processed_bases: num_processed_bases + input_length][
                            ::-1]) # flip orientation to 5' to 3'
                    num_processed_bases += input_length
                if feature_class.name == 'Acceptor':
                    feature = feature_class(
                        *string_components,
                        num_extra_threeprime=num_extra_threeprime,
                        start_index=len(self.read) - len(profiled_read) - num_processed_bases,
                        stop_index=len(self.read) - len(profiled_read))
                else:
                    feature = feature_class(
                        *string_components,
                        start_index=len(self.read) - len(profiled_read) - num_processed_bases,
                        stop_index=len(self.read) - len(profiled_read))
                # The sequence is valid if it doesn't have too many unconserved bases.
                if feature.meets_conserved_thresh:
                    if make_stem:
                        stem = stem_class(feature, threeprime_stem_seq)
                        # The stem is valid if it doesn't have too many unpaired bases.
                        if stem.meets_pair_thresh:
                            if make_arm:
                                arm = arm_class(stem, loop)
                                # The are is valid if it doesn't have too many unconserved bases.
                                if arm.meets_conserved_thresh:
                                    incremental_profile_candidates.append((
                                        unprofiled_read[: num_processed_bases][
                                            ::-1], # flip orientation to 5' to 3'
                                        [arm, stem, feature],
                                        feature.num_conserved,
                                        feature.num_unconserved,
                                        stem.num_paired,
                                        stem.num_unpaired,
                                        0)) # number of nucleotides in extrapolated 5' feature
                                    continue
                            else:
                                incremental_profile_candidates.append((
                                    unprofiled_read[: num_processed_bases][
                                        ::-1], # flip orientation to 5' to 3'
                                    [stem, feature],
                                    feature.num_conserved,
                                    feature.num_unconserved,
                                    stem.num_paired,
                                    stem.num_unpaired,
                                    0)) # number of nucleotides in extrapolated 5' feature
                                continue
                    else:
                        incremental_profile_candidates.append((
                            unprofiled_read[: num_processed_bases][
                                ::-1], # flip orientation to 5' to 3'
                            [feature],
                            feature.num_conserved,
                            feature.num_unconserved,
                            0, # number of paired nucleotides (not considering a stem)
                            0, # number of unpaired nucleotides (not considering a stem)
                            0)) # number of nucleotides in extrapolated 5' feature

                        # AVOID TESTING PARTIAL ACCEPTOR SEQUENCES IF THE FULL ONE WAS JUST FOUND.
                        if feature.name == 'Acceptor':
                            if feature.string == 'CCA':
                                break

                        continue

                # Only consider extra 3' nucleotides and a full acceptor, CCA,
                # not a partial acceptor, CC or C.
                elif feature.name == 'Acceptor':
                    if summed_input_length == 3:
                        if num_extra_threeprime > 0:
                            break

            # AVOID TESTING PARTIAL ACCEPTOR SEQUENCES IF EXTRA 3' NUCLEOTIDES ARE CONSIDERED.
            if feature_class.name == 'Acceptor':
                if num_extra_threeprime > 0:
                    break


        if not incremental_profile_candidates: # the feature didn't pass muster
            # Try adding an extra 3' "nontemplated" nucleotide, up to 2 extra.
            if feature_class.name == 'Acceptor':
                if num_extra_threeprime < 2:
                    # Since the first feature, the acceptor, was not found,
                    # return the result from adding extra nucleotides,
                    # rather than comparing to the current (nonexistent) result.
                    return self.get_profile(
                        unprofiled_read[1: ],
                        unprofiled_read[0] + profiled_read, # 5' to 3' orientation
                        [], # extra 3' nucleotides are not counted as a feature
                        num_extra_threeprime=num_extra_threeprime + 1)

            if len(unprofiled_read) > (
                self.fiveprime_max_lengths[feature_index] + self.long_read_unprofiled_thresh):
                is_long_read = True
            else:
                is_long_read = False
            return (
                profiled_read,
                features,
                num_conserved,
                num_unconserved,
                num_paired,
                num_unpaired,
                num_extra_threeprime,
                0, # number of nucleotides in extrapolated 5' feature
                is_mature,
                is_long_read)


        # Sort candidates by
        # 1. number of features identified (at most, sequence + stem + arm) (descending),
        # 2. number of unconserved + unpaired nucleotides (ascending),
        # 3. incompleteness of the last (most 5') feature (ascending).
        # This sort also happens later for full sequence profiles,
        # but this first sort is useful for seeking out and returning "flawless" mature tRNA.
        incremental_profile_candidates.sort(key=lambda p: (-len(p[1]), p[3] + p[5], p[6]))
        # Continue finding features in reads that have not been fully profiled --
        # do not recurse profile candidates in which the final feature was extrapolated.
        for ipc in incremental_profile_candidates:
            if ipc[6] == 0: # 5' feature was NOT extrapolated from unprofiled read
                if is_mature or feature_class.name == '5\' Acceptor Stem Sequence':
                    profile_candidate = self.get_profile( # RECURSE
                        unprofiled_read[len(ipc[0]): ],
                        ipc[0] + profiled_read,
                        ipc[1] + features,
                        num_conserved=ipc[2] + num_conserved,
                        num_unconserved=ipc[3] + num_unconserved,
                        num_paired=ipc[4] + num_paired,
                        num_unpaired=ipc[5] + num_unpaired,
                        num_extra_threeprime=num_extra_threeprime,
                        feature_index=feature_index + len(ipc[1]),
                        is_mature=True)
                else:
                    profile_candidate = self.get_profile( # RECURSE
                        unprofiled_read[len(ipc[0]): ],
                        ipc[0] + profiled_read,
                        ipc[1] + features,
                        num_conserved=ipc[2] + num_conserved,
                        num_unconserved=ipc[3] + num_unconserved,
                        num_paired=ipc[4] + num_paired,
                        num_unpaired=ipc[5] + num_unpaired,
                        num_extra_threeprime=num_extra_threeprime,
                        feature_index=feature_index + len(ipc[1]),
                        is_mature=False)
                if (profile_candidate[8] # is mature (full-length tRNA)
                    and profile_candidate[3] == 0 # no unconserved
                    and profile_candidate[5] == 0): # no unpaired
                    return profile_candidate
                else:
                    profile_candidates.append(profile_candidate)
            else: # 5' feature was extrapolated from unprofiled read
                # Extrapolated features are grounded
                # in conserved nucleotides and/or paired nucleotides in stems.
                profile_candidates.append((
                    ipc[0] + profiled_read,
                    ipc[1] + features,
                    ipc[2] + num_conserved,
                    ipc[3] + num_unconserved,
                    ipc[4] + num_paired,
                    ipc[5] + num_unpaired,
                    num_extra_threeprime,
                    ipc[6], # number of nucleotides in extrapolated 5' feature
                    False, # not a full-length mature tRNA
                    False)) # read is not longer than full-length tRNA

        # Try adding an extra 3' "nontemplated" nucleotide, up to 2 extra.
        if feature_class.name == 'Acceptor':
            if num_extra_threeprime < 2:
                profile_candidates.append(self.get_profile(
                    unprofiled_read[1: ],
                    unprofiled_read[0] + profiled_read, # 5' to 3' orientation
                    [], # extra 3' nucleotides are not counted as a feature
                    num_extra_threeprime=num_extra_threeprime + 1))

        # Do not add 5' features to profiled 3' features if the additional features
        # have no grounding in conserved nucleotides or paired nucleotides in stems.
        profile_candidates = [
            p for p in profile_candidates if p[2] + p[4] > num_conserved + num_paired]
        if profile_candidates:
            profile_candidates.sort(key=lambda p: (-len(p[1]), p[3] + p[5], p[7]))
            return profile_candidates[0]
        else:
            if len(unprofiled_read) > (
                self.fiveprime_max_lengths[feature_index] + self.long_read_unprofiled_thresh):
                is_long_read = True
            else:
                is_long_read = False
            return (
                profiled_read,
                features,
                num_conserved,
                num_unconserved,
                num_paired,
                num_unpaired,
                num_extra_threeprime,
                0, # number of nucleotides in extrapolated 5' feature
                is_mature,
                is_long_read)


    def get_unconserved_positions(self):
        unconserved_info = []
        for feature in self.features:
            # Only _Nucleotide and _Sequence subclasses have the attribute.
            if hasattr(feature, 'conserved_status'):
                component_start_index = feature.start_index
                # Conserved nucleotides are indexed within the string "component" (substring).
                for string_component_statuses, string_component in zip(
                    feature.conserved_status, feature.string_components):
                    for (
                        nucleotide_index, is_conserved, observed_nucleotide, expected_nucleotides
                        ) in string_component_statuses:
                        # Avoid N padding in an extrapolated 5' feature.
                        if not is_conserved and observed_nucleotide != 'N':
                            unconserved_info.append(
                                (component_start_index + nucleotide_index,
                                 observed_nucleotide,
                                 ','.join(expected_nucleotides)))
                    component_start_index += len(string_component)
        return unconserved_info


    def get_unpaired_positions(self):
        unpaired_info = []
        for feature in self.features:
            # Only the _Stem class has the attribute.
            if hasattr(feature, 'paired_status'):
                for nucleotide_index, (
                    is_paired, fiveprime_nucleotide, threeprime_nucleotide
                    ) in enumerate(feature.paired_status):
                    # Avoid N padding in an extrapolated 5' feature.
                    if not is_paired and fiveprime_nucleotide != 'N':
                        unpaired_info.append((
                            feature.fiveprime_seq.start_index + nucleotide_index,
                            feature.threeprime_seq.stop_index - nucleotide_index - 1,
                            fiveprime_nucleotide,
                            threeprime_nucleotide))
        return unpaired_info

##### TEST CASES #####

# E. coli tRNA-Ala-GGC-1-1
# forward = 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Gln-CTG-1-1
# forward = 'TGGGGTATCGCCAAGCGGTAAGGCACCGGATTCTGATTCCGGCATTCCGAGGTTCGAATCCTCGTACCCCAGCCA'
# E. coli tRNA-Glu-TTC-1-1
# forward = 'GTCCCCTTCGTCTAGAGGCCCAGGACACCGCCCTTTCACGGCGGTAACAGGGGTTCGAATCCCCTAGGGGACGCCA'
# E. coli tRNA-Leu-CAG-1-1
# forward = 'GCGAAGGTGGCGGAATTGGTAGACGCGCTAGCTTCAGGTGTTAGTGTCCTTACGGACGTGGGGGTTCAAGTCCCCCCCCTCGCACCA'
# E. coli tRNA-Leu-TAA-1-1
# forward = 'GCCCGGATGGTGGAATCGGTAGACACAAGGGATTTAAAATCCCTCGGCGTTCGCGCTGTGCGGGTTCAAGTCCCGCTCCGGGTACCA'
# E. coli tRNA-Ile-GAT-1-1
# forward = 'AGGCTTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCAAGTCCACTCAGGCCTACCA'
# E. coli tRNA-Leu-CAG-2-1
# forward = 'GCGAAGGTGGCGGAATTGGTAGACGCGCTAGCTTCAGGTGTTAGTGTTCTTACGGACGTGGGGGTTCAAGTCCCCCCCCTCGCACCA'
# E. coli tRNA-fMet-CAT-1-1
# forward = 'CGCGGGGTGGAGCAGCCTGGTAGCTCGTCGGGCTCATAACCCGAAGGTCGTCGGTTCAAATCCGGCCCCCGCAACCA'
# E. coli tRNA-Leu-TAA-1-1
# forward = 'GCCCGGATGGTGGAATCGGTAGACACAAGGGATTTAAAATCCCTCGGCGTTCGCGCTGTGCGGGTTCAAGTCCCGCTCCGGGTACCA'
# E. coli tRNA-His-GTG-1-1: includes the 5' G
# forward = 'GGTGGCTATAGCTCAGTTGGTAGAGCCCTGGATTGTGATTCCAGTTGTCGTGGGTTCGAATCCCATTAGCCACCCCA'
# H. sapiens tRNA-SeC-TCA-1-1
# forward = 'GCCCGGATGATCCTCAGTGGTCTGGGGTGCAGGCTTCAAACCTGTAGCTGTCTAGCGACAGAGTGGTTCAATTCCACCTTTCGGGCGCCA'
# C. elegans tRNA-SeC-TCA-1-1
# forward = 'GCCCGGATGAACCATGGCGGTCTGTGGTGCAGACTTCAAATCTGTAGGCGGTTAGCGCCGCAGTGGTTCGACTCCACCTTTCGGGTGCCA'
# E. coli tRNA-SeC-TCA-1-1
# forward = 'GGAAGATCGTCGTCTCCGGTGAGGCGGCTGGACTTCAAATCCAGTTGGGGCCGCCAGCGGTCCCGGGCAGGTTCGACTCCTGTGATCTTCCGCCA'
# A. fulgidus DSM 4304 tRNA-Glu-TTC-1-1: the gene has introns, so the mature tRNA may be weird
# forward = 'GCUCCGGUGGUGUAGCCCGGCCAAUCAUUCCGGCCUUUCGAGCCGGCGACCCGGGUUCAAAUCCCGGCCGGAGCACCA'.replace('U', 'T')

# Truncated
# E. coli tRNA-Ala-GGC-1-1: start at second nucleotide
# forward = 'GGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start at second position of 5' strand of acceptor stem
# forward = 'ATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start at second position of 5' strand of D stem
# forward = 'CTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start at second position of 5' strand of D stem, "mutate" "pos 37" (index 26) A -> T
# forward = 'CTCAGCTGGGAGAGCGCTTGCATGGCTTGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start at second position of 5' strand of D stem, "mutate" "pos 37" (index 26) A -> T, "mutate" "pos 18" (index 7) G -> T
# forward = 'CTCAGCTTGGAGAGCGCTTGCATGGCTTGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start at second position of 5' strand of D stem, "mutate" "pos 31" (index 20) C -> A, which "unpairs" with "pos 39" (index 28) G
# forward = 'CTCAGCTGGGAGAGCGCTTGAATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start at first position of 3' strand of anticodon stem
# forward = 'GCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'
# E. coli tRNA-Ala-GGC-1-1: start first position of 3' strand of D stem
# forward = 'GAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCACCA'

# profile = Profile(forward)
# print(profile.profiled_seq)
# print(profile.features)
# print(profile.num_unconserved)
# print(profile.num_unpaired)
# print(profile.num_in_extrapolated_fiveprime_feature)
# print(profile.is_mature)
# print(profile.is_fully_profiled)
# print(profile.get_unconserved_positions())
# print(profile.get_unpaired_positions())
# print(profile.alpha_start)
# print(profile.beta_start)
