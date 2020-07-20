# -*- coding: utf-8
# pylint: disable=line-too-long
"""tRNA identification from input sequence"""

import itertools

from anvio.constants import WC_plus_wobble_base_pairs as WC_PLUS_WOBBLE_BASE_PAIRS
from anvio.constants import anticodon_to_AA as ANTICODON_TO_AA
from anvio.constants_package.trnaseq import LONGEST_KNOWN_TRNA_LENGTH
from anvio.errors import TRNAIdentifierError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


THREEPRIME_VARIANTS = ['CCA',
                       'C', 'CC',
                       'CCC', 'CCG', 'CCT',
                       'CAA', 'CGA', 'CTA',
                       'ACA', 'GCA', 'TCA',
                       'CCAA', 'CCAC', 'CCAG', 'CCAT',
                       'CCAAA', 'CCAAC', 'CCAAG', 'CCAAT',
                       'CCACA', 'CCACC', 'CCACG', 'CCACT',
                       'CCAGA', 'CCAGC', 'CCAGG', 'CCAGT',
                       'CCATA', 'CCATC', 'CCATG', 'CCATT']


class _TRNAFeature:
    conserved_nucleotides = ({}, )

    def __init__(self,
                 string_components, # ex. ('CCA', )
                 num_allowed_unconserved=-1, # ex. 0
                 cautious=False):

        if cautious:
            if type(string_components) != tuple:
                raise TRNAIdentifierError("`string_components` must be in the form of a tuple, e.g., ('ACTGG', 'CCAGT'). "
                                          "Your `string_components` were %s" % (string_components, ))
        self.string_components = string_components

        self.num_nucleotides = sum(map(len, string_components))

        # By default, base conservation is not enforced.
        if num_allowed_unconserved == -1:
            self.num_allowed_unconserved = sum(len(d) for d in self.conserved_nucleotides)
        else:
            self.num_allowed_unconserved = num_allowed_unconserved


    def check_conserved_nucleotides(self):
        num_conserved = 0
        num_unconserved = 0 # can include N "padding" in extrapolated 5' feature
        conserved_status = []
        for substring, nuc_dict in zip(self.string_components, self.conserved_nucleotides):
            substring_statuses = []
            conserved_status.append(substring_statuses)
            for pos, expected_nucleotides in nuc_dict.items():
                try:
                    observed_nucleotide = substring[pos]
                except IndexError:
                    # This occurs for an Acceptor feature of CC or C rather than the canonical CCA.
                    break
                if observed_nucleotide in expected_nucleotides:
                    num_conserved += 1
                    substring_statuses.append((pos, True, observed_nucleotide, expected_nucleotides))
                else:
                    num_unconserved += 1
                    substring_statuses.append((pos, False, observed_nucleotide, expected_nucleotides))
        if num_unconserved > self.num_allowed_unconserved:
            meets_conserved_thresh = False
        else:
            meets_conserved_thresh = True
        return (meets_conserved_thresh, num_conserved, num_unconserved, conserved_status)


    @staticmethod
    def list_all_tRNA_features():
        return [TRNAHisPositionZero,
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


class _Nucleotide(_TRNAFeature):
    allowed_input_lengths = ((1, ), )
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))

    def __init__(self,
                 string, # must be a string of length 1
                 num_allowed_unconserved=-1,
                 start_index=None,
                 stop_index=None,
                 cautious=False):

        self.string = string
        self.start_index = start_index
        self.stop_index = stop_index

        super().__init__((string, ), num_allowed_unconserved=num_allowed_unconserved, cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved,
         self.conserved_status) = self.check_conserved_nucleotides()


class _Sequence(_TRNAFeature):
    def __init__(self,
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

        super().__init__(string_components, num_allowed_unconserved=num_allowed_unconserved, cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved, # can include N "padding" in extrapolated 5' feature
         self.conserved_status) = self.check_conserved_nucleotides()


class _Loop(_Sequence):
    def __init__(self,
                 substrings, # must be a string, tuple of strings, or tuple of _Nucleotide/_Sequence objects
                 num_allowed_unconserved=-1,
                 start_index=None,
                 stop_index=None,
                 cautious=False):

        super().__init__(substrings,
                         num_allowed_unconserved=num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class _Stem(_TRNAFeature):
    def __init__(self, fiveprime_seq, threeprime_seq, num_allowed_unpaired=0, num_allowed_unconserved=-1, cautious=False):

        if cautious:
            if type(fiveprime_seq) != _Sequence or type(threeprime_seq) != _Sequence:
                raise TRNAIdentifierError("You can only define a _Stem from _Sequence objects.")
        self.fiveprime_seq = fiveprime_seq
        self.threeprime_seq = threeprime_seq

        self.canonical_positions=(*self.fiveprime_seq.canonical_positions, *self.threeprime_seq.canonical_positions)
        self.conserved_nucleotides = (*self.fiveprime_seq.conserved_nucleotides, *self.threeprime_seq.conserved_nucleotides)

        if cautious:
            if tuple(map(len, self.fiveprime_seq.string_components)) != tuple(map(len, self.threeprime_seq.string_components[::-1])):
                raise TRNAIdentifierError("The two _Sequence objects, %s and %s, "
                                          "that were used to define your _Stem are not the same length."
                                          % (self.fiveprime_seq.string_components, threeprime_seq.string_components))
            length = sum(map(len, self.fiveprime_seq.string_components))
            if num_allowed_unpaired > length:
                raise TRNAIdentifierError("You tried to leave at most %d base pairs unpaired, "
                                          "but there are only %d base pairs in the stem." % (num_allowed_unpaired, length))
        self.num_allowed_unpaired = num_allowed_unpaired

        self.start_indices = (self.fiveprime_seq.start_index, self.threeprime_seq.start_index)
        self.stop_indices = (self.fiveprime_seq.stop_index, self.threeprime_seq.stop_index)

        super().__init__((*self.fiveprime_seq.string_components, *self.threeprime_seq.string_components),
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_pair_thresh,
         self.num_paired,
         self.num_unpaired, # can include N "padding" in extrapolated 5' feature
         self.paired_status) = self.check_pairs()


    def check_pairs(self):
        num_paired = 0
        num_unpaired = 0 # can include N "padding" in extrapolated 5' feature
        paired_status = []
        for fiveprime_nuc, threeprime_nuc in zip(self.fiveprime_seq.string, self.threeprime_seq.string[::-1]):
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
        return meets_pair_thresh, num_paired, num_unpaired, paired_status


class _Arm(_TRNAFeature):
    def __init__(self,
                 stem, # must be _Stem object
                 loop, # must be _Loop object
                 num_allowed_unconserved=-1,
                 cautious=False):

        if cautious:
            if type(stem) != _Stem or type(loop) != _Loop:
                raise TRNAIdentifierError("A `_Stem` and a `_Loop` are required input to create an `_Arm`.")
        self.stem = stem
        self.loop = loop

        self.canonical_positions = (*stem.fiveprime_seq.canonical_positions,
                                    *loop.canonical_positions,
                                    *stem.threeprime_seq.canonical_positions)
        if cautious:
            if (tuple(pos for component_positions in self.canonical_positions for pos in component_positions)
                != tuple(range(self.canonical_positions[0][0], (self.canonical_positions[0][0]
                                                                + len(self.canonical_positions[0])
                                                                + len(self.canonical_positions[1])
                                                                + len(self.canonical_positions[2]))))):
                raise TRNAIdentifierError("The canonical positions in an `_Arm` must be contiguous. "
                                          "These were yours: %s. "
                                          "These came from the canonical positions in _Stem, %s, "
                                          "and the canonical positions in _Loop, %s."
                                          % (self.canonical_positions,
                                             stem.canonical_positions,
                                             loop.canonical_positions))

        self.conserved_nucleotides=(*stem.fiveprime_seq.conserved_nucleotides,
                                    *loop.conserved_nucleotides,
                                    *stem.threeprime_seq.conserved_nucleotides)

        self.start_index = self.stem.start_indices[0]
        self.stop_index = self.stem.stop_indices[1]

        super().__init__((*stem.fiveprime_seq.string_components,
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


class TRNAHisPositionZero(_Nucleotide):
    name = 'tRNA-His Position 0'
    canonical_positions = ((-1, ), )
    conserved_nucleotides = ({0: 'G'}, )

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        super().__init__(string, num_allowed_unconserved=0, start_index=start_index, stop_index=stop_index, cautious=cautious)


class AcceptorStem(_Stem):
    # For our purposes, the acceptor stem only includes the base-paired nucleotides in the stem.
    name = 'Acceptor Stem'

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq, threeprime_seq, num_allowed_unpaired=1, cautious=cautious)


class FiveprimeAcceptorStemSeq(_Sequence):
    name = '5\' Acceptor Stem Sequence'
    canonical_positions = ((1, 2, 3, 4, 5, 6, 7), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({}, )
    stem_class = AcceptorStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)

        if len(self.string) != 7:
            raise TRNAIdentifierError("Your `FiveprimeAcceptorSeq` was not the required 7 bases long: %s" % self.string)


class PositionEight(_Nucleotide):
    name = 'Position 8'
    canonical_positions = ((8, ), )
    conserved_nucleotides = ({0: 'T'}, )

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        super().__init__(string, num_allowed_unconserved=1, start_index=start_index, stop_index=stop_index, cautious=cautious)


class PositionNine(_Nucleotide):
    name = 'Position 9'
    canonical_positions = ((9, ), )

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        super().__init__(string, start_index=start_index, stop_index=stop_index, cautious=cautious)


class DArm(_Arm):
    name = 'D Arm'

    def __init__(self, stem, loop, cautious=False):
        super().__init__(stem, loop, num_allowed_unconserved=2, cautious=cautious)


class DStem(_Stem):
    name = 'D Stem'
    arm_class = DArm

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq, threeprime_seq, num_allowed_unpaired=1, cautious=cautious)


class FiveprimeDStemSeq(_Sequence):
    name = '5\' D Stem Sequence'
    canonical_positions = ((10, 11, 12), (13, ))
    allowed_input_lengths = tuple(itertools.product((3, ), (0, 1)))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucleotides = ({}, {})
    arm_class = DArm
    stem_class = DStem

    def __init__(self, positions_10_to_12_string, position_13_string='', start_index=None, stop_index=None, cautious=False):

        if cautious:
            if len(positions_10_to_12_string) != 3:
                raise TRNAIdentifierError("Your `positions_10_to_12_string` was not the required 3 bases long: %s"
                                          % positions_10_to_12_string)
            if not 0 <= len(position_13_string) <= 1:
                raise TRNAIdentifierError("Your `position_13_string` was not the required 0 or 1 bases long: %s"
                                          % position_13_string)

        super().__init__((positions_10_to_12_string, position_13_string),
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=False)


class DLoop(_Loop):
    name = 'D Loop'
    canonical_positions = ((14, 15), (16, 17), (18, 19), (20, ), (21, ))
    allowed_input_lengths = tuple(itertools.product((2, ), (1, 2, 3), (2, ), (1, 2, 3), (1, )))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucleotides = ({0: 'A', 1: ('A', 'G')}, {}, {0: 'G', 1: 'G'}, {}, {0: ('A', 'G')})
    arm_class = DArm
    stem_class = DStem

    def __init__(self,
                 positions_14_to_15_string,
                 alpha_positions_string,
                 positions_18_to_19_string,
                 beta_positions_string,
                 position_21_string,
                 start_index=None,
                 stop_index=None,
                 cautious=False):

        if cautious:
            if len(positions_14_to_15_string) != 1:
                raise TRNAIdentifierError("Your `positions_14_to_15_string` was not the required 1 base long: %s" % positions_14_to_15_string)
            if not 1 <= len(alpha_positions_string) <= 3:
                raise TRNAIdentifierError("Your `alpha_positions_string` was not the required 1 to 3 bases long: %s" % alpha_positions_string)
            if len(positions_18_to_19_string) != 2:
                raise TRNAIdentifierError("Your `positions_18_to_19_string` was not the required 2 bases long: %s" % positions_18_to_19_string)
            if not 1 <= len(beta_positions_string) <= 3:
                raise TRNAIdentifierError("Your `beta_positions_string` was not the required 1 to 3 bases long: %s" % beta_positions_string)
            if len(position_21_string) != 1:
                raise TRNAIdentifierError("Your `position_21_string` was not the required 1 base long: %s" % position_21_string)

        alpha_start_index = 1
        alpha_stop_index = alpha_start_index + len(alpha_positions_string)
        self.alpha_seq = _Sequence(alpha_positions_string, start_index=alpha_start_index, stop_index=alpha_stop_index)
        beta_start_index = alpha_stop_index + 2
        beta_stop_index = beta_start_index + len(beta_positions_string)
        self.beta_seq = _Sequence(beta_positions_string, start_index=beta_start_index, stop_index=beta_stop_index)

        super().__init__((positions_14_to_15_string,
                          alpha_positions_string,
                          positions_18_to_19_string,
                          beta_positions_string,
                          position_21_string),
                         num_allowed_unconserved=2,
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

    def __init__(self, position_22_string, positions_23_to_25_string='', start_index=None, stop_index=None, cautious=False):

        if cautious:
            if not 0 <= len(position_22_string) <= 1:
                raise TRNAIdentifierError("Your `position_22_string` was not the required 1 base long: %s" % position_22_string)
            if len(positions_23_to_25_string) != 3:
                raise TRNAIdentifierError("Your `positions_23_to_25_string` was not the required 1 to 3 bases long: %s" % positions_23_to_25_string)

        super().__init__((position_22_string, positions_23_to_25_string),
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=False)


class PositionTwentySix(_Nucleotide):
    name = 'Position 26'
    canonical_positions = ((26, ), )

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        super().__init__(string, start_index=start_index, stop_index=stop_index, cautious=cautious)


class AnticodonArm(_Arm):
    name = 'Anticodon Arm'

    def __init__(self, stem, loop, cautious=False):
        self.anticodon = loop.anticodon
        super().__init__(stem, loop, num_allowed_unconserved=1, cautious=cautious)


class AnticodonStem(_Stem):
    name = 'Anticodon Stem'
    arm_class = AnticodonArm

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq, threeprime_seq, num_allowed_unpaired=1, cautious=cautious)


class FiveprimeAnticodonStemSeq(_Sequence):
    name = '5\' Anticodon Stem Sequence'
    canonical_positions = ((27, 28, 29, 30, 31), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({}, )
    stem_class = AnticodonStem
    arm_class = AnticodonArm

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)


class Anticodon(_Sequence):
    name = 'Anticodon'
    canonical_positions = ((34, 35, 36))
    allowed_input_lengths = ((3, ), )
    summed_input_lengths = (3, )
    conserved_nucleotides = ({}, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)
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

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, num_allowed_unconserved=1, start_index=start_index, stop_index=stop_index, cautious=cautious)
        if len(self.string) != 7:
            raise TRNAIdentifierError("Your `AnticodonLoop` was not the required 7 bases long: %s" % self.string)
        self.anticodon = Anticodon(self.string[2: 5])


class ThreeprimeAnticodonStemSeq(_Sequence):
    name = '3\' Anticodon Stem Sequence'
    canonical_positions = ((39, 40, 41, 42, 43), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({}, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)


class VLoop(_Loop):
    name = 'V Loop'
    canonical_positions = ((44, 45, 46, 47, 48), )
    allowed_input_lengths = tuple((l, ) for l in range(4, 24))
    summed_input_lengths = tuple(range(4, 24))
    conserved_nucleotides = ({}, )

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)
        if 4 <= len(self.string) <= 5:
            self.type = 'I'
        elif 12 <= len(self.string) <= 23:
            self.type = 'II'
        else:
            self.type = 'NA'


class TArm(_Arm):
    name = 'T Arm'
    def __init__(self, stem, loop, cautious=False):
        super().__init__(stem, loop, num_allowed_unconserved=2, cautious=cautious)


class TStem(_Stem):
    name = 'T Stem'
    arm_class = TArm

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq, threeprime_seq, num_allowed_unpaired=1, cautious=cautious)


class FiveprimeTStemSeq(_Sequence):
    name = '5\' T Stem Sequence'
    canonical_positions = ((49, 50, 51, 52, 53), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({4: 'G'}, )
    arm_class = TArm
    stem_class = TStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, num_allowed_unconserved=1, start_index=start_index, stop_index=stop_index, cautious=cautious)


class TLoop(_Loop):
    name = 'T Loop'
    canonical_positions = ((54, 55, 56, 57, 58, 59, 60), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({0: 'T', 1: 'T', 2: 'C', 3: ('A', 'G'), 4: 'A'}, )
    arm_class = TArm

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, num_allowed_unconserved=2, start_index=start_index, stop_index=stop_index, cautious=cautious)


class ThreeprimeTStemSeq(_Sequence):
    name = '3\' T Stem Sequence'
    canonical_positions = ((61, 62, 63, 64, 65), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucleotides = ({0: 'C'}, )
    arm_class = TArm
    stem_class = TStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, num_allowed_unconserved=1, start_index=start_index, stop_index=stop_index, cautious=cautious)


class ThreeprimeAcceptorStemSeq(_Sequence):
    name = '3\' Acceptor Stem Sequence'
    canonical_positions=((66, 67, 68, 69, 70, 71, 72), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucleotides = ({}, )
    stem_class = AcceptorStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)
        if len(self.string) != 7:
            raise TRNAIdentifierError("Your `ThreeprimeAcceptorSeq` was not the required 7 bases long: %s" % self.string)


class Discriminator(_Nucleotide):
    name = 'Discriminator'
    canonical_positions = ((73, ), )

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        super().__init__(string, start_index=start_index, stop_index=stop_index, cautious=cautious)


class Acceptor(_Sequence):
    name = 'Acceptor'
    canonical_positions = ((74, 75, 76), )
    allowed_input_lengths = ((3, ), (2, ), (1, )) # unprocessed tRNA can end in 3'-C or CC instead of CCA
    summed_input_lengths = (3, 2, 1)
    conserved_nucleotides = ({0: 'C', 1: 'C', 2: 'A'}, )

    def __init__(self, substrings, num_extra_threeprime=0, start_index=None, stop_index=None, cautious=False):

        # An unprocessed tRNA can have alternate 3' endings to CCA.
        # We consider what we observe to be the most common: C, CC, CCN, CNA, NCA, CCAN, and CCANN.
        # For CCAN and CCANN, num_extra_threeprime must be set to 1 and 2, respectively.
        # The N and NN of CCAN and CCANN are not treated explicitly by this class.
        # C and CC are accommodated by allowed_input_lengths being 1, 2, or 3.
        # Missing nucleotides in C (CA) and CC (A) are not recorded as unconserved nucleotides,
        # but variant nucleotides in CCN, CNA, and NCA are recorded as such.
        # Although num_allowed_unconserved is 0,
        # CCN, CNA, and NCA cause Acceptor.meets_conserved_thresh to be set to True.

        super().__init__(substrings, num_allowed_unconserved=0, start_index=start_index, stop_index=stop_index, cautious=cautious)
        if self.num_unconserved == 1:
            if self.num_nucleotides == 3:
                if num_extra_threeprime == 0:
                    self.meets_conserved_thresh = True


def profile_wrapper(input_queue, output_queue):
    ''' Iterates a Queue of FASTA records and builds a Queue of tRNA profiles '''
    while True:
        # A record is a tuple of input sequence name and sequence.
        name, input_seq = input_queue.get(True)
        output_queue.put(Profile(input_seq, name))


class Profile:
    FIVEPRIME_TO_THREEPRIME_FEATURE_CLASSES = _TRNAFeature.list_all_tRNA_features()
    THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES = FIVEPRIME_TO_THREEPRIME_FEATURE_CLASSES[::-1]
    STEM_FORMATION_TRIGGERS = [FiveprimeTStemSeq, FiveprimeAnticodonStemSeq, FiveprimeDStemSeq, FiveprimeAcceptorStemSeq]
    ARM_FORMATION_TRIGGERS = [TStem, AnticodonStem, DStem]
    T_LOOP_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(TLoop)
    T_ARM_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(TArm)
    THREEPRIME_STEM_SEQ_INDICES = {TStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(ThreeprimeTStemSeq),
                                   AnticodonStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(ThreeprimeAnticodonStemSeq),
                                   DStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(ThreeprimeDStemSeq),
                                   AcceptorStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(ThreeprimeAcceptorStemSeq)}
    D_LOOP_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(DLoop)
    ARM_LOOP_INDEX_DICT = {TArm: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(TLoop),
                           AnticodonArm: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(AnticodonLoop),
                           DArm: D_LOOP_INDEX}
    ANTICODON_LOOP_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(AnticodonLoop)
    EXTRAPOLATION_INELIGIBLE_FEATURES = [ThreeprimeAcceptorStemSeq,
                                         ThreeprimeTStemSeq,
                                         VLoop,
                                         ThreeprimeAnticodonStemSeq,
                                         ThreeprimeDStemSeq]


    def __init__(self, input_seq, name=''):
        # The input sequence is treated like a tRNA-seq read starting from the 3' end of a tRNA molecule.
        self.input_seq = input_seq
        self.name = name

        profile_dict = self.get_profile(unprofiled_seq=self.input_seq[::-1])
        self.profiled_seq = profile_dict['profiled_seq']
        self.features = profile_dict['features']
        self.num_conserved = profile_dict['num_conserved']
        self.num_unconserved = profile_dict['num_unconserved']
        self.num_paired = profile_dict['num_paired']
        self.num_unpaired = profile_dict['num_unpaired']
        self.num_extra_threeprime = profile_dict['num_extra_threeprime']
        self.num_in_extrapolated_fiveprime_feature = profile_dict['num_extrapolated']
        self.has_complete_feature_set = profile_dict['is_complete']
        self.num_extra_fiveprime = profile_dict['num_extra_fiveprime']

        self.feature_names = [f.name for f in self.features]

        if self.features:
            # The extra 3' nucleotides are not explicitly added to the Acceptor object for complex reasons.
            if self.num_extra_threeprime > 0:
                self.acceptor_variant_string = self.features[-1].string + self.profiled_seq[-self.num_extra_threeprime: ]
            else:
                self.acceptor_variant_string = self.features[-1].string

            # Explicitly record the start and stop positions within the input seq
            # of the variable-length alpha and beta regions of the D loop.
            if self.D_LOOP_INDEX < len(self.features):
                D_loop = self.features[-self.D_LOOP_INDEX - 1]
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

            if self.ANTICODON_LOOP_INDEX < len(self.features):
                anticodon = self.features[-self.ANTICODON_LOOP_INDEX - 1].anticodon
                self.anticodon_seq = anticodon.string
                self.anticodon_aa = anticodon.aa_string
            else:
                self.anticodon_seq = ''
                self.anticodon_aa = ''

            # Use a slightly less stringent standard for tRNA
            # when the acceptor sequence is CCA than C, CC, NCA, CNA, CCN, CCAN, or CCANN.
            if self.acceptor_variant_string == 'CCA':
                # The features must include the T loop.
                if len(self.features) > self.T_LOOP_INDEX:
                    self.is_predicted_trna = True
                else:
                    self.is_predicted_trna = False
            else:
                # The features must include the T arm.
                if len(self.features) > self.T_ARM_INDEX:
                    self.is_predicted_trna = True
                else:
                    self.is_predicted_trna = False
            # Furthermore, if the sequence does not have a complete feature set but is longer than tRNA should be,
            # there is a high likelihood that this is due to the sequence (tRNA-seq transcript) being a chimera
            # or another type of RNA with structural features akin to the 3' end of tRNA.
            # (Longer-than-expected sequences with a complete feature set are often pre-tRNA.)
            if not self.has_complete_feature_set:
                if len(self.input_seq) > LONGEST_KNOWN_TRNA_LENGTH:
                    self.is_predicted_trna = False
        else:
            self.acceptor_variant_string = None

            self.alpha_start = None
            self.alpha_stop = None
            self.beta_start = None
            self.beta_stop = None

            self.anticodon_seq = ''
            self.anticodon_aa = ''

            self.is_predicted_trna = False

        self.is_fully_profiled = (self.input_seq == self.profiled_seq)


    def get_profile(self,
                    unprofiled_seq='', # string in 3' to 5' direction
                    profiled_seq='', # string in 5' to 3' direction
                    features=None, # listed in the 5' to 3' direction
                    num_conserved=0,
                    num_unconserved=0,
                    num_paired=0,
                    num_unpaired=0,
                    num_extra_threeprime=0,
                    feature_index=0,
                    has_complete_feature_set=False):

        if features is None:
            features = []

        if feature_index == len(self.THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES):
            # To reach this point,
            # all tRNA features including the tRNA-His 5'-G must have been found,
            # and the input sequence must extend 5' of that.
            return {'profiled_seq': profiled_seq, # string in 5' to 3' direction
                    'features': features, # listed in the 5' to 3' direction
                    'num_conserved': num_conserved,
                    'num_unconserved': num_unconserved,
                    'num_paired': num_paired,
                    'num_unpaired': num_unpaired,
                    'num_extra_threeprime': num_extra_threeprime,
                    'num_extrapolated': 0, # number of nucleotides in an extrapolated 5' feature -- there is no extrapolated 5' feature
                    'is_complete': has_complete_feature_set,
                    'num_extra_fiveprime': len(unprofiled_seq)} # extra 5' nucleotides

        if not unprofiled_seq:
            # To reach this point,
            # the full length of the input sequence must have been profiled with tRNA features.
            return {'profiled_seq': profiled_seq,
                    'features': features, # listed in the 5' to 3' direction
                    'num_conserved': num_conserved,
                    'num_unconserved': num_unconserved,
                    'num_paired': num_paired,
                    'num_unpaired': num_unpaired,
                    'num_extra_threeprime': num_extra_threeprime,
                    'num_extrapolated': 0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                    'is_complete': has_complete_feature_set,
                    'num_extra_fiveprime': 0} # input sequence is not longer than full-length tRNA


        feature_class = self.THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES[feature_index]
        if feature_class in self.STEM_FORMATION_TRIGGERS: # 3' stem seq triggers stem formation
            # Prepare to form a stem as well as the 3' sequence.
            make_stem = True
            stem_class = feature_class.stem_class
            threeprime_stem_seq = features[-self.THREEPRIME_STEM_SEQ_INDICES[stem_class] - 1]
            if stem_class in self.ARM_FORMATION_TRIGGERS:
                # Prepare to form an arm (stem + loop) as well as the stem.
                make_arm = True
                arm_class = stem_class.arm_class
                loop = features[-self.ARM_LOOP_INDEX_DICT[arm_class] - 1]
            else:
                make_arm = False
        else:
            make_stem = False
            make_arm = False
            if feature_class == TRNAHisPositionZero:
                # Check for tRNA-His based on the anticodon sequence.
                # tRNA-His uniquely has an extra nucleotide (G) at the 5' end.
                anticodon_string = features[-self.ANTICODON_LOOP_INDEX - 1].anticodon.string
                try:
                    aa_string = ANTICODON_TO_AA[anticodon_string]
                except KeyError:
                    aa_string = 'NA'
                if aa_string != 'His':
                    # The input sequence is longer than full-length tRNA.
                    return {'profiled_seq': profiled_seq,
                            'features': features, # listed in the 5' to 3' direction
                            'num_conserved': num_conserved,
                            'num_unconserved': num_unconserved,
                            'num_paired': num_paired,
                            'num_unpaired': num_unpaired,
                            'num_extra_threeprime': num_extra_threeprime,
                            'num_extrapolated': 0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                            'is_complete': has_complete_feature_set,
                            'num_extra_fiveprime': len(unprofiled_seq)} # extra 5' nucleotides


        # This list stores the result of a recursive function call finding subsequent 5' features.
        incremental_profile_candidates = []

        # Each primary sequence feature takes (sub)sequence inputs, which can be of varying length.
        # Consider each possible combination of input lengths for the feature,
        # e.g., the D loop contains alpha and beta subsequences of variable length.
        for input_lengths, summed_input_length in zip(feature_class.allowed_input_lengths, feature_class.summed_input_lengths):

            # Strands of unequal length cannot form a stem.
            if make_stem:
                # Compare the lengths of the 5' and 3' sequence components.
                if input_lengths != tuple(map(len, threeprime_stem_seq.string_components[::-1])):
                    continue

            # Determine whether there is enough information in the remaining 5' end of the input sequence
            # to assign it to a feature despite the incompleteness of the feature sequence.
            if len(unprofiled_seq) < summed_input_length:

                # Features that lack conserved positions
                # or that do not form base pairs with a previously profiled stem sequence
                # do not contain any information for extrapolation of an incomplete sequence.
                if feature_class in self.EXTRAPOLATION_INELIGIBLE_FEATURES:
                    continue
                # The unprofiled sequence must be at least 6 nucleotides long
                # to span the 2 conserved positions in the anticodon loop (positions 33 and 37).
                # Conservation of these positions is here considered to be the minimum information
                # needed to identify the anticodon loop.
                elif feature_class == AnticodonLoop:
                    if len(unprofiled_seq) < 6:
                        continue

                string_components = [] # feature input substrings
                num_processed_bases = 0
                for input_length in input_lengths[::-1]: # create substrings from 3' to 5'
                    threeprime_to_fiveprime_string = ''
                    for _ in range(input_length):
                        if num_processed_bases >= len(unprofiled_seq):
                            threeprime_to_fiveprime_string += 'N' # pad the string at the 5' end
                        else:
                            threeprime_to_fiveprime_string += unprofiled_seq[num_processed_bases]
                        num_processed_bases += 1
                    string_components.insert(0, threeprime_to_fiveprime_string[::-1]) # 5' to 3'
                # WITH THE N PADDING, THE START INDEX (5') OF THE FEATURE IN THE INPUT SEQUENCE IS NEGATIVE.
                feature = feature_class(*string_components,
                                        start_index=len(self.input_seq) - len(profiled_seq) - num_processed_bases,
                                        stop_index=len(self.input_seq) - len(profiled_seq))

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
                                    incremental_profile_candidates.append(
                                        {'unprofiled_seq': unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                         'features': [arm, stem, feature],
                                         'num_conserved': feature.num_conserved,
                                         'num_unconserved': feature.num_unconserved,
                                         'num_paired': stem.num_paired,
                                         'num_unpaired': stem.num_unpaired,
                                         'num_extrapolated': summed_input_length - len(unprofiled_seq)}) # number of nucleotides in extrapolated 5' feature
                                    continue
                            else:
                                incremental_profile_candidates.append(
                                    {'unprofiled_seq': unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                     'features': [stem, feature],
                                     'num_conserved': feature.num_conserved,
                                     'num_unconserved': feature.num_unconserved,
                                     'num_paired': stem.num_paired,
                                     'num_unpaired': stem.num_unpaired,
                                     'num_extrapolated': summed_input_length - len(unprofiled_seq)}) # number of nucleotides in extrapolated 5' feature
                                continue
                    else:
                        incremental_profile_candidates.append(
                            {'unprofiled_seq': unprofiled_seq[::-1], # flip orientation to 5' to 3'
                             'features': [feature],
                             'num_conserved': feature.num_conserved,
                             'num_unconserved': feature.num_unconserved,
                             'num_paired': 0, # since we're not considering a stem, 0
                             'num_unpaired': 0, # since we're not considering a stem, 0
                             'num_extrapolated': summed_input_length - len(unprofiled_seq)}) # number of nucleotides in extrapolated 5' feature
                        continue

            # The procedure for assigning full-length features
            # is similar to the prior precedure for partial-length features
            # with a few efficiencies and special consideration of the acceptor (always a full-length feature).
            else:
                string_components = [] # feature input substrings
                num_processed_bases = 0
                for input_length in input_lengths[::-1]: # create substrings from 3' to 5'
                    string_components.insert(0, unprofiled_seq[num_processed_bases: num_processed_bases + input_length][::-1]) # flip orientation to 5' to 3'
                    num_processed_bases += input_length
                if feature_class.name == 'Acceptor':
                    feature = feature_class(*string_components,
                                            num_extra_threeprime=num_extra_threeprime,
                                            start_index=len(self.input_seq) - len(profiled_seq) - num_processed_bases,
                                            stop_index=len(self.input_seq) - len(profiled_seq))
                else:
                    feature = feature_class(*string_components,
                                            start_index=len(self.input_seq) - len(profiled_seq) - num_processed_bases,
                                            stop_index=len(self.input_seq) - len(profiled_seq))

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
                                    incremental_profile_candidates.append(
                                        {'unprofiled_seq': unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                         'features': [arm, stem, feature],
                                         'num_conserved': feature.num_conserved,
                                         'num_unconserved': feature.num_unconserved,
                                         'num_paired': stem.num_paired,
                                         'num_unpaired': stem.num_unpaired,
                                         'num_extrapolated': 0}) # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                                    continue
                            else:
                                incremental_profile_candidates.append(
                                    {'unprofiled_seq': unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                     'features': [stem, feature],
                                     'num_conserved': feature.num_conserved,
                                     'num_unconserved': feature.num_unconserved,
                                     'num_paired': stem.num_paired,
                                     'num_unpaired': stem.num_unpaired,
                                     'num_extrapolated': 0}) # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                                continue
                    else:
                        incremental_profile_candidates.append(
                            {'unprofiled_seq': unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                             'features': [feature],
                             'num_conserved': feature.num_conserved,
                             'num_unconserved': feature.num_unconserved,
                             'num_paired': 0, # since we're not considering a stem, 0
                             'num_unpaired': 0, # since we're not considering a stem, 0
                             'num_extrapolated': 0}) # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature

                        # Avoid testing partial acceptor sequences (C and CC) if the full one (CCA) was just found.
                        if feature.name == 'Acceptor':
                            if feature.string == 'CCA':
                                break # the following "continue" would otherwise test the length 2 acceptor sequence, CC

                        continue
                # When considering acceptor sequences with extra 3' nucleotides
                # (trying to find tRNA ending in CCAN or CCANN),
                # don't allow the CCA part to vary from CCA.
                # The "elif" block is entered if such a variant is encountered.
                # Without the "break" below, variants of CCN/CCNN and CN/CNN would then be considered.
                elif feature.name == 'Acceptor':
                    if summed_input_length == 3:
                        if num_extra_threeprime > 0:
                            break

            # Avoid testing CCN or CN if CCAN or CCANN was just found.
            if feature_class.name == 'Acceptor':
                if num_extra_threeprime > 0:
                    break


        if not incremental_profile_candidates:
            # The feature didn't pass muster.
            if feature_class.name == 'Acceptor':
                if num_extra_threeprime < 2:
                    # This will try to find an acceptor of CCAN, or if CCAN was just tested, CCANN.
                    # We "return" because we are not bothering to compare this next profile to the current non-existent profile.
                    return self.get_profile(unprofiled_seq=unprofiled_seq[1: ],
                                            profiled_seq=unprofiled_seq[0] + profiled_seq, # 5' to 3' orientation
                                            features=[], # extra 3' nucleotides are not counted as a feature
                                            num_extra_threeprime=num_extra_threeprime + 1)

            return {'profiled_seq': profiled_seq,
                    'features': features,
                    'num_conserved': num_conserved,
                    'num_unconserved': num_unconserved,
                    'num_paired': num_paired,
                    'num_unpaired': num_unpaired,
                    'num_extra_threeprime': num_extra_threeprime,
                    'num_extrapolated': 0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                    'is_complete': has_complete_feature_set,
                    'num_extra_fiveprime': 0} # input sequence is not longer than full-length tRNA


        # Sort candidates by
        # 1. number of features identified (at most, sequence + stem + arm) (descending),
        # 2. number of unconserved + unpaired nucleotides (ascending),
        # 3. incompleteness of the last (most 5') feature (ascending).
        # This sort also happens later for full sequence profiles,
        # but this first sort is useful for seeking out and returning "flawless" mature tRNA.
        incremental_profile_candidates.sort(key=lambda p: (-len(p['features']),
                                                           p['num_unconserved'] + p['num_unpaired'],
                                                           p['num_extrapolated']))
        # Continue finding features in input sequences that have not been fully profiled --
        # do not recurse profile candidates in which the final feature was extrapolated.
        profile_candidates = []
        for ipc in incremental_profile_candidates:
            if ipc['num_extrapolated'] == 0: # 5' feature was NOT extrapolated from unprofiled input sequence
                if has_complete_feature_set or feature_class.name == '5\' Acceptor Stem Sequence':
                    profile_candidate = self.get_profile(unprofiled_seq=unprofiled_seq[len(ipc['unprofiled_seq']): ], # recurse
                                                         profiled_seq=ipc['unprofiled_seq'] + profiled_seq,
                                                         features=ipc['features'] + features,
                                                         num_conserved=ipc['num_conserved'] + num_conserved,
                                                         num_unconserved=ipc['num_unconserved'] + num_unconserved,
                                                         num_paired=ipc['num_paired'] + num_paired,
                                                         num_unpaired=ipc['num_unpaired'] + num_unpaired,
                                                         num_extra_threeprime=num_extra_threeprime,
                                                         feature_index=feature_index + len(ipc['features']),
                                                         has_complete_feature_set=True)
                else:
                    profile_candidate = self.get_profile(unprofiled_seq=unprofiled_seq[len(ipc['unprofiled_seq']): ], # recurse
                                                         profiled_seq=ipc['unprofiled_seq'] + profiled_seq,
                                                         features=ipc['features'] + features,
                                                         num_conserved=ipc['num_conserved'] + num_conserved,
                                                         num_unconserved=ipc['num_unconserved'] + num_unconserved,
                                                         num_paired=ipc['num_paired'] + num_paired,
                                                         num_unpaired=ipc['num_unpaired'] + num_unpaired,
                                                         num_extra_threeprime=num_extra_threeprime,
                                                         feature_index=feature_index + len(ipc['features']),
                                                         has_complete_feature_set=False)
                if (profile_candidate['is_complete'] # has complete feature set
                    and profile_candidate['num_unconserved'] == 0
                    and profile_candidate['num_unpaired'] == 0):
                    return profile_candidate
                else:
                    profile_candidates.append(profile_candidate)
            else: # 5' feature was extrapolated from unprofiled input sequence
                profile_candidates.append(
                    {'profiled_seq': ipc['unprofiled_seq'] + profiled_seq,
                     'features': ipc['features'] + features,
                     'num_conserved': ipc['num_conserved'] + num_conserved,
                     'num_unconserved': ipc['num_unconserved'] + num_unconserved,
                     'num_paired': ipc['num_paired'] + num_paired,
                     'num_unpaired': ipc['num_unpaired'] + num_unpaired,
                     'num_extra_threeprime': num_extra_threeprime,
                     'num_extrapolated': ipc['num_extrapolated'], # number of nucleotides in extrapolated 5' feature
                     'is_complete': False, # does not have a complete feature set
                     'num_extra_fiveprime': 0}) # input sequence is not longer than full-length tRNA

        # Consider an extra 3' nucleotide in the acceptor sequence -- up to 2 extra.
        if feature_class.name == 'Acceptor':
            if num_extra_threeprime < 2:
                profile_candidates.append(self.get_profile(unprofiled_seq=unprofiled_seq[1: ],
                                                           profiled_seq=unprofiled_seq[0] + profiled_seq, # 5' to 3' orientation
                                                           features=[], # extra 3' nucleotides are not counted as a feature
                                                           num_extra_threeprime=num_extra_threeprime + 1))

        # Do not add 5' features to profiled 3' features if the additional features
        # have no grounding in conserved nucleotides or paired nucleotides in stems.
        profile_candidates = [p for p in profile_candidates
                              if p['num_conserved'] + p['num_paired'] > num_conserved + num_paired]
        if profile_candidates:
            profile_candidates.sort(key=lambda p: (-len(p['features']),
                                                   p['num_unconserved'] + p['num_unpaired'],
                                                   p['num_extrapolated']))
            return profile_candidates[0]
        else:
            return {'profiled_seq': profiled_seq,
                    'features': features,
                    'num_conserved': num_conserved,
                    'num_unconserved': num_unconserved,
                    'num_paired': num_paired,
                    'num_unpaired': num_unpaired,
                    'num_extra_threeprime': num_extra_threeprime,
                    'num_extrapolated': 0, # number of nucleotides in extrapolated 5' feature
                    'is_complete': has_complete_feature_set,
                    'num_extra_fiveprime': 0} # input sequence is not longer than full-length tRNA


    def get_unconserved_positions(self):
        unconserved_info = []
        for feature in self.features:
            # Only _Nucleotide and _Sequence subclasses have the attribute.
            if hasattr(feature, 'conserved_status'):
                component_start_index = feature.start_index
                # Conserved nucleotides are indexed within the string "component" (substring).
                for string_component_statuses, string_component in zip(
                    feature.conserved_status, feature.string_components):
                    for nucleotide_index, is_conserved, observed_nucleotide, expected_nucleotides in string_component_statuses:
                        # Avoid N padding in an extrapolated 5' feature.
                        if not is_conserved and observed_nucleotide != 'N':
                            unconserved_info.append((component_start_index + nucleotide_index,
                                                     observed_nucleotide,
                                                     ','.join(expected_nucleotides)))
                    component_start_index += len(string_component)
        return unconserved_info


    def get_unpaired_positions(self):
        unpaired_info = []
        for feature in self.features:
            # Only the _Stem class has the attribute.
            if hasattr(feature, 'paired_status'):
                for nucleotide_index, (is_paired, fiveprime_nucleotide, threeprime_nucleotide) in enumerate(feature.paired_status):
                    # Avoid N padding in an extrapolated 5' feature.
                    if not is_paired and fiveprime_nucleotide != 'N':
                        unpaired_info.append((feature.fiveprime_seq.start_index + nucleotide_index,
                                              feature.threeprime_seq.stop_index - nucleotide_index - 1,
                                              fiveprime_nucleotide,
                                              threeprime_nucleotide))
        return unpaired_info


class GeneProfile:
    def __init__(self, input_seq, name='', check_encoded_acceptor=True):
        self.unencoded_acceptor_profile = Profile(input_seq + 'CCA', name=name)
        if check_encoded_acceptor:
            self.encoded_acceptor_profile = Profile(input_seq, name=name)
            if self.unencoded_acceptor_profile.is_predicted_trna and not self.encoded_acceptor_profile.is_predicted_trna:
                self.predicted_profile = self.unencoded_acceptor_profile
                self.has_encoded_acceptor = False
            elif self.encoded_acceptor_profile.is_predicted_trna and not self.unencoded_acceptor_profile.is_predicted_trna:
                self.predicted_profile = self.encoded_acceptor_profile
                self.has_encoded_acceptor = True
            else:
                self.predicted_profile = None
                self.has_encoded_acceptor = None
        else:
            self.encoded_acceptor_profile = None
            if self.unencoded_acceptor_profile.is_predicted_trna:
                self.predicted_profile = self.unencoded_acceptor_profile
                self.has_encoded_acceptor = False
            else:
                self.predicted_profile = None
                self.has_encoded_acceptor = None


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
# print(profile.has_complete_feature_set)
# print(profile.is_fully_profiled)
# print(profile.get_unconserved_positions())
# print(profile.get_unpaired_positions())
# print(profile.alpha_start)
# print(profile.beta_start)
