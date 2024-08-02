# -*- coding: utf-8
# pylint: disable=line-too-long
"""tRNA identification from a nucleotide sequence."""

import re
import sys
import itertools
import pandas as pd

from tabulate import tabulate
from collections import OrderedDict

import anvio.filesnpaths as filesnpaths

from anvio.errors import TRNAIdentifierError
from anvio.filesnpaths import is_file_exists, is_output_file_writable
from anvio.constants import WC_BASE_PAIRS, WC_PLUS_WOBBLE_BASE_PAIRS, anticodon_to_AA as ANTICODON_TO_AA


__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


class TRNAFeature(object):
    """The most ancestral superclass of tRNA features.

    Parameters
    ==========
    string_components : tuple
        Substring components of the feature, with divisions based on "subfeatures," such as the
        variable-length alpha and beta regions of the D loop:
        e.g., ('AG', 'TT', 'GG', 'G', 'A') for the D loop of yeast tRNA-Phe-GAA

    num_allowed_unconserved : int, -1
        Number of unconserved nucleotides allowed among canonically conserved nucleotides in the
        feature. The default of -1 is equivalent to allowing all to be unconserved.

    cautious : bool, False
    """

    # Full tRNA feature parameter names found in .ini file
    INI_FEATURE_PARAMS = [
        'Conserved nucleotides',
        'Number allowed unconserved',
        'Number allowed unpaired',
        'Allowed lengths'
    ]
    # Feature parameter attributes
    ACCESSIBLE_FEATURE_PARAMS = [
        'conserved_nts',
        'num_allowed_unconserved',
        'num_allowed_unpaired',
        'allowed_input_lengths'
    ]

    def __init__(self,
                 string_components,
                 conserved_nts=None,
                 num_allowed_unconserved=-1,
                 cautious=False):
        if cautious:
            if type(string_components) != tuple:
                raise TRNAIdentifierError("`string_components` must be in the form of a tuple, e.g., ('ACTGG', 'CCAGT'). "
                                          "Your `string_components` were %s" % (string_components, ))
        self.string_components = string_components

        self.nt_count = sum(map(len, string_components))

        if conserved_nts is None:
            self.conserved_nts = ({}, )
        # By default, base conservation is not enforced.
        if num_allowed_unconserved == -1:
            self.num_allowed_unconserved = sum(len(d) for d in self.conserved_nts)
        else:
            self.num_allowed_unconserved = num_allowed_unconserved


    def check_conserved_nts(self):
        """Determine which of the canonical nucleotides in the feature are conserved.

        Returns
        =======
        meets_conserved_thresh : bool

        num_conserved : int

        num_unconserved : int

        conserved_status : list
            Nested list with the following structure:
            [[(), (), ...], [(), (), ...], ...]
            There is an inner list for each subsequence (substring component) and an inner tuple
            for each conserved nucleotide of the subsequence.
            Each inner tuple has four elements:
                1. position of conserved nucleotide in subsequence
                2. whether nucleotide is conserved (bool)
                3. observed nucleotide in subsequence (char)
                4. expected canonical nucleotide in subsequence (char)
        """
        num_conserved = 0
        num_unconserved = 0 # can include N "padding" in extrapolated 5' feature
        conserved_status = []
        for substring, nt_dict in zip(self.string_components, self.conserved_nts):
            substring_statuses = []
            conserved_status.append(substring_statuses)
            for pos, expected_nts in nt_dict.items():
                observed_nt = substring[pos]
                if observed_nt in expected_nts:
                    num_conserved += 1
                    substring_statuses.append((pos, True, observed_nt, expected_nts))
                else:
                    num_unconserved += 1
                    substring_statuses.append((pos, False, observed_nt, expected_nts))

        if num_unconserved > self.num_allowed_unconserved:
            meets_conserved_thresh = False
        else:
            meets_conserved_thresh = True

        return (meets_conserved_thresh, num_conserved, num_unconserved, conserved_status)


    @staticmethod
    def list_all_tRNA_features():
        """List all tRNA feature classes in order from 5' to 3'."""
        return [
            TRNAHisPositionZero,
            AcceptorStem,
            AcceptorStemFiveprimeStrand,
            PositionEight,
            PositionNine,
            DArm,
            DStem,
            DStemFiveprimeStrand,
            DLoop,
            DStemThreeprimeStrand,
            PositionTwentySix,
            AnticodonArm,
            AnticodonStem,
            AnticodonStemFiveprimeStrand,
            AnticodonLoop,
            AnticodonStemThreeprimeStrand,
            VLoop,
            TArm,
            TStem,
            TStemFiveprimeStrand,
            TLoop,
            TStemThreeprimeStrand,
            AcceptorStemThreeprimeStrand,
            Discriminator,
            ThreeprimeTerminus
        ]


    @staticmethod
    def list_primary_tRNA_features():
        """List tRNA feature classes not including stems and arms (secondary structures) in order
        from 5' to 3'."""
        return [
            TRNAHisPositionZero,
            AcceptorStemFiveprimeStrand,
            PositionEight,
            PositionNine,
            DStemFiveprimeStrand,
            DLoop,
            DStemThreeprimeStrand,
            PositionTwentySix,
            AnticodonStemFiveprimeStrand,
            AnticodonLoop,
            AnticodonStemThreeprimeStrand,
            VLoop,
            TStemFiveprimeStrand,
            TLoop,
            TStemThreeprimeStrand,
            AcceptorStemThreeprimeStrand,
            Discriminator,
            ThreeprimeTerminus
        ]


class Nucleotide(TRNAFeature):
    """Superclass for tRNA primary features of a single nucleotide, e.g., discriminator, position 8."""

    allowed_input_lengths = ((1, ), )
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))

    def __init__(self,
                 string, # must be a string of length 1
                 conserved_nts=None,
                 num_allowed_unconserved=-1,
                 start_pos=None,
                 stop_pos=None,
                 cautious=False):
        self.string = string
        self.start_pos = start_pos
        self.stop_pos = stop_pos

        super().__init__((string, ),
                         conserved_nts=conserved_nts,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved,
         self.conserved_status) = self.check_conserved_nts()


class Sequence(TRNAFeature):
    """Superclass for tRNA primary sequences, e.g., 5' strand of T stem."""

    def __init__(self,
                 substrings, # must be a string, tuple of strings, or tuple of Nucleotide/Sequence objects
                 conserved_nts=None,
                 num_allowed_unconserved=-1,
                 start_pos=None,
                 stop_pos=None,
                 cautious=False):
        if type(substrings) == str:
            string_components = (substrings, )
        elif all([type(s) == str for s in substrings]):
            string_components = substrings
        elif all([type(s) == Nucleotide or type(s) == Sequence for s in substrings]):
            string_components = tuple(s.string for s in substrings)
        self.string = ''.join(substrings)
        self.start_pos = start_pos
        self.stop_pos = stop_pos

        super().__init__(string_components,
                         conserved_nts=conserved_nts,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved, # can include N "padding" in extrapolated 5' feature
         self.conserved_status) = self.check_conserved_nts()


class Loop(Sequence):
    """Superclass for loops: D loop, anticodon loop, T loop."""
    def __init__(self,
                 substrings, # must be a string, tuple of strings, or tuple of Nucleotide/Sequence objects
                 conserved_nts=None,
                 num_allowed_unconserved=-1,
                 start_pos=None,
                 stop_pos=None,
                 cautious=False):
        super().__init__(substrings,
                         conserved_nts=conserved_nts,
                         num_allowed_unconserved=num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class Stem(TRNAFeature):
    """Superclass for stems: acceptor stem, D stem, anticodon stem, T stem."""

    def __init__(self,
                 fiveprime_seq, # must be Sequence object
                 threeprime_seq, # must be Sequence object
                 num_allowed_unpaired=0,
                 num_allowed_unconserved=-1,
                 cautious=False):
        if cautious:
            if type(fiveprime_seq) != Sequence or type(threeprime_seq) != Sequence:
                raise TRNAIdentifierError("You can only define a Stem from Sequence objects.")
        self.fiveprime_seq = fiveprime_seq
        self.threeprime_seq = threeprime_seq

        self.canonical_positions = (*self.fiveprime_seq.canonical_positions, *self.threeprime_seq.canonical_positions)
        self.conserved_nts = (*self.fiveprime_seq.conserved_nts, *self.threeprime_seq.conserved_nts)

        if cautious:
            if tuple(map(len, self.fiveprime_seq.string_components)) != tuple(map(len, self.threeprime_seq.string_components[::-1])):
                raise TRNAIdentifierError("The two Sequence objects, %s and %s, "
                                          "that were used to define your Stem are not the same length."
                                          % (self.fiveprime_seq.string_components, self.threeprime_seq.string_components))
            length = sum(map(len, self.fiveprime_seq.string_components))
            if num_allowed_unpaired > length:
                raise TRNAIdentifierError("You tried to leave at most %d base pairs unpaired, "
                                          "but there are only %d base pairs in the stem." % (num_allowed_unpaired, length))
        self.num_allowed_unpaired = num_allowed_unpaired

        self.start_positions = (self.fiveprime_seq.start_pos, self.threeprime_seq.start_pos)
        self.stop_positions = (self.fiveprime_seq.stop_pos, self.threeprime_seq.stop_pos)

        super().__init__((*self.fiveprime_seq.string_components, *self.threeprime_seq.string_components),
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_pair_thresh,
         self.num_paired,
         self.num_unpaired, # can include N "padding" in extrapolated 5' feature
         self.paired_status) = self.check_pairs()


    def check_pairs(self):
        """Determine base pairing in the stem.

        Returns
        =======
        meets_pair_thresh : bool

        num_paired : int

        num_unpaired : int

        paired_status : list
            List of tuples, one tuple for each nucleotide pair in the stem. Each tuple has three
            elements: whether a base pair exists, the 5' nucleotide, and the 3' nucleotide.
        """
        num_paired = 0
        num_unpaired = 0 # can include N "padding" in extrapolated 5' feature
        paired_status = []
        for fiveprime_nt, threeprime_nt in zip(self.fiveprime_seq.string, self.threeprime_seq.string[::-1]):
            if fiveprime_nt in WC_PLUS_WOBBLE_BASE_PAIRS[threeprime_nt]:
                num_paired += 1
                paired_status.append((True, fiveprime_nt, threeprime_nt))
            else:
                num_unpaired += 1
                paired_status.append((False, fiveprime_nt, threeprime_nt))

        if num_unpaired > self.num_allowed_unpaired:
            meets_pair_thresh = False
        else:
            meets_pair_thresh = True

        return meets_pair_thresh, num_paired, num_unpaired, paired_status


class Arm(TRNAFeature):
    """Superclass for arms: D arm, anticodon arm, T arm

    The number of unconserved nucleotides allowed in the arm
    can differ from the sum of the numbers of unconserved nucleotides allowed in the stem and loop.
    """

    def __init__(self,
                 stem, # must be Stem object
                 loop, # must be Loop object
                 num_allowed_unconserved=-1,
                 cautious=False):
        if cautious:
            if type(stem) != Stem or type(loop) != Loop:
                raise TRNAIdentifierError("A `Stem` and a `Loop` are required input to create an `Arm`.")
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
                raise TRNAIdentifierError("The canonical positions in an `Arm` must be contiguous. "
                                          "These were yours: %s. "
                                          "These came from the canonical positions in Stem, %s, "
                                          "and the canonical positions in Loop, %s."
                                          % (self.canonical_positions,
                                             stem.canonical_positions,
                                             loop.canonical_positions))

        self.conserved_nts=(*stem.fiveprime_seq.conserved_nts,
                            *loop.conserved_nts,
                            *stem.threeprime_seq.conserved_nts)

        self.start_pos = self.stem.start_positions[0]
        self.stop_pos = self.stem.stop_positions[1]

        super().__init__((*stem.fiveprime_seq.string_components,
                          *loop.string_components,
                          *stem.threeprime_seq.string_components),
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        if (self.stem.fiveprime_seq.num_unconserved
            + self.loop.num_unconserved
            + self.stem.threeprime_seq.num_unconserved) > self.num_allowed_unconserved:
            self.meets_conserved_thresh = False
        else:
            self.meets_conserved_thresh = True


class TRNAHisPositionZero(Nucleotide):
    """The G typically found at the 5' end of mature tRNA-His."""

    name = 'tRNA-His position 0'
    canonical_positions = ((-1, ), )
    conserved_nts = ({0: 'G'}, )
    num_allowed_unconserved = 0

    def __init__(self, string, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(string,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class AcceptorStem(Stem):
    """The base-paired nucleotides of the acceptor stem of tRNA (excludes the discriminator and
    CCA acceptor or other 3' terminus)."""

    name = 'acceptor stem'
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class AcceptorStemFiveprimeStrand(Sequence):
    """The 5' strand of the acceptor stem of tRNA."""

    name = 'acceptor stem 5\' strand'
    canonical_positions = ((1, 2, 3, 4, 5, 6, 7), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1
    stem_class = AcceptorStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class PositionEight(Nucleotide):
    """The nucleotide at canonical position 8 of tRNA, expected but not required to be T."""

    name = 'position 8'
    canonical_positions = ((8, ), )
    conserved_nts = ({0: 'T'}, )
    num_allowed_unconserved = 1

    def __init__(self, string, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(string,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class PositionNine(Nucleotide):
    """The nucleotide at canonical position 9 of tRNA."""

    name = 'position 9'
    canonical_positions = ((9, ), )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, string, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(string,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class DArm(Arm):
    """The D arm of tRNA."""

    name = 'D arm'
    num_allowed_unconserved = 1

    def __init__(self, stem, loop, cautious=False):
        super().__init__(stem,
                         loop,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         cautious=cautious)


class DStem(Stem):
    """The D stem of tRNA. Type II (long V loop) tRNAs often have D stems of length 3 rather than 4,
    but the nucleotides at canonical positions 13 and 22 are always included in the stem objects
    rather than D loop object."""

    name = 'D stem'
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1
    arm_class = DArm

    def __init__(self, fiveprime_seq, threeprime_seq, type_II_trna=False, cautious=False):
        self.type_II_trna = type_II_trna
        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


    def check_pairs(self):
        """This method overrides the one with the same name in the Stem superclass."""
        self.paired_positions_13_22_in_type_II = False
        self.unpaired_positions_13_22_in_type_II = False
        if self.type_II_trna and self.num_allowed_unpaired != 4:
            # Do not penalize type II tRNAs for having unpaired nucleotides at canonical positions
            # 13 and 22, unless the user changes the parameterization of `num_allowed_unpaired` to
            # force positions 10-13 to always pair with 25-22.
            num_paired = 0
            num_unpaired = 0 # can include N "padding" in extrapolated 5' feature
            paired_status = []
            for fiveprime_nt, threeprime_nt in zip(self.fiveprime_seq.string[: -1], self.threeprime_seq.string[::-1]):
                if fiveprime_nt in WC_PLUS_WOBBLE_BASE_PAIRS[threeprime_nt]:
                    num_paired += 1
                    paired_status.append((True, fiveprime_nt, threeprime_nt))
                else:
                    num_unpaired += 1
                    paired_status.append((False, fiveprime_nt, threeprime_nt))

            if num_unpaired > self.num_allowed_unpaired:
                meets_pair_thresh = False
            else:
                meets_pair_thresh = True

            fiveprime_nt = self.fiveprime_seq.string[-1]
            threeprime_nt = self.threeprime_seq.string[0]
            if fiveprime_nt in WC_PLUS_WOBBLE_BASE_PAIRS[threeprime_nt]:
                num_paired += 1
                paired_status.append((True, fiveprime_nt, threeprime_nt))
                self.paired_positions_13_22_in_type_II = True
            else:
                num_unpaired += 1
                paired_status.append((False, fiveprime_nt, threeprime_nt))
                self.unpaired_positions_13_22_in_type_II = True
        else:
            meets_pair_thresh, num_paired, num_unpaired, paired_status = super().check_pairs()

        return meets_pair_thresh, num_paired, num_unpaired, paired_status


class DStemFiveprimeStrand(Sequence):
    """The 5' strand of the D stem of tRNA. Type II (long V loop) tRNAs often have D stems of length
    3 rather than 4, but the nucleotides at canonical positions 13 and 22 are always included in the
    stem objects rather than D loop object."""

    name = 'D stem 5\' strand'
    canonical_positions = ((10, 11, 12, 13), )
    allowed_input_lengths = ((4, ), )
    summed_input_lengths = (4, )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1
    arm_class = DArm
    stem_class = DStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=False)


class DLoop(Loop):
    """The D loop of tRNA, allowing for variable alpha (canonical positions 16, 17, 17a, 17b) and
    beta (canonical positions 20, 20a, 20b) sections. Type II (long V loop) tRNAs often have D stems
    of length 3 rather than 4, but the nucleotides at canonical positions 13 and 22 are always
    included in the stem objects rather than D loop object."""

    name = 'D loop'
    canonical_positions = ((14, 15), (16, 17), (18, 19), (20, ), (21, ))
    allowed_section_lengths = ((2, ), (1, 2, 3), (2, ), (1, 2, 3), (1, ))
    allowed_input_lengths = tuple(itertools.product(*allowed_section_lengths))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nts = ({0: 'A', 1: ('A', 'G')}, {}, {0: 'G', 1: 'G'}, {}, {0: ('A', 'G')})
    num_allowed_unconserved = 1
    arm_class = DArm
    stem_class = DStem

    def __init__(self,
                 positions_14_to_15_string,
                 alpha_positions_string,
                 positions_18_to_19_string,
                 beta_positions_string,
                 position_21_string,
                 start_pos=None,
                 stop_pos=None,
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

        alpha_start_pos = 2
        alpha_stop_pos = alpha_start_pos + len(alpha_positions_string)
        self.alpha_seq = Sequence(alpha_positions_string, start_pos=alpha_start_pos, stop_pos=alpha_stop_pos)
        beta_start_pos = alpha_stop_pos + 2
        beta_stop_pos = beta_start_pos + len(beta_positions_string)
        self.beta_seq = Sequence(beta_positions_string, start_pos=beta_start_pos, stop_pos=beta_stop_pos)

        super().__init__((positions_14_to_15_string,
                          alpha_positions_string,
                          positions_18_to_19_string,
                          beta_positions_string,
                          position_21_string),
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class DStemThreeprimeStrand(Sequence):
    """The 3' strand of the D stem of tRNA. Type II (long V loop) tRNAs often have D stems of length
    3 rather than 4, but the nucleotides at canonical positions 13 and 22 are always included in the
    stem objects rather than D loop object."""

    name = 'D stem 3\' strand'
    canonical_positions = ((22, 23, 24, 25), )
    allowed_input_lengths = ((4, ), )
    summed_input_lengths = (4, )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1
    arm_class = DArm
    stem_class = DStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=False)


class PositionTwentySix(Nucleotide):
    """The nucleotide at canonical position 26 of tRNA."""

    name = 'position 26'
    canonical_positions = ((26, ), )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, string, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(string,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class AnticodonArm(Arm):
    """The anticodon arm of tRNA."""

    name = 'anticodon arm'
    num_allowed_unconserved = 1

    def __init__(self, stem, loop, cautious=False):
        self.anticodon = loop.anticodon
        super().__init__(stem,
                         loop,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         cautious=cautious)


class AnticodonStem(Stem):
    """The anticodon stem of tRNA."""

    name = 'anticodon stem'
    arm_class = AnticodonArm
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class AnticodonStemFiveprimeStrand(Sequence):
    """The 5' strand of the anticodon stem of tRNA."""

    name = 'anticodon stem 5\' strand'
    canonical_positions = ((27, 28, 29, 30, 31), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1
    stem_class = AnticodonStem
    arm_class = AnticodonArm

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class Anticodon(Sequence):
    """The anticodon sequence of tRNA."""

    name = 'anticodon'
    canonical_positions = ((34, 35, 36))
    allowed_input_lengths = ((3, ), )
    summed_input_lengths = (3, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings, start_pos=start_pos, stop_pos=stop_pos, cautious=cautious)
        try:
            self.aa_string = ANTICODON_TO_AA[self.string]
        except KeyError:
            self.aa_string = 'NA'


class AnticodonLoop(Loop):
    """The anticodon loop of tRNA."""

    name = 'anticodon loop'
    canonical_positions = ((32, 33, 34, 35, 36, 37, 38), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nts = ({1: 'T', 5: ('A', 'G')}, )
    num_allowed_unconserved = 1
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)
        self.anticodon = Anticodon(self.string[2: 5])


class AnticodonStemThreeprimeStrand(Sequence):
    """The 3' strand of the anticodon stem of tRNA."""

    name = 'anticodon stem 3\' strand'
    canonical_positions = ((39, 40, 41, 42, 43), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class VLoop(Loop):
    """The V arm of tRNA: no stem/loop structure and base pairing is considered, so call it a
    loop."""

    name = 'V loop'
    canonical_positions = ((44, 45, 46, 47, 48), )
    allowed_input_lengths = tuple(itertools.product(range(4, 6))) + tuple(itertools.product(range(9, 24)))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nts = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)
        if 4 <= len(self.string) <= 5:
            self.type = 'I'
        elif 12 <= len(self.string) <= 23:
            self.type = 'II'
        else:
            self.type = 'NA'


class TArm(Arm):
    """The T arm of tRNA."""

    name = 'T arm'
    num_allowed_unconserved = 1

    def __init__(self, stem, loop, cautious=False):
        super().__init__(stem,
                         loop,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         cautious=cautious)


class TStem(Stem):
    """The T stem of tRNA, with canonical positions 53 and 61 expected to be G and C."""

    name = 'T stem'
    arm_class = TArm
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class TStemFiveprimeStrand(Sequence):
    """The 5' T strand of the T stem of tRNA."""

    name = 'T stem 5\' strand'
    canonical_positions = ((49, 50, 51, 52, 53), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nts = ({4: 'G'}, )
    num_allowed_unconserved = 1
    arm_class = TArm
    stem_class = TStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class TLoop(Loop):
    """The T loop of tRNA."""

    name = 'T loop'
    canonical_positions = ((54, 55, 56, 57, 58, 59, 60), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nts = ({0: 'T', 1: 'T', 2: 'C', 3: ('A', 'G'), 4: 'A'}, )
    num_allowed_unconserved = 1
    arm_class = TArm

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class TStemThreeprimeStrand(Sequence):
    """The 3' strand of the T stem of tRNA."""

    name = 'T stem 3\' strand'
    canonical_positions = ((61, 62, 63, 64, 65), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nts = ({0: 'C'}, )
    num_allowed_unconserved = 1
    arm_class = TArm
    stem_class = TStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class AcceptorStemThreeprimeStrand(Sequence):
    """The 3' strand of the acceptor stem of tRNA."""

    name = 'acceptor stem 3\' strand'
    canonical_positions=((66, 67, 68, 69, 70, 71, 72), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1
    stem_class = AcceptorStem

    def __init__(self, substrings, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(substrings,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class Discriminator(Nucleotide):
    """The discriminator nucleotide 5' of the CCA acceptor or other 3' terminus."""

    name = 'discriminator'
    canonical_positions = ((73, ), )
    conserved_nts = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, string, start_pos=None, stop_pos=None, cautious=False):
        super().__init__(string,
                         conserved_nts=self.conserved_nts,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_pos=start_pos,
                         stop_pos=stop_pos,
                         cautious=cautious)


class ThreeprimeTerminus(object):
    """The CCA acceptor or any other sequence at the 3' end of tRNA. This feature differs from the
    others due to the variety of sequences that may be found at the 3' end of tRNA-seq reads."""

    name = '3\' terminus'
    allowed_patterns = (re.compile('CCA'), )

    def __init__(self, string, start_pos, stop_pos, priority=0):
        self.string = string
        self.start_pos = start_pos
        self.stop_pos = stop_pos
        # Priority is the index of the terminus among the possible termini. Highest priority is
        # given by a value of 0.
        self.priority = priority


class TRNAFeatureParameterizer(object):
    """Interacts with tRNA feature objects, setting and reporting their tRNA-defining attributes."""

    def __init__(self):
        # Do not include the 3' terminus among the parameterizable features.
        self.dict_mapping_feature_or_subfeature_name_to_class = {
            TRNAHisPositionZero.name: TRNAHisPositionZero,
            AcceptorStem.name: AcceptorStem,
            AcceptorStemFiveprimeStrand.name: AcceptorStemFiveprimeStrand,
            PositionEight.name: PositionEight,
            PositionNine.name: PositionNine,
            DArm.name: DArm,
            DStem.name: DStem,
            DStemFiveprimeStrand.name: DStemFiveprimeStrand,
            DLoop.name: DLoop,
            DLoop.name + "/positions 14-15": DLoop,
            DLoop.name + "/alpha positions": DLoop,
            DLoop.name + "/positions 18-19": DLoop,
            DLoop.name + "/beta positions": DLoop,
            DLoop.name + "/position 21": DLoop,
            DStemThreeprimeStrand.name: DStemThreeprimeStrand,
            PositionTwentySix.name: PositionTwentySix,
            AnticodonArm.name: AnticodonArm,
            AnticodonStem.name: AnticodonStem,
            AnticodonStemFiveprimeStrand.name: AnticodonStemFiveprimeStrand,
            AnticodonLoop.name: AnticodonLoop,
            AnticodonStemThreeprimeStrand.name: AnticodonStemThreeprimeStrand,
            VLoop.name: VLoop,
            TArm.name: TArm,
            TStem.name: TStem,
            TStemFiveprimeStrand.name: TStemFiveprimeStrand,
            TLoop.name: TLoop,
            TStemThreeprimeStrand.name: TStemThreeprimeStrand,
            AcceptorStemThreeprimeStrand.name: AcceptorStemThreeprimeStrand,
            Discriminator.name: Discriminator
        }

        self.feature_and_subfeature_names_with_accessible_lengths = [
            VLoop.name,
            DLoop.name + "/alpha positions",
            DLoop.name + "/beta positions"
        ]

        self.subfeature_section_dict = {name: 0 for name in self.dict_mapping_feature_or_subfeature_name_to_class}
        self.subfeature_section_dict[DLoop.name + "/alpha positions"] = 1
        self.subfeature_section_dict[DLoop.name + "/positions 18-19"] = 2
        self.subfeature_section_dict[DLoop.name + "/beta positions"] = 3
        self.subfeature_section_dict[DLoop.name + "/position 21"] = 4


    def write_param_file(self, feature_param_path):
        """Write the table of feature and "subfeature" parameters (class attributes)."""
        filesnpaths.is_output_file_writable(feature_param_path)
        rows = self.get_param_table_as_list()
        with open(feature_param_path, 'w') as f:
            for row in rows:
                f.write("\t".join(row) + "\n")


    def tabulate_params(self):
        """Get a nicely formatted version of the table of feature and "subfeature" parameters (class
        attributes)."""
        rows = self.get_param_table_as_list()

        return tabulate(rows, headers='firstrow', tablefmt='github')


    def list_accessible_param_tuples(self, pretty=False):
        """Get a list of tuples, each containing the name and value of an accessible parameter."""
        param_tuples = []
        if pretty:
            # Used in producing the anvi-trnaseq analysis summary file.
            param_types = TRNAFeature.INI_FEATURE_PARAMS
            for feature_info in self.get_param_table_as_list()[1: ]:
                feature_name = feature_info[0]
                for param_type, param_value in zip(param_types, feature_info[1: ]):
                    if param_value != '-':
                        if feature_name in ['tRNA-His position 0']:
                            # Avoid capitalizing certain strings.
                            param_tuples.append((feature_name + ': ' + param_type, param_value))
                        else:
                            param_tuples.append((feature_name.capitalize() + ': ' + param_type, param_value))
        else:
            # Used in setting tRNA-seq database meta-values.
            param_types = TRNAFeature.ACCESSIBLE_FEATURE_PARAMS
            for feature_info in self.get_param_table_as_list()[1: ]:
                feature_name = feature_info[0].replace(' ', '_')
                for param_type, param_value in zip(param_types, feature_info[1: ]):
                    if param_value != '-':
                        param_tuples.append((
                            (feature_name + '_' + param_type).replace('-', '_').replace('/', '_').lower(),
                            param_value)
                        )

        return param_tuples


    def get_param_table_as_list(self):
        """Get a table of user-accessible feature and "subfeature" parameters (class attributes).

        A dash indicates a parameter that cannot be set for the feature. Quotes indicate a parameter
        that is currently not set.
        """
        ACCESSIBLE_FEATURE_PARAMS = TRNAFeature.ACCESSIBLE_FEATURE_PARAMS
        feature_and_subfeature_names_with_accessible_lengths = self.feature_and_subfeature_names_with_accessible_lengths
        subfeature_section_dict = self.subfeature_section_dict

        rows = [["Feature/subfeature"] + TRNAFeature.INI_FEATURE_PARAMS]
        for feature_or_subfeature_name, feature_class in self.dict_mapping_feature_or_subfeature_name_to_class.items():
            row = [feature_or_subfeature_name]

            for param_name in ACCESSIBLE_FEATURE_PARAMS:
                if param_name == 'conserved_nts':
                    if '/' not in feature_or_subfeature_name and 'allowed_section_lengths' in feature_class.__dict__:
                        # Conserved nucleotides can only be set for individual "subfeature" sections of the D arm.
                        row.append("-")
                        continue
                elif param_name == 'num_allowed_unconserved':
                    if '/' in feature_or_subfeature_name: # dealing with a "subfeature"
                        row.append("-")
                        continue
                elif param_name == 'allowed_input_lengths':
                    if feature_or_subfeature_name not in feature_and_subfeature_names_with_accessible_lengths:
                        row.append("-")
                        continue

                try:
                    param_value = feature_class.__dict__[param_name]
                except KeyError:
                    row.append("-")
                    continue

                if param_name == 'conserved_nts':
                    conserved_nt_list = []
                    if not param_value[subfeature_section_dict[feature_or_subfeature_name]]:
                        row.append("\"\"")
                        continue
                    for pos, nt in param_value[subfeature_section_dict[feature_or_subfeature_name]].items():
                        if type(nt) == str:
                            conserved_nt_list.append(str(pos) + "," + nt)
                        elif type(nt) == tuple:
                            # The only alternatives to a single conserved nucleotide are a purine or pyrimidine.
                            if nt == ('A', 'G'):
                                conserved_nt_list.append(str(pos) + ",R")
                            elif nt == ('C', 'T'):
                                conserved_nt_list.append(str(pos) + ",Y")
                    row.append(";".join(conserved_nt_list))
                elif param_name == 'num_allowed_unconserved' or param_name == 'num_allowed_unpaired':
                    if param_value == -1:
                        row.append("\"\"")
                        continue
                    row.append(str(param_value))
                elif param_name == 'allowed_input_lengths':
                    if feature_or_subfeature_name == VLoop.name:
                        # Here is an example to show the format of VLoop.allowed_input_lengths:
                        # ((4, ), (5, ), (9, ), (10, ), ..., (23, ))
                        allowed_lengths_output = ""
                        prev_allowed_length = None
                        for allowed_length_tuple in param_value:
                            for allowed_length in allowed_length_tuple:
                                if prev_allowed_length:
                                    if allowed_length - prev_allowed_length > 1:
                                        allowed_lengths_output += "-" + str(prev_allowed_length) + "," + str(allowed_length)
                                else:
                                    allowed_lengths_output += str(allowed_length)
                                prev_allowed_length = allowed_length
                        allowed_lengths_output += "-" + str(prev_allowed_length)
                        row.append(allowed_lengths_output)
                    else:
                        # Sections of the D loop are the only other variable-length "subfeatures".
                        allowed_section_lengths = feature_class.allowed_section_lengths[subfeature_section_dict[feature_or_subfeature_name]]
                        row.append(str(allowed_section_lengths[0]) + "-" + str(allowed_section_lengths[-1]))
            rows.append(row)

        return rows


    def set_params_from_file(self, feature_param_path):
        """Set user-accessible tRNA feature parameters for de novo tRNA profiling and
        identification."""
        filesnpaths.is_file_exists(feature_param_path)
        feature_param_df = pd.read_csv(feature_param_path, sep='\t', header=0, index_col=0, keep_default_na=False)

        INI_FEATURE_PARAMS = TRNAFeature.INI_FEATURE_PARAMS
        ACCESSIBLE_FEATURE_PARAMS = TRNAFeature.ACCESSIBLE_FEATURE_PARAMS
        for feature_or_subfeature_name, row in feature_param_df.iterrows():
            for ini_param_name, param_name in zip(INI_FEATURE_PARAMS, ACCESSIBLE_FEATURE_PARAMS):
                param_value = row[ini_param_name]
                if param_value == '-':
                    continue
                self.set_param(feature_or_subfeature_name, param_name, param_value)


    def set_param(self, feature_or_subfeature_name, param_name, param_value):
        """Set tRNA feature parameter class attributes.

        Parameters
        ==========
        feature_or_subfeature_name : str
            The name of a feature or "subfeature" (in the case of a section of the D arm).

        param_name : str
            The name of a parameter that can be set.

        param_value : str
            The particularly formatted value of the parameter to be set.
        """
        if feature_or_subfeature_name not in self.dict_mapping_feature_or_subfeature_name_to_class:
            raise TRNAIdentifierError("`TRNAFeatureParameterizer.set_param` does not recognize the supplied feature "
                                      "or subfeature name, %s. Here are the recognized feature names: %s"
                                      % (feature_or_subfeature_name, ", ".join(self.dict_mapping_feature_or_subfeature_name_to_class)))

        if param_name not in TRNAFeature.ACCESSIBLE_FEATURE_PARAMS:
            raise TRNAIdentifierError("`TRNAFeatureParameterizer.set_param` does not recognize the supplied feature "
                                      "parameter name, %s. Here are the recognized parameter names: %s"
                                      % (param_name, ", ".join(TRNAFeature.ACCESSIBLE_FEATURE_PARAMS)))

        if param_name == 'conserved_nts':
            self.set_conserved_nts(feature_or_subfeature_name, param_value)
        elif param_name == 'num_allowed_unconserved':
            if not param_value:
                param_value = -1
            self.set_num_allowed_unconserved(feature_or_subfeature_name, param_value)
        elif param_name == 'num_allowed_unpaired':
            if not param_value:
                param_value = sys.maxsize
            self.set_num_allowed_unpaired(feature_or_subfeature_name, param_value)
        elif param_name == 'allowed_input_lengths':
            self.set_allowed_input_lengths(feature_or_subfeature_name, param_value)


    def set_conserved_nts(self, feature_or_subfeature_name, param_value=''):
        """Modify (update a dict in) the `conserved_nts` attribute of a tRNA feature class.

        Parameters
        ==========
        feature_or_subfeature_name : str
            The name of a feature or "subfeature" (in the case of a section of the D arm)

        param_value : str, ''
            The proper format is
            <Zero-based position relative to the 5' end of the feature>,<Conserved nucleotide symbol>;<Next index>,<Next symbol>;...
            The default empty string indicates that there are no conserved sites.
        """
        feature_class = self.dict_mapping_feature_or_subfeature_name_to_class[feature_or_subfeature_name]

        conserved_nts_dict = {}
        if param_value:
            try:
                for conserved_nt in param_value.split(';'):
                    conserved_pos, nt = conserved_nt.split(',')
                    if nt == 'R':
                        nt = ('A', 'G')
                    elif nt == 'Y':
                        nt = ('C', 'T')
                    conserved_nts_dict[int(conserved_pos)] = nt
            except:
                raise TRNAIdentifierError("The proper format of a conserved nucleotide parameter value is "
                                          "\"<Zero-based position of conserved nucleotide in feature(/subfeature) relative to 5' end of feature>,"
                                          "<Conserved nucleotide symbol (A, C, G, T, R, or Y)>;"
                                          "<Next zero-based position>,<Next conserved nucleotide symbol>;...\" "
                                          "Note that the position integer is separated from the nucleotide character by a comma. "
                                          "Note that entries for different conserved nucleotides are separated by a semicolon. "
                                          "The value provided was %s" % param_value)

            if ('conserved_nts' not in feature_class.__dict__
                or ('/' not in feature_or_subfeature_name
                    and 'allowed_section_lengths' in feature_class.__dict__)):
                #   Conserved nucleotides can only be set for individual "subfeature" sections of the D arm.
                raise TRNAIdentifierError("\"%s\" does not support assignment of conserved nucleotides. "
                                          "The conserved nucleotide input provided was %s"
                                          % (feature_or_subfeature_name, param_value))

        prev_conserved_nts_dict = feature_class.conserved_nts[self.subfeature_section_dict[feature_or_subfeature_name]]
        prev_conserved_nts_dict.clear()
        prev_conserved_nts_dict.update(conserved_nts_dict)


    def set_num_allowed_unconserved(self, feature_name, num_allowed_unconserved=-1):
        """Set the `num_allowed_unconserved` attribute of a tRNA feature class.

        Parameters
        ==========
        feature_name : str
            The name of a feature.

        num_allowed_unconserved : int, -1
            The number of conserved nucleotide positions in the feature allowed to be unconserved
            but still have positive identification of the feature. The default value of -1 means
            that an "unlimited" number of unconserved positions is allowed.
        """
        feature_class = self.dict_mapping_feature_or_subfeature_name_to_class[feature_name]

        num_allowed_unconserved = int(num_allowed_unconserved)
        if 'num_allowed_unconserved' not in feature_class.__dict__ or '/' in feature_name:
            raise TRNAIdentifierError("\"%s\" does not support assignment of an allowed number of unconserved nucleotides. "
                                      "The number of allowed unconserved nucleotides provided was %s"
                                      % (feature_name, num_allowed_unconserved))
        feature_class.num_allowed_unconserved = num_allowed_unconserved


    def set_num_allowed_unpaired(self, stem_name, num_allowed_unpaired=sys.maxsize):
        """Set the `num_allowed_unpaired` attribute of a stem class.

        Parameters
        ==========
        stem_name : str
            The name of a stem feature.

        param_value : int, sys.maxsize
            The number of positions in the stem allowed to be unpaired. Pairing means Watson-Crick
            or wobble G-T. The huge default value means that an "unlimited" number of unpaired
            positions is allowed.
        """
        feature_class = self.dict_mapping_feature_or_subfeature_name_to_class[stem_name]

        num_allowed_unpaired = int(num_allowed_unpaired)
        if 'num_allowed_unpaired' not in feature_class.__dict__:
            raise TRNAIdentifierError("\"%s\" does not support assignment of an allowed number of unpaired nucleotides. "
                                      "The number of allowed unpaired nucleotides provided was %s"
                                      % (stem_name, num_allowed_unpaired))
        feature_class.num_allowed_unpaired = num_allowed_unpaired


    def set_allowed_input_lengths(self, feature_or_subfeature_name, param_value):
        """Set attributes related to the range of allowed lengths in a variable-length tRNA feature.

        Parameters
        ==========
        feature_or_subfeature_name : str
            The name of a feature or "subfeature" (in the case of a section of the D arm).

        param_value : str
            A string representing the allowed length range,
            which can be discontinuous, with the format
            <minimum length>-<maximum length>,<next minimum length>-<next maximum length>,...
        """
        feature_class = self.dict_mapping_feature_or_subfeature_name_to_class[feature_or_subfeature_name]
        if feature_or_subfeature_name not in self.feature_and_subfeature_names_with_accessible_lengths:
            raise TRNAIdentifierError("\"%s\" does not support assignment of variable lengths. "
                                      "The length range provided was %s" % (feature_or_subfeature_name, param_value))

        allowed_lengths = tuple()
        for length_range_input in param_value.split(','):
            try:
                min_length, max_length = length_range_input.split('-')
                allowed_lengths += tuple(range(int(min_length), int(max_length) + 1))
            except:
                raise TRNAIdentifierError("The proper format of a parameter value in the allowed feature length field is "
                                          "<Minimum length integer>-<Maximum length integer>,"
                                          "<Next minimum length integer>-<Next maximum length integer>,... "
                                          "The parameter value provided was %s" % param_value)

        if feature_or_subfeature_name == VLoop.name:
            feature_class.allowed_input_lengths = tuple(itertools.product(allowed_lengths))
            # Reset the dependent class attribute.
            feature_class.summed_input_lengths = allowed_lengths
        else:
            # Beside the variable loop, subfeatures of the D loop have variable lengths.
            section = self.subfeature_section_dict[feature_or_subfeature_name]
            prev_allowed_section_lengths = feature_class.allowed_section_lengths
            feature_class.allowed_section_lengths = (prev_allowed_section_lengths[: section]
                                                     + (allowed_lengths, )
                                                     + prev_allowed_section_lengths[section + 1: ])
            # Reset the dependent class attributes.
            feature_class.allowed_input_lengths = tuple(itertools.product(*feature_class.allowed_section_lengths))
            feature_class.summed_input_lengths = tuple(map(sum, feature_class.allowed_input_lengths))


    @staticmethod
    def set_threeprime_termini(threeprime_termini):
        ThreeprimeTerminus.allowed_patterns = tuple([re.compile(threeprime_terminus.replace('N', '.')) for threeprime_terminus in threeprime_termini])


class Profile(object):
    """A profile of the tRNA features in a sequence. The function, `Profiler.profile`, creates these
    objects."""

    __slots__ = (
        'input_seq',
        'name',
        'profiled_seq',
        'features',
        'num_conserved',
        'num_unconserved',
        'num_paired',
        'num_unpaired',
        'num_in_extrapolated_fiveprime_feature',
        'has_complete_feature_set',
        'num_extra_fiveprime',
        'feature_names',
        'threeprime_terminus_seq',
        'alpha_start',
        'alpha_stop',
        'beta_start',
        'beta_stop',
        'anticodon_seq',
        'anticodon_aa',
        'is_predicted_trna',
        'trunc_profile_index',
        'unconserved_info',
        'unpaired_info',
        'is_fully_profiled'
    )

    def __init__(self, input_seq, name=''):
        self.input_seq = input_seq
        self.name = name

        # The remaining attributes are set in `Profiler.profile`.
        self.profiled_seq = None
        self.features = None
        self.num_conserved = None
        self.num_unconserved = None
        self.num_paired = None
        self.num_unpaired = None
        self.num_in_extrapolated_fiveprime_feature = None
        self.has_complete_feature_set = None
        self.num_extra_fiveprime = None
        self.feature_names = None
        self.threeprime_terminus_seq = None
        self.alpha_start = None
        self.alpha_stop = None
        self.beta_start = None
        self.beta_stop = None
        self.anticodon_seq = None
        self.anticodon_aa = None
        self.is_predicted_trna = None
        self.trunc_profile_index = None
        self.unconserved_info = None
        self.unpaired_info = None
        self.is_fully_profiled = None


    def get_unconserved_positions(self):
        """Get information on the unexpectedly unconserved positions in a tRNA feature profile.

        Returns
        =======
        unconserved_info : list
            List of tuples, one for each unconserved nucleotide in the profile.
            A tuple has three elements:
                1. position in the input sequence
                2. observed nucleotide in the input
                3. expected canonical nucleotides at the site,
                   a string with nucleotides separated by commas if more than one
        """
        unconserved_info = []
        for feature in self.features:
            # Only Nucleotide and Sequence subclasses have the attribute.
            if hasattr(feature, 'conserved_status'):
                component_start_pos = feature.start_pos
                # Conserved nucleotides are indexed within the string "component" (substring).
                for string_component_statuses, string_component in zip(feature.conserved_status, feature.string_components):
                    for nt_pos, is_conserved, observed_nt, expected_nts in string_component_statuses:
                        # Avoid N padding in an extrapolated 5' feature.
                        if not is_conserved and observed_nt != 'N':
                            unconserved_info.append((component_start_pos + nt_pos,
                                                     observed_nt,
                                                     ','.join(expected_nts)))
                    component_start_pos += len(string_component)
        return unconserved_info


    def get_unpaired_positions(self):
        """Get information on the unexpectedly unpaired bases in a tRNA feature profile.

        Returns
        =======
        unpaired_info : list
            List of tuples, one for each unpaired pair of nucleotides.
            A tuple has four elements:
                1. position of 5' nucleotide in the input sequence
                2. position of 3' nucleotide in the input sequence
                3. observed 5' nucleotide
                4. observed 3' nucleotide
        """
        unpaired_info = []
        for feature in self.features:
            # Only the Stem class has the attribute.
            if hasattr(feature, 'paired_status'):
                for nt_pos, (is_paired, fiveprime_nt, threeprime_nt) in enumerate(feature.paired_status):
                    # Avoid N padding in an extrapolated 5' feature.
                    if not is_paired and fiveprime_nt != 'N':
                        unpaired_info.append((feature.fiveprime_seq.start_pos + nt_pos,
                                              feature.threeprime_seq.stop_pos - nt_pos - 1,
                                              fiveprime_nt,
                                              threeprime_nt))
        return unpaired_info


class GeneProfile(object):
    """A profile of the tRNA features in a gene sequence.

    A gene may or may not encode the 3'-CCA acceptor sequence. The function,
    `Profiler.profile_gene`, creates these objects.
    """

    __slots__ = (
        'input_seq',
        'name',
        'unencoded_acceptor_profile',
        'encoded_acceptor_profile',
        'predicted_profile',
        'has_encoded_acceptor'
    )

    def __init__(self, input_seq, name=''):
        self.input_seq = input_seq
        self.name = name

        # The remaining attributes are set in `Profiler.profile_gene`.
        self.unencoded_acceptor_profile = None
        self.encoded_acceptor_profile = None
        self.predicted_profile = None
        self.has_encoded_acceptor = None


class Profiler(object):
    """Creates tRNA feature profiles from full-length or fragmentary transcript or gene DNA
    sequences."""

    def __init__(self):
        self.fiveprime_to_threeprime_feature_classes = TRNAFeature.list_all_tRNA_features()
        self.threeprime_to_fiveprime_feature_classes = self.fiveprime_to_threeprime_feature_classes[::-1]

        self.primary_feature_names = [feature_class.name for feature_class in TRNAFeature.list_primary_tRNA_features()]
        # The last of the "summed input lengths" possible for a feature is the maximum possible
        # length of the feature.
        self.max_length_dict = OrderedDict([(feature_class.name, feature_class.summed_input_lengths[-1])
                                            if feature_class.name != '3\' terminus'
                                            else (feature_class.name, max(map(len, [pattern.pattern for pattern in feature_class.allowed_patterns])))
                                            for feature_class in TRNAFeature.list_primary_tRNA_features()])

        self.stem_formation_triggers = [
            TStemFiveprimeStrand,
            AnticodonStemFiveprimeStrand,
            DStemFiveprimeStrand,
            AcceptorStemFiveprimeStrand
        ]
        self.arm_formation_triggers = [
            TStem,
            AnticodonStem,
            DStem
        ]

        self.d_loop_pos = self.threeprime_to_fiveprime_feature_classes.index(DLoop)
        self.d_stem_pos = self.threeprime_to_fiveprime_feature_classes.index(DStem)
        self.anticodon_loop_pos = self.threeprime_to_fiveprime_feature_classes.index(AnticodonLoop)
        self.v_loop_pos = self.threeprime_to_fiveprime_feature_classes.index(VLoop)
        self.t_arm_pos = self.threeprime_to_fiveprime_feature_classes.index(TArm)
        self.threeprime_stem_seq_positions = {
            TStem: self.threeprime_to_fiveprime_feature_classes.index(TStemThreeprimeStrand),
            AnticodonStem: self.threeprime_to_fiveprime_feature_classes.index(AnticodonStemThreeprimeStrand),
            DStem: self.threeprime_to_fiveprime_feature_classes.index(DStemThreeprimeStrand),
            AcceptorStem: self.threeprime_to_fiveprime_feature_classes.index(AcceptorStemThreeprimeStrand)
        }
        self.arm_loop_pos_dict = {
            TArm: self.threeprime_to_fiveprime_feature_classes.index(TLoop),
            AnticodonArm: self.threeprime_to_fiveprime_feature_classes.index(AnticodonLoop),
            DArm: self.d_loop_pos
        }

        self.extrapolation_ineligible_features = [
            AcceptorStemThreeprimeStrand,
            TStemThreeprimeStrand,
            VLoop,
            AnticodonStemThreeprimeStrand,
            DStemThreeprimeStrand
        ]


    def profile(self, input_seq, name=''):
        """Create a tRNA feature profile representing a tRNA transcript from a DNA sequence."""
        self.p = profile = Profile(input_seq, name=name)

        (profile.profiled_seq,
         profile.features,
         profile.num_conserved,
         profile.num_unconserved,
         profile.num_paired,
         profile.num_unpaired,
         profile.num_in_extrapolated_fiveprime_feature,
         profile.has_complete_feature_set,
         profile.num_extra_fiveprime) = self.get_profile(unprofiled_seq=profile.input_seq[::-1])

        profile.feature_names = [f.name for f in profile.features]

        if profile.features:
            if len(profile.features) > self.t_arm_pos:
                profile.is_predicted_trna = True # Liable to change below...

                profile.threeprime_terminus_seq = profile.features[-1].string

                if self.anticodon_loop_pos < len(profile.features):
                    anticodon = profile.features[-self.anticodon_loop_pos - 1].anticodon
                    profile.anticodon_seq = anticodon.string
                    profile.anticodon_aa = anticodon.aa_string

                    if self.d_loop_pos < len(profile.features):
                        # Explicitly record the start and stop positions within the input seq
                        # of the variable-length alpha and beta regions of the D loop.
                        D_loop = profile.features[-self.d_loop_pos - 1]
                        alpha_seq = D_loop.alpha_seq
                        beta_seq = D_loop.beta_seq
                        profile.alpha_start = D_loop.start_pos + alpha_seq.start_pos
                        profile.alpha_stop = D_loop.start_pos + alpha_seq.stop_pos
                        profile.beta_start = D_loop.start_pos + beta_seq.start_pos
                        profile.beta_stop = D_loop.start_pos + beta_seq.stop_pos
            else:
                profile.is_predicted_trna = False

            # tRNA profiling is complicated by the presence of "camouflaged" RNAs with structural
            # features resembling tRNA and chimeras (sequencing artifacts). The impact of
            # camouflaged RNAs is minimized by stringent limits on unconserved and unpaired bases in
            # the parameter settings. Chimeras can form between tRNA sequences and other RNA
            # sequences -- especially tRNA and, because of its abundance, rRNA. In `anvi-trnaseq`,
            # chimeras of abundant tRNA species can recruit large numbers of other sequences during
            # normalized and modified sequence formation, as chimeras are often among the longest
            # sequences, which are "favored" in these processes. In turn, chimeric normalized and
            # modified sequences become tRNA seeds (contigs) in `anvi-merge-trnaseq`. Chimeric
            # sequences containing a full-length tRNA at the 3' end -- which appear to be relatively
            # rare -- are difficult to distinguish from pre-tRNA, as both have a long 5' extension.
            # These chimeras do not impact `anvi-trnaseq` analysis, as bases 5' of the full-length
            # tRNA are trimmed off. Problematic chimeric sequences containing a fragmentary tRNA at
            # the 3' end are flagged by comparing *the number of unprofiled 5' bases* to *the
            # maximum length of the next unprofiled feature in the incomplete feature profile*; when
            # the former is greater than the latter, the sequence is labeled as not being predicted
            # tRNA, as the remaining 5' bases cannot be explained as the 3' part of a feature that
            # cannot be identified. For example, take a sequence of length 85 in which profiling
            # ends at position 27, which happens to be the 5' end of the D loop, leaving 26
            # unprofiled bases. The maximum length of the next unprofiled feature, the 5' strand of
            # the D stem, is 4 (for the sake of simplicity, this is not changed if the 3' strand of
            # the D stem was found to only have a length of 3). Therefore, the sequence is predicted
            # to not be tRNA.
            if not profile.has_complete_feature_set and profile.is_predicted_trna:
                for feature_name in profile.feature_names:
                    try:
                        primary_feature_pos = self.primary_feature_names.index(feature_name)
                    except ValueError:
                        continue
                    break
                else:
                    raise TRNAIdentifierError("To reach this point, the tRNA feature profile should contain primary sequence features, "
                                              "but that doesn't appear to be the case. "
                                              "Here are the names of the profiled features: %s" % ", ".join(profile.feature_names))
                if primary_feature_pos > 1:
                    # If the last profiled primary feature was the 5' strand of the acceptor stem,
                    # do not consider tRNA-His position 0 as the next unprofiled feature even if
                    # dealing with tRNA-His.
                    if len(input_seq) - len(profile.profiled_seq) > self.max_length_dict[self.primary_feature_names[primary_feature_pos - 1]]:
                        profile.is_predicted_trna = False
                        profile.trunc_profile_index = len(input_seq) - len(profile.profiled_seq)
        else:
            profile.anticodon_seq = ''
            profile.anticodon_aa = ''

            profile.is_predicted_trna = False

        if profile.is_predicted_trna:
            profile.unconserved_info = profile.get_unconserved_positions()
            profile.unpaired_info = profile.get_unpaired_positions()

        profile.is_fully_profiled = (profile.input_seq == profile.profiled_seq)
        return profile


    def profile_gene(self, input_seq, name='', check_encoded_acceptor=True):
        """Create a tRNA feature profile representing a tRNA gene from a DNA sequence."""
        gene_profile = GeneProfile(input_seq, name=name)
        gene_profile.unencoded_acceptor_profile = self.profile(input_seq + 'CCA', name=name)

        if check_encoded_acceptor:
            if input_seq[-3: ] == 'CCA':
                gene_profile.encoded_acceptor_profile = self.profile(input_seq, name=name)

                if gene_profile.unencoded_acceptor_profile.is_predicted_trna and not gene_profile.encoded_acceptor_profile.is_predicted_trna:
                    gene_profile.predicted_profile = gene_profile.unencoded_acceptor_profile
                    gene_profile.has_encoded_acceptor = False
                elif gene_profile.encoded_acceptor_profile.is_predicted_trna and not gene_profile.unencoded_acceptor_profile.is_predicted_trna:
                    gene_profile.predicted_profile = gene_profile.encoded_acceptor_profile
                    gene_profile.has_encoded_acceptor = True

                return gene_profile
            else:
                if gene_profile.unencoded_acceptor_profile.is_predicted_trna:
                    gene_profile.predicted_profile = gene_profile.unencoded_acceptor_profile
                    gene_profile.has_encoded_acceptor = False

        if gene_profile.unencoded_acceptor_profile.is_predicted_trna:
            gene_profile.predicted_profile = gene_profile.unencoded_acceptor_profile
            gene_profile.has_encoded_acceptor = False
        return gene_profile


    def get_profile(self,
                    unprofiled_seq='', # string in 3' to 5' direction
                    profiled_seq='', # string in 5' to 3' direction
                    features=None, # listed in the 5' to 3' direction
                    num_conserved=0,
                    num_unconserved=0,
                    num_paired=0,
                    num_unpaired=0,
                    feature_pos=0,
                    has_complete_feature_set=False):
        """A recursive function to sequentially identify tRNA features from the 3' end of a
        sequence."""
        if features is None:
            features = []

        if feature_pos == len(self.threeprime_to_fiveprime_feature_classes):
            # To reach this point, all tRNA features including the tRNA-His 5'-G must have been
            # found, and the input sequence must extend 5' of that.
            return (profiled_seq, # string in 5' to 3' direction
                    features, # listed in the 5' to 3' direction
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    0, # no extrapolated 5' nucleotides
                    has_complete_feature_set,
                    len(unprofiled_seq)) # extra 5' nucleotides

        if not unprofiled_seq:
            # To reach this point, the full length of the input sequence must have been profiled
            # with tRNA features.
            return (profiled_seq,
                    features, # listed in the 5' to 3' direction
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    0, # no extrapolated 5' nucleotides
                    has_complete_feature_set,
                    0) # input sequence is not longer than full-length tRNA


        feature_class = self.threeprime_to_fiveprime_feature_classes[feature_pos]
        if feature_class in self.stem_formation_triggers: # 3' stem seq triggers stem formation
            # Prepare to form a stem as well as the 3' sequence.
            make_stem = True
            stem_class = feature_class.stem_class
            threeprime_stem_seq = features[-self.threeprime_stem_seq_positions[stem_class] - 1]
            if stem_class in self.arm_formation_triggers:
                # Prepare to form an arm (stem + loop) as well as the stem.
                make_arm = True
                arm_class = stem_class.arm_class
                loop = features[-self.arm_loop_pos_dict[arm_class] - 1]
            else:
                make_arm = False
        else:
            make_stem = False
            make_arm = False
            if feature_class == TRNAHisPositionZero:
                # Check for tRNA-His based on the anticodon sequence. tRNA-His uniquely has an extra
                # nucleotide (G) at the 5' end.
                anticodon_string = features[-self.anticodon_loop_pos - 1].anticodon.string
                try:
                    aa_string = ANTICODON_TO_AA[anticodon_string]
                except KeyError:
                    aa_string = 'NA'
                if aa_string != 'His':
                    # The input sequence is longer than full-length tRNA.
                    return (profiled_seq,
                            features, # listed in the 5' to 3' direction
                            num_conserved,
                            num_unconserved,
                            num_paired,
                            num_unpaired,
                            0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                            has_complete_feature_set,
                            len(unprofiled_seq)) # extra 5' nucleotides


        # This list stores the result of a recursive function call finding subsequent 5' features.
        incremental_profile_candidates = []

        if feature_class == ThreeprimeTerminus:
            processed_lengths = []
            for priority, terminus_pattern in enumerate(feature_class.allowed_patterns):
                terminus_length = len(terminus_pattern.pattern)
                if terminus_length in processed_lengths:
                    continue

                fiveprime_to_threeprime_string = unprofiled_seq[: terminus_length][::-1]
                if not terminus_pattern.fullmatch(fiveprime_to_threeprime_string):
                    continue

                feature = feature_class(fiveprime_to_threeprime_string,
                                        start_pos=len(self.p.input_seq) - terminus_length,
                                        stop_pos=len(self.p.input_seq),
                                        priority=priority)
                processed_lengths.append(terminus_length)
                incremental_profile_candidates.append(
                    (unprofiled_seq[: terminus_length][::-1], # flip orientation 5' to 3'
                     [feature],
                     0, # conserved nucleotides
                     0, # unconserved nucleotides
                     0, # not a stem, so 0 paired nucleotides
                     0, # not a stem, so 0 unpaired nucleotides
                     0) # no extrapolated 5' nucleotides
                )
        else:
            # Primary sequence features take (sub)sequence inputs, which can be of varying length.
            # Consider each possible combination of input lengths for the feature, e.g., the D loop
            # contains alpha and beta subsequences of variable length.
            for input_lengths, summed_input_length in zip(feature_class.allowed_input_lengths, feature_class.summed_input_lengths):

                # Strands of unequal length cannot form a stem.
                if make_stem:
                    # Compare the lengths of the 5' and 3' sequence components.
                    if input_lengths != tuple(map(len, threeprime_stem_seq.string_components[::-1])):
                        continue

                if len(unprofiled_seq) < summed_input_length:
                    # Determine whether there is enough information in the remaining 5' end of the
                    # input sequence to assign it to a feature despite the incompleteness of the
                    # feature sequence.

                    # Features that lack conserved positions or that do not form base pairs with a
                    # previously profiled stem sequence do not contain any information for
                    # extrapolation of an incomplete sequence.
                    if feature_class in self.extrapolation_ineligible_features:
                        continue
                    # The unprofiled sequence must be at least 6 nucleotides long to span the 2
                    # conserved positions in the anticodon loop (positions 33 and 37). Conservation
                    # of these positions is here considered to be the minimum information needed to
                    # identify the anticodon loop.
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
                    # WITH THE N PADDING, THE START POSITION (5') OF THE FEATURE IN THE INPUT SEQUENCE IS NEGATIVE.
                    feature = feature_class(*string_components,
                                            start_pos=len(self.p.input_seq) - len(profiled_seq) - num_processed_bases,
                                            stop_pos=len(self.p.input_seq) - len(profiled_seq))

                    # The sequence is valid if it doesn't have too many unconserved bases.
                    if feature.meets_conserved_thresh:
                        if make_stem:
                            if stem_class == DStem:
                                stem = stem_class(feature, threeprime_stem_seq, type_II_trna=features[-self.v_loop_pos - 1].type == 'II')
                            else:
                                stem = stem_class(feature, threeprime_stem_seq)
                            # The stem is valid if it doesn't have too many unpaired bases.
                            if stem.meets_pair_thresh:
                                if make_arm:
                                    arm = arm_class(stem, loop)
                                    # The arm is valid if it doesn't have too many unconserved bases.
                                    if arm.meets_conserved_thresh:
                                        incremental_profile_candidates.append(
                                            (unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                             [arm, stem, feature],
                                             feature.num_conserved,
                                             feature.num_unconserved,
                                             stem.num_paired,
                                             stem.num_unpaired,
                                             summed_input_length - len(unprofiled_seq)) # number of nucleotides in extrapolated 5' feature
                                        )
                                        continue
                                else:
                                    incremental_profile_candidates.append(
                                        (unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                         [stem, feature],
                                         feature.num_conserved,
                                         feature.num_unconserved,
                                         stem.num_paired,
                                         stem.num_unpaired,
                                         summed_input_length - len(unprofiled_seq)) # number of nucleotides in extrapolated 5' feature
                                    )
                                    continue
                        else:
                            incremental_profile_candidates.append(
                                (unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                 [feature],
                                 feature.num_conserved,
                                 feature.num_unconserved,
                                 0, # not a stem, so 0 paired nucleotides
                                 0, # not a stem, so 0 unpaired nucleotides
                                 summed_input_length - len(unprofiled_seq)) # number of nucleotides in extrapolated 5' feature
                            )
                            continue
                # The procedure for assigning full-length features is similar to the prior precedure
                # for partial-length features with a few efficiencies.
                else:
                    string_components = [] # feature input substrings
                    num_processed_bases = 0
                    for input_length in input_lengths[::-1]: # create substrings from 3' to 5'
                        string_components.insert(0, unprofiled_seq[num_processed_bases: num_processed_bases + input_length][::-1]) # flip orientation to 5' to 3'
                        num_processed_bases += input_length
                    feature = feature_class(*string_components,
                                            start_pos=len(self.p.input_seq) - len(profiled_seq) - num_processed_bases,
                                            stop_pos=len(self.p.input_seq) - len(profiled_seq))

                    # The sequence is valid if it doesn't have too many unconserved bases.
                    if feature.meets_conserved_thresh:
                        if make_stem:
                            if stem_class == DStem:
                                stem = stem_class(feature, threeprime_stem_seq, type_II_trna=features[-self.v_loop_pos - 1].type == 'II')
                            else:
                                stem = stem_class(feature, threeprime_stem_seq)
                            # The stem is valid if it doesn't have too many unpaired bases.
                            if stem.meets_pair_thresh:
                                if make_arm:
                                    arm = arm_class(stem, loop)
                                    # The arm is valid if it doesn't have too many unconserved bases.
                                    if arm.meets_conserved_thresh:
                                        incremental_profile_candidates.append(
                                            (unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                             [arm, stem, feature],
                                             feature.num_conserved,
                                             feature.num_unconserved,
                                             stem.num_paired,
                                             stem.num_unpaired,
                                             0) # no extrapolated 5' nucleotides
                                        )
                                        continue
                                else:
                                    incremental_profile_candidates.append(
                                        (unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                         [stem, feature],
                                         feature.num_conserved,
                                         feature.num_unconserved,
                                         stem.num_paired,
                                         stem.num_unpaired,
                                         0) # no extrapolated 5' nucleotides
                                    )
                                    continue
                        else:
                            incremental_profile_candidates.append(
                                (unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                [feature],
                                feature.num_conserved,
                                feature.num_unconserved,
                                0, # not a stem, so 0 paired nucleotides
                                0, # not a stem, so 0 unpaired nucleotides
                                0) # no extrapolated 5' nucleotides
                            )
                            continue


        if not incremental_profile_candidates:
            # The feature didn't pass muster.
            if feature_class == TRNAHisPositionZero:
                # The extra 5' G added post-transcriptionally to tRNA-His was not found. Therefore,
                # unprofiled 5' nucleotides are unaccounted for and called extra 5' nucleotides.
                return (profiled_seq,
                        features,
                        num_conserved,
                        num_unconserved,
                        num_paired,
                        num_unpaired,
                        0, # no extrapolated 5' nucleotides
                        has_complete_feature_set,
                        len(unprofiled_seq))
            else:
                return (profiled_seq,
                        features,
                        num_conserved,
                        num_unconserved,
                        num_paired,
                        num_unpaired,
                        0,
                        has_complete_feature_set,
                        0) # input sequence is not longer than full-length tRNA


        # Sort candidates by
        # 1. number of features identified (at most, sequence + stem + arm) (descending),
        # 2. number of unconserved + unpaired nucleotides (ascending),
        # 3. incompleteness of the last (most 5') feature (ascending).
        # This sort also happens later for full sequence profiles, but this first sort is useful for
        # seeking out and returning "flawless" mature tRNA.
        incremental_profile_candidates.sort(key=lambda p: (-len(p[1]), p[3] + p[5], p[6]))
        # Continue finding features in input sequences that have not been fully profiled --
        # do not recurse profile candidates in which the final feature was extrapolated.
        profile_candidates = []
        for ipc in incremental_profile_candidates:
            if ipc[6] == 0: # 5' feature was NOT extrapolated from unprofiled input sequence
                if has_complete_feature_set or feature_class.name == 'acceptor stem 5\' strand':
                    profile_candidate = self.get_profile(unprofiled_seq=unprofiled_seq[len(ipc[0]): ], # recurse
                                                         profiled_seq=ipc[0] + profiled_seq,
                                                         features=ipc[1] + features,
                                                         num_conserved=ipc[2] + num_conserved,
                                                         num_unconserved=ipc[3] + num_unconserved,
                                                         num_paired=ipc[4] + num_paired,
                                                         num_unpaired=ipc[5] + num_unpaired,
                                                         feature_pos=feature_pos + len(ipc[1]),
                                                         has_complete_feature_set=True)
                else:
                    profile_candidate = self.get_profile(unprofiled_seq=unprofiled_seq[len(ipc[0]): ], # recurse
                                                         profiled_seq=ipc[0] + profiled_seq,
                                                         features=ipc[1] + features,
                                                         num_conserved=ipc[2] + num_conserved,
                                                         num_unconserved=ipc[3] + num_unconserved,
                                                         num_paired=ipc[4] + num_paired,
                                                         num_unpaired=ipc[5] + num_unpaired,
                                                         feature_pos=feature_pos + len(ipc[1]),
                                                         has_complete_feature_set=False)
                if (profile_candidate[7] # has complete feature set
                    and profile_candidate[3] == 0 # number unconserved
                    and profile_candidate[5] == 0): # number unpaired
                    return profile_candidate
                else:
                    profile_candidates.append(profile_candidate)
            else: # 5' feature was extrapolated from unprofiled input sequence
                profile_candidates.append(
                    (ipc[0] + profiled_seq,
                     ipc[1] + features,
                     ipc[2] + num_conserved,
                     ipc[3] + num_unconserved,
                     ipc[4] + num_paired,
                     ipc[5] + num_unpaired,
                     ipc[6], # number of nucleotides in extrapolated 5' feature
                     False, # does not have a complete feature set
                     0) # input sequence is not longer than full-length tRNA
                )

        # Do not add 5' features to profiled 3' features if the additional features have no
        # grounding in conserved nucleotides or paired nucleotides in stems.
        supported_profile_candidates = []
        is_D_stem_in_profile = self.d_stem_pos < len(features)
        has_profile_with_type_II_D_arm = False
        for p in profile_candidates:
            if is_D_stem_in_profile:
                # Do not favor type II tRNA profiles with paired D stem positions 13 and 22 over
                # profiles with those positions unpaired.
                if features[-self.d_stem_pos - 1].paired_positions_13_22_in_type_II:
                    if p[2] + p[4] - 1 > num_conserved + num_paired:
                        supported_profile_candidates.append(p)
                        has_profile_with_type_II_D_arm = True
                        continue
            if len(p[1]) == 1:
                # The only feature in the profile candidate is the 3' terminus.
                supported_profile_candidates.append(p)
            elif p[2] + p[4] > num_conserved + num_paired:
                supported_profile_candidates.append(p)

        if supported_profile_candidates:
            if has_profile_with_type_II_D_arm:
                # Favor type II tRNA profiles with unpaired D stem positions 13 and 22 over
                # profiles with the positions paired.
                unpaired_positions_13_22_in_type_II_profiles = [
                    1 if p[1][-self.d_stem_pos - 1].unpaired_positions_13_22_in_type_II else 0
                    for p in supported_profile_candidates
                ]
                return sorted(zip(supported_profile_candidates, unpaired_positions_13_22_in_type_II_profiles),
                              key=lambda t: (-len(t[0][1]), t[0][3] + t[0][5] - t[1], t[0][6], t[0][1][-1].priority))[0][0]
            else:
                # Favor the profile with more features, fewer unconserved and unpaired nucleotides,
                # fewer extrapolated 5' nucleotides, and greater priority of 3' terminus.
                return sorted(supported_profile_candidates, key=lambda p: (-len(p[1]), p[3] + p[5], p[6], p[1][-1].priority))[0]
        else:
            return (profiled_seq,
                    features,
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    0, # no extrapolated 5' nucleotides
                    has_complete_feature_set,
                    0) # input sequence is not longer than full-length tRNA
