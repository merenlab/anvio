# -*- coding: utf-8
# pylint: disable=line-too-long
"""tRNA identification from input sequence"""

import itertools
import pandas as pd

from tabulate import tabulate

from anvio.errors import TRNAIdentifierError
from anvio.filesnpaths import is_file_exists, is_output_file_writable
from anvio.constants import LONGEST_KNOWN_TRNA_LENGTH
from anvio.constants import anticodon_to_AA as ANTICODON_TO_AA
from anvio.constants import WC_plus_wobble_base_pairs as WC_PLUS_WOBBLE_BASE_PAIRS


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


class TRNAFeature:
    # Full tRNA feature parameter names in .ini file
    INI_FEATURE_PARAMS = ['Conserved nucleotides', 'Number allowed unconserved', 'Number allowed unpaired', 'Allowed lengths']
    # Feature parameter attributes
    ACCESSIBLE_FEATURE_PARAMS = ['conserved_nucs', 'num_allowed_unconserved', 'num_allowed_unpaired', 'allowed_input_lengths']
    # These class attributes are set by `set_feature_param_class_attributes`
    # or are set the first time `set_feature_params` or `print_default_feature_params` is called.
    # Since the attributes reference other classes, they can't be assigned yet.
    dict_mapping_feature_or_subfeature_name_to_class = None
    feature_and_subfeature_names_with_accessible_lengths = None
    subfeature_section_dict = None

    def __init__(self,
                 string_components,
                 conserved_nucs=None,
                 num_allowed_unconserved=-1,
                 cautious=False):
        """The most ancestral superclass of tRNA features

        Parameters
        ==========
        string_components : tuple
            Substring components of the feature, with divisions based on "subfeatures,"
            such as the variable-length alpha and beta regions of the D loop:
            e.g., ('AG', 'TT', 'GG', 'G', 'A') for the D loop of yeast tRNA-Phe-GAA,
                  ('CCA', ) for the acceptor sequence

        num_allowed_unconserved : int, -1
            Number of unconserved nucleotides allowed among canonically conserved nucleotides in the feature
            The default of -1 is equivalent to allowing all to be unconserved.

        cautious : bool, False
        """

        if cautious:
            if type(string_components) != tuple:
                raise TRNAIdentifierError("`string_components` must be in the form of a tuple, e.g., ('ACTGG', 'CCAGT'). "
                                          "Your `string_components` were %s" % (string_components, ))
        self.string_components = string_components

        self.nuc_count = sum(map(len, string_components))

        if conserved_nucs is None:
            self.conserved_nucs = ({}, )
        # By default, base conservation is not enforced.
        if num_allowed_unconserved == -1:
            self.num_allowed_unconserved = sum(len(d) for d in self.conserved_nucs)
        else:
            self.num_allowed_unconserved = num_allowed_unconserved


    def check_conserved_nucs(self):
        """Determine which of the canonical nucleotides in the feature are conserved.

        Returns
        =======
        meets_conserved_thresh : bool

        num_conserved : int

        num_unconserved : int

        conserved_status : list
            Nested list with the following structure:
            [[(), (), ...], [(), (), ...], ...]
            There is an inner list for each subsequence (substring component)
            and an inner tuple for each conserved nucleotide of the subsequence.
            Each inner tuple has four elements:
                1. position (index) of conserved nucleotide in subsequence
                2. whether nucleotide is conserved (bool)
                3. observed nucleotide in subsequence (char)
                4. expected canonical nucleotide in subsequence (char)
        """

        num_conserved = 0
        num_unconserved = 0 # can include N "padding" in extrapolated 5' feature
        conserved_status = []
        for substring, nuc_dict in zip(self.string_components, self.conserved_nucs):
            substring_statuses = []
            conserved_status.append(substring_statuses)
            for pos, expected_nucs in nuc_dict.items():
                try:
                    observed_nuc = substring[pos]
                except IndexError:
                    # This occurs for an Acceptor feature string that differs from the canonical "CCA", such as "CC", "C", or "CCATT".
                    break

                if observed_nuc in expected_nucs:
                    num_conserved += 1
                    substring_statuses.append((pos, True, observed_nuc, expected_nucs))
                else:
                    num_unconserved += 1
                    substring_statuses.append((pos, False, observed_nuc, expected_nucs))

        if num_unconserved > self.num_allowed_unconserved:
            meets_conserved_thresh = False
        else:
            meets_conserved_thresh = True

        return (meets_conserved_thresh, num_conserved, num_unconserved, conserved_status)


    @staticmethod
    def list_all_tRNA_features():
        """List all tRNA feature classes in order from 5' to 3'"""

        return [TRNAHisPositionZero,
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
                Acceptor]


    @staticmethod
    def set_param_refs():
        """Set class attributes that reference feature classes and thereby cannot be set during class definition."""

        TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class = {
            TRNAHisPositionZero.name: TRNAHisPositionZero,
            AcceptorStem.name: AcceptorStem,
            AcceptorStemFiveprimeStrand.name: AcceptorStemFiveprimeStrand,
            PositionEight.name: PositionEight,
            PositionNine.name: PositionNine,
            DArm.name: DArm,
            DStem.name: DStem,
            DStemFiveprimeStrand.name: DStemFiveprimeStrand,
            DStemFiveprimeStrand.name + "/positions 10-12": DStemFiveprimeStrand,
            DStemFiveprimeStrand.name + "/position 13": DStemFiveprimeStrand,
            DLoop.name: DLoop,
            DLoop.name + "/positions 14-15": DLoop,
            DLoop.name + "/alpha positions": DLoop,
            DLoop.name + "/positions 18-19": DLoop,
            DLoop.name + "/beta positions": DLoop,
            DLoop.name + "/position 21": DLoop,
            DStemThreeprimeStrand.name: DStemThreeprimeStrand,
            DStemThreeprimeStrand.name + "/position 22": DStemThreeprimeStrand,
            DStemThreeprimeStrand.name + "/positions 23-25": DStemThreeprimeStrand,
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
            Discriminator.name: Discriminator,
            Acceptor.name: Acceptor}

        TRNAFeature.feature_and_subfeature_names_with_accessible_lengths = [
            VLoop.name,
            DStemFiveprimeStrand.name + "/position 13",
            DLoop.name + "/alpha positions",
            DLoop.name + "/beta positions",
            DStemThreeprimeStrand.name + "/position 22"]

        TRNAFeature.subfeature_section_dict = {name: 0 for name
                                               in TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class}
        TRNAFeature.subfeature_section_dict[DStemFiveprimeStrand.name + "/position 13"] = 1
        TRNAFeature.subfeature_section_dict[DLoop.name + "/alpha positions"] = 1
        TRNAFeature.subfeature_section_dict[DLoop.name + "/positions 18-19"] = 2
        TRNAFeature.subfeature_section_dict[DLoop.name + "/beta positions"] = 3
        TRNAFeature.subfeature_section_dict[DLoop.name + "/position 21"] = 4
        TRNAFeature.subfeature_section_dict[DStemThreeprimeStrand.name + "/positions 23-25"] = 1


    @staticmethod
    def write_default_param_file(default_feature_param_path):
        """Write the default table of feature and "subfeature" parameters (class attributes)

        Parameters
        ==========
        default_feature_param_path : str
            Output file path
        """

        is_output_file_writable(default_feature_param_path)
        rows = TRNAFeature.get_default_param_table_as_list()
        with open(default_feature_param_path, 'w') as f:
            for row in rows:
                f.write("\t".join(row) + "\n")


    @staticmethod
    def print_default_params():
        """Print a nicely formatted version of the default table of feature and "subfeature" parameters (class attributes)"""

        rows = TRNAFeature.get_default_param_table_as_list()
        print(tabulate(rows, headers='firstrow', tablefmt='github'))


    @staticmethod
    def get_default_param_table_as_list():
        """Get a table of default user-accessible feature and "subfeature" parameters (class attributes)

        A dash indicates a parameter that cannot be set for the feature.
        Quotes indicate a parameter that is currently not set.
        """

        ACCESSIBLE_FEATURE_PARAMS = TRNAFeature.ACCESSIBLE_FEATURE_PARAMS
        if (TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class is None
            or TRNAFeature.subfeature_section_dict is None
            or TRNAFeature.feature_and_subfeature_names_with_accessible_lengths is None):
            TRNAFeature.set_param_refs()
        dict_mapping_feature_or_subfeature_name_to_class = TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class
        feature_and_subfeature_names_with_accessible_lengths = TRNAFeature.feature_and_subfeature_names_with_accessible_lengths
        subfeature_section_dict = TRNAFeature.subfeature_section_dict

        rows = [["Feature/subfeature"] + TRNAFeature.INI_FEATURE_PARAMS]
        for feature_or_subfeature_name, feature_class in dict_mapping_feature_or_subfeature_name_to_class.items():
            row = [feature_or_subfeature_name]

            if feature_or_subfeature_name == 'acceptor': # The acceptor is special.
                row += ['-' for _ in ACCESSIBLE_FEATURE_PARAMS]
                rows.append(row)
                continue

            for param_name in ACCESSIBLE_FEATURE_PARAMS:
                if param_name == 'conserved_nucs':
                    if ('/' not in feature_or_subfeature_name
                        and 'allowed_section_lengths' in feature_class.__dict__):
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

                if param_name == 'conserved_nucs':
                    conserved_nuc_list = []
                    if not param_value[subfeature_section_dict[feature_or_subfeature_name]]:
                        row.append("\"\"")
                        continue
                    for pos, nuc in param_value[subfeature_section_dict[feature_or_subfeature_name]].items():
                        if type(nuc) == str:
                            conserved_nuc_list.append(str(pos) + "," + nuc)
                        elif type(nuc) == tuple:
                            # The only alternatives to a single conserved nucleotide are a purine or pyrimidine.
                            if nuc == ('A', 'G'):
                                conserved_nuc_list.append(str(pos) + ",R")
                            elif nuc == ('C', 'T'):
                                conserved_nuc_list.append(str(pos) + ",Y")
                    row.append(";".join(conserved_nuc_list))
                elif param_name == 'num_allowed_unconserved' or param_name == 'num_allowed_unpaired':
                    if param_value == -1:
                        row.append("\"\"")
                        continue
                    row.append(str(param_value))
                elif param_name == 'allowed_input_lengths':
                    if feature_or_subfeature_name == VLoop.name:
                        # Here is an example to show the format of VLoop.allowed_input_lengths: ((4, ), (5, ), ..., (23, ))
                        row.append(str(param_value[0][0]) + "-" + str(param_value[-1][0]))
                    else:
                        # Sections of the D loop are the only other variable-length "subfeatures".
                        allowed_section_lengths = feature_class.allowed_section_lengths[subfeature_section_dict[feature_or_subfeature_name]]
                        row.append(str(allowed_section_lengths[0]) + "-" + str(allowed_section_lengths[-1]))
            rows.append(row)

        return rows


    @staticmethod
    def set_params_from_file(feature_param_path):
        """Set user-accessible tRNA feature parameters for de novo tRNA profiling and identification"""

        is_file_exists(feature_param_path)
        feature_param_df = pd.read_csv(feature_param_path, sep='\t', header=0, index_col=0, keep_default_na=False)
        for feature_or_subfeature_name, row in feature_param_df.iterrows():
            for ini_param_name, param_name in zip(TRNAFeature.INI_FEATURE_PARAMS, TRNAFeature.ACCESSIBLE_FEATURE_PARAMS):
                param_value = row[ini_param_name]
                if param_value == '-':
                    continue
                TRNAFeature.set_param(feature_or_subfeature_name, param_name, param_value)


    @staticmethod
    def set_param(feature_or_subfeature_name, param_name, param_value):
        """Set tRNA feature parameter class attributes

        Parameters
        ==========
        feature_or_subfeature_name : str
            The name of a feature or "subfeature" (in the case of a section of the D arm)

        param_name : str
            The name of a parameter that can be set

        param_value : str
            The specifically formatted value of the parameter to be set
        """

        if (TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class is None
            or TRNAFeature.feature_and_subfeature_names_with_accessible_lengths is None
            or TRNAFeature.subfeature_section_dict is None):
            TRNAFeature.set_param_refs()

        if feature_or_subfeature_name not in TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class:
            raise TRNAIdentifierError("`TRNAFeature.set_param` does not recognize the supplied feature or subfeature name, %s. "
                                      "Here are the recognized feature names: %s"
                                      % (feature_or_subfeature_name, ", ".join(TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class)))

        if param_name not in TRNAFeature.ACCESSIBLE_FEATURE_PARAMS:
            raise TRNAIdentifierError("`TRNAFeature.set_param` does not recognize the supplied feature parameter name, %s. "
                                      "Here are the recognized parameter names: %s"
                                      % (param_name, ", ".join(TRNAFeature.ACCESSIBLE_FEATURE_PARAMS)))

        if param_name == 'conserved_nucs':
            TRNAFeature.set_conserved_nucs(feature_or_subfeature_name, param_value)
        elif param_name == 'num_allowed_unconserved':
            if not param_value:
                param_value = -1
            TRNAFeature.set_num_allowed_unconserved(feature_or_subfeature_name, param_value)
        elif param_name == 'num_allowed_unpaired':
            if not param_value:
                param_value = 1000
            TRNAFeature.set_num_allowed_unpaired(feature_or_subfeature_name, param_value)
        elif param_name == 'allowed_input_lengths':
            TRNAFeature.set_allowed_input_lengths(feature_or_subfeature_name, param_value)


    @staticmethod
    def set_conserved_nucs(feature_or_subfeature_name, param_value=''):
        """Modify (update a dict in) the `conserved_nucs` attribute of a tRNA feature class

        Parameters
        ==========
        feature_or_subfeature_name : str
            The name of a feature or "subfeature" (in the case of a section of the D arm)

        param_value : str, ''
            The proper format is
            <Zero-based position relative to the 5' end of the feature>,<Conserved nucleotide symbol>;<Next index>,<Next symbol>;...
            The default empty string indicates that there are no conserved sites.
        """

        feature_class = TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class[feature_or_subfeature_name]

        conserved_nucs_dict = {}
        if param_value:
            try:
                for conserved_nuc in param_value.split(';'):
                    conserved_pos, nuc = conserved_nuc.split(',')
                    if nuc == 'R':
                        nuc = ('A', 'G')
                    elif nuc == 'Y':
                        nuc = ('C', 'T')
                    conserved_nucs_dict[int(conserved_pos)] = nuc
            except:
                raise TRNAIdentifierError("The proper format of a conserved nucleotide parameter value is "
                                          "\"<Zero-based position of conserved nucleotide in feature(/subfeature) relative to 5' end of feature>,"
                                          "<Conserved nucleotide symbol (A, C, G, T, R, or Y)>;"
                                          "<Next zero-based position>,<Next conserved nucleotide symbol>;...\" "
                                          "Note that the position integer is separated from the nucleotide character by a comma. "
                                          "Note that entries for different conserved nucleotides are separated by a semicolon. "
                                          "The value provided was %s" % param_value)

            if ('conserved_nucs' not in feature_class.__dict__
                or ('/' not in feature_or_subfeature_name
                    and 'allowed_section_lengths' in feature_class.__dict__)):
                #   Conserved nucleotides can only be set for individual "subfeature" sections of the D arm.
                raise TRNAIdentifierError("\"%s\" does not support assignment of conserved nucleotides. "
                                          "The conserved nucleotide input provided was %s"
                                          % (feature_or_subfeature_name, param_value))

        prev_conserved_nucs_dict = feature_class.conserved_nucs[TRNAFeature.subfeature_section_dict[feature_or_subfeature_name]]
        prev_conserved_nucs_dict.clear()
        prev_conserved_nucs_dict.update(conserved_nucs_dict)


    @staticmethod
    def set_num_allowed_unconserved(feature_name, num_allowed_unconserved=-1):
        """Set the `num_allowed_unconserved` attribute of a tRNA feature class

        Parameters
        ==========
        feature_name : str
            The name of a feature

        num_allowed_unconserved : int, -1
            The number of conserved nucleotide positions in the feature allowed to be unconserved
            but still have positive identification of the feature
            The default value of -1 means that an "unlimited" number of unconserved positions is allowed.
        """

        feature_class = TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class[feature_name]

        num_allowed_unconserved = int(num_allowed_unconserved)
        if ('num_allowed_unconserved' not in feature_class.__dict__
            or '/' in feature_name):
            raise TRNAIdentifierError("\"%s\" does not support assignment of an allowed number of unconserved nucleotides. "
                                      "The number of allowed unconserved nucleotides provided was %s"
                                      % (feature_name, num_allowed_unconserved))
        feature_class.num_allowed_unconserved = num_allowed_unconserved


    @staticmethod
    def set_num_allowed_unpaired(stem_name, num_allowed_unpaired=1000):
        """Set the `num_allowed_unpaired` attribute of a stem class

        Parameters
        ==========
        stem_name : str
            The name of a stem feature

        param_value : int, 1000
            The number of positions in the stem allowed to be unpaired
            Pairing means Watson-Crick or wobble G-T.
            The default value of 1000 means that an "unlimited" number of unpaired positions is allowed.
        """

        feature_class = TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class[stem_name]

        num_allowed_unpaired = int(num_allowed_unpaired)
        if 'num_allowed_unpaired' not in feature_class.__dict__:
            raise TRNAIdentifierError("\"%s\" does not support assignment of an allowed number of unpaired nucleotides. "
                                      "The number of allowed unpaired nucleotides provided was %s"
                                      % (stem_name, num_allowed_unpaired))
        feature_class.num_allowed_unpaired = num_allowed_unpaired


    @staticmethod
    def set_allowed_input_lengths(feature_or_subfeature_name, param_value):
        """Set attributes related to the range of allowed lengths in a variable-length tRNA feature

        Parameters
        ==========
        feature_or_subfeature_name : str
            The name of a feature or "subfeature" (in the case of a section of the D arm)

        param_value : str
            A string representing the allowed length range, with the format <minimum length>-<maximum length>.
        """

        feature_class = TRNAFeature.dict_mapping_feature_or_subfeature_name_to_class[feature_or_subfeature_name]

        try:
            min_length, max_length = param_value.split('-')
            allowed_lengths = tuple(range(int(min_length), int(max_length) + 1))
        except:
            raise TRNAIdentifierError("The proper format of the allowed feature length field in a parameter string is "
                                      "<Minimum length integer>-<Maximum length integer>. "
                                      "The length range provided was %s" % param_value)
        if feature_or_subfeature_name not in TRNAFeature.feature_and_subfeature_names_with_accessible_lengths:
            raise TRNAIdentifierError("\"%s\" does not support assignment of variable lengths. "
                                      "The length range provided was %s"
                                      % (feature_or_subfeature_name, param_value))

        if feature_or_subfeature_name == VLoop.name:
            feature_class.allowed_input_lengths = tuple(itertools.product(allowed_lengths))
            # Reset the dependent class attribute.
            feature_class.summed_input_lengths = allowed_lengths
        else:
            # Beside the variable loop, subfeatures of the D arm have variable lengths.
            # If the possible lengths of position 13 (on the 5' strand of the D stem)
            # or position 22 (on the 3' of the D stem) are changed
            # then the other position should be changed accordingly -- this is not enforced.
            section = TRNAFeature.subfeature_section_dict[feature_or_subfeature_name]
            prev_allowed_section_lengths = feature_class.allowed_section_lengths
            feature_class.allowed_section_lengths = (prev_allowed_section_lengths[: section]
                                                     + (allowed_lengths, )
                                                     + prev_allowed_section_lengths[section + 1: ])
            # Reset the dependent class attributes.
            feature_class.allowed_input_lengths = tuple(itertools.product(*feature_class.allowed_section_lengths))
            feature_class.summed_input_lengths = tuple(map(sum, feature_class.allowed_input_lengths))


class Nucleotide(TRNAFeature):
    allowed_input_lengths = ((1, ), )
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))

    def __init__(self,
                 string, # must be a string of length 1
                 conserved_nucs=None,
                 num_allowed_unconserved=-1,
                 start_index=None,
                 stop_index=None,
                 cautious=False):
        """Superclass for tRNA primary features of a single nucleotide, e.g., discriminator, position 8"""

        self.string = string
        self.start_index = start_index
        self.stop_index = stop_index

        super().__init__((string, ),
                         conserved_nucs=conserved_nucs,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved,
         self.conserved_status) = self.check_conserved_nucs()


class Sequence(TRNAFeature):
    def __init__(self,
                 substrings, # must be a string, tuple of strings, or tuple of Nucleotide/Sequence objects
                 conserved_nucs=None,
                 num_allowed_unconserved=-1,
                 start_index=None,
                 stop_index=None,
                 cautious=False):
        """Superclass for tRNA primary sequences, e.g., acceptor, 5' strand of T stem"""

        if type(substrings) == str:
            string_components = (substrings, )
        elif all([type(s) == str for s in substrings]):
            string_components = substrings
        elif all([type(s) == Nucleotide or type(s) == Sequence for s in substrings]):
            string_components = tuple(s.string for s in substrings)
        self.string = ''.join(substrings)
        self.start_index = start_index
        self.stop_index = stop_index

        super().__init__(string_components,
                         conserved_nucs=conserved_nucs,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_conserved_thresh,
         self.num_conserved,
         self.num_unconserved, # can include N "padding" in extrapolated 5' feature
         self.conserved_status) = self.check_conserved_nucs()


class Loop(Sequence):
    def __init__(self,
                 substrings, # must be a string, tuple of strings, or tuple of Nucleotide/Sequence objects
                 conserved_nucs=None,
                 num_allowed_unconserved=-1,
                 start_index=None,
                 stop_index=None,
                 cautious=False):
        """Superclass for loops: D loop, anticodon loop, T loop"""

        super().__init__(substrings,
                         conserved_nucs=conserved_nucs,
                         num_allowed_unconserved=num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class Stem(TRNAFeature):
    def __init__(self,
                 fiveprime_seq, # must be Sequence object
                 threeprime_seq, # must be Sequence object
                 num_allowed_unpaired=0,
                 num_allowed_unconserved=-1,
                 cautious=False):
        """Superclass for stems: acceptor stem, D stem, anticodon stem, T stem"""

        if cautious:
            if type(fiveprime_seq) != Sequence or type(threeprime_seq) != Sequence:
                raise TRNAIdentifierError("You can only define a Stem from Sequence objects.")
        self.fiveprime_seq = fiveprime_seq
        self.threeprime_seq = threeprime_seq

        self.canonical_positions = (*self.fiveprime_seq.canonical_positions, *self.threeprime_seq.canonical_positions)
        self.conserved_nucs = (*self.fiveprime_seq.conserved_nucs, *self.threeprime_seq.conserved_nucs)

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

        self.start_indices = (self.fiveprime_seq.start_index, self.threeprime_seq.start_index)
        self.stop_indices = (self.fiveprime_seq.stop_index, self.threeprime_seq.stop_index)

        super().__init__((*self.fiveprime_seq.string_components, *self.threeprime_seq.string_components),
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        (self.meets_pair_thresh,
         self.num_paired,
         self.num_unpaired, # can include N "padding" in extrapolated 5' feature
         self.paired_status) = self.check_pairs()


    def check_pairs(self):
        """Determine base pairing in the stem

        Returns
        =======
        meets_pair_thresh : bool

        num_paired : int

        num_unpaired : int

        paired_status : list
            List of tuples, one tuple for each nucleotide pair in the stem
            Each tuple has three elements: whether a base pair exists, the 5' nucleotide, and the 3' nucleotide
        """

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


class Arm(TRNAFeature):
    def __init__(self,
                 stem, # must be Stem object
                 loop, # must be Loop object
                 num_allowed_unconserved=-1,
                 cautious=False):
        """Superclass for arms: D arm, anticodon arm, T arm

        The number of unconserved nucleotides allowed in the arm
        can differ from the sum of the numbers of unconserved nucleotides allowed in the stem and loop.
        """

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

        self.conserved_nucs=(*stem.fiveprime_seq.conserved_nucs,
                             *loop.conserved_nucs,
                             *stem.threeprime_seq.conserved_nucs)

        self.start_index = self.stem.start_indices[0]
        self.stop_index = self.stem.stop_indices[1]

        super().__init__((*stem.fiveprime_seq.string_components,
                          *loop.string_components,
                          *stem.threeprime_seq.string_components),
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=num_allowed_unconserved,
                         cautious=cautious)

        if (self.stem.fiveprime_seq.num_unconserved
            + self.loop.num_unconserved
            + self.stem.threeprime_seq.num_unconserved) > self.num_allowed_unconserved:
            self.meets_conserved_thresh = False
        else:
            self.meets_conserved_thresh = True


class TRNAHisPositionZero(Nucleotide):
    name = 'tRNA-His position 0'
    canonical_positions = ((-1, ), )
    conserved_nucs = ({0: 'G'}, )
    num_allowed_unconserved = 0

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        """The G typically found at the 5' end of mature tRNA-His"""

        super().__init__(string,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class AcceptorStem(Stem):
    name = 'acceptor stem'
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        """The base-paired nucleotides of the acceptor stem of tRNA (excludes the discriminator and acceptor)"""

        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class AcceptorStemFiveprimeStrand(Sequence):
    name = 'acceptor stem 5\' strand'
    canonical_positions = ((1, 2, 3, 4, 5, 6, 7), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1
    stem_class = AcceptorStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The 5' strand of the acceptor stem of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class PositionEight(Nucleotide):
    name = 'position 8'
    canonical_positions = ((8, ), )
    conserved_nucs = ({0: 'T'}, )
    num_allowed_unconserved = 1

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        """The nucleotide at canonical position 8 of tRNA, expected but not required to be T"""

        super().__init__(string,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class PositionNine(Nucleotide):
    name = 'position 9'
    canonical_positions = ((9, ), )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        """The nucleotide at canonical position 9 of tRNA"""

        super().__init__(string,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class DArm(Arm):
    name = 'D arm'
    num_allowed_unconserved = 2

    def __init__(self, stem, loop, cautious=False):
        """The D arm of tRNA"""

        super().__init__(stem,
                         loop,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         cautious=cautious)


class DStem(Stem):
    name = 'D stem'
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1
    arm_class = DArm

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        """The D stem of tRNA, which can be 3 (Type I tRNA) or 4 (Type II) nucleotides long"""

        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class DStemFiveprimeStrand(Sequence):
    name = 'D stem 5\' strand'
    canonical_positions = ((10, 11, 12), (13, ))
    allowed_section_lengths = ((3, ), (0, 1))
    allowed_input_lengths = tuple(itertools.product(*allowed_section_lengths))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucs = ({}, {})
    num_allowed_unconserved = -1
    arm_class = DArm
    stem_class = DStem

    def __init__(self, positions_10_to_12_string, position_13_string='', start_index=None, stop_index=None, cautious=False):
        """The 5' strand of the D stem of tRNA"""

        if cautious:
            if len(positions_10_to_12_string) != 3:
                raise TRNAIdentifierError("Your `positions_10_to_12_string` was not the required 3 bases long: %s"
                                          % positions_10_to_12_string)
            if not 0 <= len(position_13_string) <= 1:
                raise TRNAIdentifierError("Your `position_13_string` was not the required 0 or 1 bases long: %s"
                                          % position_13_string)

        super().__init__((positions_10_to_12_string, position_13_string),
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=False)


class DLoop(Loop):
    name = 'D loop'
    canonical_positions = ((14, 15), (16, 17), (18, 19), (20, ), (21, ))
    allowed_section_lengths = ((2, ), (1, 2, 3), (2, ), (1, 2, 3), (1, ))
    allowed_input_lengths = tuple(itertools.product(*allowed_section_lengths))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucs = ({0: 'A', 1: ('A', 'G')}, {}, {0: 'G', 1: 'G'}, {}, {0: ('A', 'G')})
    num_allowed_unconserved = 2
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
        """The D loop of tRNA, allowing for variable alpha (canonical positions 16, 17, 17a, 17b) and beta (canonical positions 20, 20a, 20b) sections"""

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
        self.alpha_seq = Sequence(alpha_positions_string, start_index=alpha_start_index, stop_index=alpha_stop_index)
        beta_start_index = alpha_stop_index + 2
        beta_stop_index = beta_start_index + len(beta_positions_string)
        self.beta_seq = Sequence(beta_positions_string, start_index=beta_start_index, stop_index=beta_stop_index)

        super().__init__((positions_14_to_15_string,
                          alpha_positions_string,
                          positions_18_to_19_string,
                          beta_positions_string,
                          position_21_string),
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class DStemThreeprimeStrand(Sequence):
    name = 'D stem 3\' strand'
    canonical_positions = ((22, ), (23, 24, 25))
    allowed_section_lengths = ((0, 1), (1, 2, 3))
    allowed_input_lengths = tuple(itertools.product(*allowed_section_lengths))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucs = ({}, {})
    num_allowed_unconserved = -1
    arm_class = DArm
    stem_class = DStem

    def __init__(self, position_22_string, positions_23_to_25_string='', start_index=None, stop_index=None, cautious=False):
        """The 3' strand of the D stem of tRNA"""

        if cautious:
            if not 0 <= len(position_22_string) <= 1:
                raise TRNAIdentifierError("Your `position_22_string` was not the required 1 base long: %s" % position_22_string)
            if len(positions_23_to_25_string) != 3:
                raise TRNAIdentifierError("Your `positions_23_to_25_string` was not the required 1 to 3 bases long: %s" % positions_23_to_25_string)

        super().__init__((position_22_string, positions_23_to_25_string),
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=False)


class PositionTwentySix(Nucleotide):
    name = 'position 26'
    canonical_positions = ((26, ), )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        """The nucleotide at canonical position 26 of tRNA"""

        super().__init__(string,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class AnticodonArm(Arm):
    name = 'anticodon arm'
    num_allowed_unconserved = 1

    def __init__(self, stem, loop, cautious=False):
        """The anticodon arm of tRNA"""

        self.anticodon = loop.anticodon
        super().__init__(stem,
                         loop,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         cautious=cautious)


class AnticodonStem(Stem):
    name = 'anticodon stem'
    arm_class = AnticodonArm
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        """The anticodon stem of tRNA"""

        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class AnticodonStemFiveprimeStrand(Sequence):
    name = 'anticodon stem 5\' strand'
    canonical_positions = ((27, 28, 29, 30, 31), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1
    stem_class = AnticodonStem
    arm_class = AnticodonArm

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The 5' strand of the anticodon stem of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class Anticodon(Sequence):
    name = 'anticodon'
    canonical_positions = ((34, 35, 36))
    allowed_input_lengths = ((3, ), )
    summed_input_lengths = (3, )
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The anticodon sequence of tRNA"""

        super().__init__(substrings, start_index=start_index, stop_index=stop_index, cautious=cautious)
        try:
            self.aa_string = ANTICODON_TO_AA[self.string]
        except KeyError:
            self.aa_string = 'NA'


class AnticodonLoop(Loop):
    name = 'anticodon loop'
    canonical_positions = ((32, 33, 34, 35, 36, 37, 38), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucs = ({1: 'T', 5: ('A', 'G')}, )
    num_allowed_unconserved = 1
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The anticodon loop of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)
        self.anticodon = Anticodon(self.string[2: 5])


class AnticodonStemThreeprimeStrand(Sequence):
    name = 'anticodon stem 3\' strand'
    canonical_positions = ((39, 40, 41, 42, 43), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1
    arm_class = AnticodonArm
    stem_class = AnticodonStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The 3' strand of the anticodon stem of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class VLoop(Loop):
    name = 'V loop'
    canonical_positions = ((44, 45, 46, 47, 48), )
    allowed_input_lengths = tuple(itertools.product(range(4, 24)))
    summed_input_lengths = tuple(map(sum, allowed_input_lengths))
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The V loop of tRNA: no stem/loop structure and base pairing is considered"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)
        if 4 <= len(self.string) <= 5:
            self.type = 'I'
        elif 12 <= len(self.string) <= 23:
            self.type = 'II'
        else:
            self.type = 'NA'


class TArm(Arm):
    name = 'T arm'
    num_allowed_unconserved = 2

    def __init__(self, stem, loop, cautious=False):
        """The T arm of tRNA"""

        super().__init__(stem,
                         loop,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         cautious=cautious)


class TStem(Stem):
    name = 'T stem'
    arm_class = TArm
    num_allowed_unconserved = -1
    num_allowed_unpaired = 1

    def __init__(self, fiveprime_seq, threeprime_seq, cautious=False):
        """The T stem of tRNA, with canonical positions 53 and 61 expected to be G and C"""

        super().__init__(fiveprime_seq,
                         threeprime_seq,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         num_allowed_unpaired=self.num_allowed_unpaired,
                         cautious=cautious)


class TStemFiveprimeStrand(Sequence):
    name = 'T stem 5\' strand'
    canonical_positions = ((49, 50, 51, 52, 53), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucs = ({4: 'G'}, )
    num_allowed_unconserved = 1
    arm_class = TArm
    stem_class = TStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The 5' T strand of the T stem of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class TLoop(Loop):
    name = 'T loop'
    canonical_positions = ((54, 55, 56, 57, 58, 59, 60), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucs = ({0: 'T', 1: 'T', 2: 'C', 3: ('A', 'G'), 4: 'A'}, )
    num_allowed_unconserved = 2
    arm_class = TArm

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The T loop of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class TStemThreeprimeStrand(Sequence):
    name = 'T stem 3\' strand'
    canonical_positions = ((61, 62, 63, 64, 65), )
    allowed_input_lengths = ((5, ), )
    summed_input_lengths = (5, )
    conserved_nucs = ({0: 'C'}, )
    num_allowed_unconserved = 1
    arm_class = TArm
    stem_class = TStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The 3' strand of the T stem of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class AcceptorStemThreeprimeStrand(Sequence):
    name = 'acceptor stem 3\' strand'
    canonical_positions=((66, 67, 68, 69, 70, 71, 72), )
    allowed_input_lengths = ((7, ), )
    summed_input_lengths = (7, )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1
    stem_class = AcceptorStem

    def __init__(self, substrings, start_index=None, stop_index=None, cautious=False):
        """The 3' strand of the acceptor stem of tRNA"""

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class Discriminator(Nucleotide):
    name = 'discriminator'
    canonical_positions = ((73, ), )
    conserved_nucs = ({}, )
    num_allowed_unconserved = -1

    def __init__(self, string, start_index=None, stop_index=None, cautious=False):
        """The discriminator nucleotide 5' of the acceptor sequence"""

        super().__init__(string,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class Acceptor(Sequence):
    name = 'acceptor'
    canonical_positions = ((74, 75, 76), )
    allowed_input_lengths = ((3, ), (2, ), (1, )) # Input lengths of 2 and 1 are 3'-CC and 3'-C.
    summed_input_lengths = (3, 2, 1)
    conserved_nucs = ({0: 'C', 1: 'C', 2: 'A'}, ) # set in stone
    num_allowed_unconserved = 0 # set in stone
    max_extra_threeprime = 2 # set in stone -- allowing this to vary would require a lot of work in `trnaseq.py`

    def __init__(self, substrings, num_extra_threeprime=0, start_index=None, stop_index=None, cautious=False):
        """The acceptor sequence (CCA in mature tRNA) and various other 3' endings

        With default values of `self.allowed_input_lengths` and `self.num_extra_threeprime`,
        a tRNA profile will prioritize CCA > CC > C > CCAN > CCANN.
        """

        # The possibilities of tRNA ending in 3'-CC and 3'-C are accommodated by `self.allowed_input_lengths`.
        # When bases beyond the 3' end of 3'-CCA (CCAN...) in profiling,
        # `self.num_extra_threeprime` SHOULD be set accordingly (1 for CCAN, 2 for CCANN) --
        # these extra 3' bases are not treated explicitly by this class.

        super().__init__(substrings,
                         conserved_nucs=self.conserved_nucs,
                         num_allowed_unconserved=self.num_allowed_unconserved,
                         start_index=start_index,
                         stop_index=stop_index,
                         cautious=cautious)


class Profile:
    FIVEPRIME_TO_THREEPRIME_FEATURE_CLASSES = TRNAFeature.list_all_tRNA_features()
    THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES = FIVEPRIME_TO_THREEPRIME_FEATURE_CLASSES[::-1]
    STEM_FORMATION_TRIGGERS = [TStemFiveprimeStrand, AnticodonStemFiveprimeStrand, DStemFiveprimeStrand, AcceptorStemFiveprimeStrand]
    ARM_FORMATION_TRIGGERS = [TStem, AnticodonStem, DStem]
    T_ARM_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(TArm)
    THREEPRIME_STEM_SEQ_INDICES = {TStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(TStemThreeprimeStrand),
                                   AnticodonStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(AnticodonStemThreeprimeStrand),
                                   DStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(DStemThreeprimeStrand),
                                   AcceptorStem: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(AcceptorStemThreeprimeStrand)}
    D_LOOP_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(DLoop)
    ARM_LOOP_INDEX_DICT = {TArm: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(TLoop),
                           AnticodonArm: THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(AnticodonLoop),
                           DArm: D_LOOP_INDEX}
    ANTICODON_LOOP_INDEX = THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES.index(AnticodonLoop)
    EXTRAPOLATION_INELIGIBLE_FEATURES = [AcceptorStemThreeprimeStrand,
                                         TStemThreeprimeStrand,
                                         VLoop,
                                         AnticodonStemThreeprimeStrand,
                                         DStemThreeprimeStrand]

    def __init__(self, input_seq, name=''):
        """A profile of the tRNA features identified from the 3' end of the input sequence"""

        # The input sequence is treated like a tRNA-seq read starting from the 3' end of a tRNA molecule.
        self.input_seq = input_seq
        self.name = name

        (self.profiled_seq,
         self.features,
         self.num_conserved,
         self.num_unconserved,
         self.num_paired,
         self.num_unpaired,
         self.num_extra_threeprime,
         self.num_in_extrapolated_fiveprime_feature,
         self.has_complete_feature_set,
         self.num_extra_fiveprime) = self.get_profile(unprofiled_seq=self.input_seq[::-1])

        self.feature_names = [f.name for f in self.features]

        if self.features:
            # The extra 3' nucleotides are not explicitly added to the Acceptor object.
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

            if len(self.features) > self.T_ARM_INDEX:
                self.is_predicted_trna = True
            else:
                self.is_predicted_trna = False
            # If the sequence does not have a complete feature set but is longer than tRNA should be,
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
        """A recursive function to sequentially identify tRNA features from the 3' end of a sequence"""

        if features is None:
            features = []

        if feature_index == len(self.THREEPRIME_TO_FIVEPRIME_FEATURE_CLASSES):
            # To reach this point,
            # all tRNA features including the tRNA-His 5'-G must have been found,
            # and the input sequence must extend 5' of that.
            return (profiled_seq, # string in 5' to 3' direction
                    features, # listed in the 5' to 3' direction
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    num_extra_threeprime,
                    0, # number of nucleotides in an extrapolated 5' feature -- there is no extrapolated 5' feature
                    has_complete_feature_set,
                    len(unprofiled_seq)) # extra 5' nucleotides

        if not unprofiled_seq:
            # To reach this point,
            # the full length of the input sequence must have been profiled with tRNA features.
            return (profiled_seq,
                    features, # listed in the 5' to 3' direction
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    num_extra_threeprime,
                    0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                    has_complete_feature_set,
                    0) # input sequence is not longer than full-length tRNA


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
                    return (profiled_seq,
                            features, # listed in the 5' to 3' direction
                            num_conserved,
                            num_unconserved,
                            num_paired,
                            num_unpaired,
                            num_extra_threeprime,
                            0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                            has_complete_feature_set,
                            len(unprofiled_seq)) # extra 5' nucleotides


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
                                    incremental_profile_candidates.append((unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                                                           [arm, stem, feature],
                                                                           feature.num_conserved,
                                                                           feature.num_unconserved,
                                                                           stem.num_paired,
                                                                           stem.num_unpaired,
                                                                           summed_input_length - len(unprofiled_seq))) # number of nucleotides in extrapolated 5' feature
                                    continue
                            else:
                                incremental_profile_candidates.append((unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                                                       [stem, feature],
                                                                       feature.num_conserved,
                                                                       feature.num_unconserved,
                                                                       stem.num_paired,
                                                                       stem.num_unpaired,
                                                                       summed_input_length - len(unprofiled_seq))) # number of nucleotides in extrapolated 5' feature
                                continue
                    else:
                        incremental_profile_candidates.append((unprofiled_seq[::-1], # flip orientation to 5' to 3'
                                                               [feature],
                                                               feature.num_conserved,
                                                               feature.num_unconserved,
                                                               0, # since we're not considering a stem, 0
                                                               0, # since we're not considering a stem, 0
                                                               summed_input_length - len(unprofiled_seq))) # number of nucleotides in extrapolated 5' feature
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
                if feature_class.name == 'acceptor':
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
                                    incremental_profile_candidates.append((unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                                                           [arm, stem, feature],
                                                                           feature.num_conserved,
                                                                           feature.num_unconserved,
                                                                           stem.num_paired,
                                                                           stem.num_unpaired,
                                                                           0)) # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                                    continue
                            else:
                                incremental_profile_candidates.append((unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                                                       [stem, feature],
                                                                       feature.num_conserved,
                                                                       feature.num_unconserved,
                                                                       stem.num_paired,
                                                                       stem.num_unpaired,
                                                                       0)) # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                                continue
                    else:
                        incremental_profile_candidates.append((unprofiled_seq[: num_processed_bases][::-1], # flip orientation to 5' to 3'
                                                               [feature],
                                                               feature.num_conserved,
                                                               feature.num_unconserved,
                                                               0, # since we're not considering a stem, 0
                                                               0, # since we're not considering a stem, 0
                                                               0)) # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature

                        # Avoid testing 3'-CC and 3'-C if 3'-CCA, was just found.
                        if feature.name == 'acceptor':
                            if feature.string == 'CCA':
                                # The following `continue` would otherwise test 3'-CC.
                                break

                        continue
                # When considering 3'-CCAN..., don't allow the CCA part to vary from CCA.
                # The "elif" block here is entered after considering CCAN....
                # Without the "break" below, CCN..., CN..., N... would then be considered.
                elif feature.name == 'acceptor':
                    if summed_input_length == 3:
                        if num_extra_threeprime > 0:
                            break

            # Avoid testing CCN..., CN..., N... when CCAN... was just found.
            if feature_class.name == 'acceptor':
                if num_extra_threeprime > 0:
                    break


        if not incremental_profile_candidates:
            # The feature didn't pass muster.
            if feature_class.name == 'acceptor':
                if num_extra_threeprime < Acceptor.max_extra_threeprime:
                    # This will try to find 3'-CCAN..., sequentially adding extra bases.
                    # The first feature of the profile, the acceptor sequence or variant thereof, has not be found.
                    # Do not bother to compare this next profile to the current non-existent profile, so `return`.
                    return self.get_profile(unprofiled_seq=unprofiled_seq[1: ],
                                            profiled_seq=unprofiled_seq[0] + profiled_seq, # 5' to 3' orientation
                                            features=[], # extra 3' nucleotides are not counted as a feature
                                            num_extra_threeprime=num_extra_threeprime + 1)

            return (profiled_seq,
                    features,
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    num_extra_threeprime,
                    0, # number of nucleotides in extrapolated 5' feature -- there is no extrapolated 5' feature
                    has_complete_feature_set,
                    0) # input sequence is not longer than full-length tRNA


        # Sort candidates by
        # 1. number of features identified (at most, sequence + stem + arm) (descending),
        # 2. number of unconserved + unpaired nucleotides (ascending),
        # 3. incompleteness of the last (most 5') feature (ascending).
        # This sort also happens later for full sequence profiles,
        # but this first sort is useful for seeking out and returning "flawless" mature tRNA.
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
                                                         num_extra_threeprime=num_extra_threeprime,
                                                         feature_index=feature_index + len(ipc[1]),
                                                         has_complete_feature_set=True)
                else:
                    profile_candidate = self.get_profile(unprofiled_seq=unprofiled_seq[len(ipc[0]): ], # recurse
                                                         profiled_seq=ipc[0] + profiled_seq,
                                                         features=ipc[1] + features,
                                                         num_conserved=ipc[2] + num_conserved,
                                                         num_unconserved=ipc[3] + num_unconserved,
                                                         num_paired=ipc[4] + num_paired,
                                                         num_unpaired=ipc[5] + num_unpaired,
                                                         num_extra_threeprime=num_extra_threeprime,
                                                         feature_index=feature_index + len(ipc[1]),
                                                         has_complete_feature_set=False)
                if (profile_candidate[8] # has complete feature set
                    and profile_candidate[3] == 0 # number unconserved
                    and profile_candidate[5] == 0): # number unpaired
                    return profile_candidate
                else:
                    profile_candidates.append(profile_candidate)
            else: # 5' feature was extrapolated from unprofiled input sequence
                profile_candidates.append((ipc[0] + profiled_seq,
                                           ipc[1] + features,
                                           ipc[2] + num_conserved,
                                           ipc[3] + num_unconserved,
                                           ipc[4] + num_paired,
                                           ipc[5] + num_unpaired,
                                           num_extra_threeprime,
                                           ipc[6], # number of nucleotides in extrapolated 5' feature
                                           False, # does not have a complete feature set
                                           0)) # input sequence is not longer than full-length tRNA

        # Consider an extra 3' nucleotide in the acceptor sequence.
        if feature_class.name == 'acceptor':
            if num_extra_threeprime < Acceptor.max_extra_threeprime:
                profile_candidates.append(self.get_profile(unprofiled_seq=unprofiled_seq[1: ],
                                                           profiled_seq=unprofiled_seq[0] + profiled_seq, # 5' to 3' orientation
                                                           features=[], # extra 3' nucleotides are not counted as a feature
                                                           num_extra_threeprime=num_extra_threeprime + 1))

        # Do not add 5' features to profiled 3' features if the additional features
        # have no grounding in conserved nucleotides or paired nucleotides in stems.
        profile_candidates = [p for p in profile_candidates if p[2] + p[4] > num_conserved + num_paired]
        if profile_candidates:
            profile_candidates.sort(key=lambda p: (-len(p[1]), p[3] + p[5], p[7]))
            return profile_candidates[0]
        else:
            return (profiled_seq,
                    features,
                    num_conserved,
                    num_unconserved,
                    num_paired,
                    num_unpaired,
                    num_extra_threeprime,
                    0, # number of nucleotides in extrapolated 5' feature
                    has_complete_feature_set,
                    0) # input sequence is not longer than full-length tRNA


    def get_unconserved_positions(self):
        """Get information on the unexpectedly unconserved positions in a tRNA feature profile

        Returns
        =======
        unconserved_info : list
            List of tuples, one for each unconserved nucleotide in the profile
            A tuple has three elements:
                1. position (index) in the input sequence
                2. observed nucleotide in the input
                3. expected canonical nucleotides at the site,
                   a string with nucleotides separated by commas if more than one
        """

        unconserved_info = []
        for feature in self.features:
            # Only Nucleotide and Sequence subclasses have the attribute.
            if hasattr(feature, 'conserved_status'):
                component_start_index = feature.start_index
                # Conserved nucleotides are indexed within the string "component" (substring).
                for string_component_statuses, string_component in zip(
                    feature.conserved_status, feature.string_components):
                    for nuc_index, is_conserved, observed_nuc, expected_nucs in string_component_statuses:
                        # Avoid N padding in an extrapolated 5' feature.
                        if not is_conserved and observed_nuc != 'N':
                            unconserved_info.append((component_start_index + nuc_index,
                                                     observed_nuc,
                                                     ','.join(expected_nucs)))
                    component_start_index += len(string_component)

        return unconserved_info


    def get_unpaired_positions(self):
        """Get information on the unexpectedly unpaired bases in a tRNA feature profile

        Returns
        =======
        unpaired_info : list
            List of tuples, one for each unpaired pair of nucleotides
            A tuple has four elements:
                1. position (index) of 5' nucleotide in the input sequence
                2. position (index) of 3' nucleotide in the input sequence
                3. observed 5' nucleotide
                4. observed 3' nucleotide
        """

        unpaired_info = []
        for feature in self.features:
            # Only the Stem class has the attribute.
            if hasattr(feature, 'paired_status'):
                for nuc_index, (is_paired, fiveprime_nuc, threeprime_nuc) in enumerate(feature.paired_status):
                    # Avoid N padding in an extrapolated 5' feature.
                    if not is_paired and fiveprime_nuc != 'N':
                        unpaired_info.append((feature.fiveprime_seq.start_index + nuc_index,
                                              feature.threeprime_seq.stop_index - nuc_index - 1,
                                              fiveprime_nuc,
                                              threeprime_nuc))

        return unpaired_info


class GeneProfile:
    __slots__ = ('unencoded_acceptor_profile',
                 'encoded_acceptor_profile',
                 'predicted_profile',
                 'has_encoded_acceptor',
                 'predicted_genomic_seq')

    def __init__(self, input_seq, name='', check_encoded_acceptor=True):
        """Profile the tRNA features in a gene sequence

        A gene may or may not encode the 3'-CCA acceptor sequence.
        """

        self.unencoded_acceptor_profile = Profile(input_seq + 'CCA', name=name)
        if check_encoded_acceptor:
            self.encoded_acceptor_profile = Profile(input_seq, name=name)
            if self.unencoded_acceptor_profile.is_predicted_trna and not self.encoded_acceptor_profile.is_predicted_trna:
                self.predicted_profile = self.unencoded_acceptor_profile
                self.has_encoded_acceptor = False
                self.predicted_genomic_seq = input_seq
            elif self.encoded_acceptor_profile.is_predicted_trna and not self.unencoded_acceptor_profile.is_predicted_trna:
                self.predicted_profile = self.encoded_acceptor_profile
                self.has_encoded_acceptor = True
                self.predicted_genomic_seq = input_seq[: -3]
            else:
                self.predicted_profile = None
                self.has_encoded_acceptor = None
                self.predicted_genomic_seq = None
        else:
            self.encoded_acceptor_profile = None
            if self.unencoded_acceptor_profile.is_predicted_trna:
                self.predicted_profile = self.unencoded_acceptor_profile
                self.has_encoded_acceptor = False
                self.predicted_genomic_seq = input_seq
            else:
                self.predicted_profile = None
                self.has_encoded_acceptor = None
                self.predicted_genomic_seq = None
