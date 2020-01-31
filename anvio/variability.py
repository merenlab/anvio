# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of single nucleotide variation"""

import numpy as np

from collections import Counter

import anvio

from anvio.constants import nucleotides
from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class VariablityTestFactory:
    def __init__(self, params={'b': 3, 'm': 1.45, 'c': 0.05}):
        self.params = params


    def get_min_acceptable_departure_from_consensus(self, coverage):
        """Get minimum allowable departure from consensus

        Notes
        =====
        - 0 returned if self.params is None
        - https://www.desmos.com/calculator/qwocua4zi5
        - https://i.imgur.com/zd04pui.png
        """
        if self.params is None:
            if hasattr(coverage, '__len__'):
                return np.zeros(len(coverage))
            else:
                return 0

        b, m, c = self.params['b'], self.params['m'], self.params['c']

        return (1 / b) ** (coverage ** (1 / b) - m) + c


def get_competing_items(reference, items_frequency_tuples_list=[]):
    """Resolves competing nts and aas.

       This function will return None if there is no varaition and the most frequent
       item is equal to the reference. But will NOT return None if there is no
       variation AND the most frequent item is different than the reference.

       `items_frequency_tuples_list` MUST BE SORTED & should look like this:

            >>> [('Val', 69), ('Asn', 0), ('Gln', 0), ('Cys', 0), ('Glu', 0), ...]

        See issue #544 (https://github.com/merenlab/anvio/issues/544) for test data.

    """

    num_items = len(items_frequency_tuples_list)
    if not num_items:
        raise ConfigError("Wat. You sent an empty list to `get_competing_items` :(")

    # get the most frequent base
    most_frequent_item = items_frequency_tuples_list[0][0]

    if (num_items == 1 or not items_frequency_tuples_list[1][1]) and most_frequent_item == reference:
        # ^^^^^ the list has only one item, or the second item has 0 frequency ^^^^^
        # so there is no variation, and the most frequent base IS the reference.
        # nothing to see here.
        return None
    elif (num_items == 1 or not items_frequency_tuples_list[1][1]) and most_frequent_item != reference:
        # ^^^^^ the list has only one item, or the second item has 0 frequency ^^^^^
        # again, there is no variation, but the most frequent base DIFFERS from the reference.
        # much more interesting. We are not returning None here, becasue we do not want this
        # information to be confused wit hthe true case of no variation (see the case above).
        return [most_frequent_item, most_frequent_item]
    else:
        # the only other option is to have multiple bases in items_frequency_tuples_list.
        # competing nts are simply the most frequent two nucleotides in the column.
        # clearly, the `reference` nucleotide (which is the observed nucleotide in
        # the contig for this particular `pos`) may not be one of these. but here,
        # we don't care about that.

        # FIXME: CONGRATULATIONS. YOU DID FIND THE NTH SHITTIEST PIECE OF CODE IN THIS
        #        REPOSITORY. IF YOU SEND US AN E-MAIL, SOMEONE WILL RESPOND WITH A
        #        FORMAL APOLOGY FOR THIS MONSTROSITY.
        if num_items > 2 and items_frequency_tuples_list[1][1] == items_frequency_tuples_list[2][1]:
            frequency_of_the_second = items_frequency_tuples_list[1][1]
            second_most_frequent_items = sorted([tpl[0] for tpl in items_frequency_tuples_list if tpl[1] == frequency_of_the_second and tpl[0] != most_frequent_item], reverse=True)
            return sorted([most_frequent_item, second_most_frequent_items[0]])
        else:
            second_most_frequent_item = items_frequency_tuples_list[1][0]
            return sorted([most_frequent_item, second_most_frequent_item])


class ProcessAlleleCounts:
    def __init__(self, allele_counts, allele_to_array_index, sequence, min_coverage=None):
        """A class to process raw variability information for a given allele counts array

        FIXME
        """
        self.allele_counts = allele_counts
        self.allele_to_array_index = allele_to_array_index
        self.array_index_to_allele = {v: k for k, v in self.allele_to_array_index.items()}
        self.sequence = sequence
        self.min_coverage = min_coverage

        # the sequence cast as indices
        self.sequence_as_index = np.array([allele_to_array_index[item] for item in self.sequence])

        if len(self.sequence) != self.allele_counts.shape[1]:
            raise ConfigError("ProcessAlleleCounts :: allele_counts has %d positions, but sequence has %d." \
                              % (len(self.sequence), self.allele_counts.shape[1]))


    def subset_arrays_by_positions(self, positions, *arrays):
        out = []
        for array in arrays:
            if array.ndim == 1:
                out.append(array[positions])
            else:
                out.append(array[:, positions])

        return tuple(out)


    def process(self):
        coverage = self.get_coverage()

        positions_above_coverage_threshold = self.get_positions_above_coverage_threshold(coverage, self.min_coverage)
        num_positions = len(positions_above_coverage_threshold)

        if num_positions != len(self.sequence):
            self.sequence_as_index, self.allele_counts, coverage = \
                self.subset_arrays_by_positions(positions_above_coverage_threshold,
                                                self.sequence_as_index,
                                                self.allele_counts,
                                                coverage)

        reference_coverage = self.get_reference_coverage()
        departure_from_consensus = self.get_departure_from_consensus(reference_coverage, coverage)
        competing_items = self.get_competing_items(reference_coverage, coverage)




    def get_coverage(self):
        return np.sum(self.allele_counts, axis=0)


    def get_reference_coverage(self):
        return self.allele_counts[self.sequence_as_index, np.arange(self.allele_counts.shape[1])]


    def get_departure_from_consensus(self, reference_coverage=None, coverage=None):
        if reference_coverage is None:
            reference_coverage = self.get_reference_coverage()

        if coverage is None:
            coverage = self.get_coverage()

        return 1 - reference_coverage/coverage


    def get_competing_items(self, reference_coverage=None, coverage=None):
        if reference_coverage is None:
            reference_coverage = self.get_reference_coverage()

        if coverage is None:
            coverage = self.get_coverage()

        n = self.allele_counts.shape[1]

        # as a first pass, sort the row indices (-allele_counts_array is used to sort from highest -> lowest)
        competing_items_as_index = np.argsort(-self.allele_counts, axis=0)

        # take the top 2 items
        competing_items_as_index = competing_items_as_index[:2, :]

        # get the coverage of the second item
        coverage_second_item = self.allele_counts[competing_items_as_index[1, :], np.arange(n)]

        # if the coverage of the second item is 0, set the second index equal to the first
        competing_items_as_index[1, :] = np.where(coverage_second_item == 0, competing_items_as_index[0, :], competing_items_as_index[1, :])

        # sort the competing nts
        competing_items_as_index = np.sort(competing_items_as_index, axis=0)

        # make the competing nts list
        nts_1 = [self.array_index_to_allele[index_1] for index_1 in competing_items_as_index[0, :]]
        nts_2 = [self.array_index_to_allele[index_2] for index_2 in competing_items_as_index[1, :]]
        competing_items = np.fromiter((nt_1 + nt_2 for nt_1, nt_2 in zip(nts_1, nts_2)), np.dtype('<U2'), count=n)

        # If the second item is 0, and the reference is the first item, set competing_items to None.
        # This can easily be checked by seeing if reference_coverage == coverage
        competing_items = np.where(reference_coverage == coverage, None, competing_items)

        return competing_items


    def get_positions_above_coverage_threshold(self, coverage=None, threshold=None):
        if coverage is None:
            coverage = self.get_coverage()

        if threshold is None:
            if self.min_coverage:
                threshold = self.min_coverage
            else:
                # no threshold given, give all positions
                return np.arange(len(self.sequence))

        return np.where(coverage >= threshold)[0]


    def get_worth_reporting_array(self):
        return 


        # if we came all the way down here, we want this position to be reported as a variable
        # nucleotide position.
        self.profile['worth_reporting'] = True

        if test_class:
            # BUT THEN if there is a test class, we have to check whether we have enough coverage to be confident
            # to suggest that this variable position should be reported, and flip that report flag
            min_acceptable_departure_from_consensus = test_class.min_acceptable_departure_from_consensus(self.profile['coverage'])
            if departure_from_reference < min_acceptable_departure_from_consensus:
                self.profile['worth_reporting'] = False


class ColumnProfile:
    """A class to report raw variability information for a given nucleotide position"""

    def __init__(self, column, reference, coverage=None, pos=None, split_name=None, sample_id=None, test_class=None):
        # make sure we have coverage value regardless of user's input:
        coverage = coverage if coverage else len(column)

        if len(reference) != 1:
            raise ConfigError("ColumnProfile class is upset. The reference must be a single base.")

        self.profile = {'sample_id': sample_id, 'split_name': split_name, 'pos': pos, 'reference': reference,
                        'coverage': coverage, 'departure_from_reference': 0, 'competing_nts': None, 'worth_reporting': False}

        nt_counts = Counter(column)
        for nt in nucleotides:
            self.profile[nt] = nt_counts[nt]

        # sort nts based on their frequency
        nts_sorted = nt_counts.most_common()

        competing_nts = get_competing_items(reference, nts_sorted)
        if not competing_nts:
            # there is no action here, we can return without further processing.
            return

        self.profile['competing_nts'] = ''.join(competing_nts)

        # here we quantify the ratio of frequencies of non-reference-nts observed in this column
        # to the overall coverage, and store it as `departure_from_reference`:
        total_frequency_of_all_bases_but_the_reference = sum([tpl[1] for tpl in nts_sorted if tpl[0] != reference])
        departure_from_reference = total_frequency_of_all_bases_but_the_reference / coverage
        self.profile['departure_from_reference'] = departure_from_reference

        # if we came all the way down here, we want this position to be reported as a variable
        # nucleotide position.
        self.profile['worth_reporting'] = True

        if test_class:
            # BUT THEN if there is a test class, we have to check whether we have enough coverage to be confident
            # to suggest that this variable position should be reported, and flip that report flag
            min_acceptable_departure_from_consensus = test_class.min_acceptable_departure_from_consensus(self.profile['coverage'])
            if departure_from_reference < min_acceptable_departure_from_consensus:
                self.profile['worth_reporting'] = False
