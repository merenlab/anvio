# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of single nucleotide variation"""

import numpy as np
import pandas as pd

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


class ProcessAlleleCounts:
    def __init__(self, allele_counts, allele_to_array_index, sequence, min_coverage=None, test_class=None):
        """A class to process raw variability information for a given allele counts array

        FIXME
        """
        self.data = {}

        self.allele_counts = allele_counts
        self.allele_to_array_index = allele_to_array_index
        self.array_index_to_allele = {v: k for k, v in self.allele_to_array_index.items()}
        self.sequence = sequence
        self.min_coverage = min_coverage
        self.test_class = test_class
        self.positions = np.arange(len(sequence))

        # the sequence cast as indices
        self.sequence_as_index = np.array([allele_to_array_index[item] for item in self.sequence])

        if len(self.sequence) != self.allele_counts.shape[1]:
            raise ConfigError("ProcessAlleleCounts :: allele_counts has %d positions, but sequence has %d." \
                              % (len(self.sequence), self.allele_counts.shape[1]))


    def subset_by_index(self, indices, *arrays):
        out = []
        for array in arrays:
            if array.ndim == 1:
                out.append(array[indices])
            else:
                out.append(array[:, indices])

        return tuple(out)


    def process(self):
        coverage = self.get_coverage()

        indices_to_keep = self.get_indices_above_coverage_threshold(coverage, self.min_coverage)
        num_good_coverage = len(indices_to_keep)

        if num_good_coverage != len(self.positions):
            # Some positions were not well covered. Remove them
            target = (self.positions, self.sequence_as_index, self.allele_counts, coverage)

            self.positions, self.sequence_as_index, self.allele_counts, coverage = \
                self.subset_by_index(indices_to_keep, *target)

        reference_coverage = self.get_reference_coverage()
        departure_from_reference = self.get_departure_from_reference(reference_coverage, coverage)

        indices_to_keep = self.get_positions_worth_reporting(coverage, departure_from_reference)
        num_worth_reporting = len(indices_to_keep)

        if num_worth_reporting != num_good_coverage:
            # Some positions were not worth reporting. Remove them
            target = (self.positions, self.sequence_as_index, self.allele_counts, coverage, reference_coverage, departure_from_reference)

            self.positions, self.sequence_as_index, self.allele_counts, coverage, reference_coverage, departure_from_reference = \
                self.subset_by_index(indices_to_keep, *target)

        competing_items = self.get_competing_items(reference_coverage, coverage)

        self.data.update({
            'pos': self.positions,
            'reference': [self.array_index_to_allele[x] for x in self.sequence_as_index],
            'coverage': coverage,
            'departure_from_reference': departure_from_reference,
            'competing_items': competing_items,
        })

        for index, item in self.array_index_to_allele.items():
            self.data[item] = self.allele_counts[index, :]


    def get_coverage(self):
        return np.sum(self.allele_counts, axis=0)


    def get_reference_coverage(self):
        return self.allele_counts[self.sequence_as_index, np.arange(self.allele_counts.shape[1])]


    def get_departure_from_reference(self, reference_coverage=None, coverage=None):
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


    def get_indices_above_coverage_threshold(self, coverage=None, threshold=None):
        if coverage is None:
            coverage = self.get_coverage()

        if threshold is None:
            if self.min_coverage:
                threshold = self.min_coverage
            else:
                # no threshold given, give all positions
                return np.arange(len(self.sequence))

        return np.where(coverage >= threshold)[0]


    def get_positions_worth_reporting(self, coverage, departure_from_reference):
        worth_reporting = np.array([True] * len(coverage))

        if not self.test_class:
            return worth_reporting

        threshold = self.test_class.get_min_acceptable_departure_from_consensus(coverage)

        return np.where(departure_from_reference >= threshold)[0]

