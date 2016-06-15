# -*- coding: utf-8
# pylint: disable=line-too-long

"""Classes to make sense of single nucleotide variation"""

from __future__ import division
from collections import Counter

import anvio

from anvio.constants import nucleotides
from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


class VariablityTestFactory:
    """an experimental class to make sense whether the nucleotide variation in a column
       is meaningful beyond sequencing errors, given the coverage of that position."""
    def __init__(self, params={'b': 3, 'm': 1.45, 'c': 0.05}):
        self.params = params
        self.coverage_upper_limit = 500

        # for fast access
        if params:
            self.cov_var_map_dict = dict([(c, self.curve(c)) for c in range(0, self.coverage_upper_limit + 1)])
        else:
            self.cov_var_map_dict = dict([(c, 0) for c in range(0, self.coverage_upper_limit + 1)])


    def min_acceptable_departure_from_consensus(self, coverage):
        """Returns a minimum departure from consensus ratio for a given coverage
           based on test factory settings."""
        if coverage >= self.coverage_upper_limit:
            coverage = self.coverage_upper_limit

        return self.cov_var_map_dict[coverage]


    def curve(self, coverage):
        # https://www.desmos.com/calculator/qwocua4zi5
        # and/or https://i.imgur.com/zd04pui.png
        b, m, c = self.params['b'], self.params['m'], self.params['c']
        y = ((1 / b) ** ((coverage ** (1 / b)) - m)) + c
        return y


class ColumnProfile:
    """A class to report raw variability information for a given nucleotide position"""

    def __init__(self, column, reference, coverage=None, pos=None, split_name=None, sample_id=None, test_class=None):
        # make sure we have coverage value regardless of user's input:
        coverage = coverage if coverage else len(column)

        if len(reference) != 1:
            raise ConfigError, "ColumnProfile class is upset. The reference must be a single base."

        self.profile = {'sample_id': sample_id, 'split_name': split_name, 'pos': pos, 'reference': reference,
                        'coverage': coverage, 'departure_from_reference': 0, 'competing_nts': None, 'worth_reporting': False}


        nt_counts = Counter(column)
        for nt in nucleotides:
            self.profile[nt] = nt_counts[nt]

        # sort nts based on their frequency
        nts_sorted = nt_counts.most_common()

        # get the most frequent base
        most_frequent_base = nts_sorted[0][0]

        if len(nts_sorted) == 1 and most_frequent_base == reference:
            # there is no variation, and the most frequent base is the reference.
            # nothing to see here.
            return
        elif len(nts_sorted) == 1 and most_frequent_base != reference:
            # there is no variation, but the most frequent base differs from the reference.
            # much more interesting.
            self.profile['competing_nts'] = most_frequent_base + most_frequent_base
        else:
            # the only other option is to have multiple bases in nts_sorted.
            # competing nts are simply the most frequent two nucleotides in the column.
            # clearly, the `reference` nucleotide (which is the observed nucleotide in
            # the contig for this particular `pos`) may not be one of these. but here,
            # we don't care about that.
            second_most_frequent_base = nts_sorted[1][0]
            self.profile['competing_nts'] = ''.join(sorted(most_frequent_base + second_most_frequent_base))

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
