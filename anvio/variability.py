# -*- coding: utf-8

"""Module to make sense of variability (SNPs) across nucleotide positions"""

from __future__ import division
from collections import Counter

import anvio

from anvio.constants import nucleotides


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class VariablityTestFactory:
    # an experimental class to make sense whether the nucleotide variation in a column
    # is meaningful beyond sequencing errors, given the coverage of that position.
    def __init__(self, params = {'b': 3, 'm': 1.45, 'c': 0.05}):
        self.params = params
        self.coverage_upper_limit = 500

        # for fast access
        if params:
            self.cov_var_map_dict = dict([(c, self.curve(c)) for c in range(0, self.coverage_upper_limit + 1)])
        else:
            self.cov_var_map_dict = dict([(c, 0) for c in range(0, self.coverage_upper_limit + 1)])


    def min_acceptable_ratio_given_coverage(self, coverage):
        if coverage >= self.coverage_upper_limit:
            coverage = self.coverage_upper_limit

        return self.cov_var_map_dict[coverage]


    def curve(self, coverage):
        # https://www.desmos.com/calculator/qwocua4zi5
        # and/or https://i.imgur.com/zd04pui.png
        b, m, c = self.params['b'], self.params['m'], self.params['c']
        y = ((1 / b) ** ((coverage ** (1/b)) - m)) + c
        return y


class ColumnProfile:
    def __init__(self, column, coverage=None, pos=None, split_name=None, sample_id=None, test_class=None):
        self.profile = {'sample_id': sample_id, 'split_name': split_name, 'pos': pos, 'consensus': None,
                        'coverage': coverage if coverage else len(column),
                        'n2n1ratio': 0, 'competing_nts': None}

        nt_counts = Counter(column)
        for nt in nucleotides:
            self.profile[nt] = nt_counts[nt] if nt_counts.has_key(nt) else 0

        competing_two = nt_counts.most_common(2)
        if len(competing_two) == 1:
            # no variation.
            self.profile['consensus'] = competing_two[0][0]
            return

        n1_tuple, n2_tuple = competing_two
        competing_nts = n1_tuple[0] + n2_tuple[0]

        if n1_tuple[0] == 'N':
            return

        self.profile['consensus'] = n1_tuple[0]

        if n2_tuple[0] == 'N':
            return

        n2n1ratio = n2_tuple[1] / n1_tuple[1]

        if test_class:
            if n2n1ratio > test_class.min_acceptable_ratio_given_coverage(self.profile['coverage']):
                self.profile['competing_nts'] = competing_nts
                self.profile['n2n1ratio'] = n2n1ratio
        else:
            # if there is no test class, just report everything.
            self.profile['competing_nts'] = competing_nts
            self.profile['n2n1ratio'] = n2n1ratio
