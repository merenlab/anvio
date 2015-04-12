# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

from __future__ import division
from collections import Counter

from PaPi.constants import nucleotides


class VariablityTestFactory:
    # an experimental class to make sense whether the nucleotide variation in a column
    # is meaningful beyond sequencing errors, given the coverage of that position.
    def __init__(self, params = {'b': 3, 'm': 1.45, 'c': 0.05}):
        self.params = params
        self.coverage_upper_limit = 500

        # for fast access
        self.cov_var_map_dict = dict([(c, self.curve(c)) for c in range(0, self.coverage_upper_limit + 1)])


    def min_acceptable_ratio_given_coverage(self, coverage):
        if coverage >= self.coverage_upper_limit:
            coverage = self.coverage_upper_limit

        return self.cov_var_map_dict[coverage]


    def curve(self, coverage, b=3, m=1.45, c=0.05):
        # https://www.desmos.com/calculator/qwocua4zi5
        # and/or https://i.imgur.com/zd04pui.png
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
