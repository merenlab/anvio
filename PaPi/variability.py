# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import operator
from scipy import log2 as log
from numpy import sqrt


class ColumnProfile:
    def __init__(self, column, pos):
        self.pos = pos
        self.coverage = len(column)
        self.nucleotide_counts = {}
        self.consensus_nucleotide = None
        self.n2n1ratio = 0.0
        self.competing_nucleotides = ''

        nucleotides = list(set(column))
        self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])
        nucleotides_sorted_by_occurence = [x[0] for x in sorted(self.nucleotide_counts.iteritems(), key=operator.itemgetter(1), reverse=True)]

        if len(nucleotides_sorted_by_occurence) == 1:
            # no variation.
            self.consensus_nucleotide = column[0]
            return

        nucleotides = list(set(column))
        denominator = float(len(column))
        self.nucleotide_counts = dict([(n, column.count(n)) for n in nucleotides])

        n1 = nucleotides_sorted_by_occurence[0]
        n2 = nucleotides_sorted_by_occurence[1]

        self.consensus_nucleotide = n1
        self.competing_nucleotides = n1 + n2

        if 'N' in self.competing_nucleotides or 'n' in self.competing_nucleotides:
            return

        if self.nucleotide_counts[n1] * 1.0 / self.nucleotide_counts[n2] > 10:
            return

        if len(column) < 4:
            return

        self.n2n1ratio = 1.0 * self.nucleotide_counts[n2] / self.nucleotide_counts[n1]