# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import PaPi.utils

kmers = PaPi.utils.KMers()


class Composition:
    def __init__(self, reference, auxiliary):
        """Gets the contig auxiliary dict, generaetes TNF and other compositional stats"""
        self.reference = reference
        self.auxiliary = auxiliary
        self.composition = {}

    def report(self):
        sequence = self.auxiliary['rep_seq']
        raw_length = len(sequence)
        
        A = sequence.count('A')
        T = sequence.count('T')
        C = sequence.count('C')
        G = sequence.count('G')
        UNKNOWN = raw_length - (A + T + C + G)

        self.composition['A'] = A
        self.composition['T'] = T
        self.composition['C'] = C
        self.composition['G'] = G
        self.composition['N'] = UNKNOWN

        length = raw_length - UNKNOWN

        self.composition['GC_content'] = (G + C) * 1.0 / length
        self.composition['tnf'] = kmers.get_kmer_frequency(sequence)

        return self.composition