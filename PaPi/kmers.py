# -*- coding: utf-8
#
# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import string
import itertools



complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV',\
                               'tgcayrkmvhdbTGCAYRKMVHDB')

def rev_comp(seq):
    return seq.translate(complements)[::-1]

class KMers:
    def __init__(self, k = 4):
        self.kmers = {}
        self.k = k
        
        self.get_kmers()

    def get_kmers(self):
        k = self.k
        arg = ['ATCG'] * k
        kmers = set()
        
        for item in itertools.product(*arg):
            kmer = ''.join(item)
            if rev_comp(kmer) not in kmers:
                kmers.add(kmer)
        
        self.kmers[k] = kmers


    def get_kmer_frequency(self, sequence, dist_metric_safe = True):
        k = self.k
        sequence = sequence.upper()

        if len(sequence) < k:
            return None

        if not self.kmers.has_key(k):
            self.get_kmers(k)
        
        kmers = self.kmers[k]
        frequencies = dict(zip(kmers, [0] * len(kmers)))
        
        for i in range(0, len(sequence) - (k - 1)):
            kmer = sequence[i:i + k]
            
            # FIXME: this can be faster/better
            if len([n for n in kmer if n not in 'ATCG']):
                continue

            if frequencies.has_key(kmer):
                frequencies[kmer] += 1
            else:
                frequencies[rev_comp(kmer)] += 1

        if dist_metric_safe:
            # we don't want all kmer freq values to be zero. so the distance
            # metrics wouldn't go crazy. instead we fill it with 1. which
            # doesn't affect relative distances.
            if sum(frequencies.values()) == 0:
                words = self.kmers[self.k]
                frequencies = dict(zip(words, [1] * len(words)))

        return frequencies
