# -*- coding: utf-8
"""Simple KMers class to compute kmer-nucleotide frequecies"""

import itertools

import anvio

from anvio.constants import complements


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


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
