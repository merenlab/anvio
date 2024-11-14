# -*- coding: utf-8
# pylint: disable=line-too-long
"""Simple KMers class to compute kmer frequecies

   This module should not be used for k > 4.
"""

import itertools

from collections import Counter

import anvio

from anvio.constants import complements


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


def rev_comp(seq):
    return seq.translate(complements)[::-1]


class KMers:
    def __init__(self, k=4, alphabet="ATCG", consider_rev_comps=True):
        self.kmers = {}
        self.alphabet = alphabet
        self.consider_rev_comps = consider_rev_comps
        self.k = k

        self.get_kmers()

    def get_kmers(self):
        k = self.k
        arg = [self.alphabet] * k
        kmers = set()

        for item in itertools.product(*arg):
            kmer = "".join(item)
            if self.consider_rev_comps:
                if rev_comp(kmer) not in kmers:
                    kmers.add(kmer)
            else:
                kmers.add(kmer)

        self.kmers[k] = kmers

    def get_kmer_frequency(self, sequence, dist_metric_safe=False):
        k = self.k
        sequence = sequence.upper()

        if len(sequence) < k:
            return None

        if k not in self.kmers:
            self.get_kmers(k)

        kmers = self.kmers[k]

        frequencies = Counter({})
        for i in range(0, len(sequence) - (k - 1)):
            kmer = sequence[i : i + k]

            if self.consider_rev_comps:
                if kmer in kmers:
                    frequencies[kmer] += 1
                else:
                    frequencies[rev_comp(kmer)] += 1
            else:
                frequencies[kmer] += 1

        if dist_metric_safe:
            # we don't want all kmer freq values to be zero. so the distance
            # metrics wouldn't go crazy. instead we fill it with 1. which
            # doesn't affect relative distances.
            if not len(frequencies):
                frequencies = dict(list(zip(kmers, [1] * len(kmers))))

        return frequencies
