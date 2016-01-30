# -*- coding: utf-8

'''Primitive classes for basic DNA sequence properties.'''

import numpy
import collections

import anvio

__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


class Composition:
    def __init__(self, sequence):
        self.sequence = sequence
        self.GC_content = 0.0

        self.report()


    def report(self):
        s = self.sequence
        raw_length = len(s)
        
        self.A = s.count('A')
        self.T = s.count('T')
        self.C = s.count('C')
        self.G = s.count('G')
        self.N = raw_length - (self.A + self.T + self.C + self.G)
        length = raw_length - self.N
    
        if not length:
            # sequence is composed of only N's
            self.GC_content = 0.0
        else:
            self.GC_content = (self.G + self.C) * 1.0 / length


class Coverage:
    def __init__(self):
        self.c = []
        self.min = 0
        self.max = 0
        self.std = 0.0
        self.mean = 0.0
        self.median = 0.0
        self.portion_covered = 0.0
        self.mean_Q1Q3 = 0.0


    def run(self, bam, split):
        coverage_profile = {}
        for pileupcolumn in bam.pileup(split.parent, split.start, split.end):
            if pileupcolumn.pos < split.start or pileupcolumn.pos >= split.end:
                continue

            coverage_profile[pileupcolumn.pos] = pileupcolumn.n

        for i in range(split.start, split.end):
            if coverage_profile.has_key(i):
                self.c.append(coverage_profile[i])
            else:
                self.c.append(0)

        if self.c:
            split.explicit_length = len(self.c)
            self.process_c(self.c)

    def process_c(self, c):
        c = numpy.asarray(c)
        self.min = numpy.amin(c)
        self.max = numpy.amax(c)
        self.median = numpy.median(c)
        self.mean = numpy.mean(c)
        self.std = numpy.std(c)
        self.portion_covered = 1 - (float(collections.Counter(c)[0]) / len(c))

        if c.size < 4:
            self.mean_Q1Q3 = self.mean
        else:
            c.sort()
            Q = int(c.size * 0.25)
            Q1Q3 = c[Q:-Q]
            self.mean_Q1Q3 = numpy.mean(Q1Q3)
