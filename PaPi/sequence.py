# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

'''Primitive classes for basic DNA sequence properties.'''


import copy
import numpy


class Composition:
    def __init__(self, sequence):
        self.sequence = sequence
        self.A = 0
        self.T = 0
        self.C = 0
        self.G = 0
        self.N = 0
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
        self.normalized = 0.0


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
        self.min = numpy.min(c)
        self.max = numpy.max(c)
        self.median = numpy.median(c)
        self.mean = numpy.mean(c)
        self.std = numpy.std(c)
        self.portion_covered = 1 - (float(c.count(0)) / len(c))
        self.normalized = self.mean * self.portion_covered