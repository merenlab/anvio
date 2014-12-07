# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import copy
import numpy

from PaPi.kmers import KMers
from PaPi.variability import ColumnProfile

kmers = KMers()


class Contig:
    def __init__(self, name):
        self.name = name
        self.splits = []
        self.length = 0
        self.mean_coverage = 0.0
        self.tnf = {}


    def analyze_coverage(self, bam, progress):
        contig_coverage = []
        for split in self.splits:
            progress.update('Coverage (split: %d of %d)' % (split.order, len(self.splits)))
            split.coverage = Coverage(split, bam)
            contig_coverage.extend(split.coverage.c)

        self.mean_coverage = numpy.mean(contig_coverage)
        self.std_coverage = numpy.std(contig_coverage)


    def analyze_auxiliary(self, bam, progress):
        for split in self.splits:
            progress.update('Auxiliary stats (split: %d of %d) CMC: %.1f :: SMC: %.1f'\
                                 % (split.order, len(self.splits), self.mean_coverage, split.coverage.mean))
            split.auxiliary = Auxiliary(split, bam)


    def analyze_composition(self, bam, progress):
        for split in self.splits:
            progress.update('Composition (split: %d of %d)' % (split.order, len(self.splits)))
            split.composition = Composition(split.auxiliary.rep_seq)


    def get_rep_seq(self):
        return ''.join([s.auxiliary.rep_seq for s in self.splits])


    def analyze_tnf(self, progress):
        progress.update('TNF')
        rep_seq = self.get_rep_seq()
        self.tnf = kmers.get_kmer_frequency(rep_seq)


class Split:
    def __init__(self, parent, order, start, end):
        self.name = '_'.join([parent, 'split', '%05d' % order])
        self.parent = parent
        self.end = end
        self.order = order
        self.start = start
        self.length = end - start
        self.explicit_length = 0


class Auxiliary:
    def __init__(self, split, bam):
        self.rep_seq = ''
        self.split = split
        self.variability_score = 0.0
        self.v = []
        self.competing_nucleotides = {}

        self.run(bam)


    def run(self, bam):
        column_profile = {}
        ratios = []
        for pileupcolumn in bam.pileup(self.split.parent, self.split.start, self.split.end):
            if pileupcolumn.pos < self.split.start or pileupcolumn.pos >= self.split.end:
                continue

            column = ''.join([pileupread.alignment.seq[pileupread.qpos] for pileupread in pileupcolumn.pileups])

            column_profile[pileupcolumn.pos] = ColumnProfile(column, pileupcolumn.pos)

            c = column_profile[pileupcolumn.pos]
            ratios.append((c.n2n1ratio, c.coverage), )

        # take top 100 based on n2n1ratio, then take top 50 of those with highest coverage:
        variable_positions_with_high_cov = sorted([(x[1], x[0]) for x in sorted(ratios, reverse=True)[0:100]], reverse=True)[0:50]
        self.variability_score = sum([x[1] for x in variable_positions_with_high_cov])

        for i in range(self.split.start, self.split.end):
            if column_profile.has_key(i):
                self.rep_seq += column_profile[i].consensus_nucleotide
                self.v.append(column_profile[i].n2n1ratio)
                if column_profile[i].n2n1ratio > 0:
                    # here populating the dict with i - self.split start, instead if i, because I want to
                    # have a record of the relative position of the competing nucleotide withing the context
                    # of the split. i itself holds the position for the entire contig
                    self.competing_nucleotides[i - self.split.start] = column_profile[i].competing_nucleotides
            else:
                self.rep_seq += 'N'
                self.v.append(0)


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

        #self.composition['tnf'] = kmers.get_kmer_frequency(sequence)
    

class Coverage:
    def __init__(self, split, bam):
        self.split = split
        self.c = []
        self.min = 0
        self.max = 0
        self.std = 0.0
        self.mean = 0.0
        self.median = 0.0

        self.run(bam)


    def run(self, bam):
        coverage_profile = {}
        for pileupcolumn in bam.pileup(self.split.parent, self.split.start, self.split.end):
            if pileupcolumn.pos < self.split.start or pileupcolumn.pos >= self.split.end:
                continue

            coverage_profile[pileupcolumn.pos] = pileupcolumn.n

        for i in range(self.split.start, self.split.end):
            if coverage_profile.has_key(i):
                self.c.append(coverage_profile[i])
            else:
                self.c.append(0)

        if self.c:
            self.split.explicit_length = len(self.c) 
            self.min = numpy.min(self.c)
            self.max = numpy.max(self.c)
            self.median = numpy.median(self.c)
            self.mean = numpy.mean(self.c)
            self.std = numpy.std(self.c)
