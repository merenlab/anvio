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


from PaPi.variability import ColumnProfile, VariablityTestFactory
from PaPi.sequence import Coverage, Composition


variability_test_class = VariablityTestFactory()

def set_contigs_abundance(contigs):
    """takes a list of contigs (of Contig class) and sets abundance values. a better way to do this is to implement
       a Contigs wrapper .. maybe later."""

    # first calculate the mean coverage
    overall_mean_coverage = sum([c.coverage.mean * c.length for c in contigs.values()]) / sum(c.length for c in contigs.values())

    # set normalized abundance factor for each contig
    for contig in contigs:
        contigs[contig].abundance = contigs[contig].coverage.mean / overall_mean_coverage
        for split in contigs[contig].splits:
            split.abundance = split.coverage.mean / overall_mean_coverage


def gen_split_name(parent_name, order):
    return '_'.join([parent_name, 'split', '%05d' % (order + 1)])


class Contig:
    def __init__(self, name):
        self.name = name
        self.parent = None
        self.splits = []
        self.length = 0
        self.abundance = 0.0
        self.coverage = Coverage()
        self.composition = None


    def get_metadata_dict(self):
        d = {'length': self.length,
             'GC_content': self.composition.GC_content,
             'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'normalized_coverage': self.coverage.normalized,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'portion_covered': self.coverage.portion_covered,
             'abundance': self.abundance,
             'variability': sum(s.auxiliary.variability_score for s in self.splits),
             '__parent__': None}

        return d


    def analyze_coverage(self, bam, progress):
        contig_coverage = []
        for split in self.splits:
            progress.update('Coverage (split: %d of %d)' % (split.order, len(self.splits)))
            split.coverage = Coverage()
            split.coverage.run(bam, split)
            contig_coverage.extend(split.coverage.c)

        self.coverage.process_c(contig_coverage)


    def analyze_auxiliary(self, bam, progress):
        for split in self.splits:
            progress.update('Auxiliary stats (split: %d of %d) CMC: %.1f :: SMC: %.1f'\
                                 % (split.order, len(self.splits), self.coverage.mean, split.coverage.mean))
            split.auxiliary = Auxiliary(split, bam)


    def analyze_composition(self, bam, progress):
        for split in self.splits:
            progress.update('Composition (split: %d of %d)' % (split.order, len(self.splits)))
            split.composition = Composition(split.auxiliary.rep_seq)
        progress.update('Composition (split: %d of %d)' % (split.order, len(self.splits)))
        self.composition = Composition(self.get_rep_seq())


    def get_rep_seq(self):
        return ''.join([s.auxiliary.rep_seq for s in self.splits])



class Split:
    def __init__(self, name, parent, order, start = 0, end = 0):
        self.name = name
        self.parent = parent
        self.end = end
        self.order = order
        self.start = start
        self.length = end - start
        self.explicit_length = 0
        self.abundance = 0.0
        self.column_profiles = {}

    def get_metadata_dict(self):
        d = {'length': self.length,
             'GC_content': self.composition.GC_content,
             'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'normalized_coverage': self.coverage.normalized,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'portion_covered': self.coverage.portion_covered,
             'abundance': self.abundance,
             'variability': self.auxiliary.variability_score,
             '__parent__': self.parent}

        return d


class Auxiliary:
    def __init__(self, split, bam, min_coverage = 10):
        self.rep_seq = ''
        self.min_coverage = min_coverage
        self.split = split
        self.column_profile = self.split.column_profiles
        self.variability_score = 0.0
        self.v = []
        self.competing_nucleotides = {}

        self.run(bam)


    def run(self, bam):
        ratios = []

        for pileupcolumn in bam.pileup(self.split.parent, self.split.start, self.split.end):
            if pileupcolumn.pos < self.split.start or pileupcolumn.pos >= self.split.end:
                continue

            coverage = pileupcolumn.n
            if coverage < self.min_coverage:
                continue

            column = ''.join([pileupread.alignment.seq[pileupread.qpos] for pileupread in pileupcolumn.pileups])

            cp = ColumnProfile(column,
                               coverage = coverage,
                               split_name = self.split.name,
                               pos = pileupcolumn.pos - self.split.start,
                               test_class = variability_test_class).profile

            if cp['n2n1ratio']:
                ratios.append((cp['n2n1ratio'], cp['coverage']), )
                self.column_profile[pileupcolumn.pos] = cp

        # take top 50 based on n2n1ratio, then take top 25 of those with highest coverage:
        variable_positions_with_high_cov = sorted([(x[1], x[0]) for x in sorted(ratios, reverse=True)[0:50]], reverse=True)[0:25]
        self.variability_score = sum([x[1] * x[0] for x in variable_positions_with_high_cov])

        for i in range(self.split.start, self.split.end):
            if self.column_profile.has_key(i):
                self.rep_seq += self.column_profile[i]['consensus']
                self.v.append(self.column_profile[i]['n2n1ratio'])
                self.competing_nucleotides[self.column_profile[i]['pos']] = self.column_profile[i]['competing_nts']
            else:
                self.rep_seq += 'N'
                self.v.append(0)


