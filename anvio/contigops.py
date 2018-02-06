# -*- coding: utf-8
# pylint: disable=line-too-long
"""Classes and functions for handling, storing, and retrieving atomic data from contigs and splits"""

import anvio

from anvio.sequence import Coverage
from anvio.terminal import Run, Progress
from anvio.variability import VariablityTestFactory

import anvio.tables as t

run = Run()
progress = Progress()
progress.verbose = False

from anvio.variability import ColumnProfile


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = ["Faruk Uzun"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


variability_test_class_default = VariablityTestFactory(params={'b': 2, 'm': 1.45, 'c': 0.05})
variability_test_class_null = VariablityTestFactory(params=None) # get everything for every coverage level


def gen_split_name(parent_name, order):
    return '_'.join([parent_name, 'split', '%05d' % (order + 1)])


def get_atomic_data_dicts(sample_id, contigs):
    """Takes a list of contigops.Contig objects, and returns contigs and splits atomic data
       dictionaries"""
    atomic_data_contigs = {}
    atomic_data_splits = {}

    # this loop will get atomic_data information from Contig instanes and store them into the db
    # at once. this was broken down into about 10 functions, but this structure seems to be the most efficient
    # although it looks crappy:
    for contig in contigs:
        contig_atomic_data = contig.get_atomic_data_dict()

        for split in contig.splits:
            atomic_data_contigs[split.name] = {'contig': contig.name}
            for atomic_data_field in t.atomic_data_table_structure[1:]:
                atomic_data_contigs[split.name][atomic_data_field] = contig_atomic_data[atomic_data_field]

        # contig is done, deal with splits in it:
        for split in contig.splits:
            split_atomic_data = split.get_atomic_data_dict()
            atomic_data_splits[split.name] = {'contig': split.name}
            for atomic_data_field in t.atomic_data_table_structure[1:]:
                atomic_data_splits[split.name][atomic_data_field] = split_atomic_data[atomic_data_field]

    return atomic_data_splits, atomic_data_contigs


class Contig:
    def __init__(self, name):
        self.name = name
        self.sequence = None
        self.parent = None
        self.splits = []
        self.length = 0
        self.abundance = 0.0
        self.coverage = Coverage()

        self.min_coverage_for_variability = 10
        self.skip_SNV_profiling = False
        self.report_variability_full = False


    def get_atomic_data_dict(self):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'mean_coverage_Q2Q3': self.coverage.mean_Q2Q3,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'detection': self.coverage.detection,
             'abundance': self.abundance,
             'variability': sum(s.auxiliary.variation_density for s in self.splits) if not self.skip_SNV_profiling else None,
             '__parent__': None}

        return d


    def analyze_coverage(self, bam):
        contig_coverage = []

        counter = 1
        for split in self.splits:
            split.coverage = Coverage()
            split.coverage.run(bam, split)
            contig_coverage.extend(split.coverage.c)

            counter += 1

        self.coverage.process_c(contig_coverage)


    def analyze_auxiliary(self, bam):
        counter = 1
        for split in self.splits:
            split.auxiliary = Auxiliary(split,
                                        bam,
                                        parent_outlier_positions=self.coverage.outlier_positions,
                                        min_coverage=self.min_coverage_for_variability,
                                        report_variability_full=self.report_variability_full)

            counter += 1


class Split:
    def __init__(self, name, sequence, parent, order, start=0, end=0):
        self.name = name
        self.sequence = sequence
        self.parent = parent
        self.end = end
        self.order = order
        self.start = start
        self.length = end - start
        self.explicit_length = 0
        self.abundance = 0.0
        self.column_profiles = {}
        self.auxiliary = None

    def get_atomic_data_dict(self):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'mean_coverage_Q2Q3': self.coverage.mean_Q2Q3,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'detection': self.coverage.detection,
             'abundance': self.abundance,
             'variability': self.auxiliary.variation_density if self.auxiliary else None,
             '__parent__': self.parent}

        return d


class Auxiliary:
    def __init__(self, split, bam, parent_outlier_positions, min_coverage=10, report_variability_full=False):
        self.v = []
        self.rep_seq = ''
        self.split = split
        self.variation_density = 0.0
        self.parent_outlier_positions = parent_outlier_positions
        self.competing_nucleotides = {}
        self.min_coverage = min_coverage
        self.column_profile = self.split.column_profiles
        self.report_variability_full = report_variability_full

        self.run(bam)


    def run(self, bam):
        ratios = []

        for pileupcolumn in bam.pileup(self.split.parent, self.split.start, self.split.end):
            pos_in_contig = pileupcolumn.pos
            if pos_in_contig < self.split.start or pos_in_contig >= self.split.end:
                continue

            valid_nts = [pileupread.alignment.seq[pileupread.query_position] for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip]

            coverage = len(valid_nts)
            if coverage < self.min_coverage:
                continue

            column = ''.join(valid_nts)

            pos_in_split = pos_in_contig - self.split.start
            base_in_contig = self.split.sequence[pos_in_split]

            cp = ColumnProfile(column,
                               reference=base_in_contig,
                               coverage=coverage,
                               split_name=self.split.name,
                               pos=pos_in_split,
                               test_class=variability_test_class_null if self.report_variability_full else variability_test_class_default).profile

            if cp['worth_reporting']:
                ratios.append((cp['departure_from_reference'], cp['coverage']), )
                cp['pos_in_contig'] = pos_in_contig
                cp['cov_outlier_in_split'] = pos_in_split in self.split.coverage.outlier_positions
                cp['cov_outlier_in_contig'] = pos_in_contig in self.parent_outlier_positions
                self.column_profile[pos_in_contig] = cp

        # variation density = number of SNVs per kb
        self.variation_density = len(ratios) * 1000.0 / self.split.length

        for i in range(self.split.start, self.split.end):
            if i in self.column_profile:
                self.rep_seq += self.column_profile[i]['reference']
                self.v.append(self.column_profile[i]['departure_from_reference'])
                self.competing_nucleotides[self.column_profile[i]['pos']] = self.column_profile[i]['competing_nts']
            else:
                self.rep_seq += 'N'
                self.v.append(0)
