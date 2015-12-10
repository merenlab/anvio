# -*- coding: utf-8
"""Classes and functions for handling, storing, and retrieving atomic data from contigs and splits"""

import anvio

from anvio.sequence import Coverage
from anvio.terminal import Run, Progress
from anvio.terminal import pretty_print as pp
from anvio.variability import VariablityTestFactory

import anvio.tables as t

run = Run()
progress = Progress()
progress.verbose = False

try:
    from anvio.columnprofile import ColumnProfile
except ImportError:
    run.info_single('C extension for ColumnProfile failed to load, falling back to the Python implementation...', mc = 'gray', nl_after = 1)
    from anvio.variability import ColumnProfile


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = ["Faruk Uzun"]
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


variability_test_class_default = VariablityTestFactory(params = {'b': 2, 'm': 1.45, 'c': 0.05})
variability_test_class_null = VariablityTestFactory(params = None) # get everything for every coverage level


def set_contigs_abundance(contigs):
    """takes a list of contigs (of Contig class) and sets abundance values. a better way to do this is to implement
       a Contigs wrapper .. maybe later."""

    # first calculate the mean coverage
    total_length_of_all_contigs = sum([c.length for c in contigs.values()])
    total_coverage_values_for_all_contigs = sum([c.coverage.mean * c.length for c in contigs.values()])
    overall_mean_coverage = total_coverage_values_for_all_contigs / total_length_of_all_contigs

    # set normalized abundance factor for each contig
    for contig in contigs:
        contigs[contig].abundance = contigs[contig].coverage.mean / overall_mean_coverage if overall_mean_coverage else 0
        for split in contigs[contig].splits:
            split.abundance = split.coverage.mean / overall_mean_coverage if overall_mean_coverage else 0


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

        self.min_coverage_for_variability = 10
        self.report_variability_full = False


    def get_atomic_data_dict(self):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'normalized_coverage': self.coverage.normalized,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'portion_covered': self.coverage.portion_covered,
             'abundance': self.abundance,
             'variability': sum(s.auxiliary.variation_density for s in self.splits),
             '__parent__': None}

        return d


    def analyze_coverage(self, bam, progress):
        contig_coverage = []
        num_splits = len(self.splits)
        counter = 1

        for split in self.splits:
            progress.update('Coverage (split: %d of %d)' % (counter, num_splits))

            split.coverage = Coverage()
            split.coverage.run(bam, split)
            contig_coverage.extend(split.coverage.c)

            counter += 1

        self.coverage.process_c(contig_coverage)


    def analyze_auxiliary(self, bam, progress):
        num_splits = len(self.splits)
        counter = 1

        for split in self.splits:
            progress.update('Auxiliary stats (split: %d of %d) CMC: %.1f :: SMC: %.1f'\
                                 % (counter, num_splits, self.coverage.mean, split.coverage.mean))

            split.auxiliary = Auxiliary(split, bam, min_coverage = self.min_coverage_for_variability,
                                        report_variability_full = self.report_variability_full)

            counter += 1


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

    def get_atomic_data_dict(self):
        d = {'std_coverage': self.coverage.std,
             'mean_coverage': self.coverage.mean,
             'normalized_coverage': self.coverage.normalized,
             'max_normalized_ratio': 1.0,
             'relative_abundance': 1.0,
             'portion_covered': self.coverage.portion_covered,
             'abundance': self.abundance,
             'variability': self.auxiliary.variation_density,
             '__parent__': self.parent}

        return d


class Auxiliary:
    def __init__(self, split, bam, min_coverage = 10, report_variability_full = False):
        self.rep_seq = ''
        self.min_coverage = min_coverage
        self.report_variability_full = report_variability_full 
        self.split = split
        self.column_profile = self.split.column_profiles
        self.variation_density = 0.0
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

            column = ''.join([pileupread.alignment.seq[pileupread.query_position] for pileupread in pileupcolumn.pileups if not pileupread.is_del and not pileupread.is_refskip])

            cp = ColumnProfile(column,
                               coverage = coverage,
                               split_name = self.split.name,
                               pos = pileupcolumn.pos - self.split.start,
                               test_class = variability_test_class_null if self.report_variability_full else variability_test_class_default).profile

            if cp['n2n1ratio']:
                ratios.append((cp['n2n1ratio'], cp['coverage']), )
                self.column_profile[pileupcolumn.pos] = cp

        # variation density = number of SNPs per kb
        self.variation_density = len(ratios) * 1000.0 / self.split.length

        for i in range(self.split.start, self.split.end):
            if self.column_profile.has_key(i):
                self.rep_seq += self.column_profile[i]['consensus']
                self.v.append(self.column_profile[i]['n2n1ratio'])
                self.competing_nucleotides[self.column_profile[i]['pos']] = self.column_profile[i]['competing_nts']
            else:
                self.rep_seq += 'N'
                self.v.append(0)


class AtomicContigSplitData:
    def __init__(self, p=progress):
        self.atomic_data_contigs = {}
        self.atomic_data_splits = {}
        self.progress = p


    def store_atomic_data_for_contigs_and_splits(self, sample_id, contigs, db):
        self.progress.new('Storing atomic_data')

        num_contigs = pp(len(contigs))
        cur_contig = 1

        # this loop will get atomic_data information from Contig instanes and store them into the db
        # at once. this was broken down into about 10 functions, but this structure seems to be the most efficient
        # although it looks crappy:
        for contig_name in contigs:
            self.progress.update("Processing contig %s of %s" % (pp(cur_contig), num_contigs))
            contig = contigs[contig_name]
            contig_atomic_data = contig.get_atomic_data_dict()

            self.atomic_data_contigs[contig.name] = {'contig': contig.name}
            for atomic_data_field in t.atomic_data_table_structure[1:]:
                self.atomic_data_contigs[contig.name][atomic_data_field] = contig_atomic_data[atomic_data_field]

            # contig is done, deal with splits in it:
            for split in contig.splits:
                split_atomic_data = split.get_atomic_data_dict()
                self.atomic_data_splits[split.name] = {'contig': split.name}
                for atomic_data_field in t.atomic_data_table_structure[1:]:
                    self.atomic_data_splits[split.name][atomic_data_field] = split_atomic_data[atomic_data_field]


        self.progress.update("Generating tables ...")
        gen_atomic_data_tables_for_contigs_and_splits(self.atomic_data_splits, self.atomic_data_contigs, db)
        self.progress.end()



def gen_atomic_data_tables_for_contigs_and_splits(atomic_data_splits, atomic_data_contigs, db):
    # all objects are ready, creating tables next.
    db.create_table('atomic_data_splits', t.atomic_data_table_structure, t.atomic_data_table_types)
    db_entries = [tuple([split] + [atomic_data_splits[split][h] for h in t.atomic_data_table_structure[1:]]) for split in atomic_data_splits]
    db._exec_many('''INSERT INTO atomic_data_splits VALUES (?,?,?,?,?,?,?,?,?,?)''', db_entries)

    db.create_table('atomic_data_contigs', t.atomic_data_table_structure, t.atomic_data_table_types)
    db_entries = [tuple([split] + [atomic_data_contigs[atomic_data_splits[split]['__parent__']][h] for h in t.atomic_data_table_structure[1:]]) for split in atomic_data_splits]
    db._exec_many('''INSERT INTO atomic_data_contigs VALUES (?,?,?,?,?,?,?,?,?,?)''', db_entries)

    db.commit()
