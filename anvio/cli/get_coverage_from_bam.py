#!/usr/bin/env python
# -*- coding: utf-8

import sys
import numpy as np
import pandas as pd

import anvio
import anvio.bamops as bamops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['ekiefl']
__provides__ = ["coverages-txt"]
__requires__ = ["bam-file", "collection-txt",]
__description__ = ("Get nucleotide-level, contig-level, or bin-level coverage values from a BAM file "
                   "very rapidly. For other anvi'o programs that are designed to profile BAM files, "
                   "see `anvi-profile` and `anvi-profile-blitz`")


class ComputeCoverage(object):
    def __init__(self, args, run=terminal.Run(), progress=terminal.Progress()):
        self.args = args
        self.run = run
        self.progress = progress

        self.progress.new('Initializing')

        self.sanity_check()
        self.f = self.init_output()

        self.bam = bamops.BAMFileObject(self.args.bam_file, 'rb')

        # lengths is {contig_name: contig_length} of each contig in the bam
        self.lengths = dict(zip(self.bam.references, self.bam.lengths))

        self.progress.end()


    def sanity_check(self):
        has_arg = lambda x: True if self.args.__dict__.get(x) else False

        if sum([has_arg(x) for x in ('contig_name', 'contigs_of_interest', 'collection_txt')]) != 1:
            raise ConfigError("Provide input for exactly one of --contig-name, --contigs-of-interest, or --collections-txt")

        if self.args.contig_name:
            self.input_type = 'contig_name'
        elif self.args.contigs_of_interest:
            self.input_type = 'contigs_of_interest'
        else:
            self.input_type = 'collection_txt'

        filesnpaths.is_file_exists(self.args.bam_file)

        self.run.info('BAM file', self.args.bam_file)
        self.run.info('Input type', self.input_type)
        self.run.info('Coverage will be reported for each', self.args.method)

        if (self.input_type == 'contig_name' or self.input_type == 'contigs_of_interest') and self.args.method == 'bin':
            raise ConfigError("Choosing --method bin is only available with --collection-txt")


    def get_column_names(self):
        if self.input_type == 'contig_name':
            if self.args.method == 'pos':
                colnames = ('pos', 'coverage')
            elif self.args.method == 'contig':
                colnames = ('contig', 'coverage')

        elif self.input_type == 'contigs_of_interest':
            if self.args.method == 'pos':
                colnames = ('contig', 'pos', 'coverage')
            elif self.args.method == 'contig':
                colnames = ('contig', 'coverage')

        elif self.input_type == 'collection_txt':
            if self.args.method == 'pos':
                colnames = ('bin', 'contig', 'pos', 'coverage')
            elif self.args.method == 'contig':
                colnames = ('bin', 'contig', 'coverage')
            elif self.args.method == 'bin':
                colnames = ('bin', 'coverage')

        return colnames


    def is_contigs_in_bam(self, contigs):
        """Raises error if any contigs are not found in the bam"""

        for contig_name in contigs:
            if contig_name not in self.bam.references:
                raise ConfigError("%s is not in your bam file. Make sure you're using contig names "
                                  "and not split names." % contig_name)


    def init_output(self):
        # write the headers
        f = open(self.args.output, 'w')
        f.write('\t'.join(self.get_column_names()) + '\n')
        f.close()

        # reopen file in append mode
        f = open(self.args.output, 'ab')

        return f


    def write_cache_to_output(self, array):
        np.savetxt(self.f, array, delimiter='\t', fmt="%s")


    def get_nt_coverage_of_contig(self, contig_name):
        """Returns nucleotide array of coverage values

        Parameters
        ==========
        contig_name : str
        """

        self.progress.increment()
        self.progress.update('%s (%d/%d)' % (contig_name, self.progress.progress_current_item, self.progress.progress_total_items))

        cov = bamops.Coverage()
        cov.run(self.bam, contig_name, read_iterator='fetch', skip_coverage_stats=True)

        return cov.c


    def run_pos_method_for_contig(self, contig_name, columns):
        nt_cov = self.get_nt_coverage_of_contig(contig_name)

        array_tuple = tuple()
        for column in columns:
            array_tuple += (np.array([column] * len(nt_cov)), )

        data_cache = np.vstack(
            array_tuple + (np.arange(len(nt_cov)), nt_cov)
        ).T

        self.write_cache_to_output(data_cache)


    def run_contig_method_for_contig(self, contig_name, columns):
        nt_cov = self.get_nt_coverage_of_contig(contig_name)

        data_cache = np.array([
            columns + (np.mean(nt_cov), )
        ])

        self.write_cache_to_output(data_cache)


    def process(self):
        if self.input_type == 'contig_name':
            self.process_contig()
        elif self.input_type == 'contigs_of_interest':
            self.process_contigs_of_interest()
        elif self.input_type == 'collection_txt':
            self.process_collection_txt()

        self.progress.end()
        self.f.close()
        self.run.info_single('Output file %s has been written' % self.args.output, nl_before=1)


    def process_contig(self):
        m = self.args.method
        contig_name = self.args.contig_name

        self.progress.new('Processing', progress_total_items=1)

        if m == 'contig':
            self.run_contig_method_for_contig(contig_name, (contig_name, ))
        else: # m == 'pos'
            self.run.warning("Warning, your file will contain %d lines" % self.lengths[contig_name])

            self.run_pos_method_for_contig(contig_name, tuple())


    def process_contigs_of_interest(self):
        m = self.args.method
        contigs_of_interest = self.setup_contigs_of_interest()

        self.progress.new('Processing', progress_total_items=len(contigs_of_interest))

        if m == 'contig':
            for contig_name in contigs_of_interest:
                self.run_contig_method_for_contig(contig_name, (contig_name, ))

        else: # m == 'pos'
            for contig_name in contigs_of_interest:
                self.run_pos_method_for_contig(contig_name, (contig_name, ))


    def setup_contigs_of_interest(self):
        """Check contigs of interest file and then return list of contigs"""
        filesnpaths.is_file_exists(self.args.contigs_of_interest)
        contigs_of_interest = [x.strip() for x in open(self.args.contigs_of_interest).readlines()]

        if not self.args.skip_contigs_check:
            self.is_contigs_in_bam(contigs_of_interest)

        return contigs_of_interest


    def process_collection_txt(self):
        """Compute average when --collection-txt is chosen

        Because each bin could contain an extremely large number of nucleotide position, which could
        be too much information to allocate in memory, if the chosen method is 'bin', i.e the user wants
        the mean coverage for each bin, bin averages are actually weighted averages of mean contig
        coverages, where each each is the contig length.
        """

        m = self.args.method
        collection = self.setup_collection_txt()

        self.progress.new('Processing', progress_total_items=len(collection['contig_name']))

        if m == 'bin':
            for bin_name, contigs_in_bin in collection.groupby('bin')['contig_name']:
                # allocate empty array to store each contig's coverage
                contig_coverages = np.zeros(len(contigs_in_bin))

                # get a contig length array for weighting
                contig_lengths = np.array([self.lengths[c] for c in contigs_in_bin])

                for i, contig_name in enumerate(contigs_in_bin):
                    contig_coverages[i] = np.mean(self.get_nt_coverage_of_contig(contig_name))

                data_cache = np.array([[bin_name, np.average(contig_coverages, weights=contig_lengths/np.sum(contig_lengths))]])

                self.write_cache_to_output(data_cache)

        elif m == 'contig':
            for bin_name, contigs_in_bin in collection.groupby('bin')['contig_name']:
                for contig_name in contigs_in_bin:
                    self.run_contig_method_for_contig(contig_name, (bin_name, contig_name))

        else: # m == 'pos'
            for bin_name, contigs_in_bin in collection.groupby('bin')['contig_name']:
                for i, contig_name in enumerate(contigs_in_bin):
                    self.run_pos_method_for_contig(contig_name, (bin_name, contig_name))


    def setup_collection_txt(self):
        filesnpaths.is_file_exists(self.args.collection_txt)
        collection = pd.read_csv(self.args.collection_txt, sep='\t', names=('contig_name', 'bin'))

        if not self.args.skip_contigs_check:
            self.is_contigs_in_bam(collection['contig_name'])

        return collection


@terminal.time_program
def main():
    args = get_args()

    try:
        c = ComputeCoverage(args, terminal.Run(), terminal.Progress())
        c.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('REQUIRED', 'Declare your BAM file here')
    groupA.add_argument('-b', '--bam-file', required=True, help="Sorted and indexed BAM file to analyze.")

    groupB = parser.add_argument_group('OPTION #1', "This is the first and simplest option. Provide a contig name")
    groupB.add_argument('-c', '--contig-name', default=None, help="The name of a single contig")

    groupC = parser.add_argument_group('OPTION #2', "Use this to characterize coverage for a list of contigs")
    groupC.add_argument('-l', '--contigs-of-interest', default=None, help="Provide here a file where each line is a contig name.")

    groupD = parser.add_argument_group('OPTION #3', "Use this to characterize coverage for a collection of contig sets (bins)")
    groupD.add_argument('-C', '--collection-txt', default=None, help="Provide a collection text file. The first "
        "column should be contig names and the second column should be the bin to which the contig belongs. If "
        "you have a collection from a profile database, you can export it in this format with anvi-export-collection.")

    groupE = parser.add_argument_group('METHOD',  "Do you want to report coverage at a nucleotide level? Contig averages? Bin averages? Pick the method here.")
    groupE.add_argument('-m', '--method', choices=('pos','contig','bin'), required=True, help="If pos, each nucleotide position "
        "will be reported (valid for OPTION #1, #2, #3). If contig, report contains contig averages (valid for OPTION #2, #3). "
        "If bin, report contains bin averages (valid for OPTION #3).")

    groupF = parser.add_argument_group('OUTPUT',  "Your output file is decided here. Keep in mind if you use --method pos, "
        "this file will contain as many lines as there are nucleotides defined by your input option")
    groupF.add_argument('-o', '--output', required=True, help="Output tab-delimited file path. Will overwrite existing files.")

    groupG = parser.add_argument_group('EXTRAS', "All the misfits")
    groupG.add_argument('--skip-contigs-check', action='store_true', help="Checking to see that your collection text or "
        "contigs of interest file has correct names can take a really long time if you have a large enough number of "
        "contigs. Use this flag to forego checking, and find out the hard way.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
