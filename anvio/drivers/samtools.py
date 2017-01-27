# coding: utf-8
"""Interface to Samtools."""

import os

import anvio
import numpy
import subprocess
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "Özcan Esen"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


progress = terminal.Progress()
pp = terminal.pretty_print


class Samtools:
    def __init__(self, contig_lengths, contig_names, skip_SNV_profiling, progress=progress):
        utils.is_program_exists('samtools')
        self.progress = progress
        self.skip_SNV_profiling = skip_SNV_profiling
        self.num_contigs = len(contig_names)

        self.contig_name_to_length = {}
        for i in range(self.num_contigs):
            self.contig_name_to_length[contig_names[i]] = contig_lengths[i]

        # we will store output in dictionaries below
        self.coverages = {}
        self.column_nucleotide_counts = {}

    def run(self, bam_file):
        self.progress.new('Reading coverage information from samtools')
        process = subprocess.Popen(["samtools", "mpileup", "-Q", "0", bam_file], 
            stdout=subprocess.PIPE, stderr=open(os.devnull, 'w'))

        while True:
            output = process.stdout.readline().decode()
            if output == '' and process.poll() is not None:
                break
            if output:
                output = output.split("\t")
                # note about output columns
                # 0 -> contig_name
                # 1 -> pos (index starts from 1, we subtract 1)
                # 3 -> coverage 
                # 4 -> column
                self.process(output[0], int(output[1]) - 1, int(output[3]), output[4])

        self.progress.end()

    def process(self, contig_name, pos, coverage, column):
        length = self.contig_name_to_length[contig_name]

        if not contig_name in self.coverages:
            self.coverages[contig_name] = numpy.zeros((length,), dtype=numpy.uint16)

            if not self.skip_SNV_profiling:
                self.column_nucleotide_counts[contig_name] = numpy.zeros((length,5), dtype=numpy.uint16)

            self.progress.update('Received %d of %d.' % (len(self.coverages), self.num_contigs))

        if (pos < 0 or pos >= length):
            return

        self.coverages[contig_name][pos] = coverage

        if not self.skip_SNV_profiling:
            for nucleotide in column:
                if nucleotide == 'A' or nucleotide == 'a':
                    self.column_nucleotide_counts[contig_name][pos][0] += 1
                elif nucleotide == 'T' or nucleotide == 't':
                    self.column_nucleotide_counts[contig_name][pos][1] += 1
                elif nucleotide == 'G' or nucleotide == 'g':
                    self.column_nucleotide_counts[contig_name][pos][2] += 1
                elif nucleotide == 'C' or nucleotide == 'c':
                    self.column_nucleotide_counts[contig_name][pos][3] += 1
                elif nucleotide == 'N' or nucleotide == 'n':
                    self.column_nucleotide_counts[contig_name][pos][4] += 1

