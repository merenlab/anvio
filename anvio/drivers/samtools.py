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


class Samtools:
    def __init__(self, bam_file, contig_name, contig_length, skip_SNV_profiling):
        utils.is_program_exists('samtools')
        self.skip_SNV_profiling = skip_SNV_profiling
        self.contig_name = contig_name
        self.bam_file = bam_file
        self.contig_length = contig_length

        # we will store output in arrays below
        self.coverages = numpy.zeros((self.contig_length, ), dtype=numpy.uint16)
        
        if not self.skip_SNV_profiling:
            self.column_nucleotide_counts = numpy.zeros((self.contig_length, 5), dtype=numpy.uint16)

    def run(self):
        process = subprocess.Popen(["samtools", "mpileup", "-Q", "0", "-r", self.contig_name, self.bam_file], 
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
                self.process(int(output[1]) - 1, int(output[3]), output[4])

    def process(self, pos, coverage, column):
        if (pos < 0 or pos >= self.contig_length):
            return

        self.coverages[pos] = coverage

        if not self.skip_SNV_profiling:
            for nucleotide in column:
                if nucleotide == 'A' or nucleotide == 'a':
                    self.column_nucleotide_counts[pos][0] += 1
                elif nucleotide == 'T' or nucleotide == 't':
                    self.column_nucleotide_counts[pos][1] += 1
                elif nucleotide == 'G' or nucleotide == 'g':
                    self.column_nucleotide_counts[pos][2] += 1
                elif nucleotide == 'C' or nucleotide == 'c':
                    self.column_nucleotide_counts[pos][3] += 1
                elif nucleotide == 'N' or nucleotide == 'n':
                    self.column_nucleotide_counts[pos][4] += 1

