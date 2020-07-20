# coding: utf-8
"""Interface to Bowtie2."""

import os

import anvio
import anvio.terminal as terminal
import anvio.utils as utils


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


class Bowtie2:
    def __init__(self, query_fasta=None, target_fasta=None, run=terminal.Run(), progress=terminal.Progress(), num_threads=1):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads

        utils.is_program_exists('bowtie2')
        utils.is_program_exists('bowtie2-build')

        self.k = None
        self.L = None
        self.N = None
        self.score_min = None
        self.rdg = None
        self.rfg = None
        self.options = [self.k, self.L, self.N, self.score_min, self.rdg, self.rfg]

        # At the moment, this only allows a FASTA query file rather than paired FASTQ files.
        self.query_fasta = query_fasta
        self.target_fasta = target_fasta

        if not self.target_fasta:
            self.target_fasta = self.query_fasta

        self.index = os.path.splitext(self.target_fasta)[0]
        self.sam = os.path.splitext(self.query_fasta)[0] + '.sam'

        if not self.run.log_file_path:
            self.run.log_file_path = 'bowtie2_log_file.txt'


    def build_index(self, output_index_files_prefix=None):
        self.progress.new("Bowtie2")
        self.progress.update("Building an index")

        if output_index_files_prefix:
            self.index = output_index_files_prefix
        output_index_files_prefix = self.index

        cmd_line = ['bowtie2-build', '-f', self.target_fasta, output_index_files_prefix]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('bowtie2-build cmd', ' '.join([str(x) for x in cmd_line]), quiet=(not anvio.DEBUG))
        self.run.info('Bowtie2 index file prefix', output_index_files_prefix, quiet=(not anvio.DEBUG))


    def align(self, output_sam_path=None):
        self.progress.new("Bowtie2")
        self.progress.update("Mapping using %d threads" % self.num_threads)

        single_dash_params = ['k', 'L', 'N']
        double_dash_params = ['score_min']
        optional_args = []
        for param in single_dash_params:
            arg = self.__getattribute__(param)
            if arg:
                optional_args.append('-' + param)
                optional_args.append(arg)
        for param in double_dash_params:
            arg = self.__getattribute__(param)
            if arg:
                optional_args.append('--' + param.replace('_', '-'))
                optional_args.append(arg)

        if output_sam_path:
            self.sam = output_sam_path
        output_sam_path = self.sam

        cmd_line = ['bowtie2',
                    '-x', self.index,
                    '-f', self.query_fasta,
                    '-S', output_sam_path,
                    '-p', self.num_threads] + optional_args

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('Bowtie2 bowtie2-align cmd', ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))
        self.run.info('Bowtie2 mapping result', output_sam_path, quiet=(not anvio.DEBUG))