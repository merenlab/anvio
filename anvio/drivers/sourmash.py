# coding: utf-8
"""Interface to sourmash"""

import os
import shutil

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Mahmoud Yousef"
__email__ = "mahmoudyousef@uchicago.edu"


class Sourmash:
    def __init__(self, args={}, run=terminal.Run(), progress=terminal.Progress(),
                 program_name='sourmash'):
        self.run = run
        self.progress = progress
        self.program_name = program_name
        self.check_program()

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.log_file_path = os.path.abspath(A('log_file') or filesnpaths.get_temp_file_path())
        self.num_threads = A('num_threads') or 1
        self.kmer_size = A('kmer_size') or 51
        self.scale = A('scale') or 1000

        self.run.warning("Anvi'o will use 'sourmash' by Brown et al. (DOI: 10.21105/joss.00027) to compute kmer sequences and determine mash distances.\
                         If you publish your findings, please do not forget to properly credit their work",
                         lc='green', header="CITATION")

        if self.num_threads != 1:
            self.run.warning("Anvi'o speaking: sourmash currently doesn't support multithreading.\
                             Anvi'o will have to reduce your number of threads to one :(")
            self.num_threads = 1

        self.run.info('[sourmash] Kmer size', self.kmer_size)
        self.run.info('[sourmash] Compression ratio', self.scale)
        self.run.info('[sourmash] Log file path', self.log_file_path, nl_after=1)


    def check_program(self):
        utils.is_program_exists(self.program_name)


    def process(self, input_path, fasta_files):
        # backup the old working directory before changing the directory
        old_wd = os.getcwd()
        os.chdir(input_path)

        self.progress.new('Sourmash')
        self.progress.update('Computing fasta signatures...')

        scale = '--scaled=%i' % self.scale
        compute_command = [self.program_name, 'compute',
                           '-k', self.kmer_size,
                           '-f', scale]
        compute_command.extend(fasta_files)

        exit_code = utils.run_command(compute_command, self.log_file_path)
        if int(exit_code):
            self.progress.end()
            raise ConfigError("sourmash returned with non-zero exit code, there may be some errors.\
                              Please check the log file `%s` for details. Offending command: \
                              `%s` ..." % (self.log_file_path, ' '.join([str(x) for x in compute_command[:7]])))

        self.progress.update('Computing distance matrix...')
        compare_command = [self.program_name, 'compare',
                           '-k', self.kmer_size,
                           '--csv', 'mash_distance.txt']
        for f in fasta_files:
            compare_command.append(f + ".sig")

        exit_code = utils.run_command(compare_command, self.log_file_path)
        if int(exit_code):
            self.progress.end()
            raise ConfigError("sourmash returned with non-zero exit code, there may be some errors.\
                              Please check the log file `%s` for details. Offending command: \
                              `%s` ..." % (self.log_file_path, ' '.join([str(x) for x in compute_command[:7]])))

        matrix = utils.get_TAB_delimited_file_as_dictionary('mash_distance.txt', indexing_field=-1, separator=',')

        self.progress.end()

        # restore old working directory
        os.chdir(old_wd)

        return matrix
