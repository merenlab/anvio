# coding: utf-8
"""Interface to fastANI."""

import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class FastANIDriver:
    """ Super class for all fastANI usage """
    def __init__(self, args={}, run=terminal.Run(), progress=terminal.Progress(), program_name='fastANI'):
        self.run = run
        self.progress = progress
        self.program_name = program_name

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.kmer_size = A('fastani_kmer_size') or 16
        self.fragment_length = A('fragment_length') or 3000
        self.min_num_fragments = A('min_num_fragments') or 50
        self.num_threads = A('num_threads') or 1
        self.log_file_path = os.path.abspath(A('log_file') or filesnpaths.get_temp_file_path())
        self.quiet = A('quiet')

        self.check_programs()

        self.run.warning("Anvi'o will use 'fastANI' by Jain et al. (DOI: 10.1038/s41467-018-07641-9) to compute ANI. \
                          If you publish your findings, please do not forget to properly credit their work.", lc='green', header="CITATION")

        self.run.info('[fastANI] Kmer size', self.kmer_size)
        self.run.info('[fastANI] Fragment length', self.fragment_length)
        self.run.info('[fastANI] Min num of fragments', self.min_num_fragments)
        self.run.info('[fastANI] Num threads to use', self.num_threads)
        self.run.info('[fastANI] Log file path', self.log_file_path, nl_after=1)


    def check_programs(self):
        utils.is_program_exists(self.program_name)


    def run_command(self, input_path, *args):
        """ Run the command

        Parameters
        ==========
        input_path : string or Path-like
            The directory in which the binary will be executed
        *args :
            all arguments needed to generate the commmand, i.e.
            arguments passed to get_command(). See Notes.

        Notes
        =====
        - This command assumes all child classes have a method called `get_command(...)`.
          As an example, refer to ManyToMany.get_command.
        """
        # backup the old working directory before changing the directory
        old_wd = os.getcwd()
        os.chdir(input_path)

        self.progress.new('fastANI')
        self.progress.update('Running ...')
        exit_code = utils.run_command(self.get_command(*args), self.log_file_path)
        self.progress.end()

        if int(exit_code):
            raise ConfigError("fastANI returned with non-zero exit code, there may be some errors. \
                               please check the log file for details.")

        # restore old working directory
        os.chdir(old_wd)

        return


class ManyToMany(FastANIDriver):
    """ Compare many genomes to many other genomes """
    def __init__(self, args={}, run=terminal.Run(), progress=terminal.Progress(), program_name='fastANI'):
        FastANIDriver.__init__(self, args, run, progress, program_name)


    def get_command(self, fastANI_input_file):
        return [self.program_name,
                '--ql', fastANI_input_file,
                '--rl', fastANI_input_file,
                '-k', self.kmer_size,
                '--fragLen', self.fragment_length,
                '--minFrag', self.min_num_fragments,
                '-t', self.num_threads,
                '-o', 'many_to_many_output.txt']
