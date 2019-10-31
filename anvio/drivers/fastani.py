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
    """ Parent class for all fastANI usage

    Parameters
    ==========
    program_name : str, fastANI
        What is the fastANI binary called?
    """

    def __init__(self, program_name='fastANI', args={}, run=terminal.Run(), progress=terminal.Progress()):
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


    def check_programs(self):
        utils.is_program_exists(self.program_name)


    def add_run_info(self):
        self.run.info('[fastANI] Kmer size', self.kmer_size)
        self.run.info('[fastANI] Fragment length', self.fragment_length)
        self.run.info('[fastANI] Min num of fragments', self.min_num_fragments)
        self.run.info('[fastANI] Num threads to use', self.num_threads)
        self.run.info('[fastANI] Log file path', self.log_file_path, nl_after=1)



class ManyToMany(FastANIDriver):
    """ Compare many genomes (the queries) to many other genomes (the references)

    If you're confused by the concept of many queries and many references, it's worth mentioning
    that if you want to do ALL the pairwise comparisons of a set of genomes will call A, then both the
    queries and the references are A
    """

    def __init__(self, program_name='fastANI', args={}, run=terminal.Run(), progress=terminal.Progress()):
        FastANIDriver.__init__(self, program_name, args, run, progress)


    def gen_results_dict(self):
        d = {}


    def run_command(self, query_targets, reference_targets, output_path, run_dir=os.getcwd()):
        """ Run the command

        Parameters
        ==========
        query_targets : string or Path-like
            The query set of genomes (--ql). It should be a list of filepaths, one per line
        reference_targets : string or Path-like
            The reference set of genomes (--rl). It should be a list of filepaths, one per line
        output_path : string or Path-like
            Where should the raw fastANI output file be created? Relative to current working
            directory, not `run_dir`
        run_dir : string or Path-like, os.getcwd()
            Where should the command be run? The current directory is the default

        Returns
        =======
        results : dictionary
            results dictionary
        """

        self.add_run_info()

        command = [self.program_name,
                   '--ql', query_targets,
                   '--rl', reference_targets,
                   '-k', self.kmer_size,
                   '--fragLen', self.fragment_length,
                   '--minFrag', self.min_num_fragments,
                   '-t', self.num_threads,
                   '-o', output_path]

        self.progress.new('fastANI')
        self.progress.update('Many to Many ...')

        exit_code = utils.run_command(command, self.log_file_path, run_dir=run_dir)

        self.progress.end()

        if int(exit_code):
            raise ConfigError("fastANI returned with non-zero exit code, there may be some errors. \
                               please check the log file for details.")

        return self.gen_results_dict(output_path)



