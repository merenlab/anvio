# coding: utf-8
"""Interface to fastANI."""

import os
import pandas as pd

from itertools import product

import anvio
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.utils.commandline import run_command
from anvio.utils.files import store_dataframe_as_TAB_delimited_file
from anvio.utils.system import is_program_exists


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
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
        self.min_fraction = A('min_fraction') or 0.25
        self.num_threads = A('num_threads') or 1
        self.log_file_path = os.path.abspath(A('log_file') or filesnpaths.get_temp_file_path())
        self.quiet = A('quiet')

        self.check_programs()

        self.run.warning("Anvi'o will use 'fastANI' by Jain et al. (DOI: 10.1038/s41467-018-07641-9) to compute ANI. "
                         "If you publish your findings, please do not forget to properly credit their work.", lc='green', header="CITATION")


    def check_programs(self):
        is_program_exists(self.program_name)


    def add_run_info(self):
        self.run.info('[fastANI] Kmer size', self.kmer_size)
        self.run.info('[fastANI] Fragment length', self.fragment_length)
        self.run.info('[fastANI] Min fraction of alignment', self.min_fraction)
        self.run.info('[fastANI] Num threads to use', self.num_threads)
        self.run.info('[fastANI] Log file path', self.log_file_path, nl_after=1)


    def gen_results_dict(self):
        results = {
            'ani': {},
            'alignment_fraction': {},
            'mapping_fragments': {},
            'total_fragments': {},
        }

        for _, row in self.fastANI_output.iterrows():
            query, reference = row['query'], row['reference']
            for result_name in results:
                if query not in results[result_name]:
                    results[result_name][query] = {}
                results[result_name][query][reference] = row[result_name]

        return results


class ManyToMany(FastANIDriver):
    """ Compare many genomes (the queries) to many other genomes (the references)

    If you're confused by the concept of many queries and many references, it's worth mentioning
    that if you want to do ALL the pairwise comparisons of a set of genomes will call A, then both the
    queries and the references are A
    """

    def __init__(self, program_name='fastANI', args={}, run=terminal.Run(), progress=terminal.Progress()):
        FastANIDriver.__init__(self, program_name, args, run, progress)


    def fill_missing_data(self, fastANI_output):
        """If alignment is insufficient, fastANI removes the results from the output. This puts them back"""

        all_query_reference_combinations = set(product(self.query_names, self.reference_names))
        present_query_reference_combinations = set([tuple(r) for r in fastANI_output[['query', 'reference']].values])
        missing_query_reference_combinations = all_query_reference_combinations - present_query_reference_combinations

        def missing_data_template(x, y):
            return {
                'query': x,
                'reference': y,
                'ani': 0,
                'mapping_fragments': 0,
                'total_fragments': 0,
                'alignment_fraction': 0,
            }

        null_rows = pd.DataFrame([missing_data_template(x, y) for x, y in missing_query_reference_combinations])
        fastANI_output = pd.concat([fastANI_output, null_rows], ignore_index=True, verify_integrity=False)
        fastANI_output = fastANI_output.sort_values(by=['query', 'reference'])
        fastANI_output = fastANI_output[['query', 'reference', 'ani', 'mapping_fragments', 'total_fragments', 'alignment_fraction']]
        return fastANI_output


    def load_output_as_dataframe(self, output_path, name_conversion_dict=None):
        names = ('query', 'reference', 'ani', 'mapping_fragments', 'total_fragments')
        fastANI_output = pd.read_csv(output_path, sep='\t', header=None, names=names)
        fastANI_output['alignment_fraction'] = fastANI_output['mapping_fragments'] / fastANI_output['total_fragments']
        fastANI_output['ani'] /= 100

        fastANI_output = self.fill_missing_data(fastANI_output)

        if name_conversion_dict:
            for target in ['query', 'reference']:
                fastANI_output[target] = fastANI_output[target].map(name_conversion_dict)

        return fastANI_output


    def get_all_query_and_reference_names(self, query_targets, reference_targets):
        query_names = set([x.strip() for x in open(query_targets).readlines()])
        reference_names = set([x.strip() for x in open(reference_targets).readlines()])
        return query_names, reference_names


    def run_command(self, query_targets, reference_targets, output_path, run_dir=os.getcwd(), name_conversion_dict=None):
        """Run the command

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
        name_conversion_dict : dict, None
            The keys of `results` are by default the file paths of the genomes, since that's what
            fastANI outputs. Pass an optional dictionary with <path>:<name> to convert the output.
            Note: this effects both the raw output in `output_path` and `results`

        Returns
        =======
        results : dictionary
            results dictionary
        """

        self.add_run_info()

        self.query_names, self.reference_names = self.get_all_query_and_reference_names(query_targets, reference_targets)

        command = [self.program_name,
                   '--ql', query_targets,
                   '--rl', reference_targets,
                   '-k', self.kmer_size,
                   '--fragLen', self.fragment_length,
                   '--minFraction', self.min_fraction,
                   '-t', self.num_threads,
                   '-o', output_path]

        self.progress.new('fastANI')
        self.progress.update('Many to Many ...')

        with utils.RunInDirectory(run_dir):
            exit_code = run_command(command, self.log_file_path)

        self.progress.end()

        if int(exit_code):
            raise ConfigError("fastANI returned with non-zero exit code, there may be some errors. "
                   "please check the log file for details.")

        self.fastANI_output = self.load_output_as_dataframe(output_path, name_conversion_dict)
        store_dataframe_as_TAB_delimited_file(self.fastANI_output, output_path)

        self.results = self.gen_results_dict()
        return self.results


