# coding: utf-8
"""Interface to sourmash"""

import os
import numpy as np
import pandas as pd
import shutil

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from scipy.stats import entropy, skew, kurtosis

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Mahmoud Yousef"
__email__ = "mahmoudyousef@uchicago.edu"


class Sourmash:
    """This calculates a single kmer signature, and computes similarities.

    Feel free to buff this to suit your needs
    """

    def __init__(self, args={}, run=terminal.Run(), progress=terminal.Progress(), program_name='sourmash'):
        self.run = run
        self.progress = progress
        self.program_name = program_name
        self.check_program()

        self.results = {}

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.log_file_path = os.path.abspath(A('log_file') or filesnpaths.get_temp_file_path())
        self.num_threads = A('num_threads') or 1
        self.kmer_size = A('kmer_size') or 51
        self.scale = A('scale') or 1000

        self.run.warning("Anvi'o will use 'sourmash' by Brown et al. (DOI: 10.21105/joss.00027) to compute kmer sequences and determine mash distances.\
                         If you publish your findings, please do not forget to properly credit their work",
                         lc='green', header="CITATION")

        if self.num_threads != 1:
            self.num_threads = 1
            self.run.warning("Anvi'o speaking: sourmash currently doesn't support multithreading.\
                             Anvi'o will have to reduce your number of threads to one :(")

        self.run.info('[sourmash] Log file path', self.log_file_path, nl_after=1)


    def check_program(self):
        utils.is_program_exists(self.program_name)


    def process(self, input_path, fasta_files):
        self.run.info('[sourmash] Kmer size', self.kmer_size, nl_before=1)
        self.run.info('[sourmash] Compression ratio', self.scale)

        report_name = 'kmer_%d_mash_similarity' % self.kmer_size

        # backup the old working directory before changing the directory
        old_wd = os.getcwd()
        os.chdir(input_path)
        if not os.path.exists('output'):
            os.mkdir('output')
        else:
            pass

        self.progress.new('Sourmash')
        self.progress.update('Computing fasta signatures for kmer=%d, scale=%d' % (self.kmer_size, self.scale))

        scale = '--scaled=%i' % self.scale
        compute_command = [self.program_name, 'compute',
                           '-k', self.kmer_size,
                           '-f', scale]
        compute_command.extend(fasta_files)

        exit_code = utils.run_command(compute_command, self.log_file_path, remove_log_file_if_exists=False)
        if int(exit_code):
            self.progress.end()
            raise ConfigError("sourmash returned with non-zero exit code, there may be some errors.\
                              Please check the log file `%s` for details. Offending command: \
                              `%s` ..." % (self.log_file_path, ' '.join([str(x) for x in compute_command[:7]])))

        self.progress.update('Computing similarity matrix for kmer=%d, scale=%d' % (self.kmer_size, self.scale))
        compare_command = [self.program_name, 'compare',
                           '-k', self.kmer_size,
                           '--csv', os.path.join('output', report_name + '.txt')]
        for f in fasta_files:
            compare_command.append(f + ".sig")

        exit_code = utils.run_command(compare_command, self.log_file_path, remove_log_file_if_exists=False)
        if int(exit_code):
            self.progress.end()
            raise ConfigError("sourmash returned with non-zero exit code, there may be some errors.\
                              Please check the log file `%s` for details. Offending command: \
                              `%s` ..." % (self.log_file_path, ' '.join([str(x) for x in compute_command[:7]])))

        self.results[report_name] = utils.get_TAB_delimited_file_as_dictionary(os.path.join('output', report_name + '.txt'),
                                                                               indexing_field=-1,
                                                                               separator=',')

        self.progress.end()

        # restore old working directory
        os.chdir(old_wd)

        return self.results


class IterateKmerSourmash(Sourmash):
    """ This class runs sourmash iteratively until it finds the 'best' kmer based on entropy maximization.

    Parameters
    ==========
    method : str, "scan"
        The method to determine maximal entropy. "scan" simply tries all of the kmer values
        within a range and picks the highest.

    Notes
    =====
    - Why on earth are you using this class?
    """

    def __init__(self, args, method='scan'):
        self.args = args
        self.method = method

        # Rename Sourmash.process to Sourmash_process
        self.Sourmash_process = Sourmash.process
        Sourmash.__init__(self, args=self.args)

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.comprehensive = A('comprehensive')
        self.store_all = A('store_all')
        self.just_do_it = A('just_do_it')

        # various methods by which the maximum is found. 'scan' rudimentarily calculates all values
        # in a range and picks the highest value. method_name: (method, method_kwargs)
        self.available_methods = {
            'scan': (self.scan_for_max, {
                'step': A('step') or 1,
                'lower_bound': A('lower_bound') or 5,
                'upper_bound': A('upper_bound') or 35,
            }),
        }

        if self.comprehensive:
            if not self.just_do_it:
                raise ConfigError("Don't use --developer.")
            self.setup_developer_environment()


    def setup_developer_environment(self):
        """ Setup class so that `comprehensive_stats.txt` is output at the end of the run. """
        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        self.method = 'scan'

        self.comprehensive_stats = {
            'kmer': [],
            'min_distance_from_0_or_1': [],
            'mean': [],
            'variance': [],
            'skewness': [],
            'kurtosis': [],
        }

        self.entropy_binning_methods = [
            'array_len_div_4',
            'array_len_div_2',
            'array_len',
            'array_len_times_2',
            'array_len_times_4',
            'array_len_times_8',
            'auto',
            'fd',
            'doane',
            'scott',
            'rice',
            'sturges',
            'sqrt',
        ]

        for binning_method in self.entropy_binning_methods:
            self.comprehensive_stats['entropy_' + binning_method] = []


    def process(self, input_path, fasta_files):
        kmer_search, kmer_search_args = self.determine_kmer_search_method()
        best_kmer = kmer_search(input_path, fasta_files, **kmer_search_args)
        best_report_name = 'kmer_%d_mash_similarity' % best_kmer

        self.run.info_single('Optimal kmer value was %d' % best_kmer, nl_before=1, nl_after=1)

        if not self.store_all:
            self.results = {best_report_name: self.results[best_report_name]}

        if self.comprehensive:
            self.report_comprehensive()

        return self.results


    def report_comprehensive(self):
        self.comprehensive_stats = pd.DataFrame(self.comprehensive_stats)
        self.comprehensive_stats.to_csv(self.comprehensive, sep='\t', index=False)


    def scan_for_max(self, input_path, fasta_files, step=1, lower_bound=5, upper_bound=51):
        entropy_dict = {} # kmer: entropy

        for kmer in range(lower_bound, upper_bound + 1, step):
            self.kmer_size = kmer
            report_name = 'kmer_%d_mash_similarity' % self.kmer_size

            self.Sourmash_process(self, input_path, fasta_files)

            upper_triangular = self.get_upper_triangular(self.results[report_name])

            entropy_dict[kmer] = calculate_entropy(upper_triangular)
            self.run.info('[sourmash] Similarity entropy', entropy_dict[kmer])

            if self.comprehensive:
                self.calculate_comprehensive_stats(upper_triangular)

        return max(entropy_dict, key=entropy_dict.get)


    def get_upper_triangular(self, matrix):
        upper_triangular = pd.DataFrame(matrix).values.astype(float)
        upper_triangular = upper_triangular[np.triu_indices(upper_triangular.shape[0], k=1)]
        return upper_triangular


    def determine_kmer_search_method(self):
        self.run.info('[sourmash] Method to find optimal kmer', self.method)

        try:
            return self.available_methods[self.method]
        except KeyError as e:
            raise ConfigError("IterateKmerSourmash :: .%s is not a valid method for entropy maximization." % self.method)


    def calculate_comprehensive_stats(self, upper_triangular):
        self.comprehensive_stats['kmer'].append(self.kmer_size)
        self.comprehensive_stats['min_distance_from_0_or_1'].append(calculate_average_min_distance_from_0_or_1(upper_triangular))
        self.comprehensive_stats['mean'].append(calculate_mean(upper_triangular))
        self.comprehensive_stats['variance'].append(calculate_variance(upper_triangular))
        self.comprehensive_stats['skewness'].append(calculate_skewness(upper_triangular))
        self.comprehensive_stats['kurtosis'].append(calculate_kurtosis(upper_triangular))

        for bin_method in self.entropy_binning_methods:
            self.comprehensive_stats['entropy_' + bin_method].append(calculate_entropy(upper_triangular, bin_method=bin_method))

        for metric in [metric for metric in self.comprehensive_stats.keys() if metric != 'kmer']:
            self.run.info('[sourmash] ' + metric, self.comprehensive_stats[metric][-1])


def calculate_entropy(array, bin_method='auto'):
    """Calculates Shannon entropy of array.

    This method is very sensitive to how many bins are used.

    Parameters
    ==========
    array : array-like
        Input array
    bin_method : str, optional
        How are the bins selected? Default is numpy's 'auto'

    Returns
    =======
    entropy : float
        Shannon entropy
    """
    if bin_method == 'array_len_div_4':
        bins = np.linspace(0, 1, int((len(array) + 1)/4))
    elif bin_method == 'array_len_div_2':
        bins = np.linspace(0, 1, int((len(array) + 1)/2))
    elif bin_method == 'array_len':
        bins = np.linspace(0, 1, len(array) + 1)
    elif bin_method == 'array_len_times_2':
        bins = np.linspace(0, 1, 2*(len(array) + 1))
    elif bin_method == 'array_len_times_4':
        bins = np.linspace(0, 1, 4*(len(array) + 1))
    elif bin_method == 'array_len_times_8':
        bins = np.linspace(0, 1, 8*(len(array) + 1))
    else:
        bins = bin_method

    return entropy(np.histogram(array, bins=bins, range=(0,1))[0])


def calculate_average_min_distance_from_0_or_1(array):
    """Calculates the element-wise distance to 0 or 1, whichever is closer. Takes the average.

    Parameters
    ==========
    array : array-like
        Input array

    Returns
    =======
    output : type
        The mean of the minimum distance to 0 or 1
    """

    return np.mean(np.min(np.stack((array, 1 - array)), axis=0))


def calculate_mean(array):
    return np.mean(array)


def calculate_variance(array):
    return np.var(array)


def calculate_skewness(array):
    return skew(array)


def calculate_kurtosis(array):
    return kurtosis(array)
