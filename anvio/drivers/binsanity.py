# coding: utf-8
"""Interface to BinSanity."""
import os
import glob
import shutil

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class BinSanity:
    arguments = {
        'preference': (
                ['--preference'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Specify a preference (default is -3)\
                           Note: decreasing the preference leads to more lumping,\
                           increasing will lead to more splitting. If your range\
                           of coverages are low you will want to decrease the preference,\
                           if you have 10 or less replicates increasing the preference could\
                           benefit you."}
                    ),
        'maxiter': (
                ['--maxiter'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Specify a max number of iterations [default is 2000]"}
                    ),
        'conviter': (
                ['--conviter'],
                {'metavar': "INT",
                 'required': False,
                 'help': " Specify the convergence iteration number (default is 200)\
                           e.g Number of iterations with no change in the number\
                           of estimated clusters that stops the convergence."}
                    ),
        'damp': (
                ['--damp'],
                {'metavar': "FLOAT",
                 'required': False,
                 'help': "Specify a damping factor between 0.5 and 1, default is 0.9"}
                    ),
        'contigsize': (
                ['--contigsize'],
                {'metavar': "INT",
                 'required': False,
                 'help': "TNF probability cutoff for building TNF graph. Use it to skip the\
                                    preparation step. (0: auto)."}
                    ),
    }
    citation = "Graham ED, Heidelberg JF, Tully BJ. (2017) BinSanity: unsupervised \
                clustering of environmental microbial assemblies using coverage and \
                affinity propagation. PeerJ 5:e3035 https://doi.org/10.7717/peerj.3035"

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = 'Binsanity'

        utils.is_program_exists(self.program_name)


    def cluster(self, input_files, args, threads=1, splits_mode=False):
        self.temp_path = filesnpaths.get_temp_directory_path()
        self.run.info_single("If you publish results from this workflow, \
                               please do not forget to cite \n%s" % BinSanity.citation,
                               nl_before=1, nl_after=1, mc='green')

        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)

        log_path = os.path.join(self.temp_path, 'logs.txt')

        translation = {
            'preference': 'p',
            'maxiter': 'm',
            'conviter': 'v',
            'damp': 'd',
            'contigsize': 'x'
        }

        cmd_line = [self.program_name,
            '-c', input_files.coverage, 
            '-f', os.path.dirname(input_files.fasta),
            '-l', os.path.basename(input_files.fasta),
            '-o', self.temp_path,
            *utils.serialize_args(args, single_dash=True, translate=translation)]

        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        clusters = {}
        bin_count = 0
        for bin_file in glob.glob(os.path.join(self.temp_path, '*.fna')):
            bin_count += 1
            with open(bin_file, 'r') as f:
                pretty_bin_name = os.path.basename(bin_file)
                pretty_bin_name = pretty_bin_name.replace('sequence_', '')
                pretty_bin_name = pretty_bin_name.replace('.fna', '')
                pretty_bin_name = pretty_bin_name.replace('-', '_')

                clusters[pretty_bin_name] = [line.strip().replace('>', '') for line in f if line.startswith('>')]

        self.run.info('Bins formed', bin_count)

        if not anvio.DEBUG:
            shutil.rmtree(self.temp_path)

        return clusters
        