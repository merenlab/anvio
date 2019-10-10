# coding: utf-8
"""Interface to MaxBin2."""
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


class MaxBin2:
    arguments = {
        'min_contig_length': (
                ['--min-contig-length'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum contig length. Default: 1000."}
                    ),
        'max_iteration': (
                ['--max-iteration'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Maximum Expectation-Maximization algorithm iteration number. Default 50."}
                    ),
        'prob_threshold': (
                ['--prob-threshold'],
                {'metavar': "FLOAT",
                 'required': False,
                 'help': "Probability threshold for EM final classification. Default 0.9."}
                    ),
        'markerset': (
                ['--merkerset'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum contig length. Default: 1000"}
                    ),
    }
    citation = "Yu-Wei Wu, Blake A. Simmons, Steven W. Singer, MaxBin 2.0: \
                an automated binning algorithm to recover genomes from multiple \
                metagenomic datasets, Bioinformatics, Volume 32, Issue 4, 15 February 2016,\
                Pages 605–607, https://doi.org/10.1093/bioinformatics/btv638"


    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = 'run_MaxBin.pl'

        utils.is_program_exists(self.program_name)


    def cluster(self, input_files, args, threads=1, splits_mode=False):
        self.run.info_single("If you publish results from this workflow, \
                               please do not forget to cite \n%s" % MaxBin2.citation,
                               nl_before=1, nl_after=1, mc='green')

        self.temp_path = filesnpaths.get_temp_directory_path()

        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)

        bin_prefix = os.path.join(self.temp_path, 'Bin')
        log_path = os.path.join(self.temp_path, 'logs.txt')

        cmd_line = [self.program_name,
            '-contig', input_files.fasta,
            '-abund', input_files.coverage,
            '-out', bin_prefix,
            '-thread', str(threads),
            *utils.serialize_args(args, single_dash=True, use_underscore=True)]


        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        clusters = {}
        bin_count = 0

        for bin_file in glob.glob(bin_prefix + '*.fasta'):
            bin_count += 1
            with open(bin_file, 'r') as f:
                bin_name = os.path.basename(bin_file).replace('.fasta', '')
                bin_name = bin_name.replace('.', '_')

                clusters[bin_name] = []

                for line in f.readlines():
                    if line.startswith('>'):
                        clusters[bin_name].append(line[1:].strip())

        self.run.info('Bins formed', bin_count)

        if not anvio.DEBUG:
            shutil.rmtree(self.temp_path)

        return clusters
