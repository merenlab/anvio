# -*- coding: utf-8
# pylint: disable=line-too-long
"""Anvi'o - CONCOCT interface for unsupervised clustering"""
import os

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


class CONCOCT:
    arguments = {
        'clusters': (
                ['--clusters'],
                {'metavar': "INT",
                 'required': False,
                 'help': "specify maximal number of clusters for VGMM, default\
                        400"}
                    ),
        'kmer_length': (
                ['--kmer-length'],
                {'metavar': "INT",
                 'required': False,
                 'help': "pecify kmer length, default 4"}
                    ),
        'length_threshold': (
                ['--length-threshold'],
                {'metavar': "INT",
                 'required': False,
                 'help': "specify the sequence length threshold, contigs shorter\
                        than this value will not be included. Defaults to\
                        1000"}
                    ),
        'read_length': (
                ['--read-length'],
                {'metavar': "INT",
                 'help': "specify read length for coverage, default 100"}
                    ),
        'no_cov_normalization': (
                ['--no-cov-normalization'],
                {'required': False,
                 'action': 'store_true',
                 'help': "By default the coverage is normalized with regards to\
                        samples, then normalized with regards of contigs and\
                        finally log transformed. By setting this flag you skip\
                        the normalization and only do log transorm of the\
                        coverage."}
                    ),
        'total_percentage_pca': (
                ['--total-percentage-pca'],
                {'metavar': "INT",
                 'required': False,
                 'help': "The percentage of variance explained by the principal\
                        components for the combined data."}
                    ),
        'epsilon': (
                ['--epsilon'],
                {'metavar': "FLOAT",
                 'type': float,
                 'required': False,
                 'help': "Specify the epsilon for VBGMM. Default value is 1.0e-6"}
                    ),
        'iterations': (
                ['--iterations'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Specify maximum number of iterations for the VBGMM.\
                        Default value is 500"}
                    ),
        'seed': (
                ['--seed'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Specify an integer to use as seed for clustering. 0\
                        gives a random seed, 1 is the default seed and any\
                        other positive integer can be used. Other values give\
                        ArgumentTypeError."}
                    )
    }
    citation = "Johannes Alneberg, Brynjar Smári Bjarnason, Ino de Bruijn, \
                Melanie Schirmer, Joshua Quick, Umer Z Ijaz, Leo Lahti,\
                Nicholas J Loman, Anders F Andersson & Christopher Quince. \
                2014. Binning metagenomic contigs by coverage and composition. \
                Nature Methods, doi: 10.1038/nmeth.3103"

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.program_name = 'concoct'
        utils.is_program_exists(self.program_name)


    def cluster(self, input_files, args, threads=1, splits_mode=False):
        self.run.info_single("If you publish results from this workflow, \
                               please do not forget to cite \n%s" % CONCOCT.citation,
                               nl_before=1, nl_after=1, mc='green')

        self.temp_path = filesnpaths.get_temp_directory_path()
        log_path = os.path.join(self.temp_path, 'logs.txt')

        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)

        cmd_line = [self.program_name,
            '--coverage_file', input_files.coverage,
            '--composition_file', input_files.fasta,
            '--basename', self.temp_path,
            '--threads', threads,
             *utils.serialize_args(args, use_underscore=True)]

        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        clusters = {}
        threshold = args.length_threshold or '1000'

        with open(os.path.join(self.temp_path, 'clustering_gt%s.csv' % threshold), 'r') as f:
            lines = f.readlines()[1:]

            for entry in lines:
                contig, bin_name = map(str.strip, entry.split(','))

                pretty_bin_name = 'Bin_' + bin_name

                if pretty_bin_name not in clusters:
                    clusters[pretty_bin_name] = []

                clusters[pretty_bin_name].append(contig)

        return clusters
