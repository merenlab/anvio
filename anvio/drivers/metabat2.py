# coding: utf-8
"""Interface to MetaBAT2."""
import os
import glob
import shutil

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


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


class MetaBAT2:
    arguments = {
        'minContig': (
                ['-m', '--minContig'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum size of a contig for binning (should be >=1500)"}
                    ),
        'maxP': (
                ['--maxP'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Percentage of 'good' contigs considered for binning decided by connection\
                                    among contigs. The greater, the more sensitive."}
                    ),
        'minS': (
                ['--minS'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum score of a edge for binning (should be between 1 and 99). The\
                                    greater, the more specific."}
                    ),
        'maxEdges': (
                ['--maxEdges'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Maximum number of edges per node. The greater, the more sensitive."}
                    ),
        'pTNF': (
                ['--pTNF'],
                {'metavar': "INT",
                 'required': False,
                 'help': "TNF probability cutoff for building TNF graph. Use it to skip the\
                                    preparation step. (0: auto)."}
                    ),
        'noAdd': (
                ['--noAdd'],
                {'action': 'store_true',
                 'default': False,
                 'help': "Turning off additional binning for lost or small contigs."}
                    ),
        'minCV': (
                ['--minCV'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum mean coverage of a contig in each library for binning."}
                    ),
        'minCVSum': (
                ['--minCVSum'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum size of a bin as the output."}
                    ),
        'minClsSize': (
                ['--minClsSize'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Minimum size of a bin as the output."}
                    ),
        'seed': (
                ['--seed'],
                {'metavar': "INT",
                 'required': False,
                 'help': "For exact reproducibility. (0: use random seed)."}
                    ),
    }

    citation = "Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. 2019. \
                MetaBAT 2: an adaptive binning algorithm for robust and efficient \
                genome reconstruction from metagenome assemblies. \
                PeerJ 7:e7359 https://doi.org/10.7717/peerj.7359"

    cluster_type = 'contig'


    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = 'metabat2'

        utils.is_program_exists(self.program_name)


    def cluster(self, input_files, args, work_dir, threads=1):
        J = lambda p: os.path.join(work_dir, p)

        bin_prefix = J('METABAT_')
        log_path = J('logs.txt')

        cmd_line = [self.program_name,
            '-i', input_files.contigs_fasta,
            '-a', input_files.contig_coverages,
            '-o', bin_prefix,
            '--cvExt',
            '-l',
            *utils.serialize_args(args)]


        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        output_file_paths = glob.glob(J(bin_prefix + '*'))
        if not len(output_file_paths):
            raise ConfigError("Some critical output files are missing. Please take a look at the "
                              "log file: %s" % (log_path))

        clusters = {}
        bin_count = 0
        for bin_file in output_file_paths:
            bin_count += 1
            with open(bin_file, 'r') as f:
                pretty_bin_name = os.path.basename(bin_file).replace('.', '_')
                clusters[pretty_bin_name] = list(map(str.strip, f.readlines()))

        return clusters
