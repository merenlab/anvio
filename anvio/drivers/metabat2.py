# coding: utf-8
"""Interface to MetaBAT2."""
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

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = 'metabat2'

        utils.is_program_exists(self.program_name)


    def cluster(self, input_files, args, threads=1, splits_mode=False):
        self.temp_path = filesnpaths.get_temp_directory_path()
        self.run.info_single("If you publish results from this workflow, \
                               please do not forget to cite \n%s" % MetaBAT2.citation,
                               nl_before=1, nl_after=1, mc='green')

        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)

        bin_prefix = os.path.join(self.temp_path, 'Bin')
        log_path = os.path.join(self.temp_path, 'logs.txt')

        cmd_line = [self.program_name,
            '-i', input_files.fasta,
            '-a', input_files.coverage,
            '-o', bin_prefix,
            '--cvExt',
            '-l',
            *utils.serialize_args(args)]


        self.progress.new(self.program_name)
        self.progress.update('Running using %d threads...' % threads)
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        clusters = {}
        bin_count = 0
        for bin_file in glob.glob(bin_prefix + '*'):
            bin_count += 1
            with open(bin_file, 'r') as f:
                pretty_bin_name = os.path.basename(bin_file).replace('.', '_')
                clusters[pretty_bin_name] = list(map(str.strip, f.readlines()))

        self.run.info('Bins formed', bin_count)

        if not anvio.DEBUG:
            shutil.rmtree(self.temp_path)

        return clusters
