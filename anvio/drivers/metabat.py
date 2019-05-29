# coding: utf-8
"""Interface to MetaBAT."""
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


class MetaBAT:
    arguments = {
        "infile"
        'seed': (
                ['--seed'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Seed for random numbers"}
                    ),
        'threads': (
                ['-T', '--threads'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Number of threads"}
                    ),
    }
    citation = "Citation here"

    def __init__(self, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = 'metabat2'

        utils.is_program_exists(self.program_name)
        self.temp_path = filesnpaths.get_temp_directory_path()


    def cluster(self, input_files, args, threads=1):
        if anvio.DEBUG:
            self.run.info('Working directory', self.temp_path)

        bin_prefix = os.path.join(self.temp_path, 'Bin')
        log_path = os.path.join(self.temp_path, 'logs.txt')

        cmd_line = [self.program_name,
            '-i', input_files.fasta, 
            '-a', input_files.coverage, 
            '-o', bin_prefix,
            '--cvExt',
            '-l']

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
        