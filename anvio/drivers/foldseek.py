# coding: utf-8
""" Foldseek to Pangenome"""

import os
import pandas as pd
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Metehan Sever"
__email__ = "metehaansever@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class Foldseek():
    
    def __init__(self, query_fasta=None, run=run, progress=progress, num_threads=1, overwrite_output_destinations=False, output_file_path=None, weight=None):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations
        self.output_file_path = output_file_path
        self.weight = weight

        utils.is_program_exists('foldseek')

        self.query_fasta = query_fasta

        if not self.run.log_file_path:
            self.run.log_file_path = 'foldseek-log-file.txt'

        self.additional_params_for_blastp = ""

    def create_db(self):
        self.run.warning(None, header="FOLDSEEK CREATEDB", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['foldseek',
                    'createdb',
                    self.query_fasta,
                    self.output_file_path + 'db',
                    '--prostt5-model', self.weight, # Where should the weight of Prostt5 be placed?
                    '--threads', str(self.num_threads)
                    ]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        expected_output = self.output_file_path + 'db.ss'

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Foldseek search DB', expected_output)

    def search(self, query_db, target_db, output_dir):
        self.run.warning(None, header="FOLDSEEK SEARCH", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('Running search using Foldseek ...')

        cmd_line = [
            'foldseek',
            'search',
            query_db,
            target_db,
            output_dir + '/results',
            '--threads', self.num_threads
        ]

        utils.run_command(cmd_line, self.run.log_file_path)
        self.progress.end()