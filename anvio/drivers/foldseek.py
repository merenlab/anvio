#!/usr/bin/env python
# coding: utf-8
""" Foldseek Driver ,"""

import os
import argparse
import pandas as pd
import tempfile

import anvio
import anvio.fastalib as f
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.constants as constants

from collections import defaultdict
from anvio.drivers.mcl import MCL
from anvio.errors import ConfigError
from anvio.filesnpaths import AppendableFile


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
    
    def __init__(self, query_fasta=None, run=run, progress=progress, num_threads=1, weight_dir=None, overwrite_output_destinations=False):
        self.run = run
        self.progress = progress

        utils.is_program_exists('foldseek')

        self.query_fasta = query_fasta
        self.num_threads = num_threads
        self.weight_dir = weight_dir or constants.default_foldseek_weight_path
        self.overwrite_output_destinations = overwrite_output_destinations
        self.tmp_dir = tempfile.gettempdir()

        filesnpaths.is_file_exists(self.weight_dir)

        self.output_file = 'foldseek-search-results'

        if not self.run.log_file_path:
            self.run.log_file_path = filesnpaths.get_temp_file_path()

        self.names_dict = None

    def create_db(self):
        self.run.warning(None, header="FOLDSEEK CREATEDB", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        expected_output_dir = os.path.join(self.output_file, "db")
        expected_output_file = os.path.join(expected_output_dir, "search_db")

        filesnpaths.gen_output_directory(expected_output_dir, delete_if_exists=False)

        cmd_line = ['foldseek',
                    'createdb',
                    self.query_fasta,
                    expected_output_file,
                    '--prostt5-model', self.weight_dir,
                    '--threads', self.num_threads
                    ]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()
        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Foldseek search DB', expected_output_file)

    def search(self, query_db, target_db):
        self.run.warning(None, header="FOLDSEEK EASY SEARCH", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('Running search using Foldseek ...')

        query_db = os.path.join(query_db, 'search_db')
        target_db = os.path.join(target_db, 'search_db')

        result_file_dir = os.path.join(self.output_file, 'result')

        cmd_line = [
            'foldseek',
            'easy-search',
            query_db,
            target_db,
            result_file_dir,
            self.tmp_dir,
            '--threads', self.num_threads
        ]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Foldseek search Result', result_file_dir)

    def process(self, query_db, target_db):

        self.create_db()
        self.search(query_db, target_db)

    def get_foldseek_results(self):
        """ Return result.m8 file """
        force_makedb, force_search = False, False

        result_dir = os.path.join(self.output_file, 'result')

        return result_dir
