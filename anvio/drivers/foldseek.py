#!/usr/bin/env python
# coding: utf-8
""" Foldseek Driver """

import os
import shutil
import argparse
import tempfile

import anvio
import anvio.fastalib as f
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.constants as constants

from collections import defaultdict
from anvio.drivers.mcl import MCL
from anvio.errors import ConfigError, FilesNPathsError
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
        self.weight_dir = weight_dir or constants.default_prostt5_weight_path
        self.overwrite_output_destinations = overwrite_output_destinations
        self.tmp_dir = tempfile.gettempdir()

        try: 
            filesnpaths.is_file_exists(self.weight_dir)
        except FilesNPathsError:
            run.warning("Anvi'o requires to have ProstT5 to run --pan-mode structure."
                        " You can easily download that with the command down below.",
                        header="⚠️  YOUR ATTENTION PLEASE ⚠️", overwrite_verbose=True, lc='yellow')
            run.info_single("anvi-setup-prostt5", level=0, overwrite_verbose=True)
            raise ConfigError("It seems like you forgot to download the ProstT5 model."
                            " Please run 'anvi-setup-prostt5' and try again.")

        if not self.run.log_file_path:
            self.run.log_file_path = filesnpaths.get_temp_file_path()

        self.names_dict = None

    def create_db(self, output_file):
        self.run.warning(None, header="FOLDSEEK CREATEDB", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        expected_output_dir = os.path.join(output_file, "db")
        expected_output_file = os.path.join(expected_output_dir, "search_db")

        filesnpaths.gen_output_directory(expected_output_dir, delete_if_exists=False)

        cmd_line = ['foldseek',
                    'createdb',
                    self.query_fasta,
                    expected_output_file,
                    '--prostt5-model', self.weight_dir,
                    '--threads', self.num_threads
                    ]

        try:
            utils.run_command(cmd_line, self.run.log_file_path)
        except FilesNPathsError:
            run.warning("Opss! CREATEDB not working. Probably you are giving wrong file path :/")

        self.progress.end()
        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        #self.run.info('Foldseek search DB', expected_output_file)

    def search(self, query_db, target_db):
        self.run.warning(None, header="FOLDSEEK EASY SEARCH", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('Running search using Foldseek ...')

        query_db = os.path.join(query_db, 'db', 'search_db')
        target_db = os.path.join(target_db, 'db', 'search_db')

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

    def process(self, output_file):

        self.create_db(output_file)
        self.search(output_file, output_file)

    def get_foldseek_results(self):
        """ Return result.m8 file """
        force_makedb, force_search = False, False

        result_dir = os.path.join(self.output_file, 'result')

        return result_dir


class Prostt5SetupWeight:
    """A class to download and setup the weights of PROSTT5"""
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        utils.is_program_exists('foldseek')

        self.weight_dir = args.prostt5_data_dir

        if self.weight_dir and args.reset:
            raise ConfigError("You are attempting to run ProstT5 setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --weight-dir. If you really want to get rid of this "
                              "directory and regenerate it with InteracDome data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.weight_dir, self.weight_dir))

        if not self.weight_dir:
            self.weight_dir = constants.default_prostt5_weight_path

        self.run.warning('', header='Setting up ProstT5 Weights', lc='yellow')
        self.run.info('Data directory', self.weight_dir)
        self.run.info('Reset contents', args.reset)

        filesnpaths.gen_output_directory(self.weight_dir, delete_if_exists=args.reset)

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.weight_dir)))

        if not args.reset and not anvio.DEBUG:
            self.is_weight_exists()

        if not self.run.log_file_path:
            self.run.log_file_path = os.path.join(self.weight_dir, 'foldseek-setup-log-file.txt')


    def is_weight_exists(self):
        """Raise ConfigError if weight exists"""

        if os.path.exists(self.weight_dir) and os.listdir(self.weight_dir):
            raise ConfigError("It seems you already have the ProstT5 Weights downloaded in '%s', please "
                              "use --reset flag if you want to re-download it." % self.weight_dir)

    def setup(self):
        """ Sets up the ProstT5 Weights directory for Foldseek """

        self.run.warning('', header='Downloading Weight Model', lc='yellow')
        self.download_foldseek_weight()

    def download_foldseek_weight(self):
        """Download the weights of ProstT5 models and clean up temporary files"""

        self.progress.new('FOLDSEEK')
        self.progress.update('Downloading ...')

        cmd_line = [
            'foldseek',
            'databases',
            'ProstT5',
            self.weight_dir,
            os.path.join(self.weight_dir, 'tmp')
        ]

        result = utils.run_command(cmd_line, self.run.log_file_path)

        # Successful condition
        if result == 0:
            self.progress.end()

            self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
            self.run.info('Log file', self.run.log_file_path)

            # Remove tmp folder after download model
            tmp_dir = os.path.join(self.weight_dir, 'tmp')
            if os.path.exists(tmp_dir):
                try:
                    shutil.rmtree(tmp_dir)
                    self.run.info('Temporary folder', f"'{tmp_dir}' was successfully deleted.")
                except Exception as e:
                    self.run.warning('Temporary folder cleanup failed', str(e))
        else:
            self.run.warning('Download failed', 'Please check the log file for more details.')