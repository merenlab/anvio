#!/usr/bin/env python
# coding: utf-8
""" Foldseek to Pangenome"""

import os
import shutil
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
from anvio.drivers.foldseek import Foldseek
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

class Prostt5SetupWeight:
    """A class to download and setup the weights of PROSTT5"""
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        utils.is_program_exists('foldseek')

        self.weight_dir = args.prostt5_weight_dir

        if self.weight_dir and args.reset:
            raise ConfigError("You are attempting to run PROSTT5 setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --weight-dir. If you really want to get rid of this "
                              "directory and regenerate it with InteracDome data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.weight_dir, self.weight_dir))

        if not self.weight_dir:
            self.weight_dir = constants.default_prostt5_weight_path

        self.run.warning('', header='Setting up PROSTT5 Weights', lc='yellow')
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
            raise ConfigError("It seems you already have the PROSTT5 Weights downloaded in '%s', please "
                              "use --reset flag if you want to re-download it." % self.weight_dir)

    def setup(self):
        """ Sets up the PROSTT5 Weights directory for Foldseek """

        self.run.warning('', header='Downloading Weight Model', lc='yellow')
        self.download_foldseek_weight()

    def download_foldseek_weight(self):
        """Download the weights of PROSTT5 models and clean up temporary files"""

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