#!/usr/bin/env python
# -*- coding: utf-8
"""This file contains CAZyme related classes."""

import os

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew Schechter"
__email__ = "mschechter@uchicago.edu"

run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class CAZymeSetup(object):
    def __init__(self, args, run=run, progress=progress):
        """Setup a CAZyme database for anvi'o

        http://www.cazy.org/

        Parameters
        ==========
        args : argparse.Namespace
            See `bin/anvi-setup-cazymes` for available arguments
            - cazyme_data_dir : str, optional
                The directory where the CAZyme data should be stored. If not provided, the data will be stored in the anvi'o data directory.
            - reset : bool, optional
                If True, the data directory will be deleted and recreated. Defaults to False.
        run : terminal.Run, optional
            An object for printing messages to the console.
        progress : terminal.Progress, optional
            An object for printing progress bars to the console.
        """

        self.args = args
        self.run = run
        self.progress = progress
        self.cazyme_data_dir = args.cazyme_data_dir

        filesnpaths.is_program_exists('hmmpress')
        
        if self.cazyme_data_dir and args.reset:
            raise ConfigError("You are attempting to run CAZyme setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --cazyme-data-dir. If you really want to get rid of this "
                              "directory and regenerate it with CAZyme data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.cazyme_data_dir, self.cazyme_data_dir))

        if not self.cazyme_data_dir:
            self.cazyme_data_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/CAZyme')

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.cazyme_data_dir)))

        self.resolve_database_url()

        if not args.reset and not anvio.DEBUG:
            self.is_database_exists()

        if args.reset:
            filesnpaths.gen_output_directory(self.cazyme_data_dir, delete_if_exists=True, dont_warn=True)
        else:
            filesnpaths.gen_output_directory(self.cazyme_data_dir)

    def resolve_database_url(self):
        """Create path to CAZyme ftp

        Added self values
        ================= 
        - self.page_index : string
            version of CAZyme database

        """
        if self.args.cazyme_version:
            self.page_index = self.args.cazyme_version 
            self.run.info('Attempting to use version', self.args.cazyme_version)
        else:
            self.page_index = 'V11'
            self.run.info_single('No CAZyme version specified. Using current release.')

        self.database_url = os.path.join("https://bcb.unl.edu/dbCAN2/download/Databases", f"{self.page_index}", f"dbCAN-HMMdb-{self.page_index}.txt") 

    def is_database_exists(self):
        """Determine if CAZyme database has already been downloaded"""
        if os.path.exists(os.path.join(self.cazyme_data_dir, f"dbCAN-HMMdb-{self.page_index}.txt")):
            raise ConfigError(f"It seems you already have CAZyme database installed in {self.cazyme_data_dir}, please use --reset flag if you want to re-download it.")

    def download(self, hmmpress_files=True):
        """Download CAZyme database and compress with hmmpress"""
        self.run.info("Database URL", self.database_url)

        utils.download_file(self.database_url, os.path.join(self.cazyme_data_dir, os.path.basename(self.database_url)) , progress=self.progress, run=self.run)

        if hmmpress_files:
            self.hmmpress_files()

    def hmmpress_files(self):
        """Runs hmmpress on CAZyme HMM profiles."""

        file_path = os.path.join(self.cazyme_data_dir, os.path.basename(self.database_url))
        cmd_line = ['hmmpress', file_path]
        log_file_path = os.path.join(self.cazyme_data_dir, '00_hmmpress_log.txt')
        ret_val = utils.run_command(cmd_line, log_file_path)

        if ret_val:
            raise ConfigError("Hmm. There was an error while running `hmmpress` on the Pfam HMM profiles. "
                                "Check out the log file ('%s') to see what went wrong." % (log_file_path))
        else:
            # getting rid of the log file because hmmpress was successful
            os.remove(log_file_path)