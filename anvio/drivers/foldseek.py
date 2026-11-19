#!/usr/bin/env python
""" Foldseek Driver """

import os
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

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

# Output column contract for foldseek easy-search. Default 12 columns plus
# qtmscore and ttmscore. Downstream consumers (e.g. anvio.panops) unpack the
# result file using this exact column order.
FOLDSEEK_OUTPUT_COLUMNS = (
    "query,target,fident,alnlen,mismatch,gapopen,"
    "qstart,qend,tstart,tend,evalue,bits,qtmscore,ttmscore"
)


class Foldseek():

    def __init__(self, structure_dir=None, run=run, progress=progress,
                 num_threads=1, overwrite_output_destinations=False):
        self.run = run
        self.progress = progress

        utils.is_program_exists('foldseek')

        if structure_dir is None:
            raise ConfigError("The Foldseek driver requires a `structure_dir`: a directory of structure files.")

        if not os.path.isdir(structure_dir):
            raise ConfigError(f"`structure_dir` must be an existing directory; got '{structure_dir}'.")

        self.structure_dir = structure_dir
        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations
        self.tmp_dir = tempfile.gettempdir()

        self.output_file = None
        self.result_file_path = None

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

        # Passing a directory rather than splatting every file keeps us safely under
        # ARG_MAX no matter how many structures the user has. Foldseek's createdb walks
        # the directory and picks up every structure file in it.
        cmd_line = ['foldseek',
                    'createdb',
                    self.structure_dir,
                    expected_output_file,
                    '--threads', self.num_threads
                    ]

        try:
            utils.run_command(cmd_line, self.run.log_file_path)
        except ConfigError:
            self.progress.end()
            raise ConfigError(f"Foldseek createdb failed. Please check the log file at '{self.run.log_file_path}' "
                              f"for more details. This could be due to malformed structure files in "
                              f"'{self.structure_dir}'.")

        self.progress.end()
        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Foldseek search DB', expected_output_file)

    def search(self, query_db, target_db):
        self.run.warning(None, header="FOLDSEEK EASY SEARCH", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('Running search using Foldseek ...')

        query_db = os.path.join(query_db, 'db', 'search_db')
        target_db = os.path.join(target_db, 'db', 'search_db')

        self.result_file_path = os.path.join(self.output_file, 'result')

        cmd_line = [
            'foldseek',
            'easy-search',
            query_db,
            target_db,
            self.result_file_path,
            self.tmp_dir,
            '--format-output', FOLDSEEK_OUTPUT_COLUMNS,
            '--threads', self.num_threads
        ]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Foldseek search result', self.result_file_path)


    def process(self, output_file):

        self.output_file = output_file

        self.create_db()
        self.search(output_file, output_file)


    def get_foldseek_results(self):
        """Returns the path to the foldseek search results file."""

        if not self.result_file_path or not os.path.exists(self.result_file_path):
            raise ConfigError(f"Foldseek search results file was not found at '{self.result_file_path}'. "
                              f"This could mean that the foldseek search did not complete successfully. "
                              f"Please check the log file at '{self.run.log_file_path}' for more details.")

        if os.path.getsize(self.result_file_path) == 0:
            raise ConfigError(f"The foldseek search results file at '{self.result_file_path}' is empty. "
                              f"This could mean that foldseek did not find any structural similarities "
                              f"between your gene cluster representatives. Please check the log file at "
                              f"'{self.run.log_file_path}' for more details.")

        return self.result_file_path
