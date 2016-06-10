# coding: utf-8
"""Interface to Diamond."""

import os
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2016, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class Diamond:
    def __init__(self, query_fasta, run=run, progress=progress, num_threads=1, overwrite_output_destinations=False):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations

        utils.is_program_exists('diamond')

        self.tmp_dir = tempfile.gettempdir()

        self.query_fasta = query_fasta
        self.target_db_path = 'diamond-target'
        self.search_output_path = 'diamond-search-resuults'
        self.tabular_output_path = 'diamond-search-results.txt'

        if not self.run.log_file_path:
            self.run.log_file_path = 'diamond-log-file.txt'

        self.sensitive = False

        # if names_dict is None, all fine. if not, the query_fasta is assumed to be uniqued, and names_dict is
        # the dictionary that connects the ids in the fasta file, to ids that were identical to it.
        self.names_dict = None


    def get_blastall_results(self):
        force_makedb, force_blastp, force_view = False, False, False

        if self.overwrite_output_destinations:
            force_makedb = True

        if os.path.exists(self.target_db_path + '.dmnd') and not force_makedb:
            self.run.warning("Notice: A diamond database is found in the output directory, and will be used!")
        else:
            self.makedb()
            force_blastp, forrce_view = True, True

        if os.path.exists(self.search_output_path + '.daa') and not force_blastp:
            self.run.warning("Notice: A DIAMOND search result is found in the output directory: skipping BLASTP!")
        else:
            self.blastp()
            force_view = True

        if os.path.exists(self.tabular_output_path) and not force_view:
            self.run.warning("Notice: A DIAMOND tabular output is found in the output directory. Anvi'o will not generate another one!")
        else:
            self.view()

        return self.tabular_output_path


    def check_output(self, expected_output, process='diamond'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError, "Pfft. Something probably went wrong with Diamond's '%s' since one of the expected output files are missing.\
                                Please check the log file here: '%s." % (process, self.run.log_file_path)


    def makedb(self):
        self.progress.new('DIAMOND')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('diamond makedb --in %s -d %s -p %d >> "%s" 2>&1' % (self.query_fasta,
                                                                         self.target_db_path,
                                                                         self.num_threads,
                                                                         self.run.log_file_path))

        self.run.info('diamond makedb cmd', cmd_line, quiet=True)

        utils.run_command(cmd_line)

        self.progress.end()

        expected_output = self.target_db_path + '.dmnd'
        self.check_output(expected_output, 'makedb')

        self.run.info('Diamond temp search db', expected_output)


    def blastp(self):
        self.run.info('DIAMOND is set to be', 'Sensitive' if self.sensitive else 'Fast')

        self.progress.new('DIAMOND')
        self.progress.update('running blastp (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('diamond blastp -q %s -d %s -a %s -t %s -p %d %s -k 1000000 >> "%s" 2>&1' % (self.query_fasta,
                                                                                   self.target_db_path,
                                                                                   self.search_output_path,
                                                                                   self.tmp_dir,
                                                                                   self.num_threads,
                                                                                   '--sensitive' if self.sensitive else '',
                                                                                   self.run.log_file_path))

        self.run.info('diamond blastp cmd', cmd_line, quiet=True)

        utils.run_command(cmd_line)

        self.progress.end()

        expected_output = self.search_output_path + '.daa'
        self.check_output(expected_output, 'blastp')

        self.run.info('Diamond blastp results', expected_output)


    def view(self):
        self.progress.new('DIAMOND')
        self.progress.update('generating tabular output (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('diamond view -a %s -o %s -p %d -k 1000000 >> "%s" 2>&1' % (self.search_output_path + '.daa',
                                                                     self.tabular_output_path,
                                                                     self.num_threads,
                                                                     self.run.log_file_path))

        self.run.info('diamond view cmd', cmd_line, quiet=True)

        utils.run_command(cmd_line)

        self.check_output(self.tabular_output_path, 'view')

        if self.names_dict:
            self.run.info('self.names_dict is found', 'Un-uniqueing the tabular output', quiet=True)
            self.progress.update('Un-uniqueing the tabular output ...')
            # if we are here, this means the self.tabular_output_path contains information only about unique
            # entries. We will expand it here so downstream analyses do not need to pay attention to this
            # detail.
            utils.ununique_BLAST_tabular_output(self.tabular_output_path, self.names_dict)

        self.progress.end()

        self.run.info('Diamond %stabular output file' % ('un-uniqued' if len(self.names_dict) else ''), self.tabular_output_path)
