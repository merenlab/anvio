# coding: utf-8
"""Interface to NCBI's BLAST."""

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


class BLAST:
    def __init__(self, query_fasta, run=run, progress=progress, num_threads=1, overwrite_output_destinations=False):
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.evalue = 1e-05
        self.overwrite_output_destinations = overwrite_output_destinations

        utils.is_program_exists('makeblastdb')
        utils.is_program_exists('blastp')

        self.tmp_dir = tempfile.gettempdir()

        self.query_fasta = query_fasta
        self.target_db_path = 'blast-target'
        self.search_output_path = 'blast-search-results.txt'
        self.max_target_seqs = None


        if not self.run.log_file_path:
            self.run.log_file_path = 'blast-log-file.txt'


        # if names_dict is None, all fine. if not, the query_fasta is assumed to be uniqued, and names_dict is
        # the dictionary that connects the ids in the fasta file, to ids that were identical to it.
        self.names_dict = None


    def get_blastall_results(self):
        force_makedb, force_blastp = False, False

        if self.overwrite_output_destinations:
            force_makedb = True

        if os.path.exists(self.target_db_path + '.phr') and os.path.exists(self.target_db_path + '.pin')\
                                                and os.path.exists(self.target_db_path + '.psq') and not force_makedb:
            self.run.warning("Notice: A BLAST database is found in the output directory, and will be used!")
        else:
            self.makedb()
            force_blastp = True

        if os.path.exists(self.search_output_path) and not force_blastp:
            self.run.warning("Notice: A BLAST search result is found in the output directory: skipping BLASTP!")
        else:
            self.blastp()

        return self.search_output_path


    def check_output(self, expected_output, process='blast'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError, "Pfft. Something probably went wrong with '%s' process since one of the expected output files are missing.\
                                Please check the log file here: '%s." % (process, self.run.log_file_path)


    def makedb(self):
        self.progress.new('BLAST')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('makeblastdb -in %s -dbtype prot -out %s' % (self.query_fasta,
                                                                 self.target_db_path))

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        expected_output = self.target_db_path + '.phr'
        self.check_output(expected_output, 'makeblastdb')

        self.run.info('blast makeblast cmd', cmd_line, quiet=True)
        self.run.info('BLAST search db', self.target_db_path)


    def blastp(self):
        self.progress.new('BLASTP')
        self.progress.update('running blastp (using %d thread(s)) ...' % self.num_threads)
        cmd_line = ('blastp -query %s -db %s -evalue %f -outfmt 6 -out %s -num_threads %d' % (self.query_fasta,
                                                                                              self.target_db_path,
                                                                                              self.evalue,
                                                                                              self.search_output_path,
                                                                                              self.num_threads))

        if self.max_target_seqs:
            cmd_line += ' -max_target_seqs %d' % self.max_target_seqs

        self.run.info('blast blastp cmd', cmd_line, quiet=True)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.check_output(self.search_output_path, 'blastp')

        if self.names_dict:
            self.ununique_search_results()

        self.run.info('BLASTP results', self.search_output_path)


    def ununique_search_results(self):
        self.run.info('self.names_dict is found', 'Un-uniqueing the tabular output.', quiet=True)

        self.progress.new('BLAST')
        self.progress.update('Un-uniqueing the tabular output ...')

        utils.ununique_BLAST_tabular_output(self.search_output_path, self.names_dict)

        self.progress.end()
