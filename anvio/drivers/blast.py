# coding: utf-8
"""Interface to NCBI's BLAST."""

import os
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class BLAST:
    def __init__(self, query_fasta, target_fasta=None, search_program='blastp', run=run, progress=progress, num_threads=1, overwrite_output_destinations=False):
        """BLAST driver.

           We generate target database from the `target_fasta`. If `target_fasta` is None,
           `query_fasta` is treated as `target_fasta`. If you don't have a FASTA file, but
           all you have are X.phr, X.pin, and x.psq files, you can set `target_fasta` to
           '/path/to/X' and it will still be OK. Calling the target FASTA creates some
           confusion, but we hope if you are reading these lines you have the potential to
           survive anything, so we are not that concerned really.
        """
        self.run = run
        self.progress = progress

        self.num_threads = num_threads
        self.evalue = 1e-05
        self.overwrite_output_destinations = overwrite_output_destinations

        self.tmp_dir = tempfile.gettempdir()

        self.query_fasta = query_fasta
        self.target_fasta = target_fasta

        if not self.target_fasta:
            self.target_fasta = self.query_fasta
        elif self.target_db_path:
            self.target_fasta = self.target_db_path

        self.search_program = search_program
        self.search_output_path = 'blast-search-results.txt'
        self.max_target_seqs = None

        utils.is_program_exists('makeblastdb')
        utils.is_program_exists(self.search_program)

        if not self.run.log_file_path:
            self.run.log_file_path = 'blast-log-file.txt'

        # if names_dict is None, all fine. if not, the query_fasta is assumed to be uniqued, and names_dict is
        # the dictionary that connects the ids in the fasta file, to ids that were identical to it.
        self.names_dict = None


    def get_blast_results(self):
        force_makedb, force_blast = False, False

        if self.overwrite_output_destinations:
            force_makedb = True

        if os.path.exists(self.target_fasta + '.phr') and os.path.exists(self.target_fasta + '.pin')\
                                                and os.path.exists(self.target_fasta + '.psq') and not force_makedb:
            self.run.warning("Notice: A BLAST database is found in the output directory, and will be used!")
        else:
            self.makedb()
            force_blast = True

        if os.path.exists(self.search_output_path) and not force_blast:
            self.run.warning("Notice: A BLAST search result is found in the output directory: skipping BLASTP!")
        else:
            self.blast()

        return self.search_output_path


    def check_output(self, expected_output, process='blast'):
        if not os.path.exists(expected_output):
            self.progress.end()
            raise ConfigError("Pfft. Something probably went wrong with '%s' process since one of the expected output files are missing.\
                               Please check the log file here: '%s'" % (process, self.run.log_file_path))


    def makedb(self, output_db_path=None):
        self.progress.new('BLAST')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['makeblastdb',
                    '-in', self.target_fasta,
                    '-dbtype', 'prot',
                    '-out', output_db_path or self.target_fasta]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        expected_output = (output_db_path or self.target_fasta) + '.phr'
        self.check_output(expected_output, 'makeblastdb')

        self.run.info('blast makeblast cmd', cmd_line, quiet=True)
        self.run.info('BLAST search db', self.target_fasta)


    def blast(self):
        cmd_line = [self.search_program,
                    '-query', self.query_fasta,
                    '-db', self.target_fasta,
                    '-evalue', self.evalue,
                    '-outfmt', '6',
                    '-out', self.search_output_path,
                    '-num_threads', self.num_threads]

        if self.max_target_seqs:
            cmd_line += ['-max_target_seqs', self.max_target_seqs]

        self.run.info('NCBI %s cmd' % self.search_program, ' '.join([str(p) for p in cmd_line]), quiet=(not anvio.DEBUG))

        self.progress.new('BLAST')
        self.progress.update('running search (using %s with %d thread(s)) ...' % (self.search_program, self.num_threads))

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.check_output(self.search_output_path, self.search_program)

        if self.names_dict:
            self.ununique_search_results()

        self.run.info('BLAST results', self.search_output_path)


    def ununique_search_results(self):
        self.run.info('self.names_dict is found', 'Un-uniqueing the tabular output.', quiet=True)

        self.progress.new('BLAST')
        self.progress.update('Un-uniqueing the tabular output ...')

        utils.ununique_BLAST_tabular_output(self.search_output_path, self.names_dict)

        self.progress.end()
