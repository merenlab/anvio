# coding: utf-8
""" Foldseek to Pangenome"""

import os
import pandas as pd
import tempfile

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.constants as constants

from anvio.errors import ConfigError
from filesnpaths import AppendableFile


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
    
    def __init__(self, query_fasta=None, run=run, progress=progress, num_threads=1, output_file_path=None, weight_dir=None):
        self.run = run
        self.progress = progress

        utils.is_program_exists('foldseek')

        self.query_fasta = query_fasta
        self.num_threads = num_threads
        self.weight_dir = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/PROSTT5/weights')

        filesnpaths.is_file_exists(self.weight_dir)

        if output_file_path and filesnpaths.check_output_directory(output_file_path):
            self.output_file_path = output_file_path.rstrip('/')
        else:
            raise ConfigError("Oopss. Something probably went wrong with your output file path's '%s'" % (output_file_path))

        if not self.run.log_file_path:
            self.run.log_file_path = 'foldseek-log-file.txt'

        self.additional_params_for_blastp = ""

    def create_db(self):
        self.run.warning(None, header="FOLDSEEK CREATEDB", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('creating the search database (using %d thread(s)) ...' % self.num_threads)

        expected_output_dir = os.path.join(self.output_file_path, "db")
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

    def search(self, query_db, target_db, tmp):
        self.run.warning(None, header="FOLDSEEK EASY SEARCH", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('Running search using Foldseek ...')

        result_file_dir = os.path.join(self.output_file_path, 'result')
        tmp_dir = os.path.join(self.output_file_path, tmp)

        cmd_line = [
            'foldseek',
            'easy-search',
            query_db,
            target_db,
            result_file_dir,
            tmp_dir,
            '--threads', self.num_threads
        ]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Foldseek search Result', result_file_dir)

    def process_result(self, result):
        self.run.warning(None, header="FOLDSEEK PROCESS RESULT", lc="green")
        self.progress.new('FOLDSEEK')
        self.progress.update('Processing Foldseek result...')

        m8_file = os.path.join(self.output_file_path, 'result')
        fasta_file = self.query_fasta
        output_file = os.path.join(self.output_file_path, 'processed_output.m8')

        fasta_lengths = self.read_fasta_lengths(fasta_file)

        self_bit_scores = self.read_self_bit_scores(m8_file)

        selected_rows = []

        # FIXME Use here AppendableFile to control open/write operation
        with open(m8_file, 'r') as file:
            for line in file:
                fields = line.strip().split('\t')
                query = fields[0]
                target = fields[1]
                fident = float(fields[2])
                qstart = int(fields[6])
                qend = int(fields[7])
                bits = float(fields[11])

                query_length = fasta_lengths.get(query, 0)
                coverage = (qend - qstart + 1) / query_length if query_length != 0 else 0

                if query in self_bit_scores and target in self_bit_scores:
                    minbit = bits / min(self_bit_scores[query], self_bit_scores[target])
                else:
                    continue

                if minbit < 0.5 or fident < 0.0:
                    continue

                score = coverage * fident

                selected_rows.append([query, target, str(fident)])

        with open(output_file, 'w') as out_file:
            header = ["Query", "Target", "Identity"]
            out_file.write('\t'.join(header) + '\n')
            for row in selected_rows:
                out_file.write('\t'.join(row) + '\n')

        self.progress.end()
        self.run.info('Processed Foldseek Result', output_file)


    def read_fasta_lengths(self, fasta_file):
        lengths = {}
        fasta = f.SequenceSource(fasta_file)
        while next(fasta):
            lengths[fasta.id] = len(fasta.seq)
        return lengths

    def read_self_bit_scores(self, foldseek_results_file):
        self_bit_scores = {}
        try:
            with open(foldseek_results_file, 'r') as file:
                for line in file:
                    fields = line.strip().split('\t')
                    query = fields[0]
                    target = fields[1]
                    bits = float(fields[11])

                    if query == target:
                        self_bit_scores[query] = bits
        except Exception as e:
            print(f"Error reading self bit scores: {e}")
        return self_bit_scores


class FoldseekSetupWeight:
    """A class to download and setup the weights of PROSTT5"""
    def __init__(self, args=Args(), weight_dir=None, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.weight_dir = args.weight_dir

        if self.weight_dir and args.reset:
            raise ConfigError("You are attempting to run Foldseek setup on a non-default data directory (%s) using the --reset flag. "
                              "To avoid automatically deleting a directory that may be important to you, anvi'o refuses to reset "
                              "directories that have been specified with --weight-dir. If you really want to get rid of this "
                              "directory and regenerate it with InteracDome data inside, then please remove the directory yourself using "
                              "a command like `rm -r %s`. We are sorry to make you go through this extra trouble, but it really is "
                              "the safest way to handle things." % (self.weight_dir, self.weight_dir))

        if not self.weight_dir:
            self.weight_dir = constants.default_foldseek_weight_path

        self.run.warning('', header='Setting up Foldseek Weight', lc='yellow')
        self.run.info('Data directory', self.weight_dir)
        self.run.info('Reset contents', args.reset)

        filesnpaths.is_output_dir_writable(os.path.dirname(os.path.abspath(self.weight_dir)))

        if not args.reset and not anvio.DEBUG:
            self.is_weight_exists()



    def is_weight_exists(self):
        """Raise ConfigError if weight exists"""

        # WE have readme file in that dir rn! you should delete that first
        if (os.path.exists(self.weight_dir)):
            raise ConfigError("It seems you already have the Foldseek Weight downloaded in '%s', please "
                              "use --reset flag if you want to re-download it." % self.weight_dir)

    def setup(self):
        """ Sets up the Foldseek Weight directory for PROSTT5 usage """

        self.run.warning('', header='Downloading Weight Model', lc='yellow')
        self.download_foldseek_weight()

    def download_foldseek_weight(self):
        """ Download the Weights of PROSTT5 models weight """ 

        cmd_line = [
            'foldseek',
            'databases',
            'PROSTT5',
            self.weight_dir,
            'tmp'
        ]

        utils.run_command(cmd_line, self.run.log_file_path)