# coding: utf-8
"""Interface to MetaBAT."""
import os
import glob
from subprocess import Popen, PIPE

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Özcan Esen"
__email__ = "ozcanesen@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class MetaBAT:
    arguments = {
        "infile"
        'seed': (
                ['--seed'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Seed for random numbers"}
                    ),
        'threads': (
                ['-T', '--threads'],
                {'metavar': "INT",
                 'required': False,
                 'help': "Number of threads"}
                    ),
    }
    citation = "asdads"

    def __init__(self, profile_db=None, contigs_db=None, run=run, progress=progress):
        self.run = run
        self.progress = progress
        self.program_name = 'metabat2'

        self.profile_db = profile_db
        self.contigs_db = contigs_db

        utils.is_program_exists(self.program_name)

        self.temp_path = filesnpaths.get_temp_directory_path()


    def cluster(self, input_files, args):
        input_files.coverage

        self.run.info('Working directory', self.temp_path)

        P = lambda x: os.path.join(self.temp_path, x)
        coverage_path = P('coverages.txt')
        fasta_path = P('sequence.fa')
        log_path = P('logs.txt')
        bin_prefix = P('Bin')

        coverages = self.profile_db.db.get_table_as_dict('mean_coverage_contigs')

        utils.export_sequences_from_contigs_db(self.contigs_db.db_path, fasta_path, seq_names_to_export=sorted(coverages.keys()), splits_mode=True)
        utils.store_dict_as_TAB_delimited_file(coverages, coverage_path, ['contig'] + sorted(list(self.profile_db.samples)))

        cmd_line = [self.program_name, '-i', fasta_path, '-a', coverage_path, '-o', bin_prefix ,'--cvExt', '-l', '--verbose', '--debug']

        self.progress.new(self.program_name)
        self.progress.update('Running...')
        utils.run_command(cmd_line, log_path)
        self.progress.end()

        clusters = {}
        bin_count = 0
        for bin_file in glob.glob(bin_prefix + '*'):
            bin_count += 1
            with open(bin_file, 'r') as f:
                bin_name = os.path.basename(bin_file).replace('.', '_')
                clusters[bin_name] = list(map(str.strip, f.readlines()))

        self.run.info('Bins formed', bin_count)

        return clusters
        