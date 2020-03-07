"""Interface to mmseqs2."""

import anvio

import anvio.filesnpaths as filesnpaths
import anvio.terminal as terminal
import anvio.utils as utils

import os
import shutil

from multiprocessing import cpu_count

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"

pp = terminal.pretty_print

class MMseqs2:
    def __init__(
        self,
        output_dir_path,
        num_threads=cpu_count(),
        verbosity_level=3,
        progress=terminal.Progress(),
        run=terminal.Run()):
        """A class to streamline MMseqs2 commands."""
        self.top_output_dir = filesnpaths.gen_output_directory(output_dir_path)
        self.num_threads = num_threads
        self.verbosity_level = verbosity_level
        self.run = run
        self.progress = progress
        # ULTIMATLEY CHANGE NAMED DIRS TO RANDOMLY GENERATED TEMP DIRS

        # utils.is_program_exists('mmseqs')


    def run_cluster_pipeline(self, fasta_file_paths):
        if type(fasta_file_paths) == str:
            fasta_file_paths = [fasta_file_paths]

        self.run.warning("Anvi'o will use 'MMseqs2' "
                         "by Steinegger et al. (DOI: 10.1038/nbt.3988) to cluster. "
                         "If you publish your findings, "
                         "please do not forget to properly credit their work.",
                         lc='green', header="CITATION")

        seq_db_file = self.run_createdb(fasta_file_paths)
        cluster_db_file = self.run_cluster(seq_db_file)
        cluster_tsv_file = self.run_createtsv(seq_db_file, seq_db_file, cluster_db_file)


    def run_createdb(self, fasta_file_paths, compressed=0):
        seq_db_dir = os.path.join(self.top_output_dir, 'DB')
        shutil.rmtree(seq_db_dir, ignore_errors=True)
        filesnpaths.gen_output_directory(seq_db_dir)
        seq_db = os.path.join(seq_db_dir, 'DB')

        cmd_line = ['mmseqs',
                    'createdb',
                    *fasta_file_paths,
                    seq_db,
                    '--compressed', compressed,
                    '-v', self.verbosity_level]
        log_file = os.path.join(seq_db_dir, '00_log.txt')

        self.run.info('[MMseqs2 createdb] Input FASTA file paths', ' ; '.join(fasta_file_paths))
        self.run.info('[MMseqs2 createdb] Output sequence database directory', seq_db_dir)
        self.run.info('[MMseqs2 createdb] Compressed', compressed)
        self.run.info('[MMSeqs2 createdb] Log file path', log_file)

        utils.run_command(cmd_line, log_file)

        return seq_db


    def run_cluster(self, seq_db, cov_thresh=0.8, cov_mode=2, min_seq_id=0.98):
        cluster_db_dir = os.path.join(self.top_output_dir, 'CLUSTERDB')
        shutil.rmtree(cluster_db_dir, ignore_errors=True)
        filesnpaths.gen_output_directory(cluster_db_dir)
        cluster_db_tmp_dir = filesnpaths.gen_output_directory(
            os.path.join(cluster_db_dir, 'TMP'))
        cluster_db = os.path.join(cluster_db_dir, 'CLUSTERDB')

        cmd_line = ['mmseqs',
                    'cluster',
                    seq_db,
                    cluster_db,
                    cluster_db_tmp_dir,
                    '-c', cov_thresh,
                    '--cov-mode', cov_mode,
                    '-a',
                    '--min-seq-id', min_seq_id,
                    '--threads', self.num_threads,
                    '-v', self.verbosity_level]
        log_file = os.path.join(cluster_db_dir, '00_log.txt')

        self.run.info('[MMseqs2 cluster] Input sequence database directory',
                      os.path.dirname(seq_db))
        self.run.info('[MMseqs2 cluster] Output cluster database directory', cluster_db_dir)
        self.run.info('[MMseqs2 cluster] Output cluster database temp directory',
                      cluster_db_tmp_dir)
        self.run.info('[MMseqs2 cluster] Coverage threshold', cov_thresh)
        self.run.info('[MMseqs2 cluster] Coverage mode', cov_mode)
        self.run.info('[MMseqs2 cluster] Sequence identity threshold', min_seq_id)
        self.run.info('[MMseqs2 cluster] Number of threads', self.num_threads)
        self.run.info('[MMseqs2 cluster] Log file path', log_file)

        utils.run_command(cmd_line, log_file)

        return cluster_db


    def run_createtsv(self, query_db, target_db, result_db):
        tsv_dir = os.path.join(os.path.dirname(result_db), 'TSV')
        shutil.rmtree(tsv_dir, ignore_errors=True)
        filesnpaths.gen_output_directory(tsv_dir)
        tsv_file = os.path.join(tsv_dir, os.path.basename(result_db) + '.tsv')

        cmd_line = ['mmseqs',
                    'createtsv',
                    query_db,
                    target_db,
                    result_db,
                    tsv_file,
                    '--threads', self.num_threads,
                    '-v', self.verbosity_level]
        log_file = os.path.join(tsv_dir, '00_log.txt')

        self.run.info('[MMseqs2 createtsv] Query database directory', os.path.dirname(query_db))
        self.run.info('[MMseqs2 createtsv] Target database directory', os.path.dirname(target_db))
        self.run.info('[MMseqs2 createtsv] Result database directory', os.path.dirname(result_db))
        self.run.info('[MMseqs2 createtsv] Output tsv file', tsv_file)
        self.run.info('[MMseqs2 createtsv] Number of threads', self.num_threads)
        self.run.info('[MMseqs2 createtsv] Log file path', log_file)

        utils.run_command(cmd_line, log_file)

        return tsv_file

mmseqs2 = MMseqs2('/Users/sammiller/Documents/work/mmseqs2_test/ANVIO_OUTPUT')
mmseqs2.run_cluster_pipeline('/Users/sammiller/Documents/work/mmseqs2_test/FASTA/test.fasta')