"""Interface to MMseqs2 (currently only `easy-cluster`)."""

import os
import shutil

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


class MMseqs2:
    """Thin wrapper around `mmseqs easy-cluster` returning {cluster_id: [members]}.

    Parameters can be set after instantiation (e.g., `m.min_seq_id = 0.9`).
    Call `get_clusters_dict()` to run clustering and parse the output.
    """

    def __init__(self, input_fasta_path, run=run, progress=progress, num_threads=1):
        self.run = run
        self.progress = progress

        self.input_fasta_path = input_fasta_path
        self.num_threads = num_threads

        # clustering parameters (sensible defaults; callers should override as needed)
        self.min_seq_id = 0.9
        self.coverage = 0.8
        self.cov_mode = 0
        self.cluster_mode = 0
        self.additional_params = ''

        utils.is_program_exists('mmseqs')

        filesnpaths.is_file_fasta_formatted(self.input_fasta_path)

        # output goes into a working directory that we manage internally so we can
        # clean up after ourselves. callers can override after instantiation if they
        # want to keep the files around.
        self.work_dir = filesnpaths.get_temp_directory_path()
        self.output_prefix = os.path.join(self.work_dir, 'mmseqs')
        self.tmp_dir = os.path.join(self.work_dir, 'tmp')
        self.cluster_tsv_path = self.output_prefix + '_cluster.tsv'
        self.rep_fasta_path = self.output_prefix + '_rep_seq.fasta'

        if not self.run.log_file_path:
            self.run.log_file_path = filesnpaths.get_temp_file_path()


    def check_output(self):
        if not os.path.exists(self.cluster_tsv_path):
            self.progress.end()
            raise ConfigError("MMseqs2 did not produce a cluster TSV at the expected path '%s'. "
                              "Please check the log file here: '%s'." % (self.cluster_tsv_path, self.run.log_file_path))


    def get_clusters_dict(self, name_prefix='C'):
        self.cluster()

        # mmseqs `_cluster.tsv` is a two-column file: <representative>\t<member>
        # (the representative is listed as a member of its own cluster).
        clusters = {}
        with open(self.cluster_tsv_path) as f:
            for line in f:
                rep, member = line.rstrip('\n').split('\t')
                clusters.setdefault(rep, []).append(member)

        # assign deterministic, ordered cluster IDs
        clusters_dict = {}
        representatives = {}
        for i, rep in enumerate(sorted(clusters.keys()), start=1):
            cluster_id = f"{name_prefix}_{i:08d}"
            clusters_dict[cluster_id] = clusters[rep]
            representatives[cluster_id] = rep

        self.run.info('Number of MMseqs2 clusters', '%s' % pp(len(clusters_dict)))

        # expose the representative id for each cluster so callers can pick out
        # the representative sequence later if they want.
        self.representatives = representatives

        return clusters_dict


    def cluster(self):
        self.run.warning(None, header="MMseqs2 easy-cluster", lc="green")
        self.run.info('Min sequence identity', self.min_seq_id)
        self.run.info('Coverage threshold', self.coverage)
        self.run.info('Coverage mode', self.cov_mode)

        self.progress.new('MMseqs2')
        self.progress.update('clustering (using %d thread(s)) ...' % self.num_threads)

        cmd_line = ['mmseqs', 'easy-cluster',
                    self.input_fasta_path,
                    self.output_prefix,
                    self.tmp_dir,
                    '--threads', self.num_threads,
                    '--min-seq-id', self.min_seq_id,
                    '-c', self.coverage,
                    '--cov-mode', self.cov_mode,
                    '--cluster-mode', self.cluster_mode]

        if self.additional_params:
            cmd_line += self.additional_params.split()

        self.run.info('mmseqs cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.check_output()

        self.run.info('MMseqs2 cluster TSV', self.cluster_tsv_path)


    def get_representative_sequences(self):
        """Return {representative_id: sequence} parsed from the rep_seq.fasta."""
        reps = {}
        if not os.path.exists(self.rep_fasta_path):
            return reps

        current_id = None
        current_seq = []
        with open(self.rep_fasta_path) as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_id is not None:
                        reps[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id is not None:
                reps[current_id] = ''.join(current_seq)

        return reps


    def cleanup(self):
        """Remove the temporary working directory."""
        if os.path.isdir(self.work_dir):
            shutil.rmtree(self.work_dir, ignore_errors=True)
