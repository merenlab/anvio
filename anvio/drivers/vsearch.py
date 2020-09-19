# coding: utf-8
"""Interface to VSEARCH."""

import pandas as pd
import os.path
import tempfile

import anvio
import anvio.terminal as terminal
import anvio.utils as utils

from anvio.errors import ConfigError


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"


class Vsearch:
    def __init__(
        self,
        input_fasta_path=None,
        num_threads=1,
        overwrite_output_destinations=False,
        temp_dir=None,
        run=None,
        progress=None):

        self.input_fasta_path = input_fasta_path
        if run:
            self.run = run
        else:
            terminal.Run()
        if progress:
            self.progress = progress
        else:
            terminal.Progress()
        self.num_threads = num_threads
        self.overwrite_output_destinations = overwrite_output_destinations
        self.temp_dir = temp_dir if temp_dir else tempfile.gettempdir()
        self.temp_file_path_list = []

        if not self.run.log_file_path:
            self.run.log_file_path = 'vsearch-log-file.txt'


    def get_output_path(self, output_path=None, file_name='output_file'):
        if output_path:
            output_path = output_path
        else:
            output_path = os.path.join(self.temp_dir, file_name)

        if os.path.exists(output_path):
            if self.overwrite_output_destinations:
                os.remove(output_path)
            else:
                raise ConfigError(
                    "VSEARCH was going to create an output file, %s, but it already exists. "
                    "You can set the class attribute, overwrite_output_destinations, to True "
                    "to automatically remove conflicting files." % output_path)

        return output_path


    def remove_temp_files(self):
        for temp_file in self.temp_file_path_list:
            os.remove(temp_file)


class Dereplication(Vsearch):
    def __init__(
        self,
        input_fasta_path=None,
        output_fasta_path=None,
        output_cluster_path=None,
        min_seq_length=32, # The VSEARCH default is 32
        run=None,
        progress=None,
        num_threads=1,
        overwrite_output_destinations=False):

        super().__init__(
            input_fasta_path=input_fasta_path,
            run=run,
            progress=progress,
            num_threads=1,
            overwrite_output_destinations=False)

        if output_fasta_path:
            self.output_fasta_path = output_fasta_path
        else:
            self.temp_file_path_list.append(self.get_output_path(file_name='rep_seqs.fasta'))
            self.output_fasta_path = self.temp_file_path_list[-1]

        # USEARCH cluster (UC) formatted tab-delimited file
        if output_cluster_path:
            self.output_cluster_path = output_cluster_path
        else:
            self.temp_file_path_list.append(self.get_output_path(file_name='clusters_uc.txt'))
            self.output_cluster_path = self.temp_file_path_list[-1]

        self.min_seq_length = min_seq_length


    def get_cluster_membership_dict(self):
        ''' Get a dict of cluster membership from the USEARCH-formatted cluster (UC) file '''

        # Files with this format are produced by USEARCH and VSEARCH.
        # The table contains mixed information in two parts.
        # The first part (rows that come first) contains a line for each input sequence.
        # Clusters are reported in blocks in order of cluster size,
        # with the seed sequence being the first line (record type "S").
        # The seed sequence line has the sequence name as the query label and an asterisk as the target label.
        # Member sequence lines (record type "H") have the seed sequence name as the target label.
        # The "size" column in this part means sequence length.
        # The second part contains a line for each cluster (record type "C").
        # The "size" column in this part means cluster size.

        cluster_df = pd.read_csv(
            self.output_cluster_path,
            sep='\t',
            names=['record_type', 'clust_num', 'size', 'pc_ident', 'strand', 'obsolete1', 'obsolete2', 'align', 'query_label', 'target_label'])
        cluster_df = cluster_df[(cluster_df['record_type'] == 'S') | (cluster_df['record_type'] == 'H')]

        cluster_membership_dict = {}
        for query_label, target_label in zip(
            cluster_df['query_label'].tolist(), cluster_df['target_label'].tolist()):
            if target_label == '*':
                cluster_membership_dict[query_label] = [query_label]
            else:
                cluster_membership_dict[target_label].append(query_label)

        return cluster_membership_dict


class FullLengthDereplication(Dereplication):
    def __init__(
        self,
        input_fasta_path=None,
        output_fasta_path=None,
        output_cluster_path=None,
        min_seq_length=32, # The VSEARCH default is 32
        run=None,
        progress=None,
        num_threads=1,
        overwrite_output_destinations=False):

        super().__init__(
            input_fasta_path=input_fasta_path,
            output_fasta_path=output_fasta_path,
            output_cluster_path=output_cluster_path,
            min_seq_length=min_seq_length, # The VSEARCH default is 32
            run=run,
            progress=progress,
            num_threads=1,
            overwrite_output_destinations=overwrite_output_destinations)

        if not progress:
            self.progress.new("VSEARCH full-length dereplication")


    def dereplicate(self):

        cmd_line = [
            'vsearch',
            '--derep_fulllength', self.input_fasta_path,
            '--output', self.output_fasta_path,
            '--uc', self.output_cluster_path,
            '--minseqlength', self.min_seq_length,
            '--fasta_width', 0, # A value of 0 means input fasta lines are not wrapped
            '--threads', self.num_threads]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.run.info('vsearch full-length derep cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('vsearch full-length derep fasta', self.output_fasta_path)
        self.run.info('vsearch full-length derep uc clusters', self.output_cluster_path)


class PrefixDereplication(Dereplication):
    def __init__(
        self,
        input_fasta_path=None,
        output_fasta_path=None,
        output_cluster_path=None,
        min_seq_length=32, # The VSEARCH default is 32
        run=None,
        progress=None,
        num_threads=1,
        overwrite_output_destinations=False):

        super().__init__(
            input_fasta_path=input_fasta_path,
            output_fasta_path=output_fasta_path,
            output_cluster_path=output_cluster_path,
            min_seq_length=min_seq_length, # The VSEARCH default is 32
            run=run,
            progress=progress,
            num_threads=1,
            overwrite_output_destinations=overwrite_output_destinations)

        if not progress:
            self.progress.new("VSEARCH prefix dereplication")


    def dereplicate(self):

        cmd_line = [
            'vsearch',
            '--derep_prefix', self.input_fasta_path,
            '--output', self.output_fasta_path,
            '--uc', self.output_cluster_path,
            '--minseqlength', self.min_seq_length,
            '--fasta_width', 0, # A value of 0 means input fasta lines are not wrapped
            '--threads', self.num_threads]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.run.info('vsearch prefix derep cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('vsearch prefix derep fasta', self.output_fasta_path)
        self.run.info('vsearch prefix derep uc clusters', self.output_cluster_path)


class ReverseComplementation(Vsearch):
    def __init__(
        self,
        input_fasta_path=None,
        output_fasta_path=None,
        run=None,
        progress=None,
        num_threads=1,
        overwrite_output_destinations=False):

        super().__init__(
            input_fasta_path=input_fasta_path,
            run=run,
            progress=progress,
            num_threads=1,
            overwrite_output_destinations=overwrite_output_destinations)

        if output_fasta_path:
            self.output_fasta_path = output_fasta_path
        else:
            self.temp_file_path_list.append(self.get_output_path(file_name='rev_comp.fasta'))
            self.output_fasta_path = self.temp_file_path_list[-1]

        if not progress:
            self.progress.new("VSEARCH reverse complementation")


    def reverse_complement(self):

        cmd_line = [
            'vsearch',
            '--fastx_revcomp', self.input_fasta_path,
            '--output', self.output_fasta_path,
            '--minseqlength', self.min_seq_length,
            '--fasta_width', 0, # A value of 0 means input fasta lines are not wrapped
            '--threads', self.num_threads]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.run.info('vsearch revcomp cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('vsearch revcomp fasta', self.output_fasta_path)