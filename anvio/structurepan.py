#!/usr/bin/env python
# coding: utf-8
""" Foldseek to Pangenome"""

import os
import argparse
import pandas as pd
import tempfile

import anvio
import anvio.fastalib as f
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.constants as constants

from collections import defaultdict
from anvio.drivers.mcl import MCL
from anvio.drivers.foldseek import Foldseek
from anvio.errors import ConfigError
from anvio.filesnpaths import AppendableFile


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Metehan Sever"
__email__ = "metehaansever@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print

class StructurePan(object):
    
    def __init__(self, query_fasta=None, run=run, progress=progress, num_threads=1, output_file_path=None, result=None):
        self.run = run
        self.progress = progress

        utils.is_program_exists('foldseek')

        if not query_fasta:
            raise ConfigError("Oopss. Something probably went wrong with your query fasta file path's '%s'" % (query_fasta))

        self.query_fasta = query_fasta

        if result and filesnpaths.is_file_exists(result):
            self.result = result
        else:
            raise ConfigError("Oopss. Something probably went wrong with your result file path's '%s'" % (result))

        self.num_threads = num_threads

        if output_file_path:
            self.output_file_path = output_file_path.rstrip('/')
        else:
            raise ConfigError("You must provide a output directory path '%s'" % (output_file_path))

        if not self.run.log_file_path:
            self.run.log_file_path = 'structurepan-log-file.txt'

        
    def process_result(self):
        self.run.warning(None, header="FOLDSEEK PROCESS RESULT", lc="green")
        self.progress.new('FOLDSEEK PROCESS RESULT')
        self.progress.update('Processing Foldseek result...')

        m8_file = os.path.join(self.output_file_path, self.result)
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

    def mcl_network(self, processed_result):
        """ Turn to Search Results to Network with MCL"""

        self.run.warning(None, header="FOLDSEEK MCL NETWORK", lc="green")
        self.progress.new('FOLDSEEK MCL NETWORK')
        self.progress.update('Processing Foldseek mcl network...')

        # Convert M8 to ABC format
        abc_file = os.path.join(self.output_file_path, 'mcl_input.abc')
        self.convert_m8_to_abc(processed_result, abc_file)

        # Create an instance of MCL
        mcl_instance = MCL(mcl_input_file_path=abc_file, num_threads=self.num_threads)

        # Set the inflation parameter (MCL-specific)
        # FIXME default value should be 2 and should be taken parametrically from the user 
        mcl_instance.inflation = 2

        # Run the clustering algorithm
        clusters_dict = mcl_instance.get_clusters_dict(name_prefix='Cluster')

        output_clusters_file = os.path.join(self.output_file_path, 'clusters.txt')
        tsv_output_file = os.path.join(self.output_file_path, 'clusters_output.tsv')

        # Output clusters to a text file
        # FIXME Use here AppendableFile to control open/write operation 
        with open(output_clusters_file, 'w') as outfile:
            for cluster_name, members in clusters_dict.items():
                outfile.write(f"{cluster_name}\t{','.join(members)}\n")

        print(f"Clusters have been written to {output_clusters_file}")

        # Create a set of all unique identifiers
        unique_identifiers = set()
        for members in clusters_dict.values():
            for member in members:
                unique_identifiers.add(self.extract_unique_identifier(member))

        # Create a DataFrame for the clusters
        cluster_data = defaultdict(lambda: {uid: 0 for uid in unique_identifiers})

        for cluster_name, members in clusters_dict.items():
            for member in members:
                uid = self.extract_unique_identifier(member)
                cluster_data[cluster_name][uid] += 1

        # Convert to DataFrame
        cluster_df = pd.DataFrame(cluster_data).T.reset_index().rename(columns={"index": "Cluster"})

        # Reorder columns
        columns = ["Cluster"] + sorted(unique_identifiers)
        cluster_df = cluster_df[columns]

        # Write to TSV (TAB-delimited) file
        cluster_df.to_csv(tsv_output_file, sep='\t', index=False)
        print(f"Clusters have been written to {tsv_output_file}")

        self.progress.end()
        self.run.info('MCL Network', tsv_output_file)

    def convert_m8_to_abc(self, m8_file, abc_file):
        with open(m8_file, 'r') as infile, open(abc_file, 'w') as outfile:
            next(infile)  # Skip header
            for line in infile:
                fields = line.strip().split('\t')
                query = fields[0]
                target = fields[1]
                bit_score = fields[2]
                outfile.write(f"{query}\t{target}\t{bit_score}\n")

    def extract_unique_identifier(self, member):
        return member.split('_')[0]


class FoldseekSetupWeight:
    """A class to download and setup the weights of PROSTT5"""
    def __init__(self, args, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.weight_dir = args.foldseek_weight_dir

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

        if not self.run.log_file_path:
            self.run.log_file_path = os.path.join(self.weight_dir, 'foldseek-setup-log-file.txt')

        filesnpaths.gen_output_directory(self.weight_dir, delete_if_exists=args.reset)

    def is_weight_exists(self):
        """Raise ConfigError if weight exists"""

        if os.path.exists(self.weight_dir) and os.listdir(self.weight_dir):
            raise ConfigError("It seems you already have the Foldseek Weight downloaded in '%s', please "
                              "use --reset flag if you want to re-download it." % self.weight_dir)

    def setup(self):
        """ Sets up the Foldseek Weight directory for PROSTT5 usage """

        self.run.warning('', header='Downloading Weight Model', lc='yellow')
        self.download_foldseek_weight()

    def download_foldseek_weight(self):
        """ Download the Weights of PROSTT5 models weight """ 

        self.progress.new('FOLDSEEK')
        self.progress.update('Downloading ...')

        cmd_line = [
            'foldseek',
            'databases',
            'ProstT5',
            self.weight_dir,
            os.path.join(self.weight_dir, 'tmp')
        ]

        utils.run_command(cmd_line, self.run.log_file_path)

        self.progress.end()

        self.run.info('Command line', ' '.join([str(x) for x in cmd_line]), quiet=True)
        self.run.info('Log file', self.run.log_file_path)