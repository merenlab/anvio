# -*- coding: utf-8
"""
    Classes for gene calling.
"""

import os
import sys
import shutil

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.fastalib as fastalib
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

run = terminal.Run()
progress = terminal.Progress()


class Prodigal:
    def __init__(self, progress = progress, run = run):
        self.progress = progress
        self.run = run

        self.gene_calls_dict = {} # each entry must contain {'contig', 'start', stop, 'direcion', 'partial'} items.
        self.protein_sequences_dict = {}


    def process(self, fasta_file_path, output_dir):
        self.genes_in_contigs = os.path.join(output_dir, 'contigs.genes')
        self.proteins_in_contigs = os.path.join(output_dir, 'contigs.proteins')

        log_file_path = os.path.join(output_dir, '00_log.txt')

        self.run.warning('', header = 'Finding ORFs in contigs', lc = 'green')
        self.run.info('Genes', self.genes_in_contigs)
        self.run.info('Proteins', self.proteins_in_contigs)
        self.run.info('Log file', log_file_path)

        self.progress.new('Processing')
        self.progress.update('Identifying ORFs in contigs ...')
        cmd_line = ('prodigal -i "%s" -o "%s" -a "%s" -p meta >> "%s" 2>&1' % (fasta_file_path,
                                                                               self.genes_in_contigs,
                                                                               self.proteins_in_contigs,
                                                                               log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)

        if not os.path.exists(self.proteins_in_contigs):
            self.progress.end()
            raise ConfigError, "Something went wrong with prodigal, and it failed to generate the\
                                expected output :/ Fortunately, this log file should tell you what\
                                might be the problem: '%s'. Please do not forget to include this\
                                file if you were to ask for help." % log_file_path

        self.progress.update('Processing gene calls ...')

        fasta = fastalib.SequenceSource(self.proteins_in_contigs)
        while fasta.next():
            print fasta.id
        fasta.close()

        self.progress.end()


        sys.exit()


class GeneCaller:
    def __init__(self, fasta_file_path, gene_caller = 'prodigal', progress = progress, run = run, debug = False):
        filesnpaths.is_file_exists(fasta_file_path)
        filesnpaths.is_file_fasta_formatted(fasta_file_path)

        self.fasta_file_path = fasta_file_path

        self.run = run
        self.progress = progress

        self.debug = debug
        self.tmp_dirs = []

        self.gene_callers = {'prodigal': Prodigal}

        self.gene_caller = gene_caller

        if self.gene_caller not in self.gene_callers:
            raise ConfigError, "The gene caller you requested ('%s') is not available at this point.\
                                here is a list of what we have: %s." % (', '.join(self.gene_callers))


    def process(self):
        output_dir = filesnpaths.get_temp_directory_path()
        self.tmp_dirs.append(output_dir)
        gene_caller = self.gene_callers[self.gene_caller]()

        gene_calls_dict, protein_sequences_dict = gene_caller.process(self.fasta_file_path, output_dir)

        if not self.debug:
            self.clean_tmp_dirs()
 
        return gene_calls_dict, protein_sequences_dict


    def clean_tmp_dirs(self):
        for tmp_dir in self.tmp_dirs:
            shutil.rmtree(tmp_dir)
