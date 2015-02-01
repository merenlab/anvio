# -*- coding: utf-8

# Copyright (C) 2014, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

"""
    Here described a class that takes care of the single copy gene analysis given the genes.txt and genes.hmm.gz files,
    and returns a dictionary back.
"""

import os
import gzip
import shutil

import PaPi.utils as utils
import PaPi.terminal as terminal
import PaPi.filesnpaths as filesnpaths

run = terminal.Run()
progress = terminal.Progress()

class SingleCopyGenes:
    def __init__(self, fasta, genes, hmm, ref = None, progress = progress, run = run):
        self.fasta = fasta
        self.single_copy_genes = genes
        self.hmm = hmm
        self.ref = ref
        self.progress = progress
        self.run = run

    def get_results_dict(self):
        self.run.info('Single-copy gene analysis', '', header=True)
        self.run.info('Reference', open(self.ref).readline().strip() if self.ref else 'unknown')
        self.run.info('Pfam model', self.hmm)
        self.run.info('Number of genes', len(self.single_copy_genes))

        results_dict = {}

        tmp_dir = filesnpaths.get_temp_directory_path()
        self.run.info('Temporary work dir', tmp_dir)

        log_file_path = os.path.join(tmp_dir, '00_log.txt')
        self.run.info('Log file', log_file_path)

        self.progress.new('Unpacking the model into temporary work directory')
        self.progress.update('...')
        hmm_file_path = os.path.join(tmp_dir, 'hmm.txt')
        hmm_file = open(hmm_file_path, 'w')
        hmm_file.write(gzip.open(self.hmm, 'rb').read())
        hmm_file.close()
        self.progress.end()

        self.progress.new('Copying contigs into temporary work directory')
        self.progress.update('...')
        fasta_file_path = os.path.join(tmp_dir, 'contigs.fa')
        shutil.copyfile(self.fasta, fasta_file_path)
        self.progress.end()

        self.progress.new('Processing')
        self.progress.update('Compressing the pfam model')
        cmd_line = ('hmmpress "%s" >> "%s" 2>&1' % (hmm_file_path, log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)
        self.progress.end()

        genes_in_contigs = os.path.join(tmp_dir, 'contigs.genes')
        proteins_in_contigs = os.path.join(tmp_dir, 'contigs.proteins')
        self.progress.new('Processing')
        self.progress.update('Identifying open reading frames in contigs ...')
        cmd_line = ('prodigal -i "%s" -o "%s" -a "%s" -p meta >> "%s" 2>&1' % (fasta_file_path,
                                                                               genes_in_contigs,
                                                                               proteins_in_contigs,
                                                                               log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)
        self.progress.end()

        hmm_scan_output = os.path.join(tmp_dir, 'hmm.output')
        hmm_scan_hits = os.path.join(tmp_dir, 'hmm.hits')
        self.progress.new('Processing')
        self.progress.update('Searching for proteins ...')
        cmd_line = ('hmmscan -o "%s" --tblout "%s" "%s" "%s" >> "%s" 2>&1' % (hmm_scan_output,
                                                                               hmm_scan_hits,
                                                                               hmm_file_path,
                                                                               proteins_in_contigs,
                                                                               log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)
        self.progress.end()


        shutil.rmtree(tmp_dir)
        return results_dict