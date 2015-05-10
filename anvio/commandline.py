# -*- coding: utf-8
"""
    A place for classes for binary calls.

    At this point there is one class described which takes care of search using HMM profiles.
    It simply takes genes.txt and genes.hmm.gz files as input and returns a dictionary back
    with results. See anvio/data/hmm directory for examples.
"""

import os
import gzip
import shutil

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths


run = terminal.Run()
progress = terminal.Progress()


class HMMSearch:
    def __init__(self, progress = progress, run = run):
        self.progress = progress
        self.run = run

        # hmm_scan_hits and proteins_in_contigs are the files to access later on for
        # parsing:
        self.hmm_scan_output = None
        self.hmm_scan_hits = None
        self.genes_in_contigs = None
        self.proteins_in_contigs = None

        self.tmp_dirs = []


    def run_prodigal(self, fasta_file_path):
        tmp_dir = filesnpaths.get_temp_directory_path()
        self.tmp_dirs.append(tmp_dir)

        self.genes_in_contigs = os.path.join(tmp_dir, 'contigs.genes')
        self.proteins_in_contigs = os.path.join(tmp_dir, 'contigs.proteins')

        log_file_path = os.path.join(tmp_dir, '00_log.txt')

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
        self.progress.end()

        return self.proteins_in_contigs


    def run_hmmscan(self, source, genes_in_model, hmm, ref, cut_off_flag = "--cut_ga"):
        self.run.warning('', header = 'HMM Profiling for %s' % source, lc = 'green')
        self.run.info('Reference', ref if ref else 'unknown')
        self.run.info('Pfam model', hmm)
        self.run.info('Number of genes', len(genes_in_model))

        tmp_dir = filesnpaths.get_temp_directory_path()
        self.tmp_dirs.append(tmp_dir)

        self.hmm_scan_output = os.path.join(tmp_dir, 'hmm.output')
        self.hmm_scan_hits = os.path.join(tmp_dir, 'hmm.hits')
        self.hmm_scan_hits_shitty = os.path.join(tmp_dir, 'hmm.hits.shitty')
        log_file_path = os.path.join(tmp_dir, '00_log.txt')

        self.run.info('Temporary work dir', tmp_dir)
        self.run.info('HMM scan output', self.hmm_scan_output)
        self.run.info('HMM scan hits', self.hmm_scan_hits)
        self.run.info('Log file', log_file_path)

        self.progress.new('Unpacking the model into temporary work directory')
        self.progress.update('...')
        hmm_file_path = os.path.join(tmp_dir, 'hmm.txt')
        hmm_file = open(hmm_file_path, 'w')
        hmm_file.write(gzip.open(hmm, 'rb').read())
        hmm_file.close()
        self.progress.end()

        self.progress.new('Processing')
        self.progress.update('Compressing the pfam model')
        cmd_line = ('hmmpress "%s" >> "%s" 2>&1' % (hmm_file_path, log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)
        self.progress.end()

        self.progress.new('Processing')
        self.progress.update('Performing HMM scan ...')
        cmd_line = ('hmmscan -o "%s" %s --tblout "%s" "%s" "%s" >> "%s" 2>&1' % (self.hmm_scan_output,
                                                                              cut_off_flag,
                                                                              self.hmm_scan_hits_shitty,
                                                                              hmm_file_path,
                                                                              self.proteins_in_contigs,
                                                                              log_file_path))
        with open(log_file_path, "a") as myfile: myfile.write('CMD: ' + cmd_line + '\n')
        utils.run_command(cmd_line)
        self.progress.end()

        # this shit is here because hmmscan does not courteous enough to generate a TAB-delimited
        # file......
        parseable_output = open(self.hmm_scan_hits, 'w')
        for line in open(self.hmm_scan_hits_shitty).readlines():
            if line.startswith('#'):
                continue
            parseable_output.write('\t'.join(line.split()[0:18]) + '\n')
        parseable_output.close()

        num_raw_hits = filesnpaths.get_num_lines_in_file(self.hmm_scan_hits)
        self.run.info('Number of raw hits', num_raw_hits)

        return self.hmm_scan_hits if num_raw_hits else None


    def clean_tmp_dirs(self):
        for tmp_dir in self.tmp_dirs:
            shutil.rmtree(tmp_dir)
