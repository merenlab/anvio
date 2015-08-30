# -*- coding: utf-8
"""
    Classes for HMM related operations.

    * HMMSearch takes care of searches using HMM profiles. It simply takes genes.txt and
    genes.hmm.gz files as input and returns a dictionary back with results. See anvio/data/hmm
    directory for examples.
"""

import os
import gzip
import shutil
import textwrap

import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError

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

        if not os.path.exists(self.proteins_in_contigs):
            raise ConfigError, "Something went wrong with prodigal, and it failed to generate the\
                                expected output :/ Fortunately, this log file should tell you what\
                                might be the problem: '%s'. Please do not forget to include this\
                                file if you were to ask for help." % log_file_path

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
        ret_val = utils.run_command(cmd_line)
        if ret_val:
            raise ConfigError, "The last call did not work quite well. Most probably the version of HMMER\
                                you have installed is not up-to-date enough. Just to make sure what went\
                                wrong please take a look at the log file ('%s'). Please visit %s to see what\
                                is the latest version availalbe. You can learn which version of HMMER you have\
                                on your system by typing 'hmmpress -h'"\
                                        % (log_file_path, 'http://hmmer.janelia.org/download.html')
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

        if not os.path.exists(self.hmm_scan_hits_shitty):
            raise ConfigError, "Something went wrong with hmmscan, and it failed to generate the\
                                expected output :/ Fortunately, this log file should tell you what\
                                might be the problem: '%s'. Please do not forget to include this\
                                file if you were to ask for help." % log_file_path

        self.progress.end()

        # thank you, hmmscan, for not generating a simple TAB-delimited, because we programmers
        # love to write little hacks like this into our code:
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


class SequencesForHMMHits:
    def __init__(self, contigs_db_path, sources = set([]), run = run, progress = progress):
        if type(sources) != type(set([])):
            raise ConfigError, "'sources' variable has to be a set instance."

        self.sources = set([s for s in sources if s])

        # take care of contigs db related stuff and move on:
        contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)
        self.search_info_table = contigs_db.get_table_as_dict(t.hmm_hits_info_table_name)
        self.search_table_splits = contigs_db.get_table_as_dict(t.hmm_hits_splits_table_name)
        self.search_table_contigs = contigs_db.get_table_as_dict(t.hmm_hits_contigs_table_name)
        self.contig_sequences = contigs_db.get_table_as_dict(t.contig_sequences_table_name, string_the_key = True)
        contigs_db.disconnect()

        missing_sources = [s for s in self.sources if s not in self.search_info_table]
        if len(missing_sources):
            raise ConfigError, 'Some of the requested sources were not found in the contigs database :/\
                                Here is a list of the ones that are missing: %s' % ', '.join(missing_sources)

        if len(self.sources):
            self.search_table_splits = utils.get_filtered_dict(self.search_table_splits, 'source', self.sources)
            self.search_table_contigs = utils.get_filtered_dict(self.search_table_contigs, 'source', self.sources)
        else:
            self.sources = self.search_info_table.keys()

        # create a map of all unique gene ids to contigs table entry ids for fast access:
        self.unique_id_to_contig_entry_id = {}
        for entry_id in self.search_table_contigs:
            unique_id = self.search_table_contigs[entry_id]['gene_unique_identifier']
            self.unique_id_to_contig_entry_id[unique_id] = entry_id


    def get_hmm_sequences_dict_for_splits(self, splits_dict):
        """splits dict is what you get from ccollections.GetSplitNamesInBins(args).get_dict(), and
           its struture goes like this:

                {
                    'bin_x': set['split_a, split_b, ...'],
                    'bin_y': set['split_c, split_d, ...'],
                    ...
                }
        """

        split_names = set([])
        for s in splits_dict.values():
            split_names.update(s)

        hits = utils.get_filtered_dict(self.search_table_splits, 'split', split_names)

        split_name_to_bin_id = {}
        for bin_id in splits_dict:
            for split_name in splits_dict[bin_id]:
                split_name_to_bin_id[split_name] = bin_id

        if not hits:
            return {}

        hmm_sequences_dict_for_splits = {}

        unique_ids_taken_care_of = set([])
        for entry in hits.values():
            split_name = entry['split']
            source = entry['source']
            gene_name = entry['gene_name']
            e_value = entry['e_value']
            gene_unique_id = entry['gene_unique_identifier']

            if gene_unique_id in unique_ids_taken_care_of:
                continue
            else:
                unique_ids_taken_care_of.add(gene_unique_id)

            contig_entry = self.search_table_contigs[self.unique_id_to_contig_entry_id[gene_unique_id]]
            contig_name = contig_entry['contig']
            start, stop = contig_entry['start'], contig_entry['stop']
            sequence = self.contig_sequences[contig_name]['sequence'][start:stop]

            hmm_sequences_dict_for_splits[gene_unique_id] = {'sequence': sequence,
                                                             'source': source,
                                                             'bin_id': split_name_to_bin_id[split_name],
                                                             'gene_name': gene_name,
                                                             'e_value': e_value,
                                                             'contig': contig_entry['contig'],
                                                             'start': start,
                                                             'stop': stop,
                                                             'length': stop - start}

        return hmm_sequences_dict_for_splits


    def get_FASTA_header_and_sequence_for_gene_unique_id(self, hmm_sequences_dict_for_splits, gene_unique_id):
        entry = hmm_sequences_dict_for_splits[gene_unique_id]
        header = '%s___%s|' % (entry['gene_name'], gene_unique_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['bin_id', 'source', 'e_value', 'contig', 'start', 'stop', 'length']])
        sequence = hmm_sequences_dict_for_splits[gene_unique_id]['sequence']
        return (header, sequence)


    def store_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap = 200):
        filesnpaths.is_output_file_writable(output_file_path)

        if type(wrap) != int:
            raise ConfigError, '"wrap" has to be an integer instance'

        f = open(output_file_path, 'w')

        for gene_unique_id in hmm_sequences_dict_for_splits:
            header, sequence = self.get_FASTA_header_and_sequence_for_gene_unique_id(hmm_sequences_dict_for_splits, gene_unique_id)

            if wrap:
                sequence = textwrap.fill(sequence, wrap, break_on_hyphens = False)

            f.write('>%s\n' % header)
            f.write('%s\n' % sequence)


