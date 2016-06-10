# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for HMM related operations.

    * HMMSearch takes care of searches using HMM profiles. It simply takes genes.txt and
    genes.hmm.gz files as input and returns a dictionary back with results. See anvio/data/hmm
    directory for examples.
"""

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


class SequencesForHMMHits:
    def __init__(self, contigs_db_path, sources=set([]), run=run, progress=progress):
        if not isinstance(sources, type(set([]))):
            raise ConfigError, "'sources' variable has to be a set instance."

        self.sources = set([s for s in sources if s])

        # take care of contigs db related stuff and move on:
        contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)
        self.hmm_hits = contigs_db.get_table_as_dict(t.hmm_hits_table_name)
        self.hmm_hits_info = contigs_db.get_table_as_dict(t.hmm_hits_info_table_name)
        self.hmm_hits_splits = contigs_db.get_table_as_dict(t.hmm_hits_splits_table_name)
        self.contig_sequences = contigs_db.get_table_as_dict(t.contig_sequences_table_name, string_the_key=True)
        self.genes_in_contigs = contigs_db.get_table_as_dict(t.genes_in_contigs_table_name)
        contigs_db.disconnect()

        missing_sources = [s for s in self.sources if s not in self.hmm_hits_info]
        if len(missing_sources):
            raise ConfigError, 'Some of the requested sources were not found in the contigs database :/\
                                Here is a list of the ones that are missing: %s' % ', '.join(missing_sources)

        if len(self.sources):
            self.hmm_hits_splits = utils.get_filtered_dict(self.hmm_hits_splits, 'source', self.sources)
            self.hmm_hits = utils.get_filtered_dict(self.hmm_hits, 'source', self.sources)
        else:
            self.sources = self.hmm_hits_info.keys()


    def get_hmm_hits_in_splits(self, splits_dict):
        split_names = set([])
        for s in splits_dict.values():
            split_names.update(s)

        hits_in_splits = utils.get_filtered_dict(self.hmm_hits_splits, 'split', split_names)

        split_name_to_bin_id = {}
        for bin_id in splits_dict:
            for split_name in splits_dict[bin_id]:
                split_name_to_bin_id[split_name] = bin_id

        return hits_in_splits, split_name_to_bin_id


    def get_hmm_hits_per_bin(self, splits_dict, source):
        hits_in_splits, split_name_to_bin_id = self.get_hmm_hits_in_splits(splits_dict)

        hmm_hits_per_bin = {}
        for bin_name in splits_dict.keys():
            hmm_hits_per_bin[bin_name] = {}

        unique_ids_taken_care_of = set([])
        for split_entry in hits_in_splits.values():
            hmm_hit = self.hmm_hits[split_entry['hmm_hit_entry_id']]

            hit_source = hmm_hit['source']

            if source and hit_source != source:
                continue

            split_name = split_entry['split']
            gene_name = hmm_hit['gene_name']
            gene_unique_id = hmm_hit['gene_unique_identifier']

            if gene_unique_id in unique_ids_taken_care_of:
                continue
            else:
                unique_ids_taken_care_of.add(gene_unique_id)

            bin_id = split_name_to_bin_id[split_name]

            if gene_name not in hmm_hits_per_bin[bin_id]:
                hmm_hits_per_bin[bin_id][gene_name] = 1
            else:
                hmm_hits_per_bin[bin_id][gene_name] += 1

        return hmm_hits_per_bin


    def get_hmm_sequences_dict_for_splits(self, splits_dict):
        """splits dict is what you get from ccollections.GetSplitNamesInBins(args).get_dict(), and
           its struture goes like this:

                {
                    'bin_x': set['split_a, split_b, ...'],
                    'bin_y': set['split_c, split_d, ...'],
                    ...
                }
        """

        hits_in_splits, split_name_to_bin_id = self.get_hmm_hits_in_splits(splits_dict)

        hmm_sequences_dict_for_splits = {}

        unique_ids_taken_care_of = set([])
        for split_entry in hits_in_splits.values():
            hmm_hit = self.hmm_hits[split_entry['hmm_hit_entry_id']]

            split_name = split_entry['split']
            source = hmm_hit['source']
            gene_name = hmm_hit['gene_name']
            e_value = hmm_hit['e_value']
            gene_unique_id = hmm_hit['gene_unique_identifier']

            if gene_unique_id in unique_ids_taken_care_of:
                continue
            else:
                unique_ids_taken_care_of.add(gene_unique_id)

            gene_call = self.genes_in_contigs[hmm_hit['gene_callers_id']]

            contig_name = gene_call['contig']
            start, stop = gene_call['start'], gene_call['stop']
            sequence = self.contig_sequences[contig_name]['sequence'][start:stop]

            hmm_sequences_dict_for_splits[gene_unique_id] = {'sequence': sequence,
                                                             'source': source,
                                                             'bin_id': split_name_to_bin_id[split_name],
                                                             'gene_name': gene_name,
                                                             'e_value': e_value,
                                                             'contig': contig_name,
                                                             'start': start,
                                                             'stop': stop,
                                                             'length': stop - start}

        return hmm_sequences_dict_for_splits


    def get_FASTA_header_and_sequence_for_gene_unique_id(self, hmm_sequences_dict_for_splits, gene_unique_id):
        entry = hmm_sequences_dict_for_splits[gene_unique_id]
        header = '%s___%s|' % (entry['gene_name'], gene_unique_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['bin_id', 'source', 'e_value', 'contig', 'start', 'stop', 'length']])
        sequence = hmm_sequences_dict_for_splits[gene_unique_id]['sequence']
        return (header, sequence)


    def store_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120):
        filesnpaths.is_output_file_writable(output_file_path)

        if not isinstance(wrap, int):
            raise ConfigError, '"wrap" has to be an integer instance'

        f = open(output_file_path, 'w')

        for gene_unique_id in hmm_sequences_dict_for_splits:
            header, sequence = self.get_FASTA_header_and_sequence_for_gene_unique_id(hmm_sequences_dict_for_splits, gene_unique_id)

            if wrap:
                sequence = textwrap.fill(sequence, wrap, break_on_hyphens=False)

            f.write('>%s\n' % header)
            f.write('%s\n' % sequence)
