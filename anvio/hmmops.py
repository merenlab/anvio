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
from anvio.drivers.muscle import Muscle

run = terminal.Run()
progress = terminal.Progress()


class SequencesForHMMHits:
    def __init__(self, contigs_db_path, sources=set([]), run=run, progress=progress):
        if not isinstance(sources, type(set([]))):
            raise ConfigError("'sources' variable has to be a set instance.")

        self.sources = set([s for s in sources if s])

        # take care of contigs db related stuff and move on:
        contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)
        self.hmm_hits = contigs_db.get_table_as_dict(t.hmm_hits_table_name)
        self.hmm_hits_info = contigs_db.get_table_as_dict(t.hmm_hits_info_table_name)
        self.hmm_hits_splits = contigs_db.get_table_as_dict(t.hmm_hits_splits_table_name)
        self.contig_sequences = contigs_db.get_table_as_dict(t.contig_sequences_table_name, string_the_key=True)
        self.aa_sequences = contigs_db.get_table_as_dict(t.gene_protein_sequences_table_name)
        self.genes_in_contigs = contigs_db.get_table_as_dict(t.genes_in_contigs_table_name)
        contigs_db.disconnect()

        missing_sources = [s for s in self.sources if s not in self.hmm_hits_info]
        if len(missing_sources):
            raise ConfigError('Some of the requested sources were not found in the contigs database :/\
                                Here is a list of the ones that are missing: %s' % ', '.join(missing_sources))

        if len(self.sources):
            self.hmm_hits_splits = utils.get_filtered_dict(self.hmm_hits_splits, 'source', self.sources)
            self.hmm_hits = utils.get_filtered_dict(self.hmm_hits, 'source', self.sources)
        else:
            self.sources = list(self.hmm_hits_info.keys())


    def get_hmm_hits_in_splits(self, splits_dict):
        split_names = set([])
        for s in list(splits_dict.values()):
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
        for bin_name in list(splits_dict.keys()):
            hmm_hits_per_bin[bin_name] = {}

        unique_ids_taken_care_of = set([])
        for split_entry in list(hits_in_splits.values()):
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


    def get_sequences_dict_for_hmm_hits_in_splits(self, splits_dict, return_amino_acid_sequences=False, return_best_hits=False):
        """splits dict is what you get from ccollections.GetSplitNamesInBins(args).get_dict(), and
           its struture goes like this:

                {
                    'bin_x': set['split_a, split_b, ...'],
                    'bin_y': set['split_c, split_d, ...'],
                    ...
                }

            This function will return DNA seqeunces by default. If `return_amino_acid_sequences` parameter
            is True, it will return AA sequences instead.

            `return_best_hit=True` will filter the resulting dictionary to remove weak hits if there are more
            than one hit for a given gene name in a bin for a given hmm source.
        """

        hits_in_splits, split_name_to_bin_id = self.get_hmm_hits_in_splits(splits_dict)

        hmm_sequences_dict_for_splits = {}

        unique_hits_taken_care_of = set([])
        for split_entry in list(hits_in_splits.values()):
            hmm_hit = self.hmm_hits[split_entry['hmm_hit_entry_id']]

            split_name = split_entry['split']
            source = hmm_hit['source']
            gene_name = hmm_hit['gene_name']
            e_value = hmm_hit['e_value']
            hit_unique_id = '___'.join([source, hmm_hit['gene_unique_identifier']])

            if hit_unique_id in unique_hits_taken_care_of:
                continue
            else:
                unique_hits_taken_care_of.add(hit_unique_id)

            gene_callers_id = hmm_hit['gene_callers_id']
            gene_call = self.genes_in_contigs[gene_callers_id]

            contig_name = gene_call['contig']
            start, stop = gene_call['start'], gene_call['stop']

            if return_amino_acid_sequences:
                sequence = self.aa_sequences[gene_callers_id]['sequence']
            else:
                sequence = self.contig_sequences[contig_name]['sequence'][start:stop]

            hmm_sequences_dict_for_splits[hit_unique_id] = {'sequence': sequence,
                                                            'source': source,
                                                            'bin_id': split_name_to_bin_id[split_name],
                                                            'gene_name': gene_name,
                                                            'e_value': e_value,
                                                            'contig': contig_name,
                                                            'start': start,
                                                            'stop': stop,
                                                            'gene_callers_id': gene_callers_id,
                                                            'length': stop - start}

        if return_best_hits:
            return self.filter_hmm_sequences_dict_for_splits_to_keep_only_best_hits(hmm_sequences_dict_for_splits)
        else:
            return hmm_sequences_dict_for_splits


    def filter_hmm_sequences_dict_for_splits_to_keep_only_best_hits(self, hmm_sequences_dict_for_splits):
        """This takes the output of `get_sequences_dict_for_hmm_hits_in_splits`, and goes through every hit\
           to identify for each bin_id hits with the same gene name and source. If there are multiple gene\
           names and source, removes every other except the one with the smallest e-value.

           Say, if there are multiple RecA hits in a bin based on Campbell_et_al, this will keep only the most\
           significant one.
        """

        bin_names, hmm_sources = set([]), set([])
        for v in list(hmm_sequences_dict_for_splits.values()):
            bin_names.add(v['bin_id'])
            hmm_sources.add(v['source'])

        # this dictionary will keep track of the occurrence of each gene name in each hmm source and bin:
        d = {}

        for bin_name in bin_names:
            d[bin_name] = {}
            for hmm_source in hmm_sources:
                d[bin_name][hmm_source] = {}

        # fill in gene_names
        for v in list(hmm_sequences_dict_for_splits.values()):
            d[v['bin_id']][v['source']][v['gene_name']] = []

        # add genes into lists
        for hit_unique_id in hmm_sequences_dict_for_splits:
            h = hmm_sequences_dict_for_splits[hit_unique_id]
            d[h['bin_id']][h['source']][h['gene_name']].append((h['e_value'], hit_unique_id), )

        # find the ones that occur twice:
        hit_unique_ids_to_remove = set([])
        for bin_name in d:
            for hmm_source in d[bin_name]:
                for gene_name in d[bin_name][hmm_source]:
                    if len(d[bin_name][hmm_source][gene_name]) > 1:
                        # here `d[bin_name][hmm_source][gene_name]` should look like this for a single gene name from
                        # a single source in a single bin:
                        #
                        #   [(6.2e-19, 'Rinke_et_al___7f25'), (1.8e-26, 'Rinke_et_al___57b0'), (7.7e-30, 'Rinke_et_al___4e43')]
                        #
                        # so we will sort from small e_value to big, and add unique ids to the list of shit ids to
                        # remove them from the dictionary we got.
                        hit_unique_ids_to_remove.update([t[1] for t in sorted(d[bin_name][hmm_source][gene_name])[1:]])

        for hit_unique_id in hit_unique_ids_to_remove:
            hmm_sequences_dict_for_splits.pop(hit_unique_id)


        return hmm_sequences_dict_for_splits


    def get_FASTA_header_and_sequence_for_gene_unique_id(self, hmm_sequences_dict_for_splits, gene_unique_id):
        entry = hmm_sequences_dict_for_splits[gene_unique_id]
        header = '%s___%s|' % (entry['gene_name'], gene_unique_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['bin_id', 'source', 'e_value', 'contig', 'gene_callers_id', 'start', 'stop', 'length']])
        sequence = hmm_sequences_dict_for_splits[gene_unique_id]['sequence']
        return (header, sequence)


    def store_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120, concatenate_genes=False, separator = 'XXX', genes_order=None):
        filesnpaths.is_output_file_writable(output_file_path)

        if not isinstance(wrap, int):
            raise ConfigError('"wrap" has to be an integer instance')

        if concatenate_genes:
            if len(self.sources) != 1:
                raise ConfigError("If you want your genes to be concatenated, you should be requesting a single HMM source. Why?\
                                   In fact we are not exactly sure why. But when we think of it, we couldn't come up with a \
                                   scenario where the user might truly be interested in concatenating genes from multiple HMM\
                                   sources, and we wanted to add a control in case they are making a mistake w/o realizing. If you\
                                   are sure this is what you must do for the question you are interested in, please send an\
                                   e-mail to the anvi'o discussion group, and convince us .. or you can just delete this if block\
                                   to avoid this check if you are not in the mood. We know the feeling.")

            hmm_source = self.sources.pop()
            gene_names_in_source = [g.strip() for g in self.hmm_hits_info[hmm_source]['genes'].split(',')]

            # the user wants to play rough. FINE. we will concatenate genes for phylogenomic analyses.
            gene_names = None

            # lets learn about what we have in this dictionary first.
            bin_names_in_dict = list(set([x['bin_id'] for x in hmm_sequences_dict_for_splits.values()]))
            gene_names_in_dict = sorted(list(set([x['gene_name'] for x in hmm_sequences_dict_for_splits.values()])))

            # if the function is called with a particular set and order of genes, use those, otherwise
            # stick with the gene names / order we found in the dictionary.
            if genes_order:
                genes_in_genes_order_but_missing_in_hmm_source = [g for g in genes_order if g not in gene_names_in_source]
                if len(genes_in_genes_order_but_missing_in_hmm_source):
                    raise ConfigError("One or more gene names in the genes order list does seem to appear among the genes described\
                                       by the HMM source %s (which translates to 'terrible news'). Here are the genes that cause this\
                                       issue if you want to fix this: '%s'" \
                                                  % (hmm_source, ', '.join(genes_in_genes_order_but_missing_in_hmm_source)))
                gene_names = genes_order
            else:
                gene_names = gene_names_in_dict

            # gene lenghts are especially important to accommodate missing genes with proper number of
            # gap characters
            gene_lengths = {}

            # buld a simpler dict that keeps genes sequences for each bin for a given gene name
            genes_in_bins_dict = {}
            for entry in hmm_sequences_dict_for_splits.values():
                gene_name = entry['gene_name']
                bin_name = entry['bin_id']
                sequence = entry['sequence']
                if gene_name in genes_in_bins_dict:
                    genes_in_bins_dict[gene_name][bin_name] = sequence
                else:
                    genes_in_bins_dict[gene_name] = {bin_name: sequence}


             # align homolog sequences across bins
            m = Muscle(run=terminal.Run(verbose=False))
            for gene_name in genes_in_bins_dict:
                genes_list = [(bin_name, genes_in_bins_dict[gene_name][bin_name]) \
                                                            for bin_name in genes_in_bins_dict[gene_name] \
                                                                               if bin_name in genes_in_bins_dict[gene_name]]
                genes_in_bins_dict[gene_name] = m.run_muscle_stdin(genes_list)
                gene_lengths[gene_name] = len(list(genes_in_bins_dict[gene_name].values())[0])

            # concatenate all of them and write them in a file
            f = open(output_file_path, 'w')
            gene_names_missing_from_everywhere = []
            for bin_name in bin_names_in_dict:
                sequences_list = []

                for gene_name in gene_names:
                    if gene_name in genes_in_bins_dict:
                        if bin_name in genes_in_bins_dict[gene_name]:
                            sequences_list.append(genes_in_bins_dict[gene_name][bin_name])
                        else:
                            sequences_list.append('-' * gene_lengths[gene_name])
                    else:
                        # if we are here, it means this is a gene that has been missing form the hmm hits dict, since it
                        # was not in any of the bins the dict described, but the user requested to have it in the
                        # alignment anyway. This can happen when the user wants to concatanate genes from one or more
                        # low-completion bins. We will keep track of them, and tell the user.
                        sequences_list.append('-' * 42)
                        gene_names_missing_from_everywhere.append(gene_name)

                sequence = separator.join(sequences_list)

                if wrap:
                    sequence = textwrap.fill(sequence, wrap, break_on_hyphens=False)

                f.write('>%s|genes:%s|separator:%s\n' % (bin_name, ','.join(gene_names), separator))
                f.write('%s\n' % sequence)

            if len(gene_names_missing_from_everywhere):
                run.warning("You asked for some genes that were missing from all bins this class had in the\
                HMM hits dictionary (here is a list of them: '%s'). Not knowing what to do with this werid\
                situation, anvi'o put gap characters for all of them and retained your order. Here are those\
                genes that missed the party: '%s'" % \
                    (', '.join(bin_names_in_dict), ', '.join(gene_names_missing_from_everywhere)))

            f.close()
        else:
            f = open(output_file_path, 'w')

            for gene_unique_id in hmm_sequences_dict_for_splits:
                header, sequence = self.get_FASTA_header_and_sequence_for_gene_unique_id(hmm_sequences_dict_for_splits, gene_unique_id)

                if wrap:
                    sequence = textwrap.fill(sequence, wrap, break_on_hyphens=False)

                f.write('>%s\n' % header)
                f.write('%s\n' % sequence)

            f.close()
