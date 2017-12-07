# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    HMM related operations.
"""

import textwrap
from scipy import stats

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
    def __init__(self, contigs_db_path, sources=set([]), init=True, run=run, progress=progress):
        self.run = run
        self.progress = progress

        if not isinstance(sources, type(set([]))):
            raise ConfigError("'sources' variable has to be a set instance.")

        self.sources = set([s for s in sources if s])
        self.hmm_hits = {}
        self.hmm_hits_info ={}
        self.hmm_hits_splits = {}
        self.contig_sequences = {}
        self.aa_sequences = {}
        self.genes_in_contigs = {}
        self.splits_in_contigs = {}

        if contigs_db_path:
            self.init_dicts(contigs_db_path)


    def init_dicts(self, contigs_db_path):
        # take care of contigs db related stuff and move on:
        contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__)
        self.hmm_hits = contigs_db.get_table_as_dict(t.hmm_hits_table_name)
        self.hmm_hits_info = contigs_db.get_table_as_dict(t.hmm_hits_info_table_name)
        self.hmm_hits_splits = contigs_db.get_table_as_dict(t.hmm_hits_splits_table_name)
        self.contig_sequences = contigs_db.get_table_as_dict(t.contig_sequences_table_name, string_the_key=True)
        self.aa_sequences = contigs_db.get_table_as_dict(t.gene_amino_acid_sequences_table_name)
        self.genes_in_contigs = contigs_db.get_table_as_dict(t.genes_in_contigs_table_name)
        self.splits_in_contigs = list(contigs_db.get_table_as_dict(t.splits_info_table_name).keys())
        contigs_db.disconnect()

        missing_sources = [s for s in self.sources if s not in self.hmm_hits_info]
        if len(missing_sources):
            raise ConfigError('Some of the requested sources were not found in the contigs database :/\
                                Here is a list of the ones that are missing: %s' % ', '.join(missing_sources))


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


    def get_gene_hit_counts_per_hmm_source(self, sources=None):
        if not sources:
            sources = [source for source in self.hmm_hits_info]
        else:
            if not isinstance(sources, list):
                raise ConfigError("get_gene_hit_counts_per_hmm_source speaking: `sources` variable must be of type `list`.")

            missing_sources = [source for source in sources if source not in self.hmm_hits_info]
            if len(missing_sources):
                self.progress.end()
                raise ConfigError("Anvi'o was trying to generate information regarding all the hits per HMM source stored\
                                   in its databases, but some of the sources you requested do not seem to be found anywhere :/\
                                   Here is the list of those that failed you: '%s'." % (','.join(sources)))

        gene_hit_counts = {}
        for source in sources:
            gene_hit_counts[source] = {}

            for gene_name in self.hmm_hits_info[source]['genes'].split(','):
                gene_hit_counts[source][gene_name.strip()] = 0

        for entry in list(self.hmm_hits.values()):
            source    = entry['source']
            gene_name = entry['gene_name']

            if source in sources:
                gene_hit_counts[source][gene_name.strip()] += 1

        return gene_hit_counts


    def get_num_genomes_from_SCG_sources_dict(self):
        SCG_sources = [key for key in self.hmm_hits_info if self.hmm_hits_info[key]['search_type'] == 'singlecopy']

        if not len(SCG_sources):
            return {}

        gene_hit_counts_per_hmm_source = self.get_gene_hit_counts_per_hmm_source(SCG_sources)

        num_genomes_per_SCG_source = {}
        for SCG_source in SCG_sources:
            l = list(gene_hit_counts_per_hmm_source[SCG_source].values())
            num_genomes_per_SCG_source[SCG_source] = {'num_genomes': int(stats.mode(l).mode[0]),
                                                      'domain': self.hmm_hits_info[SCG_source]['domain']}

        return num_genomes_per_SCG_source


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

        # trim hmm hits if sources
        if len(self.sources):
            self.hmm_hits_splits = utils.get_filtered_dict(self.hmm_hits_splits, 'source', self.sources)
            self.hmm_hits = utils.get_filtered_dict(self.hmm_hits, 'source', self.sources)
        else:
            self.sources = list(self.hmm_hits_info.keys())

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
            start, stop, forward = gene_call['start'], gene_call['stop'], gene_call['direction'] == 'f'

            if return_amino_acid_sequences:
                sequence = self.aa_sequences[gene_callers_id]['sequence']
            else:
                sequence = self.contig_sequences[contig_name]['sequence'][start:stop]
                if not forward:
                    sequence = utils.rev_comp(sequence)

            hmm_sequences_dict_for_splits[hit_unique_id] = {'sequence': sequence,
                                                            'source': source,
                                                            'bin_id': split_name_to_bin_id[split_name],
                                                            'gene_name': gene_name,
                                                            'e_value': e_value,
                                                            'contig': contig_name,
                                                            'start': start,
                                                            'stop': stop,
                                                            'gene_callers_id': gene_callers_id,
                                                            'rev_comped': (not forward),
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


    def get_gene_num_occurrences_across_bins(self, hmm_sequences_dict_for_splits):
        """Get a dictionary of gene names and their number of occurrences across all bins"""

        all_bins = set([])
        all_genes = set([])

        for entry in hmm_sequences_dict_for_splits.values():
            all_bins.add(entry['bin_id'])
            all_genes.add(entry['gene_name'])

        gene_num_occurrences_across_bins = dict([(gene_name, set([])) for gene_name in all_genes])

        for entry in hmm_sequences_dict_for_splits.values():
            gene_num_occurrences_across_bins[entry['gene_name']].add(entry['bin_id'])

        for gene_name in gene_num_occurrences_across_bins:
            gene_num_occurrences_across_bins[gene_name] = len(gene_num_occurrences_across_bins[gene_name])

        return gene_num_occurrences_across_bins


    def get_num_genes_missing_per_bin_dict(self, hmm_sequences_dict_for_splits, gene_names):
        """Get a dictionary of how many genes each bin is missing from a list of `gene_names`"""

        all_bins = set([])

        for entry in hmm_sequences_dict_for_splits.values():
            all_bins.add(entry['bin_id'])

        genes_in_bins_dict = dict([(bin_name, set([])) for bin_name in all_bins])
        num_genes_missing_per_bin = dict([(bin_name, 0) for bin_name in all_bins])

        for entry in hmm_sequences_dict_for_splits.values():
            genes_in_bins_dict[entry['bin_id']].add(entry['gene_name'])

        for bin_name in all_bins:
            for gene_name in gene_names:
                if gene_name not in genes_in_bins_dict[bin_name]:
                    num_genes_missing_per_bin[bin_name] += 1

        return num_genes_missing_per_bin


    def filter_hmm_sequences_dict_from_genes_that_occur_in_less_than_N_bins(self, hmm_sequences_dict_for_splits, min_num_bins_gene_occurs=None):
        """This takes in your `hmm_sequences_dict_for_splits`, and removes genes that rarely occurs across bins.

           The `min_num_bins_gene_occurs` parameter defines what is the minimum number of bins you want a gene to
           be present. It removes all the genes that do not fit into that criterion."""

        if not isinstance(min_num_bins_gene_occurs, int):
            raise ConfigError("Funny. Someone called the function to filter gene names from HMM sequences dictionary if they occur in less than\
                               a certain amount. But they didn't sen an integer for that amount :/")

        if min_num_bins_gene_occurs < 0:
            raise ConfigError("But the minimum number of bins a gene is expected to be found can't be a negative value now. Right? :/")

        all_bins = set([])

        for entry in hmm_sequences_dict_for_splits.values():
            all_bins.add(entry['bin_id'])

        if min_num_bins_gene_occurs > len(all_bins):
            raise ConfigError("OK. Well. This is awkward. You have like %d bins, eh? And you are asking anvi'o to remove any\
                               that occurs in less than %d bins. Do you see the problem here? Maybe it is time to take a break\
                               from work :(" % (len(all_bins), min_num_bins_gene_occurs))

        gene_occurrences_accross_bins = self.get_gene_num_occurrences_across_bins(hmm_sequences_dict_for_splits)

        genes_to_remove = set([])
        all_genes = set(list(gene_occurrences_accross_bins.keys()))
        for gene_name in all_genes:
            if gene_occurrences_accross_bins[gene_name] < min_num_bins_gene_occurs:
                genes_to_remove.add(gene_name)

        genes_to_keep = all_genes.difference(genes_to_remove)

        self.run.info_single("Hi! The anvi'o funciton that was supposed to remove genes that were occurring in\
                              less than X number of bins due to the use of `--min-num-bins-gene-occurs` is \
                              speaking. What follows is a report of what happened after anvi'o tried to remove\
                              genes that were occurring in at least %d of the %d bins you had at this point." \
                                    % (min_num_bins_gene_occurs, len(all_bins)), nl_before=1, nl_after=1)

        self.run.info('All genes (%d)' % len(all_genes), ', '.join(all_genes), nl_after=1)
        self.run.info('Genes occurred in at least %d of %d bins (%d)' % (min_num_bins_gene_occurs, len(all_bins), len(genes_to_keep)), ', '.join(genes_to_keep), nl_after=1, mc='green')
        self.run.info('Genes that are no more in the analysis (%d)' % (len(genes_to_remove)), ', '.join(genes_to_remove) if genes_to_remove else 'None.', nl_after=1, mc='red')

        if len(genes_to_remove):
            return (utils.get_filtered_dict(hmm_sequences_dict_for_splits, 'gene_name', genes_to_keep), genes_to_remove)
        else:
            return (hmm_sequences_dict_for_splits, set([]))


    def filter_hmm_sequences_dict_for_bins_that_lack_more_than_N_genes(self, hmm_sequences_dict_for_splits, gene_names, max_num_genes_missing=0):
        """This takes the output of `get_sequences_dict_for_hmm_hits_in_splits`, and goes through every bin\
           to identify bins or genomes that have lack more than `max_num_genes_missing` from a list of genes.

           Note that it returns a filtered dictionary, AND the bins that are removed."""

        num_genes_missing_per_bin = self.get_num_genes_missing_per_bin_dict(hmm_sequences_dict_for_splits, gene_names)

        bins_to_remove = set([])
        all_bins = set(list(num_genes_missing_per_bin.keys()))
        for bin_name in num_genes_missing_per_bin:
            if num_genes_missing_per_bin[bin_name] > max_num_genes_missing:
                bins_to_remove.add(bin_name)

        bins_to_keep = all_bins.difference(bins_to_remove)

        self.run.info_single("Hi there! The anvi'o function that kills bins is speaking (we are here because you used\
                              the --max-num-genes-missing-from-bin parameter to remove bins that are not good enough for\
                              your analysis becasue they are missing lots of genes. What follows is a report of what \
                              happened.", nl_before=1, nl_after=1)

        self.run.info('All bins (%d)' % len(all_bins), ', '.join(all_bins), nl_after=1)
        self.run.info('Bins that missed at most %d of %d genes (%d)' % (max_num_genes_missing, len(gene_names), len(bins_to_keep)), ', '.join(bins_to_keep), nl_after=1, mc='green')
        self.run.info('Bins that are no more in the analysis (%d)' % (len(bins_to_remove)), ', '.join(bins_to_remove) if bins_to_remove else 'None. Lovely.', nl_after=1, mc='red')


        if len(bins_to_remove):
            return (utils.get_filtered_dict(hmm_sequences_dict_for_splits, 'bin_id', bins_to_keep), bins_to_remove)
        else:
            return (hmm_sequences_dict_for_splits, set([]))


    def get_FASTA_header_and_sequence_for_gene_unique_id(self, hmm_sequences_dict_for_splits, gene_unique_id):
        entry = hmm_sequences_dict_for_splits[gene_unique_id]
        header = '%s___%s ' % (entry['gene_name'], gene_unique_id) + '|'.join(['%s:%s' % (k, str(entry[k])) for k in ['bin_id', 'source', 'e_value', 'contig', 'gene_callers_id', 'start', 'stop', 'length']])
        sequence = hmm_sequences_dict_for_splits[gene_unique_id]['sequence']
        return (header, sequence)


    def get_aligner(self, align_with=None):
        """Return an instance of an aligner"""

        from anvio.drivers import Aligners

        return Aligners().select(align_with)


    def __store_concatenated_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120, concatenate_genes=False, separator = 'XXX', genes_order=None, align_with=None):
        """Generates concatenated sequences from `hmm_sequences_dict_for_splits` dict.

           Please do NOT directly access to this function, and use `store_hmm_sequences_into_FASTA`
           instead.
        """

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

        # let's get an instance of the aligner early on so we learn about issues before its too late.
        aligner = self.get_aligner(align_with)

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
            self.run.warning("You did not define any gene names. Bold move. Now anvi'o will attempt to report a file with all\
                              genes defined in the HMM source '%s'." % hmm_source)

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
        self.progress.new('Aligning homolog gene sequences pre-concatenation')
        all_gene_names = list(genes_in_bins_dict.keys())
        num_genes = len(all_gene_names)
        for i in range(0, num_genes):
            gene_name = all_gene_names[i]
            self.progress.update('working on %s (%d of %d) ...' % (gene_name, i + 1, num_genes))
            genes_list = [(bin_name, genes_in_bins_dict[gene_name][bin_name]) \
                                                        for bin_name in genes_in_bins_dict[gene_name] \
                                                                           if bin_name in genes_in_bins_dict[gene_name]]
            genes_in_bins_dict[gene_name] = aligner(run=terminal.Run(verbose=False)).run_stdin(genes_list)
            gene_lengths[gene_name] = len(list(genes_in_bins_dict[gene_name].values())[0])
        self.progress.end()

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

            f.write('>%s num_genes:%d|genes:%s|separator:%s\n' % (bin_name, len(gene_names), ','.join(gene_names), separator))
            f.write('%s\n' % sequence)

        if len(gene_names_missing_from_everywhere):
            run.warning("You asked for some genes that were missing from all bins this class had in the\
            HMM hits dictionary (here is a list of them: '%s'). Not knowing what to do with this werid\
            situation, anvi'o put gap characters for all of them and retained your order. Here are those\
            genes that missed the party: '%s'" % \
                (', '.join(bin_names_in_dict), ', '.join(gene_names_missing_from_everywhere)))

        f.close()


    def __store_individual_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120, concatenate_genes=False, separator = 'XXX', genes_order=None, align_with=None):
        """Stores every sequence in hmm_sequences_dict_for_splits into the `output_file_path`.

           Please do NOT directly access to this function, and use `store_hmm_sequences_into_FASTA`
           instead.
        """

        # if the user wants alignment, lets update the input dictionary with aligned sequences
        if align_with:
            self.run.info('Sequence aligner', align_with)
            self.run.warning("Anvi'o will align your sequences since you explicitly asked for an aligner. However, you are not\
                              concatenating your genes. If you are working with multiple gene names, your multiple sequence alignment\
                              may contain genes that are evolutionarily very distant from each other, and the resulting multiple\
                              sequence alignment may be irrelevant to answer any biologically relevant questions. So here anvi'o will\
                              do what you want, assuming that you know what you are doing.")

            aligner = self.get_aligner(align_with)
            genes_list = [(gene_id, hmm_sequences_dict_for_splits[gene_id]['sequence']) for gene_id in hmm_sequences_dict_for_splits]

            self.progress.new('Alignment')
            self.progress.update('Working on %d sequences ...' % (len(genes_list)))
            genes_aligned = aligner(run=terminal.Run(verbose=False)).run_stdin(genes_list)
            self.progress.end()

            for gene_id in genes_aligned:
                hmm_sequences_dict_for_splits[gene_id]['sequence'] = genes_aligned[gene_id]

        f = open(output_file_path, 'w')

        for gene_unique_id in hmm_sequences_dict_for_splits:
            header, sequence = self.get_FASTA_header_and_sequence_for_gene_unique_id(hmm_sequences_dict_for_splits, gene_unique_id)

            if wrap:
                sequence = textwrap.fill(sequence, wrap, break_on_hyphens=False)

            f.write('>%s\n' % header)
            f.write('%s\n' % sequence)

        f.close()


    def store_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120, concatenate_genes=False, separator=None, genes_order=None, align_with=None):
        """Stores HMM sequences into a FASTA file."""

        filesnpaths.is_output_file_writable(output_file_path)

        if not isinstance(wrap, int):
            raise ConfigError('"wrap" has to be an integer instance')

        if concatenate_genes:
            self.__store_concatenated_hmm_sequences_into_FASTA(hmm_sequences_dict_for_splits, output_file_path, wrap, concatenate_genes, separator, genes_order, align_with)
        else:
            self.__store_individual_hmm_sequences_into_FASTA(hmm_sequences_dict_for_splits, output_file_path, wrap, concatenate_genes, separator, genes_order, align_with)
