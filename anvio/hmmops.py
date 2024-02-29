# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    HMM related operations.
"""

import sys
import textwrap
from collections import Counter

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
P = terminal.pluralize


class SequencesForHMMHits:
    def __init__(self, contigs_db_path, sources=set([]), split_names_of_interest=set([]), init=True, run=run, progress=progress, bin_name=None):
        self.run = run
        self.progress = progress

        if not isinstance(sources, type(set([]))):
            raise ConfigError("'sources' variable has to be a set instance.")

        if not isinstance(split_names_of_interest, type(set([]))):
            raise ConfigError("'split_names_of_interest' variable has to be a set instance.")

        self.bin_name = bin_name # this is only relevant for some output messages
        self.sources = set([s for s in sources if s])
        self.hmm_hits = {}
        self.hmm_hits_info ={}
        self.hmm_hits_splits = {}
        self.contig_sequences = {}
        self.aa_sequences = {}
        self.genes_in_contigs = {}
        self.splits_in_contigs = {}

        if contigs_db_path:
            self.init_dicts(contigs_db_path, split_names_of_interest)
            self.initialized = True
        else:
            self.initialized = False


    def init_dicts(self, contigs_db_path, split_names_of_interest=set([])):
        """Initialize essential data for HMM stuff.

           This function will do its best to not load any data that will not
           be used later for HMM related operations. For instance, it will
           learn which gene caller ids are of interest based on HMM sources,
           and only recover data for splits and contigs based on that information,
           not accessing a large fraction of a given contigs database.
        """

        utils.is_contigs_db(contigs_db_path)
        contigs_db = db.DB(contigs_db_path, anvio.__contigs__version__, run=self.run, progress=self.progress)
        self.hmm_hits_info = contigs_db.get_table_as_dict(t.hmm_hits_info_table_name)

        missing_sources = [s for s in self.sources if s not in self.hmm_hits_info]
        if len(missing_sources):
            contigs_db.disconnect()
            progress.reset()
            raise ConfigError("Bad news, Houston :/ The contigs database '%s' is missing one or more HMM sources "
                              "that you wished it didn't: '%s'." % (contigs_db_path, ', '.join(missing_sources)))

        if not self.sources:
            self.sources = set(list(self.hmm_hits_info.keys()))

        if not self.sources:
            # there is nothing to initialize..
            return

        self.progress.new("Recovering sequences for HMM Hits")
        self.progress.update('...')

        # get data from HMM tables based on sources of interest
        self.progress.update('Getting data from HMM tables for %d source(s)' % len(self.sources))
        where_clause_for_sources = "source in (%s)" % ', '.join(['"%s"' % s for s in self.sources])
        self.hmm_hits = contigs_db.get_some_rows_from_table_as_dict(t.hmm_hits_table_name,
                                                                    where_clause=where_clause_for_sources,
                                                                    error_if_no_data=False)
        self.hmm_hits_splits = contigs_db.get_some_rows_from_table_as_dict(t.hmm_hits_splits_table_name,
                                                                    where_clause=where_clause_for_sources,
                                                                    error_if_no_data=False)

        # if the user sent a split names of interest, it means they are interested in hits that only occur
        # in a specific set of split names. NOTE: this really makes it very difficult to deal with HMM hits
        # that span through multiple splits. here we will mark such HMM hits in `hmm_hits_splits_entry_ids_to_remove`
        # for removal, but we will also keep track of those guys and let the user know what happened.
        if len(split_names_of_interest):
            total_num_split_names = len(set(self.hmm_hits_splits[entry_id]['split'] for entry_id in self.hmm_hits_splits))
            hmm_hits_splits_entry_ids_to_remove = set([])
            hmm_hits_entry_ids_to_remove = set([])
            hmm_hits_entry_ids_associated_with_fragmented_hmm_hits = set([])
            hmm_sources_associated_with_fragmented_hmm_hits = set([])
            for entry_id in self.hmm_hits_splits:
                if self.hmm_hits_splits[entry_id]['split'] not in split_names_of_interest:
                    hmm_hits_splits_entry_ids_to_remove.add(entry_id)
                    hmm_hits_entry_ids_to_remove.add(self.hmm_hits_splits[entry_id]['hmm_hit_entry_id'])

                    if not self.hmm_hits_splits[entry_id]['percentage_in_split'] == 100:
                        # this is important. if we are here, there is a bit more to do since it means that
                        # the split name associated with self.hmm_hits_splits[entry_id] is not in
                        # `split_names_of_interest`. but since the `percentage_in_split` for this HMM hit is NOT
                        # 100%, it could be the case that other splits that contain pieces of this HMM hit may
                        # still be in `split_names_of_interest`. But here we are setting the stage for this
                        # HMM hit to be removed from `self.hmm_hits` altogether. To make things right, we must
                        # go through `hmm_hits_splits`, and remove remaining entries there that is associated with
                        # this HMM hit later using the contents of this variable:
                        hmm_hits_entry_ids_associated_with_fragmented_hmm_hits.add(self.hmm_hits_splits[entry_id]['hmm_hit_entry_id'])
                        hmm_sources_associated_with_fragmented_hmm_hits.add(self.hmm_hits_splits[entry_id]['source'])

            if len(hmm_hits_entry_ids_associated_with_fragmented_hmm_hits):
                # if we are here, we will have to update `hmm_hits_splits_entry_ids_to_remove` carefully:
                additional_entry_ids_to_be_removed = set([e for e in self.hmm_hits_splits if self.hmm_hits_splits[e]['hmm_hit_entry_id'] in hmm_hits_entry_ids_associated_with_fragmented_hmm_hits])
                hmm_hits_splits_entry_ids_to_remove.update(additional_entry_ids_to_be_removed)

                # let's warn the user while we're at it so they panic, too.
                progress.reset()
                if self.bin_name:
                    header = f"A WARNING RELATED TO HMMs IN '{self.bin_name}'"
                else:
                    header = "WARNING"

                self.run.warning(f"While anvi'o was trying to finalize HMM hits associated with splits of interest, "
                                 f"it realized that there were one or more HMM hits that spanned through multiple splits "
                                 f"yet not all of those splits were among the splits of interest. This can happen if you "
                                 f"refined a contig by excluding some of its splits either manually during binning or "
                                 f"automatically during whatever black magic you were engaged in. Anvi'o does not judge. But "
                                 f"In these situations anvi'o excludes the entire HMM hit from being reported to be on the "
                                 f"safe side. Some HMM hits coming from {P('HMM source', len(hmm_sources_associated_with_fragmented_hmm_hits))} "
                                 f"(\"{', '.join(hmm_sources_associated_with_fragmented_hmm_hits)}\"), will not appear in your "
                                 f"downstream analyses (including in your `anvi-summarize` outputs). If you are really really "
                                 f"interested in those partial HMM hits and sequences associated with them, there are multiple "
                                 f"ways to recover them (one of the best way to do it involves the use of `anvi-split` and "
                                 f"re-running HMMs in your final bins). Please feel free to reach out to the anvi'o community "
                                 f"for ideas.", header=header)


            if len(hmm_hits_splits_entry_ids_to_remove):
                for entry_id in hmm_hits_splits_entry_ids_to_remove:
                    self.hmm_hits_splits.pop(entry_id)

            if len(hmm_hits_entry_ids_to_remove):
                for entry_id in hmm_hits_entry_ids_to_remove:
                    self.hmm_hits.pop(entry_id)

            filtered_num_split_names = len(set(self.hmm_hits_splits[entry_id]['split'] for entry_id in self.hmm_hits_splits))

            if anvio.DEBUG:
                self.progress.end()

                self.run.warning(None, header="SequencesForHMMHits info")
                self.run.info_single('%d split names of interest are found' % len(split_names_of_interest))
                self.run.info('Total split names w/HMM hits', total_num_split_names)
                self.run.info('Final split names w/HMM hits', filtered_num_split_names, nl_after=1)

                self.progress.new("Recovering sequences for HMM Hits")
                self.progress.update('...')

        if not len(self.hmm_hits):
            # there are HMMs but no hits. FINE.
            self.progress.end()
            contigs_db.disconnect()
            progress.reset()
            self.run.warning(f"SequencesForHMMHits class here. The current database (at {contigs_db_path}) "
                             f"contains 0 HMM hits, at least within the HMM sources or splits that were "
                             f"requested. It might not be a problem for your case, but we just thought you should "
                             f"know, in case it is. So there you have it.")
            return

        gene_caller_ids_of_interest = set([e['gene_callers_id'] for e in self.hmm_hits.values()])
        where_clause_for_genes = "gene_callers_id in (%s)" % ', '.join(['%d' % g for g in gene_caller_ids_of_interest])

        self.progress.update('Recovering split and contig names for %d genes' % (len(gene_caller_ids_of_interest)))
        split_names_of_interest, contig_names_of_interest = utils.get_split_and_contig_names_of_interest(contigs_db_path, gene_caller_ids_of_interest)

        self.progress.update('Recovering contig seqs for %d genes' % (len(gene_caller_ids_of_interest)))
        where_clause_for_contigs = "contig in (%s)" % ', '.join(['"%s"' % s for s in contig_names_of_interest])
        self.contig_sequences = contigs_db.get_some_rows_from_table_as_dict(t.contig_sequences_table_name, \
                                                                            string_the_key=True, \
                                                                            where_clause=where_clause_for_contigs)

        self.progress.update('Recovering amino acid seqs for %d genes' % (len(gene_caller_ids_of_interest)))
        self.aa_sequences = contigs_db.get_some_rows_from_table_as_dict(t.gene_amino_acid_sequences_table_name, \
                                                                        where_clause=where_clause_for_genes)

        self.genes_in_contigs = contigs_db.get_some_rows_from_table_as_dict(t.genes_in_contigs_table_name, \
                                                                            where_clause=where_clause_for_genes)

        self.splits_in_contigs = list(split_names_of_interest)

        self.progress.end()
        contigs_db.disconnect()


    def check_init(self):
        if not self.initialized:
            raise ConfigError("This SequencesForHMMHits instance cannot do what its asked for "
                              "because it is not yet initialized :/ The programmer could call "
                              "`init_dicts` first, but clearly they didn't care.")


    def list_available_gene_names(self, sources=[], dont_quit=False):
        self.check_init()

        if sources:
            if not isinstance(sources, list):
                raise ConfigError('HMM sources must be of type `list`.')

            missing_sources = [s for s in sources if s not in self.hmm_hits_info]
            if len(missing_sources):
                raise ConfigError("Some of the HMM sources you are looking for are not known "
                                  "to this database: %s." % (', '.join(['"%s"' % s for s in missing_sources])))
            hmm_sources = sources
        else:
            hmm_sources = self.sources

        for source in hmm_sources:
            t = self.hmm_hits_info[source]
            run.info_single('%s [type: %s]: %s' % (source, t['search_type'], ', '.join(sorted(t['genes'].split(',')))), nl_after = 2)

        if dont_quit:
            return

        sys.exit(0)


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


    def get_gene_hit_counts_per_hmm_source(self, sources=None, dont_include_models_with_multiple_domain_hits=False):
        """This function returns a 2-level dictionary mapping genes to their number of HMM hits.
        
        The outer dictionary is keyed by HMM source and the inner dictionary is keyed by gene name.

        PARAMETERS
        ==========
        sources : list of str
            A list of HMM sources to count hits for. If not provided, will use all possible sources in the HMM hits table
        dont_include_models_with_multiple_domain_hits : Boolean
            A flag variable to control whether we return counts for models that belong to multple HMM sources. This is 
            relevant to the NumGenomesEstimator class, in which case we need to remove any genes that have hits from multiple 
            single-copy core gene domains to avoid double-counting. See https://github.com/merenlab/anvio/issues/2231 for details
        """

        if not sources:
            sources = [source for source in self.hmm_hits_info]
        else:
            if not isinstance(sources, list):
                raise ConfigError("get_gene_hit_counts_per_hmm_source speaking: `sources` variable must be of type `list`.")

            missing_sources = [source for source in sources if source not in self.hmm_hits_info]
            if len(missing_sources):
                self.progress.end()
                raise ConfigError("Anvi'o was trying to generate information regarding all the hits per HMM source stored "
                                  "in its databases, but some of the sources you requested do not seem to be found anywhere :/ "
                                  "Here is the list of those that failed you: '%s'." % (','.join(sources)))

        hmm_hits = list(self.hmm_hits.values())
        models_to_remove = {}
        if dont_include_models_with_multiple_domain_hits:
            observed_gcids_to_sources = {}
            # first get list of HMM hits, HMM model names and sources for each gene callers id
            for index, entry in enumerate(hmm_hits):
                source = entry['source']
                gcid = entry['gene_callers_id']
                model = entry['gene_name']
                if source in sources:
                    if gcid in observed_gcids_to_sources:
                        observed_gcids_to_sources[gcid]['hit_indices'].append(index)
                        if source not in observed_gcids_to_sources[gcid]['hmm_models']:
                            observed_gcids_to_sources[gcid]['hmm_models'][source] = set([])
                        observed_gcids_to_sources[gcid]['hmm_models'][source].add(model)
                    else:
                        observed_gcids_to_sources[gcid] = {'hit_indices': [index],
                                                           'hmm_models': {source: set([model])}}

            # then we identify the models with multi-domain hits
            for g, info in observed_gcids_to_sources.items():
                s_list = info['hmm_models'].keys()
                if len(s_list) > 1:
                    for s, m_set in info['hmm_models'].items():
                        if s in models_to_remove:
                            models_to_remove[s].update(m_set)
                        else:
                            models_to_remove[s] = m_set
            
            # inform the user what is going on
            num_models_affected = 0
            num_model_strs = []
            for s, m_set in models_to_remove.items():
                num = len(m_set)
                num_models_affected += num
                num_model_strs.append(f"{num} from {s}")
            self.run.warning(f"Hello there from the SequencesForHMMHits.get_gene_hit_counts_per_hmm_source() function. "
                             f"Just so you know, someone asked for SCG HMMs that belong to multiple sources *not* to be "
                             f"counted, and this will result in {num_models_affected} models "
                             f"to be removed from our counts, more specifically: {', '.join(num_model_strs)}. You can "
                             f"run this program with the `--debug` flag if you want to see a list of the models that we "
                             f"will ignore from each HMM source.")
        
        gene_hit_counts = {}
        for source in sources:
            gene_hit_counts[source] = {}

            if anvio.DEBUG and models_to_remove[source]:
                self.run.info_single(f"Models to be ignored for source {source}: {', '.join(models_to_remove[source])}")

            for gene_name in self.hmm_hits_info[source]['genes'].split(','):
        for entry in hmm_hits:
            source    = entry['source']
            gene_name = entry['gene_name']

            if source in sources:
                gene_hit_counts[source][gene_name.strip()] += 1

        return gene_hit_counts


    def get_num_genomes_from_SCG_sources_dict(self):
        SCG_sources = [key for key in self.hmm_hits_info if self.hmm_hits_info[key]['search_type'] == 'singlecopy']

        if not len(SCG_sources):
            return {}

        gene_hit_counts_per_hmm_source = self.get_gene_hit_counts_per_hmm_source(SCG_sources, dont_include_models_with_multiple_domain_hits=True)

        num_genomes_per_SCG_source = {}
        for SCG_source in SCG_sources:
            l = list(gene_hit_counts_per_hmm_source[SCG_source].values())
            num_genomes_per_SCG_source[SCG_source] = {'num_genomes': int(stats.mode(l)[0]),
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

            This function will return DNA sequences by default. If `return_amino_acid_sequences` parameter
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


    def filter_hmm_sequences_dict_to_keep_only_unique_gene_hits(self, hmm_sequences_dict_for_splits):
        """This takes the output of `get_sequences_dict_for_hmm_hits_in_splits`, and goes through every hit
           to remove hits that resolve to the same gene, and only keeps the one with the smallest e-value.

           Say, if there are multipe HMMs hitting the same gene, this filter will only keep the one that is the
           most significant.
        """

        hits_per_gene_caller_id = {}
        for hit_id in hmm_sequences_dict_for_splits:
            hit = hmm_sequences_dict_for_splits[hit_id]

            if hit['gene_callers_id'] in hits_per_gene_caller_id:
                hits_per_gene_caller_id[hit['gene_callers_id']].add(hit_id)
            else:
                hits_per_gene_caller_id[hit['gene_callers_id']] = set([hit_id])

        for gene_callers_id in hits_per_gene_caller_id:
            if len(hits_per_gene_caller_id[gene_callers_id]) > 1:
                hit_ids_to_remove = [tpl[1] for tpl in sorted([(hmm_sequences_dict_for_splits[hit]['e_value'], hit) \
                                                                    for hit in hits_per_gene_caller_id[gene_callers_id]])[1:]]
                for hit_id_to_remove in hit_ids_to_remove:
                    hmm_sequences_dict_for_splits.pop(hit_id_to_remove)

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


    def filter_hmm_sequences_dict_for_genes_that_are_too_long(self, hmm_sequences_dict_for_splits, ignore_genes_longer_than=0):
        """This takes in your `hmm_sequences_dict_for_splits`, and removes genes that are too long"""

        if not isinstance(ignore_genes_longer_than, int):
            raise ConfigError("The `--ignore-genes-longer-than` expects an integer argument :/")

        # we will keep track of these bad bois
        gene_calls_removed = set([])

        # we identify entry ids that describe genes that are too long for removal
        entry_ids_to_remove = set([])
        for entry_id in hmm_sequences_dict_for_splits:
            if hmm_sequences_dict_for_splits[entry_id]['length'] > ignore_genes_longer_than:
                entry_ids_to_remove.add(entry_id)
                gene_calls_removed.add(hmm_sequences_dict_for_splits[entry_id]['gene_callers_id'])

        # we return early if there is nothing to be done
        if not len(entry_ids_to_remove):
            self.run.warning(f"You asked anvi'o to remove genes that are longer than {ignore_genes_longer_than} nts from your "
                             f"HMM hits. But none of the gene calls were longer than that value, so you get to keep everything.",
                             header="A MESSAGE FROM YOUR GENE LENGTH FILTER 📏")
            return (hmm_sequences_dict_for_splits, set([]), set([]))

        # if we are here, it means there are things to be gotten rid of. we will remove things
        # while keeping the user informed.
        self.run.warning(f"You asked anvi'o to remove genes that are longer than {ignore_genes_longer_than} nts from your "
                         f"HMM hits. There were a total of {len(entry_ids_to_remove)} HMM hits that matched to gene calls "
                         f"that were longer than {ignore_genes_longer_than}, and they are now removed from your analysis. "
                         f"The following lines list all these gene calls, their length, and which model they belonged.",
                         header="A MESSAGE FROM YOUR GENE LENGTH FILTER 📏", lc='yellow')

        # before we actually start removing stuff, we first learn all bin names that are in the master dict
        bin_names_in_original_dict = set([])
        for entry in hmm_sequences_dict_for_splits.values():
            bin_names_in_original_dict.add(entry['bin_id'])

        # puts in business socks
        for entry_id in entry_ids_to_remove:
            e = hmm_sequences_dict_for_splits[entry_id]
            self.run.info_single(f"Source: {e['source']} / Model: {e['gene_name']} /  Bin: {e['bin_id']} / "
                                 f"Gene call: {e['gene_callers_id']} / Length: {e['length']}",
                                 cut_after=None, level=2)
            hmm_sequences_dict_for_splits.pop(entry_id)

        # now we're done, and we will take another look at the dict to figure out remaining bins
        bin_names_in_filtered_dict = set([])
        for entry in hmm_sequences_dict_for_splits.values():
            bin_names_in_filtered_dict.add(entry['bin_id'])

        # bins we lost
        bins_removed = bin_names_in_original_dict.difference(bin_names_in_filtered_dict)

        return (hmm_sequences_dict_for_splits, gene_calls_removed, bins_removed)


    def filter_hmm_sequences_dict_from_genes_that_occur_in_less_than_N_bins(self, hmm_sequences_dict_for_splits, min_num_bins_gene_occurs=None):
        """This takes in your `hmm_sequences_dict_for_splits`, and removes genes that rarely occurs across bins.

           The `min_num_bins_gene_occurs` parameter defines what is the minimum number of bins you want a gene to
           be present. It removes all the genes that do not fit into that criterion."""

        if not isinstance(min_num_bins_gene_occurs, int):
            raise ConfigError("Funny. Someone called the function to filter gene names from HMM sequences dictionary if they occur in less than "
                              "a certain amount. But they didn't sen an integer for that amount :/")

        if min_num_bins_gene_occurs < 0:
            raise ConfigError("But the minimum number of bins a gene is expected to be found can't be a negative value now. Right? :/")

        all_bins = set([])

        for entry in hmm_sequences_dict_for_splits.values():
            all_bins.add(entry['bin_id'])

        if min_num_bins_gene_occurs > len(all_bins):
            raise ConfigError("You are asking anvi'o to remove any gene that occurs in less than %d genomes (or bins), however, it seems you have only "
                              "%s genomes. Either you set a parameter that exceeds the number of genomes you actually have, or the previous filters "
                              "applied to your set of genes have removed all genes from some or all of your genomes :/ Anvi'o cannot know here what might "
                              "have gone wrong, but it kinda believes that it is all on your at this point :/" % (min_num_bins_gene_occurs, len(all_bins)))

        gene_occurrences_accross_bins = self.get_gene_num_occurrences_across_bins(hmm_sequences_dict_for_splits)

        genes_to_remove = set([])
        all_genes = set(list(gene_occurrences_accross_bins.keys()))
        for gene_name in all_genes:
            if gene_occurrences_accross_bins[gene_name] < min_num_bins_gene_occurs:
                genes_to_remove.add(gene_name)

        genes_to_keep = all_genes.difference(genes_to_remove)

        self.run.info_single("Hi! The anvi'o function that was supposed to remove genes that were occurring in "
                             "less than X number of bins due to the use of `--min-num-bins-gene-occurs` is "
                             "speaking. What follows is a report of what happened after anvi'o tried to remove "
                             "genes that were occurring in at least %d of the %d bins you had at this point." \
                                    % (min_num_bins_gene_occurs, len(all_bins)), nl_before=1, nl_after=1)

        self.run.info('All genes (%d)' % len(all_genes), ', '.join(all_genes), nl_after=1)
        self.run.info('Genes occurred in at least %d of %d bins (%d)' % (min_num_bins_gene_occurs, len(all_bins), len(genes_to_keep)), ', '.join(genes_to_keep), nl_after=1, mc='green')
        self.run.info('Genes that are no more in the analysis (%d)' % (len(genes_to_remove)), ', '.join(genes_to_remove) if genes_to_remove else 'None.', nl_after=1, mc='red')

        if len(genes_to_remove):
            return (utils.get_filtered_dict(hmm_sequences_dict_for_splits, 'gene_name', genes_to_keep), genes_to_remove)
        else:
            return (hmm_sequences_dict_for_splits, set([]))


    def filter_hmm_sequences_dict_for_to_only_include_specific_genes(self, hmm_sequences_dict_for_splits, gene_names=[]):
        """This takes the dictionary for HMM hits, and removes all genes from it except the ones in `gene_names`.

        It is critical to keep in mind that the removal of genes can leave behind no gene at all for some of the
        genomes/bins. That's why this function tracks the genome names in the dictionary before and after to make
        sure it can report the loss of genomes for the user to consider.
        """

        # gather all bin names
        bin_names_in_original_dict = set([])
        for entry in hmm_sequences_dict_for_splits.values():
            bin_names_in_original_dict.add(entry['bin_id'])

        # filter out every gene hit except those in `gene_names`
        hmm_sequences_dict_for_splits = utils.get_filtered_dict(hmm_sequences_dict_for_splits, 'gene_name', set(gene_names))

        # gather remaining bin names in the dict
        bin_names_in_filtered_dict = set([])
        for entry in hmm_sequences_dict_for_splits.values():
            bin_names_in_filtered_dict.add(entry['bin_id'])

        bins_that_are_lost = bin_names_in_original_dict.difference(bin_names_in_filtered_dict)

        if not len(bins_that_are_lost):
            # well, all bins are still in the data structure. we're good to return everything
            return hmm_sequences_dict_for_splits, []

        # if we are still here, it means some bins were gon buh-bye. we start by letting
        # the user gently that stuff went south
        self.run.info_single("Yo yo yo! The anvi'o function that helps you focus only on a specific list of gene names "
                             "among your HMM hits is speaking (we are here most likely you used the --gene-names flag "
                             "to get rid of all the other genes in a given HMM collection). What follows is a report of "
                             "happened because ANVI'O ENDED UP LOSING SOME BINS/GENOMES FROM YOUR ANALYSIS AS THEY DID "
                             "NOT HAVE *ANY* OF THE GENES YOU LISTED (sorry for the CAPS lock here, but we wanted to "
                             "make sure you don't miss this, since this will certainly influence your downstream "
                             "analyses). If you want to keep more bins in your analysis, you can include more genes "
                             "in your `--gene-names` -- but of course it will not change the fact that some bins will "
                             "still be missing some genes, and how does this patchiness will impact your downstream "
                             "analyses (such as phylogenomics) is an important question that will require you to "
                             "consider. Pro tip: you can always use the program `anvi-script-gen-function-matrix-across-genomes` "
                             "to see the distribution of HMM hits across your bins/genomes.", nl_before=1, nl_after=1)

        self.run.info('Num bins at the beginning of this filter', len(bin_names_in_original_dict), nl_after=1)
        self.run.info(f'Num bins that lacked the {P("gene", len(gene_names))} in `--gene-names`', len(bins_that_are_lost), nl_after=1, mc='red')
        self.run.info('Bins that are no more in the analysis', ', '.join(bins_that_are_lost), nl_after=1, mc='red')

        return hmm_sequences_dict_for_splits, bins_that_are_lost


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

        self.run.info_single("Hi there! The anvi'o function that kills bins is speaking (we are here because you used "
                             "the --max-num-genes-missing-from-bin parameter to remove bins that are not good enough for "
                             "your analysis becasue they are missing lots of genes. What follows is a report of what "
                             "happened.", nl_before=1, nl_after=1)

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


    def __store_concatenated_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, partition_file_path=None, wrap=120, separator = 'XXX', genes_order=None, align_with=None, just_do_it=False):
        """Generates concatenated sequences from `hmm_sequences_dict_for_splits` dict.

           Please do NOT directly access to this function, and use `store_hmm_sequences_into_FASTA`
           instead.
        """

        if len(self.sources) != 1:
            if just_do_it:
                self.run.warning("You have asked anvi'o to not pay attention to the fact that you are asking for genes to be concatenated "
                                 "that are coming from different HMM collections. Fingers crossed. Please check the deflines of the "
                                 "resulting FASTA file carefully.")
            else:
                raise ConfigError("In theory you should be requesting a single HMM source if you want your genes to be concatenated. "
                                  "But in practice everyone has different needs, so we don't know. If this is not due to an error on "
                                  "your part, and if you think you know what you are doing, you can ask anvi'o to let you concatenate "
                                  "genes from multiple HMM sources by using the flag `--just-do-it`. In that case you will not see this "
                                  "error, but you must be extremely careful to make sure the resulting file looks like it should, and "
                                  "the information it contains makes sense. Since this not the common practice, you may run into other "
                                  "errors downstream, for which we apologize in advance.")

        # if the user did not define a single HMM source, then it will recover all genes in all HMM sources.
        gene_names_in_source = []
        for _hmm_source in self.sources:
            gene_names_in_source.extend([g.strip() for g in self.hmm_hits_info[_hmm_source]['genes'].split(',')])

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
                raise ConfigError("One or more gene names in the genes order list does seem to appear among the genes described "
                                  "by the HMM sources (which translates to 'terrible news'). Here are the genes that cause this "
                                  "issue if you want to fix this: '%s' (and here are the HMM sources you have been using for this "
                                  "operation in case it helps: '%s')." \
                                              % (', '.join(genes_in_genes_order_but_missing_in_hmm_source), ', '.join(self.sources)))
            gene_names = genes_order
        else:
            self.run.warning("You did not define any gene names. Bold move. Now anvi'o will attempt to report a file with all "
                             "genes defined in your HMM source(s). This will likely be quite ugly, so please brace yourself.")

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
            genes_in_bins_dict[gene_name] = aligner(run=terminal.Run(verbose=False)).run_default(genes_list, debug=anvio.DEBUG)
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
            run.warning("You asked for some genes that were missing from all bins this class had in the "
           "HMM hits dictionary (here is a list of them: '%s'). Not knowing what to do with this werid "
           "situation, anvi'o put gap characters for all of them and retained your order. Here are those "
           "genes that missed the party: '%s'" % \
                (', '.join(bin_names_in_dict), ', '.join(set(gene_names_missing_from_everywhere))))

        f.close()

        if partition_file_path:
            utils.gen_NEXUS_format_partition_file_for_phylogenomics(partition_file_path, [(g, gene_lengths[g]) for g in gene_names], separator, run=self.run, progress=self.progress)


    def __store_individual_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120, separator = 'XXX', genes_order=None, align_with=None):
        """Stores every sequence in hmm_sequences_dict_for_splits into the `output_file_path`.

           Please do NOT directly access to this function, and use `store_hmm_sequences_into_FASTA`
           instead.
        """

        # if the user wants alignment, lets update the input dictionary with aligned sequences
        if align_with:
            self.run.info('Sequence aligner', align_with)
            self.run.warning("Anvi'o will align your sequences since you explicitly asked for an aligner. However, you are not "
                             "concatenating your genes. If you are working with multiple gene names, your multiple sequence alignment "
                             "may contain genes that are evolutionarily very distant from each other, and the resulting multiple "
                             "sequence alignment may be irrelevant to answer any biologically relevant questions. So here anvi'o will "
                             "do what you want, assuming that you know what you are doing.")

            aligner = self.get_aligner(align_with)
            genes_list = [(gene_id, hmm_sequences_dict_for_splits[gene_id]['sequence']) for gene_id in hmm_sequences_dict_for_splits]

            self.progress.new('Alignment')
            self.progress.update('Working on %d sequences ...' % (len(genes_list)))
            genes_aligned = aligner(run=terminal.Run(verbose=False)).run_default(genes_list, debug=anvio.DEBUG)
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


    def store_hmm_sequences_into_FASTA(self, hmm_sequences_dict_for_splits, output_file_path, wrap=120, concatenate_genes=False, partition_file_path=None, separator=None, genes_order=None, align_with=None, just_do_it=False):
        """Stores HMM sequences into a FASTA file."""

        filesnpaths.is_output_file_writable(output_file_path)
        filesnpaths.is_output_file_writable(partition_file_path) if partition_file_path else None

        if wrap and not isinstance(wrap, int):
            raise ConfigError('"wrap" has to be an integer instance')

        if genes_order and concatenate_genes:
            gene_frequencies = Counter(genes_order)
            non_unique_genes = [g for g in gene_frequencies if gene_frequencies[g] > 1]
            if len(non_unique_genes):
                if just_do_it:
                    self.run.warning("Anvi'o found that some gene names occur multiple times (i.e., %s), but is letting this get away "
                                     "since the user invoked the grumpy flag." % (', '.join(non_unique_genes)), nl_before=1)
                else:
                    raise ConfigError("The list of gene names you wish to concatenate contains those that occur more than once. "
                                      "Here is the list: '%s'. While anvi'o believes it is a silly idea to have the same gene "
                                      "names multiple times, it will not care about it and will let you get away with it if you "
                                      "really want that. In which case you can use the flag `--just-do-it`, and move on with your "
                                      "very unconventional and cool analysis." % (', '.join(non_unique_genes)))

        if concatenate_genes:
            self.__store_concatenated_hmm_sequences_into_FASTA(hmm_sequences_dict_for_splits, output_file_path, partition_file_path, wrap, separator, genes_order, align_with, just_do_it)
        else:
            self.__store_individual_hmm_sequences_into_FASTA(hmm_sequences_dict_for_splits, output_file_path, wrap, separator, genes_order, align_with)


class NumGenomesEstimator(SequencesForHMMHits):
    """A simple interface to get estimated number of genomes from a contigs database.

    Notes for programmers
    =====================
    - Major changes in this class should consdier the use case example
      in `anvio/docs/artifacts/contigs-db.md`.
    """
    def __init__(self, contigs_db_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False)):
        self.run = run
        self.progress = progress
        self.contigs_db_path = contigs_db_path

        SequencesForHMMHits.__init__(self, self.contigs_db_path, run=self.run, progress=self.progress)

        self.estimates_dict = self.get_num_genomes_from_SCG_sources_dict()

        self.per_domain_totals = {}
        for entry in self.estimates_dict.values():
            if entry['domain'] not in self.per_domain_totals:
                self.per_domain_totals[entry['domain']] = 0

            self.per_domain_totals[entry['domain']] += entry['num_genomes']


    def num_genomes(self, for_domains=[]):
        """Higher-order estimate for the number of genomes

        Parameters
        ==========
        for_domains, list:
            The SCG domains you wish the get a tally for. By default, the function will return number of
            genomes found in all domains.

        Returns
        =======
            <num_genomes, domains>, tuple:
                This function will return a tuple with two items. The first is the total number of genomes
                estimated, and the second is the domains used for this estimate.
        """

        if for_domains:
            missing_domains = [d for d in for_domains if d not in self.per_domain_totals]
            if len(missing_domains):
                raise ConfigError(f"One of more of the domains you have requested do not seem to be among those that "
                                  f"anvi'o knows about \"{', '.join(missing_domains)}\". If what you are seeing is not "
                                  f"due to a typo, this can happen if you don't have the single-copy core gene set "
                                  f"to estimate number of genomes that belong to this domain. You perhaps do not have "
                                  f"them because you didn't run the program `anvi-run-hmms` on your contigs-db with all "
                                  f"the default HMM collections, or you are a dreamer who wishes some domain you are "
                                  f"working with is total news to anvi'o. If it is the latter, please get in touch with "
                                  f"us.")

        if not for_domains:
            for_domains = list(self.per_domain_totals.keys())

        return (sum([self.per_domain_totals[d] for d in for_domains]), for_domains)
