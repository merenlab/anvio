# coding: utf-8
"""Summarizes information for a collection."""

import os
import sys
import numpy
import textwrap

from collections import Counter

import anvio
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.errors import ConfigError
from anvio.dbops import DatabasesMetaclass
from anvio.hmmops import SequencesForHMMHits
from anvio.summaryhtml import SummaryHTMLOutput, humanize_n, pretty


__author__ = "A. Murat Eren"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


pp = terminal.pretty_print
run = terminal.Run()
progress = terminal.Progress()
P = lambda x, y: float(x) * 100 / float(y)


class Summarizer(DatabasesMetaclass):
    """Creates an über dictionary of 'summary'."""
    def __init__(self, args = None, r = run, p = progress):
        self.summary = {}

        self.debug = False
        self.profile_db_path = None
        self.contigs_db_path = None
        self.output_directory = None
        self.split_names_per_bin = None
        self.completeness_data_available = False
        self.gene_coverages_data_available = False
        self.non_single_copy_gene_hmm_data_available = False

        self.run = r
        self.progress = p

        DatabasesMetaclass.__init__(self, args, self.run, self.progress)

        # databases initiated, let's make sure we have gene covereges data avaialable.
        if self.gene_coverages_dict:
            self.gene_coverages_data_available = True

        self.collections = ccollections.Collections()
        self.collections.populate_sources_dict(self.contigs_db_path, anvio.__contigs__version__)
        self.collections.populate_sources_dict(self.profile_db_path, anvio.__profile__version__)

        self.collection_id = None

        if args:
            if args.list_collections:
                self.collections.list_collections()
                sys.exit()

            self.collection_id = args.collection_id
            self.output_directory = args.output_dir
            self.debug = args.debug

        self.sanity_check()

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists = True)


    def sanity_check(self):
        if not self.collection_id:
            raise ConfigError, "You must specify a collection id :/"

        if self.collection_id not in self.collections.sources_dict:
            raise ConfigError, "%s is not a valid collection ID. See a list of available ones with '--list-collections' flag" % self.collection_id

        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists = True)


    def process(self):
        # learn who you are:
        collection_dict = self.collections.get_collection_dict(self.collection_id)
        collection_colors = self.collections.get_collection_colors(self.collection_id)

        # init profile data for colletion.
        self.init_collection_profile(collection_dict)

        # load completeness information if available
        self.completeness = completeness.Completeness(self.contigs_db_path)
        if len(self.completeness.sources):
            self.completeness_data_available = True

        # load HMM sources for non-single-copy genes if available
        if self.non_singlecopy_gene_hmm_sources:
            self.init_non_singlecopy_gene_hmm_sources()
            self.non_single_copy_gene_hmm_data_available = True

        # set up the initial summary dictionary
        self.summary['meta'] = {'output_directory': self.output_directory,
                                'collection': collection_dict.keys(),
                                'num_bins': len(collection_dict.keys()),
                                'collection_id': self.collection_id,
                                'total_nts_in_collection': 0,
                                'num_contigs_in_collection': 0,
                                'anvio_version': __version__, 
                                'profile': self.p_meta,
                                'contigs': self.a_meta,
                                'gene_coverages_data_available': self.gene_coverages_data_available,
                                'completeness_data_available': self.completeness_data_available,
                                'non_single_copy_gene_hmm_data_available': self.non_single_copy_gene_hmm_data_available, 
                                'percent_contigs_nts_described_by_collection': 0.0,
                                'percent_profile_nts_described_by_collection': 0.0,
                                'percent_contigs_nts_described_by_profile': P(self.p_meta['total_length'], self.a_meta['total_length']) ,
                                'percent_contigs_contigs_described_by_profile': P(self.p_meta['num_contigs'], self.a_meta['num_contigs']) ,
                                'percent_contigs_splits_described_by_profile': P(self.p_meta['num_splits'], self.a_meta['num_splits']) ,
                                    }

        # I am not sure whether this is the best place to do this, 
        self.summary['basics_pretty'] = {'profile': [
                                                     ('Created on', self.p_meta['creation_date']),
                                                     ('Version', self.p_meta['version']),
                                                     ('Minimum conting length', pretty(self.p_meta['min_contig_length'])),
                                                     ('Number of contigs', pretty(int(self.p_meta['num_contigs']))),
                                                     ('Number of splits', pretty(int(self.p_meta['num_splits']))),
                                                     ('Total nucleotides', humanize_n(int(self.p_meta['total_length']))),
                                                    ],
                                         'contigs': [
                                                        ('Created on', self.p_meta['creation_date']),
                                                        ('Version', self.a_meta['version']),
                                                        ('Split length', pretty(int(self.a_meta['split_length']))),
                                                        ('Number of contigs', pretty(int(self.a_meta['num_contigs']))),
                                                        ('Number of splits', pretty(int(self.a_meta['num_splits']))),
                                                        ('Total nucleotides', humanize_n(int(self.a_meta['total_length']))),
                                                        ('K-mer size', self.a_meta['kmer_size']),
                                                    ],
                                        }

        self.summary['max_shown_header_items'] = 10
        self.summary['slice_header_items_tmpl'] = '0:%d' % self.summary['max_shown_header_items']
        self.summary['num_not_shown_samples'] = len(self.p_meta['samples']) - self.summary['max_shown_header_items']
        self.summary['num_not_shown_hmm_items'] = dict([(hmm_search_source, len(self.hmm_sources_info[hmm_search_source]['genes']) - self.summary['max_shown_header_items']) for hmm_search_type, hmm_search_source in self.hmm_searches_header])

        self.summary['files'] = {}
        self.summary['collection'] = {}
        self.summary['collection_profile'] = self.collection_profile # reminder; collection_profile comes from ProfileSuperclass!
        self.summary['collection_profile_items'] = self.collection_profile.values()[0].keys()

        # add hmm items for each seach type:
        if self.non_single_copy_gene_hmm_data_available:
            self.summary['meta']['hmm_items'] = dict([(hmm_search_source, self.hmm_sources_info[hmm_search_source]['genes']) for hmm_search_type, hmm_search_source in self.hmm_searches_header])

        # summarize bins:
        for bin_id in collection_dict: 
            bin = Bin(self, bin_id, collection_dict[bin_id], self.run, self.progress)
            bin.output_directory = os.path.join(self.output_directory, 'bin_by_bin', bin_id)
            bin.bin_profile = self.collection_profile[bin_id]

            self.summary['collection'][bin_id] = bin.create()
            self.summary['collection'][bin_id]['color'] = collection_colors[bin_id] or '#212121'
            self.summary['meta']['total_nts_in_collection'] += self.summary['collection'][bin_id]['total_length']
            self.summary['meta']['num_contigs_in_collection'] += self.summary['collection'][bin_id]['num_contigs'] 

        # bins are computed, add some relevant meta info:
        self.summary['meta']['percent_contigs_nts_described_by_collection'] = '%.2f' % (self.summary['meta']['total_nts_in_collection'] * 100.0 / int(self.a_meta['total_length']))
        self.summary['meta']['percent_profile_nts_described_by_collection'] = '%.2f' % (self.summary['meta']['total_nts_in_collection'] * 100.0 / int(self.p_meta['total_length']))
        self.summary['meta']['bins'] = self.get_bins_ordered_by_completeness_and_size()

        # generate a TAB-delimited text output file for bin summaries
        summary_of_bins_matrix_output = {}
        properties = ['taxon', 'total_length', 'num_contigs', 'N50', 'GC_content', 'percent_complete', 'percent_redundancy']

        for bin_name in self.summary['collection']:
            summary_of_bins_matrix_output[bin_name] = dict([(prop, self.summary['collection'][bin_name][prop]) for prop in properties])

        output_file_obj = self.get_output_file_handle(prefix = 'general_bins_summary.txt')
        utils.store_dict_as_TAB_delimited_file(summary_of_bins_matrix_output, None, headers = ['bins'] + properties, file_obj = output_file_obj)

        # save merged matrices for bins x samples
        for table_name in self.collection_profile.values()[0].keys():
            d = {}
            for bin_id in self.collection_profile:
                d[bin_id] = self.collection_profile[bin_id][table_name]

            output_file_obj = self.get_output_file_handle(sub_directory = 'bins_across_samples', prefix = '%s.txt' % table_name)
            utils.store_dict_as_TAB_delimited_file(d, None, headers = ['bins'] + sorted(self.p_meta['samples']), file_obj = output_file_obj)

        # merge and store matrices for hmm hits
        if self.non_single_copy_gene_hmm_data_available:
            for hmm_search_source in self.summary['meta']['hmm_items']:
                # this is to keep numbers per hmm item:
                d = {}

                for bin_id in self.summary['meta']['bins']:
                    d[bin_id] = self.summary['collection'][bin_id]['hmms'][hmm_search_source]

                output_file_obj = self.get_output_file_handle(sub_directory = 'bins_across_samples', prefix = '%s.txt' % hmm_search_source, within='hmms')
                utils.store_dict_as_TAB_delimited_file(d, None, headers = ['bins'] + sorted(self.summary['meta']['hmm_items'][hmm_search_source]), file_obj = output_file_obj)

            # this is to keep number of hmm hits per bin:
            n = dict([(bin_id, {}) for bin_id in self.summary['meta']['bins']])
            for hmm_search_source in self.summary['meta']['hmm_items']:
                for bin_id in self.summary['meta']['bins']:
                    n[bin_id][hmm_search_source] =  sum(self.summary['collection'][bin_id]['hmms'][hmm_search_source].values())

            output_file_obj = self.get_output_file_handle(sub_directory = 'bins_across_samples', prefix = 'hmm_hit_totals.txt')
            utils.store_dict_as_TAB_delimited_file(n, None, headers = ['bins'] + sorted(self.summary['meta']['hmm_items']), file_obj = output_file_obj)

        # store percent abundance of each bin
        self.summary['bin_percent_recruitment'] = self.bin_percent_recruitment_per_sample
        self.summary['bin_percent_abundance_items'] = sorted(self.bin_percent_recruitment_per_sample.values()[0].keys())
        output_file_obj = self.get_output_file_handle(sub_directory = 'bins_across_samples', prefix = 'bins_percent_recruitment.txt')
        utils.store_dict_as_TAB_delimited_file(self.bin_percent_recruitment_per_sample,
                                               None,
                                               headers = ['samples'] + sorted(self.collection_profile.keys()) + ['__splits_not_binned__'],
                                               file_obj = output_file_obj)


        if self.debug:
            import json
            print json.dumps(self.summary, sort_keys=True, indent=4)

        self.index_html = SummaryHTMLOutput(self.summary, r = self.run, p = self.progress).generate()


    def get_bins_ordered_by_completeness_and_size(self):
        if self.completeness_data_available:
            return [t[2] for t in sorted([(self.summary['collection'][bin]['percent_complete'], self.summary['collection'][bin]['total_length'], bin) for bin in self.summary['collection']], reverse=True)]
        else:
            return sorted(self.summary['collection'].keys())


    def get_output_file_handle(self, sub_directory = None, prefix = 'output.txt', overwrite = False, within = None):
        if sub_directory:
            output_directory = os.path.join(self.output_directory, sub_directory)
        else:
            output_directory = self.output_directory

        if not os.path.exists(output_directory):
            filesnpaths.gen_output_directory(output_directory)

        if within:
            file_path = os.path.join(output_directory, '%s_%s' % (within, prefix))
        else:
            file_path = os.path.join(output_directory, '%s' % (prefix))

        if os.path.exists(file_path) and not overwrite:
            raise ConfigError, 'get_output_file_handle: well, this file already exists: "%s"' % file_path

        key = prefix.split('.')[0].replace('-', '_')

        if within:
            if not self.summary['files'].has_key(within):
                self.summary['files'][within] = {}
            self.summary['files'][within][key] = file_path[len(self.output_directory):].strip('/')
        else:
            self.summary['files'][key] = file_path[len(self.output_directory):].strip('/')

        return open(file_path, 'w')



class Bin:
    def __init__(self, summary, bin_id, split_ids, r = run, p = progress):
        self.summary = summary
        self.bin_id = bin_id
        self.split_ids = split_ids
        self.progress = p
        self.run = r
        self.across_samples = {}
        self.bin_profile = {}

        self.bin_info_dict = {'files': {}}

        self.output_directory = None
        self.contig_lengths = []

        # make sure all split_ids in the collection is actually in the contigs database.
        # in collections stored in the contigs database, split_ids that are not in the
        # oritinal contigs used to generate contigs database *may* end up in the
        # collections table. we gotta make sure we deal with them properly:
        missing_ids = [split_id for split_id in self.split_ids if split_id not in self.summary.split_sequences]
        if len(missing_ids):
            for missing_id in missing_ids:
                self.split_ids.remove(missing_id)

            self.run.warning('%d split id(s) in bin "%s" reported by collection "%s" is not found in the\
                              contigs database and removed from the bin summary. If this does not make\
                              any sense, you may need make sure everything is in order. The thing is,\
                              sometimes external clustering results that are added to the contigs via\
                              `anvi-populate-collections-table` may include split names that are not used\
                              while the contigs database was generated.'\
                                                % (len(missing_ids), bin_id, self.summary.collection_id))


    def create(self):
        self.progress.new('[Processing "%s"]' % self.bin_id)

        self.create_bin_dir()

        self.store_sequences_for_hmm_hits()

        self.store_contigs_fasta()

        if self.summary.completeness_data_available:
            self.access_completeness_scores()

        if self.summary.non_single_copy_gene_hmm_data_available:
            self.summarize_hmm_hits()

        self.compute_basic_stats()

        if self.summary.a_meta['genes_annotation_source']:
            self.set_taxon_calls()

        if self.summary.gene_coverages_dict:
            self.store_gene_coverages_matrix()

        self.store_profile_data()

        self.progress.end()

        return self.bin_info_dict


    def create_bin_dir(self):
        self.progress.update('Creating the output directory ...')

        if not self.output_directory:
            self.progress.end()
            raise ConfigError, 'You caled Bin.create() before setting an output directory. Anvio says "nope, thanks".'

        filesnpaths.gen_output_directory(self.output_directory)


    def get_output_file_handle(self, prefix = 'output.txt', overwrite = False, key = None):
        file_path = os.path.join(self.output_directory, '%s-%s' % (self.bin_id, prefix))
        if os.path.exists(file_path) and not overwrite:
            raise ConfigError, 'get_output_file_handle: well, this file already exists: "%s"' % file_path


        if not key:
            key = prefix.split('.')[0].replace('-', '_')

        self.bin_info_dict['files'][key] = file_path[len(self.summary.output_directory):].strip('/')

        return open(file_path, 'w')


    def store_data_in_file(self, output_file_name_posfix, content):
        output_file_obj = self.get_output_file_handle(output_file_name_posfix)
        output_file_obj.write('%s\n' % content)
        output_file_obj.close()


    def access_completeness_scores(self):
        self.progress.update('Accessing completeness scores ...')

        completeness = self.summary.completeness.get_info_for_splits(set(self.split_ids))

        self.bin_info_dict['completeness'] = completeness

        num_sources = len(completeness)

        # set up for the average completeness / redundancy scores:
        for k in ['percent_redundancy', 'percent_complete']:
            self.bin_info_dict[k] = 0.0

        # go through all single-copy gene reporting sources
        for c in completeness.values():
            for k in ['percent_redundancy', 'percent_complete']:
                self.bin_info_dict[k] += c[k]

        for k in ['percent_redundancy', 'percent_complete']:
            self.bin_info_dict[k] /= num_sources
            self.store_data_in_file('%s.txt' % k, '%.4f' % self.bin_info_dict[k])


    def store_profile_data(self):
        self.progress.update('Storing profile data ...')

        for table_name in self.bin_profile:
            output_file_obj = self.get_output_file_handle('%s.txt' % table_name)
            utils.store_dict_as_TAB_delimited_file({table_name: self.bin_profile[table_name]}, None, headers = ['bin'] + self.summary.p_meta['samples'], file_obj = output_file_obj)


    def summarize_hmm_hits(self):
        """Make sense of everything there is to make sense of regarding hmm hits.
        
           Unfortunately this is *VERY* complicated. Here we try to make sense of any
           HMM collection with respect to nubmer of hits that happens to be in splits
           associated with this bin, and split - hit associations. This function fills
           all the information into self.bin_mm_profile_dict, and the process function
           up above later makes sense of all to generate files and matrices, as well as
           dictionaries to diplay part of this information in the interface.
        """

        info_dict = {}

        # lets limit our interest space into splits that are in our bin and have hmm hits from the get go:
        split_ids_with_hmm_hits = [split_id for split_id in self.split_ids if self.summary.hmm_searches_dict.has_key(split_id)]

        for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header:
            hmm_items = self.summary.hmm_sources_info[hmm_search_source]['genes']
            info_dict[hmm_search_source] = dict([(hmm_item, 0) for hmm_item in hmm_items])

            hits_in_splits = []
            # keep track of unique identifiers of hmm hits to not count a single hit that spans across multiple splits:
            unique_identifiers_seen = set([])

            for split_id in split_ids_with_hmm_hits:
                for hmm_item, unique_identifier in self.summary.hmm_searches_dict[split_id][hmm_search_source]:
                    hits_in_splits.append((split_id, hmm_item, unique_identifier),)

                    if (unique_identifier in unique_identifiers_seen):
                        continue

                    unique_identifiers_seen.add(unique_identifier)
                    info_dict[hmm_search_source][hmm_item] += 1

            output_file_obj = self.get_output_file_handle('%s-hmm-hits.txt' % hmm_search_source)
            output_file_obj.write('contigs\thmm_profile\tunique_identifier\n')
            for item in hits_in_splits:
                output_file_obj.write('%s\n' % '\t'.join(item))
            output_file_obj.close()

        self.bin_info_dict['hmms'] = info_dict


    def store_gene_coverages_matrix(self):
        self.progress.update('Storing gene coverages ...')

        info_dict = {}
        genes_dict = {}

        gene_entry_ids_in_bin = set([])
        for split_name in self.split_ids:
            gene_entry_ids_in_bin.update(self.summary.split_to_genes_in_splits_ids[split_name])

        info_dict['num_genes_found'] = len(gene_entry_ids_in_bin)

        headers = ['function', 'contig', 'start', 'stop', 'direction']
        for gene_entry_id in gene_entry_ids_in_bin:
            prot_id = self.summary.genes_in_splits[gene_entry_id]['prot']
            genes_dict[prot_id] = {}

            # first fill in sample independent information;
            for header in headers:
                genes_dict[prot_id][header] = self.summary.genes_in_contigs_dict[prot_id][header]

            # then fill in distribution across samples:
            for sample_name in self.summary.p_meta['samples']:
                genes_dict[prot_id][sample_name] = self.summary.gene_coverages_dict[prot_id][sample_name]

            # finally add the sequence:
            contig = self.summary.genes_in_contigs_dict[prot_id]['contig']
            start = self.summary.genes_in_contigs_dict[prot_id]['start']
            stop = self.summary.genes_in_contigs_dict[prot_id]['stop']
            genes_dict[prot_id]['sequence'] = self.summary.contig_sequences[contig]['sequence'][start:stop]

        output_file_obj = self.get_output_file_handle('functions.txt')
        utils.store_dict_as_TAB_delimited_file(genes_dict, None, headers = ['prot'] + headers + self.summary.p_meta['samples'] + ['sequence'], file_obj = output_file_obj)

        self.bin_info_dict['genes'] = info_dict


    def store_sequences_for_hmm_hits(self):
        s = SequencesForHMMHits(self.summary.contigs_db_path)
        hmm_sequences_dict = s.get_hmm_sequences_dict_for_splits({self.bin_id: self.split_ids})

        single_copy_gene_hmm_sources = [hmm_search_source for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header]
        non_single_copy_gene_hmm_sources = self.summary.completeness.sources

        for hmm_search_source in single_copy_gene_hmm_sources + non_single_copy_gene_hmm_sources:
            filtered_hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'source', set([hmm_search_source]))

            output_file_obj = self.get_output_file_handle('%s-hmm-sequences.txt' % hmm_search_source, key = hmm_search_source)

            for gene_unique_id in filtered_hmm_sequences_dict:
                header, sequence = s.get_FASTA_header_and_sequence_for_gene_unique_id(hmm_sequences_dict, gene_unique_id)
                output_file_obj.write('>%s\n%s\n' % (header, sequence))


    def store_contigs_fasta(self):
        """Storing contig sequences.
        
           This is not an easy problem. We split contigs into smaller sequences at the beginning. Only
           a part of a given contig may be used during the binning process. On the other hand we can't
           simply store sequences of splits, whenever possible, we must store the entire sequence of
           the contig (only if all splits are selected from a contig in to the same bin). So, this
           function first identifies all splits coming from the same parent, then identifies sequential
           blocks of splits (see `SequentialBlocks` class), then checks whether all splits of a given
           contig is included in the bin. If that is the case, it puts the contig as a single entry,
           witht he identical FASTA id to the original contigs in the assembly file. Otherwise it appends
           `_partial_X_Y` to the FASTA id, X and Y being the start and stop positions.
        """

        self.progress.update('Creating the FASTA file ...')

        # store original split names:
        self.store_data_in_file('original_split_names.txt', '\n'.join(self.split_ids))

        fasta_file = self.get_output_file_handle('contigs.fa')

        # some null values:
        self.bin_info_dict['total_length'] = 0
        self.bin_info_dict['num_contigs'] = 0

        # this dict will keep all the contig ids found in this bin with split names ordered:
        contigs_represented = utils.get_contigs_splits_dict(self.split_ids, self.summary.splits_basic_info) 

        # now it is time to go through each contig found in contigs_represented to
        # figure out what fraction of the contig is in fact in this bin
        for contig_id in contigs_represented:
            splits_order = contigs_represented[contig_id].keys()

            self.progress.update('Creating the FASTA file :: Identifying sequential blocks ...')
            # this is critical: sequential_blocks is a list of one ore more lists, where each item of this list
            # describes a range of splits that follow each other to represent a coherent
            # chunk of the parent sequence (if all splits from a contig is selected into this bin,
            # then there would be one list item that spans across the entire contig):
            sequential_blocks = ccollections.GetSequentialBlocksOfSplits(splits_order).process()

            for sequential_block in sequential_blocks:
                self.progress.update('Creating the FASTA file :: Identifying the portion of contig represented ...')
                first_split = contigs_represented[contig_id][sequential_block[0]]
                last_split = contigs_represented[contig_id][sequential_block[-1]]

                contig_sequence_start_in_splits = self.summary.splits_basic_info[first_split]['start']
                contig_sequence_end_in_splits = self.summary.splits_basic_info[last_split]['end']

                # so this much of the contig is represented by its splits:
                total_contig_length_in_splits = contig_sequence_end_in_splits - contig_sequence_start_in_splits

                # and this is is actual length:
                contig_sequence_length = self.summary.contigs_basic_info[contig_id]['length']

                if contig_sequence_length == total_contig_length_in_splits:
                    # the entireity of the contig is represented!
                    appendix = ''
                else:
                    appendix = '_partial_%d_%d' % (contig_sequence_start_in_splits, contig_sequence_end_in_splits)

                sequence = ''
                self.progress.update('Creating the FASTA file :: Reconstructing contig sequence from splits ...')
                for split_order in sequential_block:
                    sequence += self.summary.split_sequences[contigs_represented[contig_id][split_order]]

                fasta_id = contig_id + appendix

                self.progress.update('Creating the FASTA file :: Writing contig sequence into file ...')
                fasta_file.write('>%s\n' % fasta_id)
                fasta_file.write('%s\n' % textwrap.fill(sequence, 80, break_on_hyphens = False))

                # fill in basic info about contigs in bin
                len_seq = len(sequence)
                self.bin_info_dict['total_length'] += len_seq
                self.contig_lengths.append(len_seq)
                self.bin_info_dict['num_contigs'] += 1

        fasta_file.close()

        self.store_data_in_file('num_contigs.txt', '%d' % self.bin_info_dict['num_contigs'])
        self.store_data_in_file('total_length.txt', '%d' % self.bin_info_dict['total_length'])


    def set_taxon_calls(self):
        self.progress.update('Filling in taxonomy info ...')

        taxon_calls_counter = Counter()
        for split_id in self.split_ids:
            taxon_calls_counter[self.summary.genes_in_splits_summary_dict[split_id]['taxonomy']] += 1

        taxon_calls = sorted([list(tc) for tc in taxon_calls_counter.items()], key = lambda x: int(x[1]), reverse = True)

        self.bin_info_dict['taxon_calls'] = taxon_calls

        # taxon_calls = [(None, 129), ('Propionibacterium avidum', 120), ('Propionibacterium acnes', 5)]
        l = [tc for tc in taxon_calls if tc[0]]
        num_calls = sum(taxon_calls_counter.values())
        # l = [('Propionibacterium avidum', 120), ('Propionibacterium acnes', 5)]
        if l and l[0][1] > num_calls / 4.0:
            # if l[0] is associated with more than 25 percent of splits:
            self.bin_info_dict['taxon'] = l[0][0]
        else:
            self.bin_info_dict['taxon'] = 'Unknown'

        # convert to percents..
        for tc in taxon_calls:
            tc[1] = tc[1] * 100.0 / num_calls


    def compute_basic_stats(self):
        self.progress.update('Computing basic stats ...')

        self.bin_info_dict['N50'] = utils.get_N50(self.contig_lengths)
        self.bin_info_dict['GC_content'] = numpy.mean([self.summary.splits_basic_info[split_id]['gc_content'] for split_id in self.split_ids]) * 100

        self.store_data_in_file('N50.txt', '%d' % self.bin_info_dict['N50'])
        self.store_data_in_file('GC_content.txt', '%.4f' % self.bin_info_dict['GC_content'])



