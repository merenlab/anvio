# coding: utf-8
# pylint: disable=line-too-long
"""Summarizes information for a collection."""

import os
import sys
import gzip
import numpy
import shutil
import textwrap

from collections import Counter

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.sequence as sequence
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.errors import ConfigError
from anvio.dbops import DatabasesMetaclass, ContigsSuperclass, PanSuperclass
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


class ArgsTemplateForSummarizerClass:
    """Just a dummy args template for ad hoc use of summary classes.

    You can use it like this:

        >>> args = summarizer.ArgsTemplateForSummarizerClass()
        >>> args.profile_db = profile_db_path
        >>> args.contigs_db = contigs_db_path
        >>> args.collection_name = collection_name
        >>> args.output_dir = output_dir

        >>> summarizer.ProfileSummarizer(args)
    """

    def __init__(self):
        self.profile_db = None
        self.pan_db = None
        self.contigs_db = None
        self.collection_name = None
        self.taxonomic_level = 't_genus'
        self.list_collections = None
        self.debug = None
        self.quick_summary = False
        self.skip_check_collection_name = False
        self.skip_init_functions = False
        self.cog_data_dir = None
        self.output_dir = filesnpaths.get_temp_directory_path()


class SummarizerSuperClass(object):
    def __init__(self, args, r=run, p=progress):
        self.summary = {}
        self.collection_name = None

        self.collections = ccollections.Collections()

        if self.summary_type == 'pan':
            self.collections.populate_collections_dict(self.pan_db_path)
        else:
            self.collections.populate_collections_dict(self.pan_db_path if self.summary_type == 'pan' else self.profile_db_path)
            self.collections.populate_collections_dict(self.contigs_db_path) if self.contigs_db_path else None

        if args.list_collections:
            self.collections.list_collections()
            sys.exit()

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.collection_name = A('collection_name')
        self.skip_check_collection_name = A('skip_check_collection_name')
        self.skip_init_functions = A('skip_init_functions')
        self.output_directory = A('output_dir')
        self.quick = A('quick_summary')
        self.debug = A('debug')
        self.taxonomic_level = A('taxonomic_level')
        self.cog_data_dir = A('cog_data_dir')

        self.sanity_check()

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=True)


    def sanity_check(self):
        if not self.skip_check_collection_name:
            if not self.collection_name:
                raise ConfigError, "You must specify a collection id :/"

            if self.collection_name not in self.collections.collections_dict:
                raise ConfigError, "%s is not a valid collection ID. See a list of available ones with '--list-collections' flag" % self.collection_name

        self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=True)


    def get_output_file_handle(self, sub_directory=None, prefix='output.txt', overwrite=False, within=None, compress_output=False, add_project_name=False):
        if sub_directory:
            output_directory = os.path.join(self.output_directory, sub_directory)
        else:
            output_directory = self.output_directory

        if not os.path.exists(output_directory):
            filesnpaths.gen_output_directory(output_directory)

        key = prefix.split('.')[0].replace('-', '_')

        if add_project_name:
            prefix = '%s_%s' % (self.p_meta['project_name'].replace(' ', '_'), prefix)

        if within:
            file_path = os.path.join(output_directory, '%s_%s' % (within, prefix))
        else:
            file_path = os.path.join(output_directory, '%s' % (prefix))

        if compress_output:
            file_path += '.gz'

        if os.path.exists(file_path) and not overwrite:
            raise ConfigError, 'get_output_file_handle: well, this file already exists: "%s"' % file_path

        if within:
            if within not in self.summary['files']:
                self.summary['files'][within] = {}
            self.summary['files'][within][key] = file_path[len(self.output_directory):].strip('/')
        else:
            self.summary['files'][key] = file_path[len(self.output_directory):].strip('/')

        if compress_output:
            return gzip.open(file_path, 'wb')
        else:
            return open(file_path, 'w')


class PanSummarizer(PanSuperclass, SummarizerSuperClass):
    """Creates a dictionary of summary for anvi'o pan profiles"""
    def __init__(self, args=None, r=run, p=progress):
        self.summary_type = 'pan'
        self.debug = False
        self.quick = False
        self.pan_db_path = None
        self.output_directory = None
        self.genomes_storage_path = None

        PanSuperclass.__init__(self, args, run, progress)
        if not self.genomes_storage_is_available:
            raise ConfigError, "No genomes storage no summary. Yes. Very simple stuff."

        SummarizerSuperClass.__init__(self, args, self.run, self.progress)

        # init protein clusters and functins from Pan super.
        self.init_protein_clusters()

        if not self.skip_init_functions:
            self.init_protein_clusters_functions()

        # see if COG functions or categories are available
        self.cog_functions_are_called = 'COG_FUNCTION' in self.protein_clusters_function_sources
        self.cog_categories_are_called = 'COG_CATEGORY' in self.protein_clusters_function_sources


    def process(self):
        # init profile data for colletion.
        collection_dict, bins_info_dict = self.init_collection_profile(self.collection_name)

        # let bin names known to all
        bin_ids = self.collection_profile.keys()

        genome_names = ', '.join(self.protein_clusters.values()[0].keys())

        # set up the initial summary dictionary
        self.summary['meta'] = { \
                'quick': self.quick,
                'cog_functions_are_called': self.cog_functions_are_called,
                'cog_categories_are_called': self.cog_categories_are_called,
                'output_directory': self.output_directory,
                'summary_type': self.summary_type,
                'collection': bin_ids,
                'num_bins': len(bin_ids),
                'collection_name': self.collection_name,
                'total_num_genes_in_collection': 0,
                'anvio_version': __version__,
                'pan': self.p_meta,
                'genomes': {'num_genomes': self.genomes_storage.num_genomes,
                            'functions_available': True if len(self.protein_clusters_function_sources) else False,
                            'function_sources': self.protein_clusters_function_sources},
                'percent_of_genes_collection': 0.0,
                'genome_names': genome_names
        }

        # I am not sure whether this is the best place to do this,
        self.summary['basics_pretty'] = { \
                'pan': [('Created on', self.p_meta['creation_date']),
                        ('Version', anvio.__pan__version__),
                        ('Number of genes', pretty(int(self.p_meta['num_genes_in_protein_clusters']))),
                        ('Number of protein clusters', pretty(int(self.p_meta['num_protein_clusters']))),
                        ('Partial genes excluded', 'Yes' if self.p_meta['exclude_partial_gene_calls'] else 'No'),
                        ('Maxbit parameter', self.p_meta['maxbit']),
                        ('PC min occurrence parameter', pretty(int(self.p_meta['pc_min_occurrence']))),
                        ('MCL inflation parameter', self.p_meta['mcl_inflation']),
                        ('NCBI blastp or DIAMOND?', 'NCBI blastp' if self.p_meta['use_ncbi_blast'] else ('DIAMOND (and it was %s)' % ('sensitive' if self.p_meta['diamond_sensitive'] else 'not sensitive'))),
                        ('Number of genomes used', pretty(int(self.p_meta['num_genomes'])))],

                'genomes': [('Created on', 'Storage DB knows nothing :('),
                            ('Version', anvio.__genomes_storage_version__),
                            ('Number of genomes described', pretty(self.genomes_storage.num_genomes)),
                            ('Functional annotation', 'Available' if len(self.protein_clusters_function_sources) else 'Not available :/'),
                            ('Functional annotation sources', '--' if not len(self.protein_clusters_function_sources) else ', '.join(self.protein_clusters_function_sources))],
        }

        self.summary['files'] = {}
        self.summary['collection_profile'] = self.collection_profile # reminder; collection_profile comes from the superclass!

        self.generate_protein_clusters_file(collection_dict)

        if self.debug:
            import json
            print json.dumps(self.summary, sort_keys=True, indent=4)

        self.index_html = SummaryHTMLOutput(self.summary, r=self.run, p=self.progress).generate(quick=self.quick)


    def generate_protein_clusters_file(self, collection_dict, compress_output=True):
        """Generates the proteins summary file"""

        self.progress.new('Protein clusters summary file')
        self.progress.update('...')

        # generate a dict of protein cluster ~ bin id relationships
        pc_name_to_bin_name= dict(zip(self.protein_clusters_in_pan_db_but_not_binned, [None] * len(self.protein_clusters_in_pan_db_but_not_binned)))
        for bin_id in collection_dict: 
            for pc_name in collection_dict[bin_id]:
                pc_name_to_bin_name[pc_name] = bin_id

        ###############################################
        # generate an output file for protein clusters.
        ###############################################
        output_file_obj = self.get_output_file_handle(prefix='protein_clusters_summary.txt', compress_output=compress_output, add_project_name=True)

        # standard headers
        header = ['unique_id', 'protein_cluster_id', 'bin_name', 'genome_name', 'gene_callers_id']

        # extend the header with functions if there are any 
        for source in self.protein_clusters_function_sources:
            if self.quick:
                header.append(source + '_ACC')
            else:
                header.append(source + '_ACC')
                header.append(source)

        # if this is not a quick summary, have AA sequences in the output
        header.append('aa_sequence') if not self.quick else None

        # write the header
        output_file_obj.write('\t'.join(header) + '\n')

        # uber loop for the file content
        unique_id = 1
        for pc_name in self.protein_clusters:
            for genome_name in self.protein_clusters[pc_name]:
                for gene_caller_id in self.protein_clusters[pc_name][genome_name]:
                    entry = [unique_id, pc_name, pc_name_to_bin_name[pc_name], genome_name, gene_caller_id]

                    for function_source in self.protein_clusters_function_sources:
                        annotations_dict = self.protein_clusters_functions_dict[pc_name][genome_name][gene_caller_id]
                        if function_source in annotations_dict:
                            annotation_blob = self.protein_clusters_functions_dict[pc_name][genome_name][gene_caller_id][function_source]
                            accessions, annotations = [l.split('!!!') for l in annotation_blob.split("|||")]
                            entry.append('|'.join(accessions))
                            entry.append('|'.join(annotations))
                        else:
                            entry.append('')
                            entry.append('')

                    entry.append(self.genomes_storage.get_gene_sequence(genome_name, gene_caller_id)) if not self.quick else None

                    output_file_obj.write('\t'.join([str(e) if e not in [None, 'UNKNOWN'] else '' for e in entry]) + '\n')
                    unique_id += 1


        # we're done here.
        output_file_obj.close()

        self.progress.end()


class ProfileSummarizer(DatabasesMetaclass, SummarizerSuperClass):
    """Creates an über dictionary of 'summary' for anvi'o profiles."""
    def __init__(self, args=None, r=run, p=progress):
        self.summary = {}

        self.summary_type = 'profile'
        self.debug = False
        self.quick = False
        self.profile_db_path = None
        self.contigs_db_path = None
        self.output_directory = None
        self.split_names_per_bin = None
        self.completeness_data_available = False
        self.gene_coverages_data_available = False
        self.non_single_copy_gene_hmm_data_available = False

        DatabasesMetaclass.__init__(self, args, run, progress)
        SummarizerSuperClass.__init__(self, args, self.run, self.progress)

        # databases initiated, let's make sure we have gene covereges data avaialable.
        if self.gene_coverages_dict:
            self.gene_coverages_data_available = True

        self.init_splits_taxonomy(self.taxonomic_level)


    def process(self):
        # init profile data for colletion.
        collection_dict, bins_info_dict = self.init_collection_profile(self.collection_name)

        # let bin names known to all
        bin_ids = self.collection_profile.keys()

        # load completeness information if available
        self.completeness = completeness.Completeness(self.contigs_db_path)
        if len(self.completeness.sources):
            self.completeness_data_available = True

        # load HMM sources for non-single-copy genes if available
        if self.non_singlecopy_gene_hmm_sources and not self.quick:
            self.init_non_singlecopy_gene_hmm_sources()
            self.non_single_copy_gene_hmm_data_available = True

        # load gene functions from contigs db superclass
        self.init_functions()

        # set up the initial summary dictionary
        self.summary['meta'] = {'quick': self.quick,
                                'summary_type': self.summary_type,
                                'output_directory': self.output_directory,
                                'collection': bin_ids,
                                'num_bins': len(bin_ids),
                                'collection_name': self.collection_name,
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
                                'percent_contigs_nts_described_by_profile': P(self.p_meta['total_length'], self.a_meta['total_length']),
                                'percent_contigs_contigs_described_by_profile': P(self.p_meta['num_contigs'], self.a_meta['num_contigs']),
                                'percent_contigs_splits_described_by_profile': P(self.p_meta['num_splits'], self.a_meta['num_splits']),
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
        self.summary['collection_profile_items'] = [] if not len(self.collection_profile.values()) else self.collection_profile.values()[0].keys()

        # add hmm items for each seach type:
        if self.non_single_copy_gene_hmm_data_available:
            self.summary['meta']['hmm_items'] = dict([(hmm_search_source, self.hmm_sources_info[hmm_search_source]['genes']) for hmm_search_type, hmm_search_source in self.hmm_searches_header])

        # summarize bins:
        for i in range(0, len(bin_ids)):
            bin_id = bin_ids[i]
            self.progress.new('[Processing "%s" (%d of %d)]' % (bin_id, i + 1, len(bin_ids)))
            bin = Bin(self, bin_id, collection_dict[bin_id], self.run, self.progress)
            bin.output_directory = os.path.join(self.output_directory, 'bin_by_bin', bin_id)
            bin.bin_profile = self.collection_profile[bin_id]

            self.summary['collection'][bin_id] = bin.create()
            self.summary['collection'][bin_id]['color'] = bins_info_dict[bin_id]['html_color'] or '#212121'
            self.summary['collection'][bin_id]['source'] = bins_info_dict[bin_id]['source'] or 'unknown_source'
            self.summary['meta']['total_nts_in_collection'] += self.summary['collection'][bin_id]['total_length']
            self.summary['meta']['num_contigs_in_collection'] += self.summary['collection'][bin_id]['num_contigs']
            self.progress.end()

        # bins are computed, add some relevant meta info:
        self.summary['meta']['percent_contigs_nts_described_by_collection'] = '%.2f' % (self.summary['meta']['total_nts_in_collection'] * 100.0 / int(self.a_meta['total_length']))
        self.summary['meta']['percent_profile_nts_described_by_collection'] = '%.2f' % (self.summary['meta']['total_nts_in_collection'] * 100.0 / int(self.p_meta['total_length']))
        self.summary['meta']['bins'] = self.get_bins_ordered_by_completeness_and_size()

        if not self.quick:
            # generate a TAB-delimited text output file for bin summaries
            summary_of_bins_matrix_output = {}
            properties = ['taxon', 'total_length', 'num_contigs', 'N50', 'GC_content']
            if self.completeness_data_available:
                properties += ['percent_complete', 'percent_redundancy']

            for bin_name in self.summary['collection']:
                summary_of_bins_matrix_output[bin_name] = dict([(prop, self.summary['collection'][bin_name][prop]) for prop in properties])

            output_file_obj = self.get_output_file_handle(prefix='general_bins_summary.txt')
            utils.store_dict_as_TAB_delimited_file(summary_of_bins_matrix_output, None, headers=['bins'] + properties, file_obj=output_file_obj)

            # save merged matrices for bins x samples
            for table_name in self.summary['collection_profile_items']:
                d = {}
                for bin_id in self.collection_profile:
                    d[bin_id] = self.collection_profile[bin_id][table_name]

                output_file_obj = self.get_output_file_handle(sub_directory='bins_across_samples', prefix='%s.txt' % table_name)
                utils.store_dict_as_TAB_delimited_file(d, None, headers=['bins'] + sorted(self.p_meta['samples']), file_obj=output_file_obj)

            # merge and store matrices for hmm hits
            if self.non_single_copy_gene_hmm_data_available:
                for hmm_search_source in self.summary['meta']['hmm_items']:
                    # this is to keep numbers per hmm item:
                    d = {}

                    for bin_id in self.summary['meta']['bins']:
                        d[bin_id] = self.summary['collection'][bin_id]['hmms'][hmm_search_source]

                    output_file_obj = self.get_output_file_handle(sub_directory='bins_across_samples', prefix='%s.txt' % hmm_search_source, within='hmms')
                    utils.store_dict_as_TAB_delimited_file(d, None, headers=['bins'] + sorted(self.summary['meta']['hmm_items'][hmm_search_source]), file_obj=output_file_obj)

                # this is to keep number of hmm hits per bin:
                n = dict([(bin_id, {}) for bin_id in self.summary['meta']['bins']])
                for hmm_search_source in self.summary['meta']['hmm_items']:
                    for bin_id in self.summary['meta']['bins']:
                        n[bin_id][hmm_search_source] = sum(self.summary['collection'][bin_id]['hmms'][hmm_search_source].values())

                output_file_obj = self.get_output_file_handle(sub_directory='bins_across_samples', prefix='hmm_hit_totals.txt')
                utils.store_dict_as_TAB_delimited_file(n, None, headers=['bins'] + sorted(self.summary['meta']['hmm_items']), file_obj=output_file_obj)

            # store percent abundance of each bin
            if self.p_meta['blank']:
                self.summary['bin_percent_recruitment'] = None
                self.summary['bin_percent_abundance_items'] = None
            else:
                self.summary['bin_percent_recruitment'] = self.bin_percent_recruitment_per_sample
                self.summary['bin_percent_abundance_items'] = sorted(self.bin_percent_recruitment_per_sample.values()[0].keys())
                output_file_obj = self.get_output_file_handle(sub_directory='bins_across_samples', prefix='bins_percent_recruitment.txt')
                utils.store_dict_as_TAB_delimited_file(self.bin_percent_recruitment_per_sample,
                                                       None,
                                                       headers=['samples'] + sorted(self.collection_profile.keys()) + ['__splits_not_binned__'],
                                                       file_obj=output_file_obj)


        if self.debug:
            import json
            print json.dumps(self.summary, sort_keys=True, indent=4)

        self.index_html = SummaryHTMLOutput(self.summary, r=self.run, p=self.progress).generate(quick=self.quick)


    def get_bins_ordered_by_completeness_and_size(self):
        if self.completeness_data_available:
            return [t[2] for t in sorted([(self.summary['collection'][bin]['percent_complete'], self.summary['collection'][bin]['total_length'], bin) for bin in self.summary['collection']], reverse=True)]
        else:
            return sorted(self.summary['collection'].keys())


class Bin:
    def __init__(self, summary, bin_id, split_ids, r=run, p=progress):
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
                                                % (len(missing_ids), bin_id, self.summary.collection_name))


        self.gene_caller_ids = self.get_gene_caller_ids()


    def create(self):
        self.create_bin_dir()

        self.store_sequences_for_hmm_hits()

        self.store_contigs_fasta()

        if self.summary.completeness_data_available:
            self.access_completeness_scores()

        if self.summary.non_single_copy_gene_hmm_data_available:
            self.summarize_hmm_hits()

        self.compute_basic_stats()

        self.set_taxon_calls()

        if self.summary.gene_coverages_dict:
            self.store_gene_coverages_and_functions()

        self.store_profile_data()

        return self.bin_info_dict


    def create_bin_dir(self):
        self.progress.update('Creating the output directory ...')

        if not self.output_directory:
            self.progress.end()
            raise ConfigError, 'You caled Bin.create() before setting an output directory. Anvio says "nope, thanks".'

        filesnpaths.gen_output_directory(self.output_directory)


    def access_completeness_scores(self):
        self.progress.update('Accessing completeness scores ...')

        p_completion, p_redundancy, domain, domain_confidence, results_dict = self.summary.completeness.get_info_for_splits(set(self.split_ids))

        self.bin_info_dict['completeness'] = results_dict

        self.bin_info_dict['percent_redundancy'] = p_redundancy
        self.bin_info_dict['percent_complete'] = p_completion
        self.bin_info_dict['scg_domain'] = domain
        self.bin_info_dict['scg_domain_confidence'] = domain_confidence

        for k in ['percent_redundancy', 'percent_complete']:
            self.store_data_in_file('%s.txt' % k, '%.4f' % self.bin_info_dict[k])

        self.store_data_in_file('scg_domain.txt', '%s' % self.bin_info_dict['scg_domain'])
        self.store_data_in_file('scg_domain_confidence.txt', '%.2f' % self.bin_info_dict['scg_domain_confidence'])


    def store_profile_data(self):

        if self.summary.quick:
            return

        self.progress.update('Storing profile data ...')

        for table_name in self.bin_profile:
            output_file_obj = self.get_output_file_handle('%s.txt' % table_name)
            utils.store_dict_as_TAB_delimited_file({table_name: self.bin_profile[table_name]}, None, headers=['bin'] + self.summary.p_meta['samples'], file_obj=output_file_obj)


    def summarize_hmm_hits(self):
        """Make sense of everything there is to make sense of regarding hmm hits.

           Unfortunately this is *VERY* complicated. Here we try to make sense of any
           HMM collection with respect to nubmer of hits that happens to be in splits
           associated with this bin, and split - hit associations. This function fills
           all the information into self.bin_mm_profile_dict, and the process function
           up above later makes sense of all to generate files and matrices, as well as
           dictionaries to diplay part of this information in the interface.
        """

        if self.summary.quick:
            return

        info_dict = {}

        # lets limit our interest space into splits that are in our bin and have hmm hits from the get go:
        split_ids_with_hmm_hits = [split_id for split_id in self.split_ids if split_id in self.summary.hmm_searches_dict]

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


    def get_gene_caller_ids(self):
        """Returns a set of gene caller ids in the bin.

        Because splits can be cut from arbitrary locations, we may have partial hits of genes in a bin.
        we don't want genes to appear in a bin more than once due to this, or end up appearing in
        two different bins just because one bin has a fraction of a gene. here we will build the
        genes_dict, which will contain every gene hit in all splits that are found in this genome
        bin.
        """

        genes_dict = {}

        for split_name in self.split_ids:
            if split_name not in self.summary.split_name_to_gene_caller_ids_dict:
                continue

            for gene_entry_id in self.summary.split_name_to_gene_caller_ids_dict[split_name]:
                gene_call_in_split = self.summary.genes_in_splits[gene_entry_id]
                gene_callers_id = gene_call_in_split['gene_callers_id']

                if gene_callers_id in genes_dict:
                    genes_dict[gene_callers_id].append(gene_call_in_split)
                else:
                    genes_dict[gene_callers_id] = [gene_call_in_split]

        # here we have every gene hit in this bin stored in genes_dict. what we will do is to find gene
        # call ids for genes more than 90% of which apper to be in this bin (so nothing wil be reported for
        # a gene where only like 20% of it ended up in this bin).
        gene_callers_ids_for_complete_genes = set([])
        for gene_caller_id in genes_dict:
            if sum([x['percentage_in_split'] for x in genes_dict[gene_caller_id]]) > 90:
                gene_callers_ids_for_complete_genes.add(gene_caller_id)

        del genes_dict

        return gene_callers_ids_for_complete_genes


    def get_gene_coverages_across_samples_dict(self):
        d = {}

        for gene_callers_id in self.gene_caller_ids:
            d[gene_callers_id] = {}

            for sample_name in self.summary.p_meta['samples']:
                d[gene_callers_id][sample_name] = self.summary.gene_coverages_dict[gene_callers_id][sample_name]

        return d


    def store_gene_coverages_and_functions(self):
        if self.summary.quick:
            return

        self.progress.update('Sorting out gene calls ...')

        # here we get the dictionary `d`, which will be holding gene coverages
        # then fill the rest up with funcitons and others stuff in the following
        # for loop
        d = self.get_gene_coverages_across_samples_dict()

        headers = ['contig', 'start', 'stop', 'direction']
        for gene_callers_id in self.gene_caller_ids:
            # add sample independent information into `d`;
            for header in headers:
                d[gene_callers_id][header] = self.summary.genes_in_contigs_dict[gene_callers_id][header]

            self.progress.update('Sorting out functions ...')
            # add functions if there are any:
            if len(self.summary.gene_function_call_sources):
                for source in self.summary.gene_function_call_sources:
                    if gene_callers_id not in self.summary.gene_function_calls_dict:
                        # this gene did not get any functional annotation
                        d[gene_callers_id][source] = ''
                        d[gene_callers_id][source + ' (ACCESSION)'] = ''
                        continue

                    if self.summary.gene_function_calls_dict[gene_callers_id][source]:
                        d[gene_callers_id][source + ' (ACCESSION)'] = self.summary.gene_function_calls_dict[gene_callers_id][source][0]
                        d[gene_callers_id][source] = self.summary.gene_function_calls_dict[gene_callers_id][source][1]
                    else:
                        d[gene_callers_id][source + ' (ACCESSION)'] = ''
                        d[gene_callers_id][source] = ''

            # finally add the sequence:
            contig = self.summary.genes_in_contigs_dict[gene_callers_id]['contig']
            start = self.summary.genes_in_contigs_dict[gene_callers_id]['start']
            stop = self.summary.genes_in_contigs_dict[gene_callers_id]['stop']
            d[gene_callers_id]['sequence'] = self.summary.contig_sequences[contig]['sequence'][start:stop]

        output_file_obj = self.get_output_file_handle('functions.txt')

        if self.summary.gene_function_call_sources:
            sources = [[source, source + ' (ACCESSION)'] for source in self.summary.gene_function_call_sources]
            headers = ['prot'] + headers + self.summary.p_meta['samples'] + [item for sublist in sources for item in sublist] + ['sequence']
        else:
            headers = ['prot'] + headers + self.summary.p_meta['samples'] + ['sequence']

        self.progress.update('Stroing gene coverages and functions ...')
        utils.store_dict_as_TAB_delimited_file(d, None, headers=headers, file_obj=output_file_obj)

        self.bin_info_dict['genes'] = {'num_genes_found': len(self.gene_caller_ids)}


    def store_sequences_for_hmm_hits(self):
        if self.summary.quick:
            return

        s = SequencesForHMMHits(self.summary.contigs_db_path)
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits({self.bin_id: self.split_ids})

        single_copy_gene_hmm_sources = [hmm_search_source for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header]
        non_single_copy_gene_hmm_sources = self.summary.completeness.sources

        for hmm_search_source in single_copy_gene_hmm_sources + non_single_copy_gene_hmm_sources:
            filtered_hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'source', set([hmm_search_source]))

            output_file_obj = self.get_output_file_handle('%s-hmm-sequences.txt' % hmm_search_source, key=hmm_search_source)

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

        if self.summary.quick:
            self.bin_info_dict['total_length'] = sum([self.summary.splits_basic_info[split_name]['length'] for split_name in self.summary.splits_basic_info if split_name in self.split_ids])
            self.bin_info_dict['num_contigs'] = len(set([self.summary.splits_basic_info[split_name]['parent'] for split_name in self.summary.splits_basic_info if split_name in self.split_ids]))
            return

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
                fasta_file.write('%s\n' % textwrap.fill(sequence, 80, break_on_hyphens=False))

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

        self.bin_info_dict['taxon_calls'] = []
        self.bin_info_dict['taxon'] = 'Unknown'

        if not self.summary.a_meta['taxonomy_source']:
            return

        taxon_calls_counter = Counter()
        for split_id in self.split_ids:
            if split_id in self.summary.splits_taxonomy_dict:
                taxon_calls_counter[self.summary.splits_taxonomy_dict[split_id]] += 1
            else:
                taxon_calls_counter['None'] += 1

        taxon_calls = sorted([list(tc) for tc in taxon_calls_counter.items()], key=lambda x: int(x[1]), reverse=True)

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

        self.store_data_in_file('N50.txt', '%d' % self.bin_info_dict['N50'] if self.bin_info_dict['N50'] else 'NA')
        self.store_data_in_file('GC_content.txt', '%.4f' % self.bin_info_dict['GC_content'])


    def get_output_file_handle(self, prefix='output.txt', overwrite=False, key=None):
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


def get_contigs_db_info_dict(contigs_db_path, run=run, progress=progress, include_AA_counts=False, split_names=None):
    """Returns an info dict for a given contigs db"""

    class Args:
        def __init__(self):
            self.contigs_db = contigs_db_path

    args = Args()
    run = run
    progress = progress
    run.verbose = False
    progress.verbose = False
    c = ContigsSuperclass(args, r=run, p=progress)

    info_dict = {'path': contigs_db_path}

    for key in c.a_meta:
        info_dict[key] = c.a_meta[key]

    # Two different strategies here depending on whether we work with a given set if split ids or
    # everything in the contigs database.
    if split_names:
        split_names = set(split_names)
        c.init_split_sequences()
        seq = ''.join([c.split_sequences[split_name] for split_name in split_names])
        info_dict['gene_caller_ids'] = set([e['gene_callers_id'] for e in c.genes_in_splits.values() if e['split'] in split_names])
    else:
        c.init_contig_sequences()
        seq = ''.join([e['sequence'] for e in c.contig_sequences.values()])
        info_dict['gene_caller_ids'] = c.genes_in_contigs_dict.keys()

    info_dict['gc_content'] = sequence.Composition(seq).GC_content
    info_dict['total_length'] = len(seq)

    info_dict['partial_gene_calls'] = set([])
    for gene_caller_id in info_dict['gene_caller_ids']:
        if c.genes_in_contigs_dict[gene_caller_id]['partial']:
            info_dict['partial_gene_calls'].add(gene_caller_id)

    info_dict['num_genes'] = len(info_dict['gene_caller_ids'])
    info_dict['gene_lengths'] = dict([(gene_caller_id, (c.genes_in_contigs_dict[gene_caller_id]['stop'] - c.genes_in_contigs_dict[gene_caller_id]['start'])) for gene_caller_id in info_dict['gene_caller_ids']])
    info_dict['avg_gene_length'] = numpy.mean(info_dict['gene_lengths'].values())
    info_dict['num_genes_per_kb'] = info_dict['num_genes'] * 1000.0 / info_dict['total_length']

    # get completeness / contamination estimates
    p_completion, p_redundancy, domain, domain_confidence, results_dict = completeness.Completeness(contigs_db_path).get_info_for_splits(split_names if split_names else set(c.splits_basic_info.keys()))

    info_dict['percent_complete'] = p_completion
    info_dict['percent_redundancy'] = p_redundancy
    info_dict['scg_domain'] = domain
    info_dict['scg_domain_confidence'] = domain_confidence

    info_dict['hmms_for_scgs_were_run'] = True if len(results_dict) else False

    # lets get all amino acids used in all complete gene calls:
    if include_AA_counts:
        if split_names:
            AA_counts_dict = c.get_AA_counts_dict(split_names=split_names)
        else:
            AA_counts_dict = c.get_AA_counts_dict()

        info_dict['AA_counts'] = AA_counts_dict['AA_counts']
        info_dict['total_AAs'] = AA_counts_dict['total_AAs']


    missing_keys = [key for key in constants.essential_genome_info if key not in info_dict]
    if len(missing_keys):
        raise ConfigError, "We have a big problem. I am reporting from get_contigs_db_info_dict. This function must\
                            produce a dictionary that meets the requirements defined in the constants module of anvi'o\
                            for 'essential genome info'. But when I look at the resulting dictionary this funciton is\
                            about to return, I can see it is missing some stuff :/ This is not a user error, but it needs\
                            the attention of an anvi'o developer. Here are the keys that should have been in the results\
                            but missing: '%s'" % (', '.join(missing_keys))

    return info_dict


class AdHocRunGenerator:
    """From a matrix file to full-blown anvi'o interface.

       This is a class to take in a view data matrix at minimum, and create all
       necessary files for an anvi'o interactive interface call in manual mode."""

    def __init__(self, view_data_path, run=run, progress=progress):
        self.run = run
        self.progress = progress

        self.view_data_path = view_data_path

        self.tree_file_path = None
        self.matrix_data_for_clustering = None
        self.additional_view_data_file_path = None

        self.samples_info_file_path = None
        self.samples_order_file_path = None

        # for clustering
        self.distance = None
        self.linkage = None

        self.output_directory = os.path.abspath('./ad-hoc-anvio-run-directory')
        self.delete_output_directory_if_exists = False

        self.sanity_checked = False


    def sanity_check(self):
        self.distance = self.distance or constants.distance_metric_default
        self.linkage = self.linkage or constants.linkage_method_default

        clustering.is_distance_and_linkage_compatible(self.distance, self.linkage)

        filesnpaths.is_file_tab_delimited(self.view_data_path)
        if self.tree_file_path:
            filesnpaths.is_proper_newick(self.tree_file_path)

        self.check_output_directory()

        new_view_data_path = self.get_output_file_path('view_data.txt')
        shutil.copyfile(self.view_data_path, new_view_data_path)
        self.view_data_path = new_view_data_path

        if self.tree_file_path:
            new_tree_path = self.get_output_file_path('tree.txt')
            shutil.copyfile(self.tree_file_path, new_tree_path)
            self.tree_file_path = new_tree_path

        if self.additional_view_data_file_path:
            new_additional_view_data_file_path = self.get_output_file_path('additional_view_data.txt')
            shutil.copyfile(self.additional_view_data_file_path, new_additional_view_data_file_path)
            self.additional_view_data_file_path = new_additional_view_data_file_path

        if self.samples_info_file_path:
            new_samples_info_file_path = self.get_output_file_path('anvio_samples_info.txt')
            shutil.copyfile(self.samples_info_file_path, new_samples_info_file_path)
            self.samples_info_file_path = new_samples_info_file_path


        self.sanity_checked = True


    def is_good_to_go(self):
        if not self.sanity_checked:
            raise ConfigError, "AdHocRunGenerator :: You gotta be nice, and run sanity check first :/"


    def get_output_file_path(self, file_name):
        return os.path.join(self.output_directory, file_name)


    def check_output_directory(self):
        if os.path.exists(self.output_directory) and not self.delete_output_directory_if_exists:
            raise ConfigError, "AdHocRunGenerator will not work with an existing directory. Please provide a new\
                                path, or use the bool member 'delete_output_directory_if_exists' to overwrite\
                                any existing directory."

        filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=self.delete_output_directory_if_exists)


    def generate(self):
        self.sanity_check()

        if not self.tree_file_path:
            self.gen_clustering_of_view_data()

        self.gen_samples_db()

        self.run.info("Ad hoc anvi'o run files", self.output_directory)


    def gen_clustering_of_view_data(self):
        self.is_good_to_go()

        self.progress.new('Hierarchical clustering')
        self.progress.update('..')

        self.tree_file_path = self.get_output_file_path('tree.txt')

        if self.matrix_data_for_clustering:
            clustering.get_newick_tree_data(self.matrix_data_for_clustering, self.tree_file_path, distance = self.distance, linkage=self.linkage)
        else:
            clustering.get_newick_tree_data(self.view_data_path, self.tree_file_path, distance = self.distance, linkage=self.linkage)

        self.progress.end()

        self.run.info('Tree', self.tree_file_path)


    def gen_samples_order_file(self, data_file_path):
        self.progress.new('Hierarchical clustering of the (transposed) view data')
        self.progress.update('..')

        newick = clustering.get_newick_tree_data(data_file_path, transpose=True, distance = self.distance, linkage=self.linkage)

        samples_order_file_path = self.get_output_file_path('anvio-samples-order.txt')
        samples_order = open(samples_order_file_path, 'w')
        samples_order.write('attributes\tbasic\tnewick\n')
        samples_order.write('view_data\t\t%s\n' % newick)
        samples_order.close()

        self.progress.end()

        self.run.info("Anvi'o samples order", samples_order_file_path)

        return samples_order_file_path



    def gen_samples_db(self):
        if not self.samples_order_file_path:
            self.samples_order_file_path = self.gen_samples_order_file(self.view_data_path)

        samples_db_output_path = self.get_output_file_path('samples.db')
        s = dbops.SamplesInformationDatabase(samples_db_output_path, run=self.run, progress=self.progress, quiet=True)
        s.create(self.samples_info_file_path, self.samples_order_file_path)

