# coding: utf-8
# pylint: disable=line-too-long
"""Summarizes information for a collection.

It also gives access to bin data that may be useful. For instance, did you know
could totally do this to access to gene coverage and detection dicts for a given
bin in a given collection:

import anvio.summarizer as summarizer

    >>> class Args: None
    >>> args = Args()
    >>> args.profile_db = 'SAMPLES-MERGED/PROFILE.db'
    >>> args.contigs_db = 'CONTIGS.db'
    >>> args.collection_name = "CONCOCT"

    >>> summary = summarizer.ProfileSummarizer(args)
    >>> summary.init()
    >>> _bin = summarizer.Bin(summary, 'Bin_1')

    # now you have these:
    >>> _bin.gene_coverages
    >>> _bin.gene_detection

"""

import os
import sys
import gzip
import glob
import numpy
import shutil
import hashlib
import mistune
import argparse
import textwrap
import pandas as pd

from collections import Counter

import anvio
import anvio.dbops as dbops
import anvio.utils as utils
import anvio.hmmops as hmmops
import anvio.sequence as seqlib
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths
import anvio.ccollections as ccollections
import anvio.completeness as completeness

from anvio.errors import ConfigError
from anvio.dbops import DatabasesMetaclass, ContigsSuperclass, PanSuperclass
from anvio.hmmops import SequencesForHMMHits
from anvio.summaryhtml import SummaryHTMLOutput, humanize_n, pretty
from anvio.tables.miscdata import TableForLayerAdditionalData, MiscDataTableFactory


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
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
        self.list_bins = None
        self.debug = None
        self.quick_summary = False
        self.init_gene_coverages = False
        self.skip_check_collection_name = False
        self.skip_init_functions = False
        self.cog_data_dir = None
        self.output_dir = filesnpaths.get_temp_directory_path()
        self.report_aa_seqs_for_gene_calls = False


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

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        if A('list_collections'):
            self.collections.list_collections()
            sys.exit()

        self.collection_name = A('collection_name')

        if A('list_bins'):
            if not self.collection_name:
                raise ConfigError("It may come across as a surprise, but you really can't list bins in a collection\
                                   without providing a collection name. Bioinformatics. Never understands what you\
                                   need and all :/")
            self.collections.list_bins_in_collection(collection_name=self.collection_name)
            sys.exit()

        self.skip_check_collection_name = A('skip_check_collection_name')
        self.skip_init_functions = A('skip_init_functions')
        self.init_gene_coverages = A('init_gene_coverages')
        self.output_directory = A('output_dir')
        self.quick = A('quick_summary')
        self.debug = A('debug')
        self.taxonomic_level = A('taxonomic_level') or 't_genus'
        self.cog_data_dir = A('cog_data_dir')
        self.report_aa_seqs_for_gene_calls = A('report_aa_seqs_for_gene_calls')
        self.delete_output_directory_if_exists = False if A('delete_output_directory_if_exists') == None else A('delete_output_directory_if_exists')
        self.just_do_it = A('just_do_it')

        if not self.lazy_init:
            self.sanity_check()

        if self.output_directory:
            self.output_directory = filesnpaths.check_output_directory(self.output_directory, ok_if_exists=self.delete_output_directory_if_exists or self.just_do_it)
            filesnpaths.gen_output_directory(self.output_directory, delete_if_exists=self.delete_output_directory_if_exists or self.just_do_it)
        else:
            self.output_directory = "SUMMARY"


    def report_misc_data_files(self, target_table='layers'):
        if target_table == 'layer_orders':
            raise ConfigError("Report misc data files do not know how to work with layer orders yet :/")

        run_obj = terminal.Run(verbose=False)

        db_path = self.pan_db_path if self.summary_type == 'pan' else self.profile_db_path
        additional_data = MiscDataTableFactory(argparse.Namespace(pan_or_profile_db=db_path, target_data_table=target_table), r=run_obj)

        data_groups, data_dict = additional_data.get_all()
        for data_group in data_groups:
            output_file_obj = self.get_output_file_handle(sub_directory='misc_data_%s' % target_table, prefix='%s.txt' % data_group)
            output_file_path = output_file_obj.name
            output_file_obj.close()
            MiscDataTableFactory(argparse.Namespace(pan_or_profile_db=db_path, target_data_table=target_table, target_data_group=data_group), r=run_obj).export(output_file_path=output_file_path)

        if 'misc_data' not in self.summary:
            self.summary['misc_data'] = {}

        data_groups_reported = list(data_groups.keys())
        self.summary['misc_data'][target_table] = data_groups
        self.run.info('Misc data reported for %s' % target_table, ', '.join(data_groups_reported) if data_groups else 'None', nl_before=1)


    def sanity_check(self):
        if not self.skip_check_collection_name:
            if not self.collection_name:
                raise ConfigError("You must specify a collection id :/")

            if self.collection_name not in self.collections.collections_dict:
                raise ConfigError("%s is not a valid collection ID. See a list of available ones with '--list-collections' flag" % self.collection_name)


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
            raise ConfigError('get_output_file_handle: well, this file already exists: "%s"' % file_path)

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
    def __init__(self, args=None, lazy_init=False, r=run, p=progress):
        self.summary_type = 'pan'
        self.debug = False
        self.quick = False
        self.pan_db_path = None
        self.lazy_init = lazy_init
        self.output_directory = None
        self.genomes_storage_path = None

        self.run = r
        self.progress = p

        PanSuperclass.__init__(self, args, self.run, self.progress)
        if not self.genomes_storage_is_available:
            raise ConfigError("No genomes storage no summary. Yes. Very simple stuff.")

        SummarizerSuperClass.__init__(self, args, self.run, self.progress)

        # init gene clusters and functions from Pan super.
        self.init_gene_clusters()

        # init items additional data.
        self.init_items_additional_data()

        if not self.skip_init_functions:
            self.init_gene_clusters_functions()

        # see if COG functions or categories are available
        self.cog_functions_are_called = 'COG_FUNCTION' in self.gene_clusters_function_sources
        self.cog_categories_are_called = 'COG_CATEGORY' in self.gene_clusters_function_sources


    def get_occurrence_of_functions_in_pangenome(self, gene_clusters_functions_summary_dict):
        """
            For each function we will create a fake merged gc, with an occurrence vector
            which is the "or" product of all gcs that match this function.

            We also keep track of all the GCs annotated with the function.
        """
        occurrence_of_functions_in_pangenome_dict = {}

        self.progress.new('Computing presence absence for functions in genomes')
        self.progress.update('Creating a dictionary')

        for gene_cluster_id in gene_clusters_functions_summary_dict:
            gene_cluster_function = gene_clusters_functions_summary_dict[gene_cluster_id]['gene_cluster_function']
            if gene_cluster_function:
                if gene_cluster_function not in occurrence_of_functions_in_pangenome_dict:
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function] = {}
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function]['gene_clusters_ids'] = []
                    occurrence_of_functions_in_pangenome_dict[gene_cluster_function]['occurrence'] = None
                occurrence_of_functions_in_pangenome_dict[gene_cluster_function]['gene_clusters_ids'].append(gene_cluster_id)

        from anvio.dbops import PanDatabase
        pan_db = PanDatabase(self.pan_db_path)

        gene_cluster_presence_absence_dataframe = pd.DataFrame.from_dict(
                                                    pan_db.db.get_table_as_dict('gene_cluster_presence_absence'),
                                                    orient='index')

        self.progress.update('Merging presence/absence of gene clusters with the same function')

        D = {}
        for gene_cluster_function in occurrence_of_functions_in_pangenome_dict:
            v = None
            for gene_cluster_id in occurrence_of_functions_in_pangenome_dict[gene_cluster_function]['gene_clusters_ids']:
                if v is None:
                    v = gene_cluster_presence_absence_dataframe.loc[gene_cluster_id, ].astype(bool)
                else:
                    v = numpy.logical_or(v, gene_cluster_presence_absence_dataframe.loc[gene_cluster_id, ])
            D[gene_cluster_function] = {}
            for genome in v.index:
                D[gene_cluster_function][genome] = v[genome]

        self.progress.end()

        return pd.DataFrame.from_dict(D), occurrence_of_functions_in_pangenome_dict


    def functional_enrichment_stats(self):
        """
            Compute the functional enrichment stats for a pangenome

            returns the enrichment_dict, which has the form:
                value = enrichment_dict[category_name][function_name][value_name]

            If self.args.output_file_path exists then an output file is created.

            To learn more about how this works refer to the docummentation:
                anvi-script-get-enriched-functions-per-pan-group -h
        """

        A = lambda x: self.args.__dict__[x] if x in self.args.__dict__ else None
        output_file_path = A('output_file')
        category_variable = A('category_variable')
        functional_annotation_source = A('annotation_source')
        list_functional_annotation_sources = A('list_annotation_sources')
        min_function_enrichment = A('min_function_enrichment')
        core_threshold = A('core_threshold')
        fdr = A('false_detection_rate')
        min_portion_occurrence_of_function_in_group = A('min_portion_occurrence_of_function_in_group')
        functional_occurrence_table_output = A('functional_occurrence_table_output')

        if output_file_path:
            filesnpaths.is_output_file_writable(output_file_path)

        if functional_occurrence_table_output:
            filesnpaths.is_output_file_writable(functional_occurrence_table_output)

        if not self.functions_initialized:
            raise ConfigError("For some reason funtions are not initialized for this pan class instance. We\
                               can't summarize functional enrichment stats without that :/")

        if not len(self.gene_clusters_functions_dict):
            raise ConfigError("The gene clusters functions dict seems to be empty. We assume this error makes\
                               zero sense to you, and it probably will not help you to know that it also makes\
                               zero sense to anvi'o too :/ Maybe you forgot to provide a genomes storage?")

        if not category_variable:
            raise ConfigError("For this to work, you must provide a category variable .. and it better be in\
                               the misc additional layer data table, too. If you don't have any idea what is\
                               available, try `anvi-show-misc-data`.")

        if list_functional_annotation_sources:
            self.run.info('Available functional annotation sources', ', '.join(self.gene_clusters_function_sources))
            sys.exit()

        if not functional_annotation_source:
            raise ConfigError("You haven't provided a functional annotation source to make sense of functional\
                               enrichment stats as defined by the categorical variable %s. These are the functions\
                               that are available, so pick one: %s." % (category_variable, ', '.join(self.gene_clusters_function_sources)))

        if functional_annotation_source not in self.gene_clusters_function_sources:
            raise ConfigError("Your favorite functional annotation source '%s' does not seem to be among one of the sources\
                               that are available to you. Here are the ones you should choose from: %s." % (functional_annotation_source, ', '.join(self.gene_clusters_function_sources)))

        keys, categories_dict = TableForLayerAdditionalData(argparse.Namespace(pan_db=self.pan_db_path)).get(additional_data_keys_requested=[category_variable])

        values_that_are_not_none = [s for s in categories_dict.values() if s[category_variable] is not None]
        if not values_that_are_not_none:
            raise ConfigError("The variable '%s' contains only values of type None,\
                               this is probably a mistake, surely you didn't mean to provide an empty category.\
                               Do you think this is a mistake on our part? Let us know." % \
                                                                    category_variable)
        type_category_variable = type(values_that_are_not_none[0][category_variable])
        if type_category_variable != str:
            raise ConfigError("The variable '%s' does not seem to resemble anything that could be a category.\
                               Anvi'o expects these variables to be of type string, yet yours is type %s :/\
                               Do you think this is a mistake on our part? Let us know." % \
                                                                    (category_variable, type_category_variable))

        gene_clusters_functions_summary_dict = self.get_gene_clusters_functions_summary_dict(functional_annotation_source)

        self.run.info('Category', category_variable)
        self.run.info('Functional annotation source', functional_annotation_source)

        occurrence_of_functions_in_pangenome_dataframe, occurrence_of_functions_in_pangenome_dict = self.get_occurrence_of_functions_in_pangenome(gene_clusters_functions_summary_dict)

        if functional_occurrence_table_output:
            occurrence_of_functions_in_pangenome_dataframe.astype(int).transpose().to_csv(functional_occurrence_table_output, sep='\t')
            self.run.info('Presence/absence of functions summary:', functional_occurrence_table_output)

        self.progress.new('Functional enrichment analysis')
        self.progress.update('Creating a dictionary')

        # get a list of unique function names
        functions_names = set(occurrence_of_functions_in_pangenome_dataframe.columns)

        # the total occurrence of functions in all categories (it is important to this before adding the category column)
        total_occurrence_of_functions = occurrence_of_functions_in_pangenome_dataframe.sum()

        # add a category column to the dataframe
        occurrence_of_functions_in_pangenome_dataframe['category'] = occurrence_of_functions_in_pangenome_dataframe.index.map(lambda x: categories_dict[x][category_variable])

        # the sum of occurrences of each function in each category
        functions_in_categories = occurrence_of_functions_in_pangenome_dataframe.groupby('category').sum()

        # unique names of categories
        categories = set([categories_dict[g][category_variable] for g in categories_dict.keys() if\
                            categories_dict[g][category_variable] is not None])

        categories_to_genomes_dict = {}
        for c in categories:
            categories_to_genomes_dict[c] = set([genome for genome in categories_dict.keys() if categories_dict[genome][category_variable] == c])
        number_of_genomes = len(categories_dict.keys())

        enrichment_dict = {}
        z_test_p_values = {}
        for c in categories:
            self.progress.update("Working on category '%s'" % c)
            group_size = len(categories_to_genomes_dict[c])
            outgroup_size = number_of_genomes - group_size
            z_test_p_values[c] = []

            for f in functions_names:
                occurrence_in_group = functions_in_categories.loc[c, f]
                occurrence_outside_of_group = (total_occurrence_of_functions[f] - functions_in_categories.loc[c, f])
                portion_occurrence_in_group = occurrence_in_group / group_size
                portion_occurrence_outside_of_group = occurrence_outside_of_group / (number_of_genomes - group_size)
                enrichment, p_value = utils.get_two_sample_z_test_statistic(portion_occurrence_in_group, \
                                                           portion_occurrence_outside_of_group, \
                                                           group_size, \
                                                           outgroup_size)

                z_test_p_values[c].append(p_value)

                if c not in enrichment_dict:
                    enrichment_dict[c] = {}

                if f not in enrichment_dict[c]:
                    enrichment_dict[c][f] = {}

                enrichment_dict[c][f]["enrichment_score"] = enrichment
                enrichment_dict[c][f]["p_value"] = p_value
                enrichment_dict[c][f]["portion_occurrence_in_group"] = portion_occurrence_in_group
                enrichment_dict[c][f]["portion_occurrence_outside_of_group"] = portion_occurrence_outside_of_group
                enrichment_dict[c][f]["occurrence_in_group"] = occurrence_in_group
                enrichment_dict[c][f]["occurrence_outside_of_group"] = occurrence_outside_of_group
                enrichment_dict[c][f]["gene_clusters_ids"] = occurrence_of_functions_in_pangenome_dict[f]["gene_clusters_ids"]
                enrichment_dict[c][f]["core_in_group"] = False
                enrichment_dict[c][f]["core"] = False
                if enrichment_dict[c][f]["portion_occurrence_in_group"] >= core_threshold:
                    enrichment_dict[c][f]["core_in_group"] = True
                    if enrichment_dict[c][f]["portion_occurrence_outside_of_group"] >= core_threshold:
                        enrichment_dict[c][f]["core"] = True

        import statsmodels.stats.multitest as multitest
        for c in enrichment_dict:

            self.progress.update("Working on statistics for category '%s'" % c)
            # correction for multiple comparrisons
            reject, corrected_p_values_z_test, foo1, foo2 = multitest.multipletests(z_test_p_values[c], method='fdr_bh', alpha=fdr)

            i = 0
            if len(corrected_p_values_z_test) != len(enrichment_dict[c].keys()):
                raise ConfigError('This should never happen, contact Alon Shaiber now.')
            for f in enrichment_dict[c]:
                enrichment_dict[c][f]['corrected_p_value'] = corrected_p_values_z_test[i]
                i += 1

        if output_file_path:
            self.progress.update('Generating the output file')
            enrichment_data_frame = self.get_enrichment_dict_as_dataframe(enrichment_dict, functional_annotation_source)
            if min_function_enrichment > 0 or min_portion_occurrence_of_function_in_group > 0:
                enrichment_data_frame = enrichment_data_frame[enrichment_data_frame['enrichment_score'] > min_function_enrichment]
                enrichment_data_frame = enrichment_data_frame[enrichment_data_frame['portion_occurrence_in_group'] > min_portion_occurrence_of_function_in_group]

            # sort according to enrichment
            enrichment_data_frame.sort_values(by=['category', 'enrichment_score'], axis=0, ascending=False, inplace=True)

            enrichment_data_frame.to_csv(output_file_path, sep='\t', index=False, float_format='%.2f')

        self.progress.end()

        if output_file_path:
            self.run.info('Functions enrichment summary', output_file_path)

        return enrichment_dict


    def get_enrichment_dict_as_dataframe(self, enrichment_dict, functional_annotation_source):
        # convert dictionary to pandas
        # we can't use pandas from_dict because it is meant for dict of dicts (i.e. tow levels)
        # and we have a dict of dicts of dicts (three levels).
        # so we first convert it to a dict of dicts and then convert to pandas
        # because this is faster than alternatives
        i = 0
        D = {}
        for c in enrichment_dict:
            for f in enrichment_dict[c]:
                D[i] = {}
                D[i]['category'] = c
                D[i][functional_annotation_source] = f
                for key, value in enrichment_dict[c][f].items():
                    try:
                        # if there is a sequence of values
                        # merge them with commas for nicer printing
                        D[i][key] = ', '.join(iter(value))
                    except:
                        # if it is not a sequnce, it is a single value
                        D[i][key] = value
                i += 1
        enrichment_data_frame = pd.DataFrame.from_dict(D, orient='index')

        return enrichment_data_frame


    def process(self):
        # init profile data for colletion.
        collection_dict, bins_info_dict = self.init_collection_profile(self.collection_name)

        # let bin names known to all
        bin_ids = list(self.collection_profile.keys())

        genome_names = ', '.join(list(self.gene_clusters.values())[0].keys())

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
                            'functions_available': True if len(self.gene_clusters_function_sources) else False,
                            'function_sources': self.gene_clusters_function_sources},
                'percent_of_genes_collection': 0.0,
                'genome_names': genome_names
        }

        # I am not sure whether this is the best place to do this,
        self.summary['basics_pretty'] = { \
                'pan': [('Created on', self.p_meta['creation_date']),
                        ('Version', anvio.__pan__version__),
                        ('Number of genes', pretty(int(self.p_meta['num_genes_in_gene_clusters']))),
                        ('Number of gene clusters', pretty(int(self.p_meta['num_gene_clusters']))),
                        ('Partial genes excluded', 'Yes' if self.p_meta['exclude_partial_gene_calls'] else 'No'),
                        ('Minbit parameter', self.p_meta['minbit']),
                        ('Gene cluster min occurrence parameter', pretty(int(self.p_meta['gene_cluster_min_occurrence']))),
                        ('MCL inflation parameter', self.p_meta['mcl_inflation']),
                        ('NCBI blastp or DIAMOND?', 'NCBI blastp' if self.p_meta['use_ncbi_blast'] else ('DIAMOND (and it was %s)' % ('sensitive' if self.p_meta['diamond_sensitive'] else 'not sensitive'))),
                        ('Number of genomes used', pretty(int(self.p_meta['num_genomes']))),
                        ('Items aditional data keys', '--' if not self.items_additional_data_keys else ', '.join(self.items_additional_data_keys))],

                'genomes': [('Created on', 'Storage DB knows nothing :('),
                            ('Version', anvio.__genomes_storage_version__),
                            ('Number of genomes described', pretty(self.genomes_storage.num_genomes)),
                            ('Functional annotation', 'Available' if len(self.gene_clusters_function_sources) else 'Not available :/'),
                            ('Functional annotation sources', '--' if not len(self.gene_clusters_function_sources) else ', '.join(self.gene_clusters_function_sources))],
        }

        self.summary['files'] = {}
        self.summary['collection_profile'] = self.collection_profile # reminder; collection_profile comes from the superclass!

        self.generate_gene_clusters_file(collection_dict)

        self.report_misc_data_files(target_table='layers')
        self.report_misc_data_files(target_table='items')

        if self.debug:
            import json
            print(json.dumps(self.summary, sort_keys=True, indent=4))

        self.index_html = SummaryHTMLOutput(self.summary, r=self.run, p=self.progress).generate(quick=self.quick)


    def generate_gene_clusters_file(self, collection_dict, compress_output=True):
        """Generates the gene summary file"""

        # generate a dict of gene cluster ~ bin id relationships
        gene_cluster_name_to_bin_name= dict(list(zip(self.gene_clusters_in_pan_db_but_not_binned, [None] * len(self.gene_clusters_in_pan_db_but_not_binned))))
        for bin_id in collection_dict:
            for gene_cluster_name in collection_dict[bin_id]:
                gene_cluster_name_to_bin_name[gene_cluster_name] = bin_id

        ###############################################
        # generate an output file for gene clusters.
        ###############################################
        output_file_obj = self.get_output_file_handle(prefix='gene_clusters_summary.txt', compress_output=compress_output, add_project_name=True)

        # standard headers
        header = ['unique_id', 'gene_cluster_id', 'bin_name', 'genome_name', 'gene_callers_id']

        # extend the header with items additional data keys
        for items_additional_data_key in self.items_additional_data_keys:
            header.append(items_additional_data_key)

        # extend the header with functions if there are any
        for function_source in self.gene_clusters_function_sources:
            if self.quick:
                header.append(function_source + '_ACC')
            else:
                header.append(function_source + '_ACC')
                header.append(function_source)

        # if this is not a quick summary, have AA sequences in the output
        AA_sequences = None
        if not self.quick:
            header.append('aa_sequence')
            AA_sequences = self.get_sequences_for_gene_clusters(gene_cluster_names=self.gene_cluster_names)

        # write the header
        output_file_obj.write(('\t'.join(header) + '\n').encode('utf-8'))

        self.progress.new('Gene clusters summary file')
        self.progress.update('...')

        # uber loop for the file content
        unique_id = 1
        for gene_cluster_name in self.gene_clusters:
            for genome_name in self.gene_clusters[gene_cluster_name]:
                for gene_caller_id in self.gene_clusters[gene_cluster_name][genome_name]:
                    entry = [unique_id, gene_cluster_name, gene_cluster_name_to_bin_name[gene_cluster_name], genome_name, gene_caller_id]

                    # populate the entry with item aditional data
                    for items_additional_data_key in self.items_additional_data_keys:
                        if gene_cluster_name in self.items_additional_data_dict:
                            entry.append(self.items_additional_data_dict[gene_cluster_name][items_additional_data_key])
                        else:
                            entry.append('')

                    # populate the entry with functions.
                    for function_source in self.gene_clusters_function_sources:
                        annotations_dict = self.gene_clusters_functions_dict[gene_cluster_name][genome_name][gene_caller_id]
                        if function_source in annotations_dict:
                            annotation_blob = self.gene_clusters_functions_dict[gene_cluster_name][genome_name][gene_caller_id][function_source]

                            # FIXME: this is an artifact from Py2 to Py3 swtich. DBs generated in Py2 and used from Py3 will have type
                            # bytes for annotation_blob. so we will convert them to str.. If the db is generated with Py3, there is no
                            # such problem, so we can remove this extra step around July 2017.
                            if isinstance(annotation_blob, bytes):
                                annotation_blob = annotation_blob.decode('utf-8')

                            accessions, annotations = [l.split('!!!') for l in annotation_blob.split("|||")]
                            entry.append('|'.join(accessions))
                            entry.append('|'.join(annotations))
                        else:
                            entry.append('')
                            entry.append('')

                    if not self.quick:
                        entry.append(AA_sequences[gene_cluster_name][genome_name][gene_caller_id])

                    output_file_obj.write(('\t'.join([str(e) if e not in [None, 'UNKNOWN'] else '' for e in entry]) + '\n').encode('utf-8'))
                    unique_id += 1


        # we're done here.
        output_file_obj.close()

        self.progress.end()


class SAAVsAndProteinStructuresSummary:
    """A class to make sense of the SAAV Structure outputs."""

    def __init__(self, args=None, r=run, p=progress):
        self.run = run
        self.progress = progress

        self.args = args

        self.summary = {}
        self.summary_type = 'saav'

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.input_directory = A('input_dir')
        self.output_directory = A('output_dir')
        self.soft_link_images = A('soft_link_images')
        self.perspectives = A('perspectives')

        self.genes_file_path = A('genes')
        self.samples_file_path = A('samples')

        # dicts that will be recovered from input files
        self.genes = None
        self.samples = None
        self.views = None

        # this will be recovered by traversing the input directory
        self.genes_info = {}

        # dicts that will be populated by the init function
        self.by_view = {}
        self.legends = {}
        self.samples_per_view = {}

        self.initialized = False
        self.sanity_checked = False


    def sanity_check(self):
        self.input_directory = os.path.abspath(self.input_directory)
        filesnpaths.is_file_exists(self.input_directory)

        # this fallback code could have been inside init instead of sanity check function
        # but calling os.path.join can raise exception if input_directory is None, it should checked first.
        self.genes_file_path = self.genes_file_path or os.path.join(self.input_directory, '.gene_list.txt')
        self.samples_file_path = self.samples_file_path or os.path.join(self.input_directory, '.sample_groups.txt')

        if not filesnpaths.is_file_exists(self.genes_file_path, dont_raise=True):
            raise ConfigError("Anvi'o could not find gene list file '%s'. If you did not provided any as a parameter \
                               anvi'o looks for '.gene_list.txt' in input directory." % self.genes_file_path)

        if not filesnpaths.is_file_exists(self.samples_file_path, dont_raise=True):
            raise ConfigError("Anvi'o could not find sample groups file '%s'. If you did not provided any as a parameter \
                               anvi'o looks for '.sample_groups.txt' in input directory." % self.samples_file_path)

        if not self.output_directory or not self.input_directory:
            raise ConfigError("You must declare both input and output directories.")

        if self.output_directory == self.input_directory:
            raise ConfigError("The input and the output directories can't be the same.")

        if self.contigs_db_path:
            utils.is_contigs_db(self.contigs_db_path)

        self.output_directory = filesnpaths.check_output_directory(self.output_directory)
        filesnpaths.gen_output_directory(self.output_directory)
        filesnpaths.gen_output_directory(os.path.join(self.output_directory, 'images'))

        self.run.info('Contigs DB found', self.contigs_db_path is not None, mc='green' if self.contigs_db_path else 'red')
        self.run.info('Input directory', self.input_directory)
        self.run.info('Output directory', self.output_directory)

        self.sanity_checked = True


    def process_input(self):
        # FIXME: Assume its a single-column file. If it isn't, assume its a multi-column file
        try:
            self.gene_list = list(utils.get_column_data_from_TAB_delim_file(self.genes_file_path, column_indices=[0], expected_number_of_fields=1).values())[0][1:]
            self.genes = {}
            for gene in self.gene_list:
                self.genes[gene] = {}
        except:
            self.genes = utils.get_TAB_delimited_file_as_dictionary(self.genes_file_path)

        # FIXME: Assume its a single-column file. If it isn't, assume its a multi-column file
        try:
            self.sample_list = list(utils.get_column_data_from_TAB_delim_file(self.samples_file_path, column_indices=[0], expected_number_of_fields=1).values())[0][1:]
            self.samples = {}
            for sample in self.sample_list:
                self.samples[sample] = {}
        except:
            self.samples = utils.get_TAB_delimited_file_as_dictionary(self.samples_file_path)

        self.views = utils.get_columns_of_TAB_delim_file(self.samples_file_path)

        # add a samples view
        self.views.append("samples")
        for sample in self.samples:
            self.samples[sample]["samples"] = "All"


    def populate_genes_info_dict(self):
        if not self.contigs_db_path:
            # go populate yourself
            return

        contigs_db = dbops.ContigsSuperclass(self.args)

        # let's make sure all genes are here
        for gene_id in map(int, self.genes):
            if gene_id not in contigs_db.genes_in_contigs_dict:
                raise ConfigError("Gene caller id %d is not in your contigs db. You must have\
                                   provided an irrelevant contigs database to your run.")

        contigs_db.init_functions()

        if len(contigs_db.gene_function_calls_dict):
            self.run.info_single("Good news, anvi'o found functional annotations in your\
                                  contigs database and will use them to display in the HMTL\
                                  output. We hope we are not missing anything: %s." % \
                                            ', '.join(contigs_db.gene_function_call_sources),
                                            nl_before=1, nl_after=1)

        # we will either of these to show on the header:
        preferred_annotation_sources = ['Pfam', 'COG_FUNCTION', 'TIGRFAM']
        annotation_source = None
        for f in preferred_annotation_sources:
            if f in contigs_db.gene_function_call_sources:
                annotation_source = f
                break

        if annotation_source:
            self.run.info_single("%s will be used to show summaries in section headers." % \
                                            annotation_source, nl_after=1)
        else:
            self.run.info_single("Congratulations, none of the preferred functional annotation \
                                  sources (%s) were in your contigs database :( although you \
                                  have functions in the contigs database, your gene headers\
                                  will look ugly. That's OK." % \
                                            ', '.join(preferred_annotation_sources), nl_after=1)

        for gene_id in map(int, self.genes):
            f = contigs_db.gene_function_calls_dict[gene_id]
            self.genes_info[gene_id] = {}

            if annotation_source:
                if f[annotation_source]:
                    self.genes_info[gene_id]['function'] = f[annotation_source][1]
                    self.genes_info[gene_id]['accession'] = f[annotation_source][0]
                else:
                    self.genes_info[gene_id]['function'] = 'Unknown function'
                    self.genes_info[gene_id]['accession'] = None

            self.genes_info[gene_id]['functions'] = f


    def init(self):
        self.sanity_check()

        self.process_input()
        self.populate_genes_info_dict()

        perspectives_found = [os.path.basename(d) for d in glob.glob('%s/*' % self.input_directory) if os.path.isdir(d)]

        if self.perspectives:
            self.perspectives = [p.strip() for p in self.perspectives.split(',') if p.strip()]
            for perspective in self.perspectives:
                if perspective not in perspectives_found:
                    raise ConfigError("The perspectie you requested ('%s') is not one of the available perspectives in\
                                       this SAAVs structure output. These are the ones that are avilable: %s" % \
                                                                          (perspective, ', '.join(perspectives_found)))
        else:
            self.perspectives = perspectives_found

        self.run.info('Num genes', len(self.genes))
        self.run.info('Num samples', len(self.samples))
        self.run.info('Views', ', '.join(self.views))
        self.run.info('Perspectives', ', '.join(self.perspectives))
        self.run.info('Images are soft linked', self.soft_link_images, mc='red' if self.soft_link_images else 'yellow')

        if(len(self.genes)) > 50:
            self.run.warning('You seem to have a lot of genes to process. Nice. The output may be quite large,\
                              just so you know :/')

        self.summary['meta'] = {'summary_type': self.summary_type,
                                'output_directory': self.output_directory,
                                'images_soft_linked': self.soft_link_images,
                                'anvio_version': anvio.__version__}

        # populate dicts
        self.populate_samples_per_view_dict()
        self.populate_by_view_dict()

        views_and_variables = {}
        for view in self.samples_per_view:
            views_and_variables[view] = sorted(self.samples_per_view[view].keys())

        self.summary['data'] = {'gene_names': sorted(list(self.genes.keys())),
                                'samples': self.samples,
                                'by_view': self.by_view,
                                'views_and_variables': views_and_variables,
                                'views': sorted(self.samples_per_view.keys()),
                                'perspectives': sorted(self.perspectives),
                                'genes': self.genes,
                                'genes_info': self.genes_info,
                                'samples_per_view': self.samples_per_view,
                                'legends': self.legends}

        self.initialized = True


    def populate_samples_per_view_dict(self):
        self.progress.new('Populating samples per view data')
        self.progress.update('...')

        self.samples_per_view = {}

        for view in self.views:
            self.progress.update('Working on "%s" ...' % view)
            self.samples_per_view[view] = {}
            for sample in self.samples:
                r = self.samples[sample][view]
                if r not in self.samples_per_view[view]:
                    self.samples_per_view[view][r] = []

                self.samples_per_view[view][r].append(sample)

        self.progress.end()


    def copy_image_and_return_path(self, variables=None):
        image_path_template = None
        if variables['image_type'] == 'sample':
            image_path_template = "%(input_directory)s/%(perspective)s/Images/%(gene)s/%(image_type)s_%(sample)s.pse.png"
        elif variables['image_type'] == 'merged':
            image_path_template = "%(input_directory)s/%(perspective)s/Images/%(gene)s/%(image_type)s_%(view)s_%(variable)s.pse.png"

        image_path = image_path_template % variables

        # if user wants a fully populated output directory, update the image_path variable
        if not self.soft_link_images:
            #new_image_path = 'images/%s_%s_%s.png' % (str(gene), sample, hashlib.sha1(image_path.encode('utf-8')).hexdigest())
            new_image_path = 'images/%s_%s_%s.png' % (variables["gene"], variables["sample"], hashlib.sha1(image_path.encode('utf-8')).hexdigest())
            shutil.copyfile(os.path.join(self.input_directory, image_path), os.path.join(self.output_directory, new_image_path))
            image_path = new_image_path

        return image_path


    def populate_by_view_dict(self):
        """This one connects the actual data and images."""

        self.progress.new('Populating views dict')

        gene_names = sorted(self.genes.keys())
        num_genes = len(gene_names)
        for index in range(0, num_genes):
            gene = gene_names[index]
            self.progress.update("gene '%s' (%d of %d) ..." % (str(gene), index + 1, num_genes))
            self.by_view[gene] = {}
            self.legends[gene] = {}
            for view in self.samples_per_view.keys():
                self.by_view[gene][view] = {}
                for perspective in self.perspectives:
                    if perspective not in self.legends[gene]:
                        self.legends[gene][perspective] = {}

                    self.by_view[gene][view][perspective] = {}
                    for variable in sorted(self.samples_per_view[view].keys()):
                        self.by_view[gene][view][perspective][variable] = {}
                        for sample in self.samples_per_view[view][variable]:
                            image_path = self.copy_image_and_return_path(variables= {'input_directory': self.input_directory,
                                                                                            'gene': str(gene),
                                                                                            'sample': sample,
                                                                                            'perspective': perspective,
                                                                                            'image_type': 'sample',
                                                                                            'view': view,
                                                                                            'variable': variable})

                            self.by_view[gene][view][perspective][variable][sample] = image_path

                        if view != "samples":

                            image_path = self.copy_image_and_return_path(variables= {'input_directory': self.input_directory,
                                                                                            'gene': str(gene),
                                                                                            'sample': sample,
                                                                                            'perspective': perspective,
                                                                                            'image_type': 'merged',
                                                                                            'view': view,
                                                                                            'variable': variable})

                            self.by_view[gene][view][perspective][variable]['__merged__'] = image_path

                    self.legends[gene][perspective] = self.get_legend_as_dict(gene, perspective)

        self.progress.end()


    def get_legend_as_dict(self, gene, perspective):
        """
        does a global legend exist? does it have less than 15 colors? if not, load gene-specific legend for each gene

        for single legends id is a gene_id, for global legend it is string constant "global"
        """

        color_legend_path_template = "%(input_directory)s/%(perspective)s/Legends/color/%(perspective)s_%(id)s_color_legend.txt"

        global_color_legend_path = color_legend_path_template % {'input_directory' : self.input_directory,
                                                                 'id'              : 'global',
                                                                 'perspective'     : perspective}

        if os.path.isfile(global_color_legend_path):
            global_legend_content = utils.get_TAB_delimited_file_as_dictionary(global_color_legend_path)
            if len(global_legend_content) <= 15:
                return global_legend_content

        return utils.get_TAB_delimited_file_as_dictionary(color_legend_path_template %
                                                                    {'input_directory': self.input_directory,
                                                                    'id': str(gene),
                                                                    'perspective': perspective})

    def process(self):
        if not self.initialized:
            self.init()

        if not self.sanity_checked:
            self.sanity_check()

        self.index_html = SummaryHTMLOutput(self.summary, r=self.run, p=self.progress).generate(quick=False)


class ProfileSummarizer(DatabasesMetaclass, SummarizerSuperClass):
    """Creates an über dictionary of 'summary' for anvi'o profiles."""
    def __init__(self, args=None, lazy_init=False, r=run, p=progress):
        self.args = args

        self.run = r
        self.progress = p
        self.lazy_init = lazy_init

        self.summary = {}
        self.summary_type = 'profile'
        self.debug = False
        self.quick = False
        self.profile_db_path = None
        self.contigs_db_path = None
        self.output_directory = None
        self.split_names_per_bin = None
        self.completeness_data_available = False
        self.gene_level_coverage_stats_available = False
        self.non_single_copy_gene_hmm_data_available = False

        DatabasesMetaclass.__init__(self, self.args, self.run, self.progress)
        SummarizerSuperClass.__init__(self, self.args, self.run, self.progress)

        # databases initiated, let's make sure we have gene covereges data avaialable.
        if self.gene_level_coverage_stats_dict:
            self.gene_level_coverage_stats_available = True

        self.init_splits_taxonomy(self.taxonomic_level)

        self.collection_dict = {}
        self.bins_info_dict = {}
        self.initialized = False


    def init(self):
        # init profile data for colletion.
        self.collection_dict, self.bins_info_dict = self.init_collection_profile(self.collection_name)

        # let bin names known to all
        self.bin_ids = list(self.collection_profile.keys())

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
                                'collection': self.bin_ids,
                                'num_bins': len(self.bin_ids),
                                'collection_name': self.collection_name,
                                'total_nts_in_collection': 0,
                                'num_contigs_in_collection': 0,
                                'anvio_version': __version__,
                                'profile': self.p_meta,
                                'contigs': self.a_meta,
                                'gene_level_coverage_stats_available': self.gene_level_coverage_stats_available,
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
                                                     ('Minimum contig length', pretty(self.p_meta['min_contig_length'])),
                                                     ('Number of contigs', pretty(int(self.p_meta['num_contigs']))),
                                                     ('Number of splits', pretty(int(self.p_meta['num_splits']))),
                                                     ('Total nucleotides', humanize_n(int(self.p_meta['total_length']))),
                                                     ('SNVs profiled', self.p_meta['SNVs_profiled']),
                                                     ('SCVs profiled', self.p_meta['SCVs_profiled']),
                                                    ],
                                         'contigs': [
                                                        ('Created on', self.p_meta['creation_date']),
                                                        ('Version', self.a_meta['version']),
                                                        ('Split length', pretty(int(self.a_meta['split_length']))),
                                                        ('Number of contigs', pretty(int(self.a_meta['num_contigs']))),
                                                        ('Number of splits', pretty(int(self.a_meta['num_splits']))),
                                                        ('Total nucleotides', humanize_n(int(self.a_meta['total_length']))),
                                                        ('K-mer size', self.a_meta['kmer_size']),
                                                        ('Genes are called', self.a_meta['genes_are_called']),
                                                        ('Splits consider gene calls', self.a_meta['splits_consider_gene_calls']),
                                                        ('Gene function sources', ', '.join(self.gene_function_call_sources) if self.gene_function_call_sources else 'None :('),
                                                    ],
                                        'description': mistune.markdown(self.p_meta['description']),
                                        }

        self.summary['max_shown_header_items'] = 10
        self.summary['slice_header_items_tmpl'] = '0:%d' % self.summary['max_shown_header_items']
        self.summary['num_not_shown_samples'] = len(self.p_meta['samples']) - self.summary['max_shown_header_items']
        self.summary['num_not_shown_hmm_items'] = dict([(hmm_search_source, len(self.hmm_sources_info[hmm_search_source]['genes']) - self.summary['max_shown_header_items']) for hmm_search_type, hmm_search_source in self.hmm_searches_header])

        self.summary['files'] = {}
        self.summary['misc_data'] = {}
        self.summary['collection'] = {}
        self.summary['collection_profile'] = self.collection_profile # reminder; collection_profile comes from ProfileSuperclass!
        self.summary['collection_profile_items'] = [] if not len(list(self.collection_profile.values())) else list(self.collection_profile.values())[0].keys()

        # add hmm items for each seach type:
        if self.non_single_copy_gene_hmm_data_available:
            self.summary['meta']['hmm_items'] = dict([(hmm_search_source, self.hmm_sources_info[hmm_search_source]['genes']) for hmm_search_type, hmm_search_source in self.hmm_searches_header])

        # yay
        self.initialized = True


    def process(self):
        if not self.initialized:
            self.init()

        if not self.output_directory:
            raise ConfigError("It seems the summarizer class have been inherited without an `output_directory` argument :/ Show stopper\
                               mistake stopped the show. Bye!")

        # summarize bins:
        for i in range(0, len(self.bin_ids)):
            bin_id = self.bin_ids[i]
            self.progress.new('[Processing "%s" (%d of %d)]' % (bin_id, i + 1, len(self.bin_ids)))
            bin = Bin(self, bin_id, self.run, self.progress)
            bin.output_directory = os.path.join(self.output_directory, 'bin_by_bin', bin_id)
            bin.bin_profile = self.collection_profile[bin_id]

            self.summary['collection'][bin_id] = bin.create()
            self.summary['collection'][bin_id]['color'] = self.bins_info_dict[bin_id]['html_color'] or '#212121'
            self.summary['collection'][bin_id]['source'] = self.bins_info_dict[bin_id]['source'] or 'unknown_source'
            self.summary['meta']['total_nts_in_collection'] += self.summary['collection'][bin_id]['total_length']
            self.summary['meta']['num_contigs_in_collection'] += self.summary['collection'][bin_id]['num_contigs']
            self.progress.end()

        # bins are computed, add some relevant meta info:
        self.summary['meta']['percent_contigs_nts_described_by_collection'] = '%.2f' % (self.summary['meta']['total_nts_in_collection'] * 100.0 / int(self.a_meta['total_length']))
        self.summary['meta']['percent_profile_nts_described_by_collection'] = '%.2f' % (self.summary['meta']['total_nts_in_collection'] * 100.0 / int(self.p_meta['total_length']))
        self.summary['meta']['bins'] = self.get_bins_ordered_by_completeness_and_size()

        if not self.quick:
            # generate a TAB-delimited text output file for bin summaries
            summary_of_bins = {}
            properties = ['taxon', 'total_length', 'num_contigs', 'N50', 'GC_content']
            if self.completeness_data_available:
                properties += ['percent_completion', 'percent_redundancy']

            for bin_name in self.summary['collection']:
                summary_of_bins[bin_name] = dict([(prop, self.summary['collection'][bin_name][prop]) for prop in properties])

            output_file_obj = self.get_output_file_handle(prefix='bins_summary.txt')
            utils.store_dict_as_TAB_delimited_file(summary_of_bins, None, headers=['bins'] + properties, file_obj=output_file_obj)

            self.report_misc_data_files(target_table='layers')
            self.report_misc_data_files(target_table='items')

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
                self.summary['bin_percent_abundance_items'] = sorted(list(self.bin_percent_recruitment_per_sample.values())[0].keys())
                output_file_obj = self.get_output_file_handle(sub_directory='bins_across_samples', prefix='bins_percent_recruitment.txt')
                utils.store_dict_as_TAB_delimited_file(self.bin_percent_recruitment_per_sample,
                                                       None,
                                                       headers=['samples'] + sorted(self.collection_profile.keys()) + ['__splits_not_binned__'],
                                                       file_obj=output_file_obj)


        if self.debug:
            import json
            print(json.dumps(self.summary, sort_keys=True, indent=4))

        self.index_html = SummaryHTMLOutput(self.summary, r=self.run, p=self.progress).generate(quick=self.quick)


    def get_bins_ordered_by_completeness_and_size(self):
        if self.completeness_data_available:
            return [t[2] for t in sorted([(self.summary['collection'][bin]['percent_completion'], self.summary['collection'][bin]['total_length'], bin) for bin in self.summary['collection']], reverse=True)]
        else:
            return sorted(self.summary['collection'].keys())


class ContigSummarizer(SummarizerSuperClass):
    def __init__(self, contigs_db_path, run=run, progress=progress):
        self.contigs_db_path = contigs_db_path
        self.run = run
        self.progress = progress


    def get_contigs_db_info_dict(self, run=run, progress=progress, include_AA_counts=False, split_names=None, gene_caller_to_use=None):
        """Returns an info dict for a given contigs db.

           Please note that this function will only return gene calls made by `gene_caller_to_use`,
           but it will report other gene callers found in the contigs database, and how many genes
           were not reported. The client side should check for those to report to the user.
        """

        if not gene_caller_to_use:
            gene_caller_to_use = constants.default_gene_caller

        args = argparse.Namespace(contigs_db=self.contigs_db_path)

        run = terminal.Run()
        progress = terminal.Progress()
        run.verbose = False
        progress.verbose = False
        c = ContigsSuperclass(args, r=run, p=progress)

        info_dict = {'path': self.contigs_db_path,
                     'gene_caller_ids': set([]),
                     'gene_caller': gene_caller_to_use}

        for key in c.a_meta:
            info_dict[key] = c.a_meta[key]

        gene_calls_from_other_gene_callers = Counter()

        # Two different strategies here depending on whether we work with a given set if split ids or
        # everything in the contigs database.
        def process_gene_call(g):
            gene_caller = c.genes_in_contigs_dict[g]['source']
            if gene_caller == gene_caller_to_use:
                info_dict['gene_caller_ids'].add(g)
            else:
                gene_calls_from_other_gene_callers[gene_caller] += 1

        if split_names:
            split_names = set(split_names)
            c.init_split_sequences()
            seq = ''.join([c.split_sequences[split_name] for split_name in split_names])
            for e in list(c.genes_in_splits.values()):
                if e['split'] in split_names:
                    process_gene_call(e['gene_callers_id'])
        else:
            c.init_contig_sequences()
            seq = ''.join([e['sequence'] for e in list(c.contig_sequences.values())])

            for g in c.genes_in_contigs_dict:
                process_gene_call(g)

        if len(gene_calls_from_other_gene_callers):
            run.info_single('PLEASE READ CAREFULLY. Contigs db info summary will not include %d gene calls that were\
                             not identified by "%s", the default gene caller. Other gene calls found in this contigs\
                             database include, %s. If you are more interested in gene calls in any of those, you should\
                             indicate that through the `--gene-caller` parameter in your program.' \
                                                                % (sum(gene_calls_from_other_gene_callers.values()), \
                                                                   gene_caller_to_use, \
                                                                   ', '.join(['%d gene calls by %s' % (tpl[1], tpl[0]) for tpl in gene_calls_from_other_gene_callers.items()])))

        info_dict['gene_calls_from_other_gene_callers'] = gene_calls_from_other_gene_callers
        info_dict['gc_content'] = seqlib.Composition(seq).GC_content
        info_dict['total_length'] = len(seq)

        info_dict['partial_gene_calls'] = set([])
        for gene_caller_id in info_dict['gene_caller_ids']:
            if c.genes_in_contigs_dict[gene_caller_id]['partial']:
                info_dict['partial_gene_calls'].add(gene_caller_id)

        info_dict['num_genes'] = len(info_dict['gene_caller_ids'])
        if info_dict['num_genes']:
            info_dict['avg_gene_length'] = numpy.mean([c.gene_lengths[gene_caller_id] for gene_caller_id in info_dict['gene_caller_ids']])
            info_dict['num_genes_per_kb'] = info_dict['num_genes'] * 1000.0 / info_dict['total_length']
        else:
            info_dict['avg_gene_length'], info_dict['num_genes_per_kb'] = 0.0, 0

        # get completeness / contamination estimates
        p_completion, p_redundancy, domain, domain_confidence, results_dict = completeness.Completeness(self.contigs_db_path).get_info_for_splits(split_names if split_names else set(c.splits_basic_info.keys()))

        info_dict['hmm_sources_info'] = c.hmm_sources_info
        info_dict['percent_completion'] = p_completion
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
            raise ConfigError("We have a big problem. I am reporting from get_contigs_db_info_dict. This function must\
                                produce a dictionary that meets the requirements defined in the constants module of anvi'o\
                                for 'essential genome info'. But when I look at the resulting dictionary this funciton is\
                                about to return, I can see it is missing some stuff :/ This is not a user error, but it needs\
                                the attention of an anvi'o developer. Here are the keys that should have been in the results\
                                but missing: '%s'" % (', '.join(missing_keys)))

        return info_dict


    def get_summary_dict_for_assembly(self, gene_caller_to_use=None):
        """Returns a simple summary dict for a given contigs database"""
        if not gene_caller_to_use:
            gene_caller_to_use = constants.default_gene_caller

        self.progress.new('Generating contigs db summary')

        self.progress.update('Initiating contigs super for %s ...' % self.contigs_db_path)
        run, progress = terminal.Run(), terminal.Progress()
        run.verbose, progress.verbose = False, False
        c = ContigsSuperclass(argparse.Namespace(contigs_db=self.contigs_db_path), r=run, p=progress)

        self.progress.update('Recovering info about %s ...' % self.contigs_db_path)
        num_genes = len([True for v in c.genes_in_contigs_dict.values() if v['source'] == gene_caller_to_use])
        project_name = c.a_meta['project_name']
        contig_lengths = sorted([e['length'] for e in c.contigs_basic_info.values()], reverse=True)
        total_length = sum(contig_lengths)
        num_contigs = len(contig_lengths)

        self.progress.update('Figuring out HMM hits in %s ...' % self.contigs_db_path)
        hmm = hmmops.SequencesForHMMHits(self.contigs_db_path)

        self.progress.update('Summarizing %s ...' % self.contigs_db_path)
        summary = {}
        summary['project_name'] = project_name
        summary['total_length'] = total_length
        summary['num_genes'] = num_genes
        summary['gene_caller'] = gene_caller_to_use
        summary['num_contigs'] = num_contigs
        summary['n_values'] = self.calculate_N_values(contig_lengths, total_length, N=100)
        summary['contig_lengths'] = contig_lengths
        summary['gene_hit_counts_per_hmm_source'] = hmm.get_gene_hit_counts_per_hmm_source()
        summary['num_genomes_per_SCG_source_dict'] = hmm.get_num_genomes_from_SCG_sources_dict()

        self.progress.end()

        return summary


    def calculate_N_values(self, contig_lengths, total_length, N=100):
        results = []

        temp_length = 0
        contigs_index = 0
        n_index = 1

        while n_index <= N:
            if (temp_length >= int(((total_length / N) * n_index))):
                results.append({
                        'num_contigs': contigs_index,
                        'length':      contig_lengths[contigs_index - 1]
                    })
                n_index += 1
            else:
                temp_length += contig_lengths[contigs_index]
                contigs_index += 1

        return results

class Bin:
    def __init__(self, summary, bin_id, r=run, p=progress):
        self.summary = summary

        if not self.summary.initialized:
            raise ConfigError("The summary object you sent to the `Bin` class to make sense of the bin '%s' does\
                               not seem to have been initialized. Anvi'o could have taken care of it for you, but\
                                it will not (not only because anvi'o is implemented by mean people, but also it kinda\
                                likes to be explicit about this kind of stuff). Please initialize your summary object\
                                first." % (bin_id))

        if bin_id not in self.summary.bin_ids:
            raise ConfigError("Bin '%s' does not seem to be in this summary :/ These are the ones in it: %s." % (bin_id, ', '.join(self.summary.bin_ids)))

        self.bin_id = bin_id
        self.split_names = summary.collection_dict[self.bin_id]
        self.progress = p
        self.run = r

        self.across_samples = {}
        self.bin_profile = {}
        self.bin_info_dict = {'files': {}}

        # this will quickly populate self.contig_names, self.total_length, and self.num_contigs variables
        self.process_contigs(quick=True)

        self.output_directory = None
        self.contig_lengths = []

        # make sure all split_names in the collection is actually in the contigs database.
        # in collections stored in the contigs database, split_names that are not in the
        # oritinal contigs used to generate contigs database *may* end up in the
        # collections table. we gotta make sure we deal with them properly:
        missing_ids = [split_id for split_id in self.split_names if split_id not in self.summary.split_sequences]
        if len(missing_ids):
            for missing_id in missing_ids:
                self.split_names.remove(missing_id)

            self.run.warning('%d split id(s) in bin "%s" reported by collection "%s" is not found in the\
                              contigs database and removed from the bin summary. If this does not make\
                              any sense, you may need make sure everything is in order. The thing is,\
                              sometimes external clustering results that are added to the contigs via\
                              `anvi-populate-collections-table` may include split names that are not used\
                              while the contigs database was generated.'\
                                                % (len(missing_ids), bin_id, self.summary.collection_name))


        self.gene_caller_ids = self.get_gene_caller_ids()
        self.num_splits = len(self.split_names)

        # make these dicts avilable:
        self.gene_level_coverage_stats_dict = {}
        self.split_coverage_values_per_nt_dict = {}

        A = lambda x: self.summary.gene_level_coverage_stats_dict[gene_callers_id][sample_name][x]

        # populate gene coverage and detection dictionaries by subsetting them from the parent summary object
        if self.summary.gene_level_coverage_stats_dict:
            for gene_callers_id in self.gene_caller_ids:
                self.gene_level_coverage_stats_dict[gene_callers_id] = {}

                for sample_name in self.summary.p_meta['samples']:
                    self.gene_level_coverage_stats_dict[gene_callers_id][sample_name] = {'mean_coverage': A('mean_coverage'),
                                                                                         'detection': A('detection'),
                                                                                         'non_outlier_mean_coverage': A('non_outlier_mean_coverage'),
                                                                                         'non_outlier_coverage_std': A('non_outlier_coverage_std')}

                    if 'gene_coverage_values_per_nt' in self.summary.gene_level_coverage_stats_dict[gene_callers_id][sample_name]:
                        self.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['gene_coverage_values_per_nt'] = A('gene_coverage_values_per_nt')
                        self.gene_level_coverage_stats_dict[gene_callers_id][sample_name]['non_outlier_positions'] = A('non_outlier_positions')

        # populate coverage values per nucleutide for the bin.
        if self.summary.split_coverage_values_per_nt_dict:
            for split_name in self.split_names:
                self.split_coverage_values_per_nt_dict[split_name] = self.summary.split_coverage_values_per_nt_dict[split_name]


    def create(self):
        self.create_bin_dir()

        self.store_sequences_for_hmm_hits()

        self.process_contigs()

        if self.summary.completeness_data_available:
            self.access_completeness_scores()

        if self.summary.non_single_copy_gene_hmm_data_available:
            self.summarize_hmm_hits()

        self.compute_basic_stats()

        self.set_taxon_calls()

        self.store_genes_basic_info()

        self.store_gene_level_coverage_stats()

        self.store_profile_data()

        return self.bin_info_dict


    def create_bin_dir(self):
        self.progress.update('Creating the output directory ...')

        if not self.output_directory:
            self.progress.end()
            raise ConfigError('You called Bin.create() before setting an output directory. Anvio says "nope, thanks".')

        filesnpaths.gen_output_directory(self.output_directory)


    def access_completeness_scores(self):
        self.progress.update('Accessing completeness scores ...')

        p_completion, p_redundancy, domain, domain_confidence, results_dict = self.summary.completeness.get_info_for_splits(set(self.split_names))

        self.bin_info_dict['completeness'] = results_dict

        self.bin_info_dict['percent_redundancy'] = p_redundancy
        self.bin_info_dict['percent_completion'] = p_completion
        self.bin_info_dict['scg_domain'] = domain
        self.bin_info_dict['scg_domain_confidence'] = domain_confidence

        for k in ['percent_redundancy', 'percent_completion']:
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
        split_names_with_hmm_hits = [split_id for split_id in self.split_names if split_id in self.summary.hmm_searches_dict]

        for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header:
            hmm_items = self.summary.hmm_sources_info[hmm_search_source]['genes']
            info_dict[hmm_search_source] = dict([(hmm_item, 0) for hmm_item in hmm_items])

            hits_in_splits = []
            # keep track of unique identifiers of hmm hits to not count a single hit that spans across multiple splits:
            unique_identifiers_seen = set([])

            for split_id in split_names_with_hmm_hits:
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

        for split_name in self.split_names:
            if split_name not in self.summary.split_name_to_genes_in_splits_entry_ids:
                continue

            for gene_entry_id in self.summary.split_name_to_genes_in_splits_entry_ids[split_name]:
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


    def store_gene_level_coverage_stats(self):
        if self.summary.quick:
            return

        if not self.summary.gene_level_coverage_stats_dict:
            return

        self.progress.update('Storing gene coverages ...')

        headers = ['gene_callers_id'] + self.summary.p_meta['samples']

        for key, file_name in [('mean_coverage', 'gene_coverages.txt'),
                               ('detection', 'gene_detection.txt'),
                               ('non_outlier_mean_coverage', 'gene_non_outlier_coverages.txt'),
                               ('non_outlier_coverage_std', 'gene_non_outlier_coverage_stds.txt')]:
            # we will create a new dictionary here by subestting values of `key` from self.gene_level_coverage_stats_dict,
            # so we can store that information into `file_name`. magical stuff .. by us .. level 3000 wizards who can summon
            # inefficiency at most random places. SHUT UP.

            d = utils.get_values_of_gene_level_coverage_stats_as_dict(self.gene_level_coverage_stats_dict, key)

            utils.store_dict_as_TAB_delimited_file(d, None, headers=headers, file_obj=self.get_output_file_handle(file_name))


    def store_genes_basic_info(self):
        if self.summary.quick:
            return

        self.progress.update('Sorting out gene calls ...')

        d = {}

        headers = ['contig', 'start', 'stop', 'direction']
        header_items_for_gene_sequences = ['dna_sequence']
        if self.summary.report_aa_seqs_for_gene_calls:
            header_items_for_gene_sequences.append('aa_sequence')

        for gene_callers_id in self.gene_caller_ids:
            d[gene_callers_id] = {}
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

            # finally add the dna and amino acid sequence for gene calls:
            contig = self.summary.genes_in_contigs_dict[gene_callers_id]['contig']
            start = self.summary.genes_in_contigs_dict[gene_callers_id]['start']
            stop = self.summary.genes_in_contigs_dict[gene_callers_id]['stop']

            dna_sequence = self.summary.contig_sequences[contig]['sequence'][start:stop]
            if self.summary.genes_in_contigs_dict[gene_callers_id]['direction'] == 'r':
                dna_sequence = utils.rev_comp(dna_sequence)

            d[gene_callers_id]['dna_sequence'] = dna_sequence

            # if the user asked for it, report amino acid sequences as well
            if self.summary.report_aa_seqs_for_gene_calls:
                try:
                    d[gene_callers_id]['aa_sequence'] = utils.get_DNA_sequence_translated(dna_sequence, gene_callers_id)
                except:
                    d[gene_callers_id]['aa_sequence'] = ''

        output_file_obj = self.get_output_file_handle('gene_calls.txt')

        if self.summary.gene_function_call_sources:
            sources = [[source, source + ' (ACCESSION)'] for source in self.summary.gene_function_call_sources]
            headers = ['gene_callers_id'] + headers + [item for sublist in sources for item in sublist] + header_items_for_gene_sequences
        else:
            headers = ['gene_callers_id'] + headers + header_items_for_gene_sequences

        self.progress.update('Storing genes basic info ...')
        utils.store_dict_as_TAB_delimited_file(d, None, headers=headers, file_obj=output_file_obj)

        self.bin_info_dict['genes'] = {'num_genes_found': len(self.gene_caller_ids)}


    def store_sequences_for_hmm_hits(self):
        if self.summary.quick:
            return

        s = SequencesForHMMHits(self.summary.contigs_db_path)
        hmm_sequences_dict = s.get_sequences_dict_for_hmm_hits_in_splits({self.bin_id: self.split_names})

        single_copy_gene_hmm_sources = [hmm_search_source for hmm_search_type, hmm_search_source in self.summary.hmm_searches_header]
        non_single_copy_gene_hmm_sources = self.summary.completeness.sources

        for hmm_search_source in single_copy_gene_hmm_sources + non_single_copy_gene_hmm_sources:
            filtered_hmm_sequences_dict = utils.get_filtered_dict(hmm_sequences_dict, 'source', set([hmm_search_source]))

            output_file_obj = self.get_output_file_handle('%s-hmm-sequences.txt' % hmm_search_source, key=hmm_search_source)

            for gene_unique_id in filtered_hmm_sequences_dict:
                header, sequence = s.get_FASTA_header_and_sequence_for_gene_unique_id(hmm_sequences_dict, gene_unique_id)
                output_file_obj.write('>%s\n%s\n' % (header, sequence))


    def process_contigs(self, quick=False):
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

        self.contig_names = set([self.summary.splits_basic_info[split_name]['parent'] for split_name in self.summary.splits_basic_info if split_name in self.split_names])
        self.total_length = sum([self.summary.splits_basic_info[split_name]['length'] for split_name in self.summary.splits_basic_info if split_name in self.split_names])
        self.num_contigs = len(self.contig_names)

        self.bin_info_dict['total_length'] = self.total_length
        self.bin_info_dict['contig_names'] = self.contig_names
        self.bin_info_dict['num_contigs'] = len(self.contig_names)

        if self.summary.quick or quick:
            return

        self.progress.update('Creating the FASTA file ...')

        # store original split names:
        self.store_data_in_file('original_split_names.txt', '\n'.join(self.split_names))

        fasta_file = self.get_output_file_handle('contigs.fa')
        fasta_file.write(self.get_bin_sequence())
        fasta_file.close()

        self.store_data_in_file('num_contigs.txt', '%d' % self.bin_info_dict['num_contigs'])
        self.store_data_in_file('total_length.txt', '%d' % self.bin_info_dict['total_length'])


    def get_bin_sequence(self):
        output = ""

        # this dict will keep all the contig ids found in this bin with split names ordered:
        contigs_represented = utils.get_contigs_splits_dict(self.split_names, self.summary.splits_basic_info)

        # now it is time to go through each contig found in contigs_represented to
        # figure out what fraction of the contig is in fact in this bin
        for contig_id in contigs_represented:
            splits_order = list(contigs_represented[contig_id].keys())

            # this is critical: sequential_blocks is a list of one ore more lists, where each item of this list
            # describes a range of splits that follow each other to represent a coherent
            # chunk of the parent sequence (if all splits from a contig is selected into this bin,
            # then there would be one list item that spans across the entire contig):
            sequential_blocks = ccollections.GetSequentialBlocksOfSplits(splits_order).process()

            for sequential_block in sequential_blocks:
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
                for split_order in sequential_block:
                    sequence += self.summary.split_sequences[contigs_represented[contig_id][split_order]]

                fasta_id = contig_id + appendix
                self.contig_lengths.append(len(sequence))
                
                output += '>%s\n' % fasta_id
                output += '%s\n' % textwrap.fill(sequence, 80, break_on_hyphens=False)

        return output


    def set_taxon_calls(self):
        self.progress.update('Filling in taxonomy info ...')

        self.bin_info_dict['taxon_calls'] = []
        self.bin_info_dict['taxon'] = 'Unknown'

        if not self.summary.a_meta['gene_level_taxonomy_source']:
            return

        taxon_calls_counter = Counter()
        for split_id in self.split_names:
            if split_id in self.summary.splits_taxonomy_dict:
                taxon_calls_counter[self.summary.splits_taxonomy_dict[split_id]] += 1
            else:
                taxon_calls_counter['None'] += 1

        taxon_calls = sorted([list(tc) for tc in list(taxon_calls_counter.items())], key=lambda x: int(x[1]), reverse=True)

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
        self.bin_info_dict['GC_content'] = numpy.mean([self.summary.splits_basic_info[split_id]['gc_content'] for split_id in self.split_names]) * 100

        self.store_data_in_file('N50.txt', '%d' % self.bin_info_dict['N50'] if self.bin_info_dict['N50'] else 'NA')
        self.store_data_in_file('GC_content.txt', '%.4f' % self.bin_info_dict['GC_content'])


    def get_output_file_handle(self, prefix='output.txt', overwrite=False, key=None):
        file_path = os.path.join(self.output_directory, '%s-%s' % (self.bin_id, prefix))

        if os.path.exists(file_path) and not overwrite:
            raise ConfigError('get_output_file_handle: well, this file already exists: "%s"' % file_path)

        if not key:
            key = prefix.split('.')[0].replace('-', '_')

        self.bin_info_dict['files'][key] = file_path[len(self.summary.output_directory):].strip('/')

        return open(file_path, 'w')


    def store_data_in_file(self, output_file_name_posfix, content):
        output_file_obj = self.get_output_file_handle(output_file_name_posfix)
        output_file_obj.write('%s\n' % content)
        output_file_obj.close()




