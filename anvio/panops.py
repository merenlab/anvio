# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import json
import math
import copy
import argparse
import pandas as pd
import networkx as nx
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
# from scipy.stats import entropy
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, to_tree
from scipy.spatial.distance import cdist, squareform
import random

from itertools import chain

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.tables.miscdata as miscdata

from anvio.dbinfo import DBInfo
from anvio.drivers.blast import BLAST
from anvio.drivers.diamond import Diamond
from anvio.drivers.mcl import MCL
from anvio.drivers import Aligners

from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomestorage import GenomeStorage
from anvio.tables.geneclusters import TableForGeneClusters
from anvio.tables.views import TablesForViews

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
aligners = Aligners()

additional_param_sets_for_sequence_search = {'diamond'   : '--masking 0',
                                             'ncbi_blast': ''}


class Pangenome(object):
    def __init__(self, args=None, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        self.max_num_gene_clusters_for_hierarchical_clustering = constants.max_num_items_for_hierarchical_clustering

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.genome_names_to_focus = A('genome_names')
        self.genomes_storage_path = A('genomes_storage')
        self.genomes = None
        self.project_name = A('project_name')
        self.output_dir = A('output_dir')
        self.num_threads = A('num_threads')
        self.user_defined_gene_clusters = A('gene_clusters_txt')
        self.skip_alignments = A('skip_alignments')
        self.skip_homogeneity = A('skip_homogeneity')
        self.quick_homogeneity = A('quick_homogeneity')
        self.align_with = A('align_with')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = anvio.DEBUG
        self.min_percent_identity = A('min_percent_identity')
        self.gene_cluster_min_occurrence = A('min_occurrence')
        self.mcl_inflation = A('mcl_inflation')
        self.minbit = A('minbit')
        self.use_ncbi_blast = A('use_ncbi_blast')
        self.exclude_partial_gene_calls = A('exclude_partial_gene_calls')
        self.description_file_path = A('description')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')
        self.enforce_the_analysis_of_excessive_number_of_genomes = anvio.USER_KNOWS_IT_IS_NOT_A_GOOD_IDEA

        self.additional_params_for_seq_search = A('additional_params_for_seq_search')
        self.additional_params_for_seq_search_processed = False

        if not self.project_name:
            raise ConfigError("Please set a project name using --project-name or -n.")

        # when it is time to organize gene_clusters
        self.linkage = A('linkage') or constants.linkage_method_default
        self.distance = A('distance') or constants.distance_metric_default

        self.log_file_path = None

        # to be filled during init:
        self.amino_acid_sequences_dict = {}
        self.view_data = {}
        self.view_data_presence_absence = {}
        self.additional_view_data = {}
        self.aligner = None

        # we don't know what we are about
        self.description = None


    def load_genomes(self):
        # genome_name parameter can be a file or comma seperated genome names.
        if self.genome_names_to_focus:
            if filesnpaths.is_file_exists(self.genome_names_to_focus, dont_raise=True):
                self.genome_names_to_focus = utils.get_column_data_from_TAB_delim_file(self.genome_names_to_focus, column_indices=[0], expected_number_of_fields=1)[0]
            else:
                self.genome_names_to_focus = [g.strip() for g in self.genome_names_to_focus.split(',')]

            self.run.warning("A subset of genome names is found, and anvi'o will focus only on to those.")

        self.genomes_storage = GenomeStorage(self.genomes_storage_path, storage_hash=None, genome_names_to_focus=self.genome_names_to_focus)
        self.genomes = self.genomes_storage.get_genomes_dict()

        self.external_genome_names = [g for g in self.genomes if self.genomes[g]['external_genome']]
        self.internal_genome_names = [g for g in self.genomes if not self.genomes[g]['external_genome']]

        self.hash_to_genome_name = {}
        for genome_name in self.genomes:
            self.hash_to_genome_name[self.genomes[genome_name]['genome_hash']] = genome_name


    def generate_pan_db(self):
        meta_values = {'internal_genome_names': ','.join(self.internal_genome_names),
                       'external_genome_names': ','.join(self.external_genome_names),
                       'num_genomes': len(self.genomes),
                       'min_percent_identity': self.min_percent_identity,
                       'gene_cluster_min_occurrence': self.gene_cluster_min_occurrence,
                       'mcl_inflation': self.mcl_inflation,
                       'user_provided_gene_clusters_txt': True if self.user_defined_gene_clusters else False,
                       'default_view': 'gene_cluster_presence_absence',
                       'use_ncbi_blast': self.use_ncbi_blast,
                       'additional_params_for_seq_search': self.additional_params_for_seq_search,
                       'minbit': self.minbit,
                       'exclude_partial_gene_calls': self.exclude_partial_gene_calls,
                       'gene_alignments_computed': False if self.skip_alignments else True,
                       'genomes_storage_hash': self.genomes_storage.get_storage_hash(),
                       'project_name': self.project_name,
                       'items_ordered': False,
                       'reaction_network_ko_annotations_hash': None,
                       'reaction_network_kegg_database_release': None,
                       'reaction_network_modelseed_database_sha': None,
                       'reaction_network_consensus_threshold': None,
                       'reaction_network_discard_ties': None,
                       'description': self.description if self.description else '_No description is provided_',
                      }

        dbops.PanDatabase(self.pan_db_path, quiet=False).create(meta_values)

        # know thyself.
        self.args.pan_db = self.pan_db_path


    def get_output_file_path(self, file_name, delete_if_exists=False):
        output_file_path = os.path.join(self.output_dir, file_name)

        if delete_if_exists:
            if os.path.exists(output_file_path):
                os.remove(output_file_path)

        return output_file_path


    def check_programs(self):
        if self.use_ncbi_blast:
            utils.is_program_exists('blastp')
        else:
            utils.is_program_exists('diamond')

        utils.is_program_exists('mcl')


    def check_project_name(self):
        # check the project name:
        if not self.project_name:
            raise ConfigError("Please set a project name using the `--project-name` parameter, and be prepared to see "
                              "it around as anvi'o will use it for multiple things, such as setting the output directory "
                              "and naming various output files including the database file that will be generated at the "
                              "end of the process. If you set your own output directory name, you can have multiple "
                              "projects in it and all of those projects can use the same intermediate files whenever "
                              "possible.")

        utils.is_this_name_OK_for_database('pan project name', self.project_name, stringent=False)


    def check_params(self):
        # if the user did not set a specific output directory name, use the project name
        # for it:
        self.output_dir = self.output_dir if self.output_dir else self.project_name

        # deal with the output directory:
        try:
            filesnpaths.is_file_exists(self.output_dir)
        except FilesNPathsError:
            filesnpaths.gen_output_directory(self.output_dir, delete_if_exists=self.overwrite_output_destinations)

        filesnpaths.is_output_dir_writable(self.output_dir)
        self.output_dir = os.path.abspath(self.output_dir)

        if not self.log_file_path:
            self.log_file_path = self.get_output_file_path('log.txt')

        filesnpaths.is_output_file_writable(self.log_file_path)
        os.remove(self.log_file_path) if os.path.exists(self.log_file_path) else None

        if self.user_defined_gene_clusters:
            filesnpaths.is_gene_clusters_txt(self.user_defined_gene_clusters)

        if not isinstance(self.minbit, float):
            raise ConfigError("minbit value must be of type float :(")

        if self.minbit < 0 or self.minbit > 1:
            raise ConfigError("Well. minbit must be between 0 and 1. Yes. Very boring.")

        if not isinstance(self.min_percent_identity, float):
            raise ConfigError("Minimum percent identity value must be of type float :(")

        if self.min_percent_identity < 0 or self.min_percent_identity > 100:
            raise ConfigError("Minimum percent identity must be between 0%% and 100%%. Although your %.2f%% is "
                              "pretty cute, too." % self.min_percent_identity)


        if len([c for c in list(self.genomes.values()) if 'genome_hash' not in c]):
            raise ConfigError("self.genomes does not seem to be a properly formatted dictionary for "
                              "the anvi'o class Pangenome.")

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering "
                              "while also asking it to enforce it.")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'r').read()

        self.pan_db_path = self.get_output_file_path(self.project_name + '-PAN.db')


    def process_additional_params(self):
        """Process user-requested additional params for sequence search with defaults

        This is a bit complicated as there are those additional parameter sets used by
        anvi'o (which are defined in `additional_param_sets_for_sequence_search`), and
        there are those that may or may not be passed by the user. To reconcile all,
        we first need to determine which algorithm they are using (i.e., diamond or
        ncbi-blast, which are the only two options we offer at the time this function
        was written), and then see if tey passed any 'additional additional' params to
        this instance.
        """

        if self.additional_params_for_seq_search is not None:
            # the user set something (or unset everything by passing ""): we obey
            pass
        else:
            if self.use_ncbi_blast:
                self.additional_params_for_seq_search = additional_param_sets_for_sequence_search['ncbi_blast']
            else:
                self.additional_params_for_seq_search = additional_param_sets_for_sequence_search['diamond']

        self.additional_params_for_seq_search_processed = True


    def run_diamond(self, unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict):
        # in case someone called this function directly without going through
        # `self.process`:
        self.process_additional_params()

        diamond = Diamond(unique_AA_sequences_fasta_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads, overwrite_output_destinations=self.overwrite_output_destinations)

        diamond.names_dict = unique_AA_sequences_names_dict
        diamond.additional_params_for_blastp = self.additional_params_for_seq_search
        diamond.search_output_path = self.get_output_file_path('diamond-search-results')
        diamond.tabular_output_path = self.get_output_file_path('diamond-search-results.txt')

        return diamond.get_blast_results()


    def run_blast(self, unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict):
        # in case someone called this function directly without going through
        # `self.process`:
        self.process_additional_params()

        self.run.warning("You elected to use NCBI's `blastp` for amino acid sequence search. Running blastp will "
                         "be significantly slower than DIAMOND, but in some cases, slightly more sensitive. "
                         "We are unsure about whether the slight increase in sensitivity may justify significant "
                         "increase in run time, but you are the boss.", lc="cyan")

        blast = BLAST(unique_AA_sequences_fasta_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads, overwrite_output_destinations=self.overwrite_output_destinations)

        blast.names_dict = unique_AA_sequences_names_dict
        blast.additional_params_for_blast = self.additional_params_for_seq_search
        blast.log_file_path = self.log_file_path
        blast.search_output_path = self.get_output_file_path('blast-search-results.txt')

        return blast.get_blast_results()


    def run_search(self, unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict):
        if self.use_ncbi_blast:
            return self.run_blast(unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict)
        else:
            return self.run_diamond(unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict)


    def run_mcl(self, mcl_input_file_path):
        mcl = MCL(mcl_input_file_path, run=self.run, progress=self.progress, num_threads=self.num_threads)

        mcl.inflation = self.mcl_inflation
        mcl.clusters_file_path = self.get_output_file_path('mcl-clusters.txt')
        mcl.log_file_path = self.log_file_path

        return mcl.get_clusters_dict()


    def gen_mcl_input(self, blastall_results):
        self.run.warning(None, header="MCL INPUT", lc="green")

        self.progress.new('Processing search results')
        self.progress.update('...')

        all_ids = set([])

        # mapping for the fields in the blast output
        mapping = [str, str, float, int, int, int, int, int, int, int, float, float]

        # here we perform an initial pass on the blast results to fill the dict that will hold
        # the bit score for each gene when it was blasted against itself. this dictionary
        # will then be used to calculate the 'minbit' value between two genes, which I learned
        # from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8). ITEP defines minbit as
        # 'bit score between target and query / min(selfbit for query, selbit for target)'. This
        # heuristic approach provides a mean to set a cutoff to eliminate weak matches between
        # two genes. minbit value reaches to 1 for hits between two genes that are almost identical.
        self_bit_scores = {}
        line_no = 1
        self.progress.update('(initial pass of the search results to set the self bit scores ...)')
        for line in open(blastall_results):
            fields = line.strip().split('\t')

            try:
                query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                    [mapping[i](fields[i]) for i in range(0, len(mapping))]
            except Exception as e:
                self.progress.end()
                raise ConfigError("Something went wrong while processing the blastall output file in line %d. "
                                   "Here is the error from the uppoer management: '''%s'''" % (line_no, e))
            line_no += 1
            all_ids.add(query_id)
            all_ids.add(subject_id)

            if query_id == subject_id:
                self_bit_scores[query_id] = bit_score

        self.progress.end()

        ids_without_self_search = all_ids - set(self_bit_scores.keys())
        if len(ids_without_self_search):
            search_tool = 'BLAST' if self.use_ncbi_blast else 'DIAMOND'
            self.run.warning("%s did not retun search results for %d of %d the amino acid sequences in your input FASTA file. "
                             "Anvi'o will do some heuristic magic to complete the missing data in the search output to recover "
                             "from this. But since you are a scientist, here are the amino acid sequence IDs for which %s "
                             "failed to report self search results: %s." \
                                                    % (search_tool, len(ids_without_self_search), len(all_ids), \
                                                       search_tool, ', '.join(ids_without_self_search)))

        # HEURISTICS TO ADD MISSING SELF SEARCH RESULTS
        # we are here, because amino acid sequences in ids_without_self_search did not have any hits in the search output
        # although they were in the FASTA file the target database were built from. so we will make sure they are not
        # missing from self_bit_scores dict, or mcl_input (additional mcl inputs will be stored in the following dict)
        additional_mcl_input_lines = {}

        for id_without_self_search in ids_without_self_search:
            entry_hash, gene_caller_id = id_without_self_search.split('_')

            try:
                genome_name = self.hash_to_genome_name[entry_hash]
            except KeyError:
                raise ConfigError("Something horrible happened. This can only happend if you started a new analysis with "
                                   "additional genomes without cleaning the previous work directory. Sounds familiar?")

            # divide the DNA length of the gene by three to get the AA length, and multiply that by two to get an approximate
            # bit score that would have recovered from a perfect match
            gene_amino_acid_sequence_length = len(self.genomes_storage.get_gene_sequence(genome_name, int(gene_caller_id), report_DNA_sequences=False))
            self_bit_scores[id_without_self_search] = gene_amino_acid_sequence_length * 2

            # add this SOB into additional_mcl_input_lines dict.
            additional_mcl_input_lines[id_without_self_search] = '%s\t%s\t1.0\n' % (id_without_self_search, id_without_self_search)


        # CONTINUE AS IF NOTHING HAPPENED
        self.run.info('Min percent identity', self.min_percent_identity)
        self.run.info('Minbit', self.minbit)
        self.progress.new('Processing search results')

        mcl_input_file_path = self.get_output_file_path('mcl-input.txt')
        mcl_input = open(mcl_input_file_path, 'w')

        line_no = 1
        num_edges_stored = 0
        for line in open(blastall_results):
            fields = line.strip().split('\t')

            query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                [mapping[i](fields[i]) for i in range(0, len(mapping))]

            line_no += 1

            if line_no % 5000 == 0:
                self.progress.update('Lines processed %s ...' % pp(line_no))

            #
            # FILTERING BASED ON PERCENT IDENTITY
            #
            if perc_id < self.min_percent_identity:
                continue

            #
            # FILTERING BASED ON MINBIT
            #
            minbit = bit_score / min(self_bit_scores[query_id], self_bit_scores[subject_id])
            if minbit < self.minbit:
                continue

            mcl_input.write('%s\t%s\t%f\n' % (query_id, subject_id, perc_id / 100.0))
            num_edges_stored += 1

        # add additional lines if there are any:
        for line in list(additional_mcl_input_lines.values()):
            mcl_input.write(line)
            num_edges_stored += 1

        mcl_input.close()

        self.progress.end()

        self.run.info('Filtered search results', '%s edges stored' % pp(num_edges_stored))
        self.run.info('MCL input', '%s' % mcl_input_file_path)

        return mcl_input_file_path


    def process_gene_clusters(self, gene_clusters_dict):
        self.progress.new('Generating view data')
        self.progress.update('...')

        gene_clusters = list(gene_clusters_dict.keys())

        for genome_name in self.genomes:
            self.genomes[genome_name]['singleton_gene_clusters'] = 0
            self.genomes[genome_name]['num_gene_clusters_raw'] = 0

        for gene_cluster in gene_clusters:
            self.view_data[gene_cluster] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.view_data_presence_absence[gene_cluster] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.additional_view_data[gene_cluster] = {'num_genes_in_gene_cluster': 0, 'num_genomes_gene_cluster_has_hits': 0, 'SCG': 0, 'max_num_paralogs': 0}

            for gene_entry in gene_clusters_dict[gene_cluster]:
                genome_name = gene_entry['genome_name']

                self.view_data[gene_cluster][genome_name] += 1
                self.view_data_presence_absence[gene_cluster][genome_name] = 1
                self.additional_view_data[gene_cluster]['num_genes_in_gene_cluster'] += 1
                self.genomes[genome_name]['num_gene_clusters_raw'] += 1

            genomes_contributing_to_gene_cluster = [t[0] for t in self.view_data_presence_absence[gene_cluster].items() if t[1]]

            if len(genomes_contributing_to_gene_cluster) == 1:
                self.genomes[genomes_contributing_to_gene_cluster[0]]['singleton_gene_clusters'] += 1

            self.additional_view_data[gene_cluster]['SCG'] = 1 if set(self.view_data[gene_cluster].values()) == set([1]) else 0
            self.additional_view_data[gene_cluster]['max_num_paralogs'] = max(self.view_data[gene_cluster].values())

            self.additional_view_data[gene_cluster]['num_genomes_gene_cluster_has_hits'] = len([True for genome in self.view_data[gene_cluster] if self.view_data[gene_cluster][genome] > 0])

        self.progress.end()
        ########################################################################################
        #                           FILTERING BASED ON OCCURRENCE
        ########################################################################################
        gene_clusters_of_interest = set([])
        for gene_cluster in gene_clusters:
            if self.additional_view_data[gene_cluster]['num_genomes_gene_cluster_has_hits'] >= self.gene_cluster_min_occurrence:
                gene_clusters_of_interest.add(gene_cluster)

        removed_gene_clusters = 0
        for gene_cluster in gene_clusters:
            if gene_cluster not in gene_clusters_of_interest:
                self.view_data.pop(gene_cluster)
                self.view_data_presence_absence.pop(gene_cluster)
                self.additional_view_data.pop(gene_cluster)
                gene_clusters_dict.pop(gene_cluster)
                removed_gene_clusters += 1

        if self.gene_cluster_min_occurrence > 1:
            self.run.info('gene_clusters min occurrence', '%d (the filter removed %d gene_clusters)' % (self.gene_cluster_min_occurrence, removed_gene_clusters))

        ########################################################################################
        #            CAN WE CLUSTER THIS STUFF? DOES THE USER WANT US TO TRY REGARDLESS?
        ########################################################################################
        if len(gene_clusters_dict) > self.max_num_gene_clusters_for_hierarchical_clustering:
            if self.enforce_hierarchical_clustering:
                self.run.warning("You have %s gene_clusters, which exceeds the number of gene_clusters anvi'o is comfortable to cluster. But "
                                 "since you have used the flag `--enforce-hierarchical-clustering`, anvi'o will attempt "
                                 "to create a hierarchical clustering of your gene_clusters anyway. It may take a bit of "
                                 "time. Pour yourself a coffee. Or go to a nice vacation. See you in 10 mins, or next year "
                                 "or never." % pp(len(gene_clusters_dict)))
            else:
                self.run.warning("It seems you have %s gene clusters in your pangenome. This exceeds the soft limit "
                                 "of %s for anvi'o to attempt to create a hierarchical clustering of your gene clusters "
                                 "(which becomes the center tree in all anvi'o displays). If you want a hierarchical "
                                 "clustering to be done anyway, please see the flag `--enforce-hierarchical-clustering`." \
                                            % (pp(len(gene_clusters_dict)), pp(self.max_num_gene_clusters_for_hierarchical_clustering)))
                self.skip_hierarchical_clustering = True

        ########################################################################################
        #                           STORING FILTERED DATA IN THE DB
        ########################################################################################
        TablesForViews(self.pan_db_path).create_new_view(
                                        view_data=self.view_data,
                                        table_name='gene_cluster_frequencies',
                                        view_name = 'gene_cluster_frequencies',
                                        from_matrix_form=True)

        TablesForViews(self.pan_db_path).create_new_view(
                                        view_data=self.view_data_presence_absence,
                                        table_name='gene_cluster_presence_absence',
                                        view_name = 'gene_cluster_presence_absence',
                                        from_matrix_form=True)

        item_additional_data_table = miscdata.TableForItemAdditionalData(self.args, r=terminal.Run(verbose=False))
        item_additional_data_keys = ['num_genomes_gene_cluster_has_hits', 'num_genes_in_gene_cluster', 'max_num_paralogs', 'SCG']
        item_additional_data_table.add(self.additional_view_data, item_additional_data_keys, skip_check_names=True)
        #                                                                                    ^^^^^^^^^^^^^^^^^^^^^
        #                                                                                   /
        # here we say skip_check_names=True, simply because there is no gene_clusters table has not been
        # generated yet, but the check names functionality in dbops looks for the gene clsuters table to
        # be certain. it is not a big deal here, since we absoluely know what gene cluster names we are
        # working with.

        ########################################################################################
        #                   RETURN THE -LIKELY- UPDATED PROTEIN CLUSTERS DICT
        ########################################################################################
        return gene_clusters_dict


    def gen_synteny_based_ordering_of_gene_clusters(self, gene_clusters_dict):
        """Take the dictionary of gene_clusters, and order gene_clusters per genome based on synteny of genes.

           This adds more orders to the pangenomic output so the user can enforce ordering of
           gene_clusters based on the synteny of genes they contain in a given genome.

           The synteny in this context is defined by the gene caller ids. Gene caller ids
           follow a numerical order in anvi'o contigs databases for genes that are coming
           from the same contig. Of course, the synteny does not mean much for genes that
           fragmented into multiple contigs.
        """

        # yes. this is meren converting the gene_clusters_dict into a pandas data frame :/ if you are reading
        # this line and if you are not evan, don't tell evan about this. everyone else: i don't know
        # what you're talking about.
        df = pd.DataFrame(list(chain.from_iterable(list(gene_clusters_dict.values()))))
        df = df.sort_values(by=['genome_name', 'gene_caller_id'])
        df = df.reset_index(drop=True)

        # forced synteny
        for genome_name in df.genome_name.unique():
            gene_clusters_in_genome = df.loc[(df.genome_name == genome_name)].gene_cluster_id.unique()
            gene_clusters_not_described = df.loc[~df.gene_cluster_id.isin(gene_clusters_in_genome)].gene_cluster_id.unique()
            gene_clusters_order_based_on_genome_synteny = list(gene_clusters_in_genome) + list(gene_clusters_not_described)

            order_name = 'Forced synteny <> %s' % genome_name

            dbops.add_items_order_to_db(self.pan_db_path, order_name, ','.join(gene_clusters_order_based_on_genome_synteny), order_data_type_newick=False, run=terminal.Run(verbose=False))

        gene_cluster_gene_cluster_edges = []
        # network description of gene_cluster-gene_cluster relationships given the gene synteny.
        gene_ordered_list_of_gene_clusters = list(zip(df.gene_caller_id, df.gene_cluster_id))
        for index in range(1, len(gene_ordered_list_of_gene_clusters)):
            (GENE_A, gene_cluster_A), (GENE_B, gene_cluster_B) = gene_ordered_list_of_gene_clusters[index-1], gene_ordered_list_of_gene_clusters[index]
            if GENE_A == GENE_B - 1:
                gene_cluster_gene_cluster_edges.append((gene_cluster_A, gene_cluster_B), )

        # FIXME: Do something with gene_cluster_gene_cluster_edges.


    def gen_hierarchical_clustering_of_gene_clusters(self):
        """Uses a clustering configuration to add hierarchical clustering of gene clusters into the pan db

        Note how this function cheats the system to create an enchanced clustering configuration:
        We want to use the clustering configurations for pan genomomic analyses to order
        gene clusters. however, we want to add something into the clustering configuraiton
        file, which depends on the number of genomes we have. this addition is 'num_genomes_gene_cluster_has_hits'
        data, which pulls together gene clusters that are distributed across genomes similarly based
        on this extra bit of inofrmation. becasue the clustering configurations framework in anvi'o
        does not allow us to have variable information in these recipes, we are going to generate one
        on the fly to have a more capable one."""

        if self.skip_hierarchical_clustering:
            return

        updated_clustering_configs = {}

        for config_name in constants.clustering_configs['pan']:
            config_path = constants.clustering_configs['pan'][config_name]

            # now we have the config path. we first get a temporary file path:
            enhanced_config_path = filesnpaths.get_temp_file_path()

            # setup the additional section based on the number of genomes we have:
            if config_name == 'presence-absence':
                additional_config_section="""\n[AdditionalData !PAN.db::item_additional_data]\ntable_form=dataframe\ncolumns_to_use = %s\nnormalize = False\n""" \
                                        % ','.join(['num_genomes_gene_cluster_has_hits'] * (int(round(len(self.genomes) / 2))))
            elif config_name == 'frequency':
                additional_config_section="""\n[AdditionalData !PAN.db::item_additional_data]\ntable_form=dataframe\ncolumns_to_use = %s\nnormalize = False\nlog=True\n""" \
                                        % ','.join(['num_genes_in_gene_cluster'] * (int(round(math.sqrt(len(self.genomes))))))

            # write the content down in to file at the new path:
            open(enhanced_config_path, 'w').write(open(config_path).read() + additional_config_section)

            # update the clustering configs:
            updated_clustering_configs[config_name] = enhanced_config_path

            dbops.do_hierarchical_clustering_of_items(self.pan_db_path, updated_clustering_configs, database_paths={'PAN.db': self.pan_db_path},\
                                                      input_directory=self.output_dir, default_clustering_config=constants.pan_default,\
                                                      distance=self.distance, linkage=self.linkage, run=terminal.Run(verbose=False), progress=self.progress)


    def populate_gene_cluster_homogeneity_index(self, gene_clusters_dict, gene_clusters_failed_to_align=set([])):
        if self.skip_alignments:
            self.run.warning('Skipping homogeneity calculations because gene clusters are not alligned.')
            return

        if self.skip_homogeneity:
            self.run.warning("Skipping homogeneity calculations per the '--skip-homogeneity' flag.")
            return

        pan = dbops.PanSuperclass(args=self.args, r=terminal.Run(verbose=False), p=self.progress)
        gene_cluster_names = set(list(gene_clusters_dict.keys()))

        d = pan.compute_homogeneity_indices_for_gene_clusters(gene_cluster_names=gene_cluster_names,
                                                              gene_clusters_failed_to_align=gene_clusters_failed_to_align,
                                                              num_threads=self.num_threads)

        if d is None:
            self.run.warning("Anvi'o received an empty dictionary for homogeneity indices. Not good :/ Returning empty handed,\
                              without updating anything in the pan database...")
            return

        keys = ['functional_homogeneity_index', 'geometric_homogeneity_index', 'combined_homogeneity_index', 'AAI_min', 'AAI_max', 'AAI_avg']
        miscdata.TableForItemAdditionalData(self.args, r=terminal.Run(verbose=False)).add(d, keys, skip_check_names=True)


    def populate_layers_additional_data_and_orders(self):
        self.progress.new('Layers additional data and orders')
        self.progress.update('Copmputing the hierarchical clustering of the (transposed) view data')

        layer_orders_data_dict = {}
        for clustering_tuple in [('gene_cluster presence absence', self.view_data), ('gene_cluster frequencies', self.view_data_presence_absence)]:
            v, d = clustering_tuple
            newick = clustering.get_newick_tree_data_for_dict(d, transpose=True, distance = self.distance, linkage=self.linkage)
            layer_orders_data_dict[v] = {'data_type': 'newick', 'data_value': newick}

        self.progress.update('Generating layers additional data ..')

        layers_additional_data_dict = {}
        layers_additional_data_keys = ['total_length', 'gc_content']

        for h in ['percent_completion', 'percent_redundancy']:
            if h in list(self.genomes.values())[0]:
                layers_additional_data_keys.append(h)

        layers_additional_data_keys.extend(['num_genes', 'avg_gene_length', 'num_genes_per_kb',
                                            'singleton_gene_clusters'])

        if self.gene_cluster_min_occurrence > 1:
            layers_additional_data_keys.extend(['num_gene_clusters_raw'])

        for genome_name in self.genomes:
            new_dict = {}
            for key in layers_additional_data_keys:
                new_dict[key] = self.genomes[genome_name][key]

            layers_additional_data_dict[genome_name] = new_dict

        # summarize gene cluster stats across genomes
        layers_additional_data_keys.extend(['num_gene_clusters'])
        for genome_name in self.genomes:
            layers_additional_data_dict[genome_name]['num_gene_clusters'] = 0
            for gene_cluster in self.view_data_presence_absence:
                # tracking the total number of gene clusters
                if self.view_data_presence_absence[gene_cluster][genome_name]:
                    layers_additional_data_dict[genome_name]['num_gene_clusters'] += 1

        self.progress.end()

        miscdata.TableForLayerOrders(self.args, r=terminal.Run(verbose=False)).add(layer_orders_data_dict)
        miscdata.TableForLayerAdditionalData(self.args, r=terminal.Run(verbose=False)).add(layers_additional_data_dict, layers_additional_data_keys)


    def sanity_check(self):
        self.check_project_name()
        self.check_programs()

        if not isinstance(self.mcl_inflation, float):
            raise ConfigError("Well, MCL likes its inflation parameter in 'float' form...")

        if self.mcl_inflation > 100 or self.mcl_inflation < 0.1:
            raise ConfigError("MCL inflation parameter should have a reasonable value :/ Like between 0.1 and 100.0.")

        if not isinstance(self.genomes, type({})):
            raise ConfigError("self.genomes must be a dict. Anvi'o needs an adult :(")

        if len(self.genomes) < 2:
            raise ConfigError("There must be at least two genomes for this workflow to work. You have like '%d' of them :/" \
                    % len(self.genomes))

        if len(self.genomes) > 100:
            if self.enforce_the_analysis_of_excessive_number_of_genomes:
                self.run.warning(f"Because you asked for it, anvi'o WILL analyze your {pp(len(self.genomes))} genomes. Perhaps it will "
                                 f"fail, but it will fail while trying.", header="YOU HAVE USED YOUR AUTHORITY ⚔️")
            else:
                raise ConfigError(f"Sorry for making you wait this long to tell you this, but you have a total of {pp(len(self.genomes))} "
                                  f"genomes to process, and anvi'o may not be the best platform to do that. The anvi'o pangenomics workflow "
                                  f"is designed for interactive and in-depth investigations of microbial pangenomes. While one can "
                                  f"certainly perform an in-depth investigation of their {pp(len(self.genomes))} genomes by turning them "
                                  f"into a pangenome, they will unlikely manage to do it with anvi'o with more than about 100 genomes. "
                                  f"We can't stop you if you insist running a pangenomic analyses on your genomes ANYWAY. For instance, "
                                  f"you can easily force anvi'o to start the analysis using the flag `--I-know-this-is-not-a-good-idea`. "
                                  f"But by using that flag, you will be forfeiting your privilege to complain to anvi'o developers if "
                                  f"something goes wrong with your analyses :)")

        if self.skip_alignments and self.align_with:
            raise ConfigError("You are asking anvi'o to skip aligning sequences within your gene clusters, and then you "
                              "are also asking it to use '%s' for aligning sequences within your gene clusters. It is easy "
                              "to ignore this and skip the alignment, but anvi'o gets nervous when it realizes her users are "
                              "being inconsistent. Please make up your mind, and come back as the explicit person you are" \
                                                                            % self.align_with)

        self.check_params()

        self.run.log_file_path = self.log_file_path
        self.run.info('Args', (str(self.args)), quiet=True)


    def store_gene_clusters(self, gene_clusters_dict):
        self.progress.new('Storing gene clusters in the database')
        self.progress.update('...')

        table_for_gene_clusters = TableForGeneClusters(self.pan_db_path, run=self.run, progress=self.progress)

        num_genes_in_gene_clusters = 0
        for gene_cluster_name in gene_clusters_dict:
            for gene_entry in gene_clusters_dict[gene_cluster_name]:
                table_for_gene_clusters.add(gene_entry)
                num_genes_in_gene_clusters += 1

        self.progress.end()

        table_for_gene_clusters.store()

        pan_db = dbops.PanDatabase(self.pan_db_path, quiet=True)
        pan_db.db.set_meta_value('num_gene_clusters', len(gene_clusters_dict))
        pan_db.db.set_meta_value('num_genes_in_gene_clusters', num_genes_in_gene_clusters)
        pan_db.disconnect()


    def gen_gene_clusters_dict_from_mcl_clusters(self, mcl_clusters):
        self.progress.new('Generating the gene clusters dictionary from raw MCL clusters')
        self.progress.update('...')

        gene_clusters_dict = {}

        for gene_cluster in mcl_clusters:
            gene_clusters_dict[gene_cluster] = []

            for entry_hash, gene_caller_id in [e.split('_') for e in mcl_clusters[gene_cluster]]:
                try:
                    genome_name = self.hash_to_genome_name[entry_hash]
                except KeyError:
                    self.progress.end()
                    raise ConfigError("Something horrible happened. This can only happen if you started a new analysis with "
                                       "additional genomes without cleaning the previous work directory. Sounds familiar?")

                gene_clusters_dict[gene_cluster].append({'gene_caller_id': int(gene_caller_id), 'gene_cluster_id': gene_cluster, 'genome_name': genome_name, 'alignment_summary': ''})

        self.progress.end()

        return gene_clusters_dict


    def compute_alignments_for_gene_clusters(self, gene_clusters_dict):
        if self.skip_alignments:
            self.run.warning('Skipping gene alignments.')
            return gene_clusters_dict, set([])

        # we run "select aligner" to print the citation information (the actual selection is
        # done in the `alignment_worker` down below)
        aligners.select(self.align_with)

        gene_cluster_names = list(gene_clusters_dict.keys())

        # we only need to align gene clusters with more than one sequence
        non_singleton_gene_cluster_names = [g for g in gene_cluster_names if len(gene_clusters_dict[g]) > 1]
        num_non_singleton_gene_clusters = len(non_singleton_gene_cluster_names)

        self.progress.new('Aligning amino acid sequences for genes in gene clusters', progress_total_items=num_non_singleton_gene_clusters)
        self.progress.update('...')

        manager = multiprocessing.Manager()
        input_queue = manager.Queue()
        output_queue = manager.Queue()

        for gene_cluster_name in non_singleton_gene_cluster_names:
            input_queue.put(gene_cluster_name)

        workers = []
        for i in range(0, self.num_threads):
            worker = multiprocessing.Process(target=Pangenome.alignment_worker,
                args=(input_queue, output_queue, gene_clusters_dict, self.genomes_storage, self.align_with, self.run))

            workers.append(worker)
            worker.start()

        received_gene_clusters = 0
        unsuccessful_alignments = set([])
        while received_gene_clusters < num_non_singleton_gene_clusters:
            try:
                gene_clusters_item = output_queue.get()

                if not gene_clusters_item['alignment_was_successful']:
                    unsuccessful_alignments.add(gene_clusters_item['name'])

                if gene_clusters_item:
                    # worker returns None if there is nothing to align
                    # we do not need to owerwrite it to gene_clusters_dict
                    gene_clusters_dict[gene_clusters_item['name']] = gene_clusters_item['entry']

                if self.debug:
                    print()
                    print(json.dumps(gene_clusters_item, indent=2))

                received_gene_clusters += 1
                self.progress.increment()
                self.progress.update("Processed %d of %d non-singleton GCs in %d threads." %
                    (received_gene_clusters, num_non_singleton_gene_clusters, self.num_threads))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        self.progress.end()

        return gene_clusters_dict, unsuccessful_alignments


    @staticmethod
    def alignment_worker(input_queue, output_queue, gene_clusters_dict, genomes_storage, align_with, run):
        # Note for future changes, this worker should not write anything to gene_clusters_dict
        # or genome_storage, changes will not be reflected to main process or other processes.

        aligner = aligners.select(align_with, quiet=True)

        # this instance of Run is here because we don't want to create this over and over again
        # in the loop down below. there also is another run instance the worker gets to make sure
        # it can report its own messages .. don't be confused we-do-not-age-discriminate-here padawan.
        r = terminal.Run()
        r.verbose = False

        # Main process needs to kill this worker after it receives all tasks because of this infinite loop
        while True:
            gene_cluster_name = input_queue.get(True)

            if len(gene_clusters_dict[gene_cluster_name]) == 1:
                # this sequence is a singleton and does not need alignment
                output_queue.put(None)
                continue

            gene_sequences_in_gene_cluster = []

            for gene_entry in gene_clusters_dict[gene_cluster_name]:
                sequence = genomes_storage.get_gene_sequence(gene_entry['genome_name'], gene_entry['gene_caller_id'])
                gene_sequences_in_gene_cluster.append(('%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id']), sequence),)

            alignment_was_successful = False

            # sometimes alignments fail, and because pangenomic analyses can take forever,
            # everything goes into the trash bin. to prevent that, here we have a try/except
            # block with lots of warnings if something goes wrong.
            try:
                alignments = aligner(run=r).run_stdin(gene_sequences_in_gene_cluster)
                alignment_was_successful = True
            except:
                # realm of sad face. before we continue to spam the user with error messages,
                # we turn our gene sequences to alignments without alignments. this worker will
                # report raw, unaligned sequences for this gene cluster as if they were aligned
                # so things will continue working operationally, and it will be on the user to
                # make sure they went through their results carefully.
                alignments = dict(gene_sequences_in_gene_cluster)

                # constructing our #sad:
                if anvio.DEBUG:
                    temp_file_path = filesnpaths.get_temp_file_path(prefix='ANVIO_GC_%s' % (gene_cluster_name))
                    with open(temp_file_path, 'w') as reporting_file:
                        for tpl in gene_sequences_in_gene_cluster:
                            reporting_file.write('>%s\n%s\n' % (tpl[0], tpl[1]))

                    progress.reset()
                    run.warning(f"'{temp_file_path}' contains raw sequences for {len(gene_sequences_in_gene_cluster)} genes "
                                f"{aligner.__name__} couldn't align :/", header=f"FAILED ALIGNMENT FOR GENE CLUSTER '{gene_cluster_name}'",
                                lc="green", nl_before=1)

            output = {'name': gene_cluster_name, 'alignment_was_successful': alignment_was_successful, 'entry': copy.deepcopy(gene_clusters_dict[gene_cluster_name])}
            for gene_entry in output['entry']:
                gene_entry['alignment_summary'] = utils.summarize_alignment(alignments['%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id'])])

            output_queue.put(output)


    def get_gene_clusters_de_novo(self):
        """Function to compute gene clusters de novo"""

        # get all amino acid sequences:
        combined_aas_FASTA_path = self.get_output_file_path('combined-aas.fa')
        self.genomes_storage.gen_combined_aa_sequences_FASTA(combined_aas_FASTA_path,
                                                             exclude_partial_gene_calls=self.exclude_partial_gene_calls)

        # get unique amino acid sequences:
        self.progress.new('Uniquing the output FASTA file')
        self.progress.update('...')
        unique_aas_FASTA_path, unique_aas_names_file_path, unique_aas_names_dict = utils.unique_FASTA_file(combined_aas_FASTA_path, store_frequencies_in_deflines=False)
        self.progress.end()
        self.run.info('Unique AA sequences FASTA', unique_aas_FASTA_path)

        # run search
        blastall_results = self.run_search(unique_aas_FASTA_path, unique_aas_names_dict)

        # generate MCL input from filtered blastall_results
        mcl_input_file_path = self.gen_mcl_input(blastall_results)

        # get clusters from MCL
        mcl_clusters = self.run_mcl(mcl_input_file_path)

        # we have the raw gene clusters dict, but we need to re-format it for following steps
        gene_clusters_dict = self.gen_gene_clusters_dict_from_mcl_clusters(mcl_clusters)
        del mcl_clusters

        return gene_clusters_dict


    def get_gene_clusters_from_gene_clusters_txt(self):
        """Function to recover gene clusters from gene-clusters-txt"""

        gene_clusters_dict = {}
        gene_clusters_txt = utils.get_TAB_delimited_file_as_dictionary(self.user_defined_gene_clusters, indexing_field=-1)

        genomes_in_gene_clusters_txt = set()

        for v in gene_clusters_txt.values():
            if v['gene_cluster_name'] not in gene_clusters_dict:
                gene_clusters_dict[v['gene_cluster_name']] = []

            gene_clusters_dict[v['gene_cluster_name']].append({'gene_caller_id': int(v['gene_caller_id']),
                                                               'gene_cluster_id': v['gene_cluster_name'],
                                                               'genome_name': v['genome_name'],
                                                               'alignment_summary': ''})

            genomes_in_gene_clusters_txt.add(v['genome_name'])

        genomes_only_in_gene_clusters_txt = [g for g in genomes_in_gene_clusters_txt if g not in self.genomes]
        genomes_only_in_genomes_storage = [g for g in self.genomes if g not in genomes_in_gene_clusters_txt]

        if len(genomes_only_in_gene_clusters_txt):
            raise ConfigError(f"Anvi'o run into an issue while processing your gene-clusters-txt file. It seems {len(genomes_only_in_gene_clusters_txt)} of "
                              f"{len(genomes_in_gene_clusters_txt)} genomes in your gene-clusters-txt file is not in the genomes-storage-db you have "
                              f"generated from your external- and/or internal-genomes with the program `anvi-gen-genomes-storage`. Here is the list of "
                              f"genome names that cause this issue: {', '.join(genomes_only_in_gene_clusters_txt)}")

        if len(genomes_only_in_genomes_storage):
            self.run.warning("Anvi'o observed something while processing your gene-clusters-txt file. It seems {len(genomes_only_in_genomes_storage)} of "
                             "{len(self.genomes)} genomes described in your genome-storage-db does not have any gene clusters described in the "
                             "gene-clusters-txt file. This may not be the end of the world, but it is a weird situation, and may lead to some "
                             "downstream issues. Anvi'o will continue working on your data, but if your computer sets itself on fire or something "
                             "please consider that it may be because of this situation.")

        return gene_clusters_dict


    def process(self):
        # start by processing the additional params user may have passed for the blast step
        self.process_additional_params()

        # load genomes from genomes storage
        self.load_genomes()

        # check sanity
        self.sanity_check()

        # gen pan_db
        self.generate_pan_db()

        # time to get the gene clusters. by default, we compute them de novo. but we also can
        # get them from the user themselves through gene-clusters-txt
        if self.user_defined_gene_clusters:
            gene_clusters_dict = self.get_gene_clusters_from_gene_clusters_txt()
        else:
            gene_clusters_dict = self.get_gene_clusters_de_novo()

        # compute alignments for genes within each gene_cluster (or don't)
        gene_clusters_dict, unsuccessful_alignments = self.compute_alignments_for_gene_clusters(gene_clusters_dict)

        if len(unsuccessful_alignments):
            if anvio.DEBUG:
               self.run.warning(f"The alignment of sequences failed for {len(unsuccessful_alignments)} of the total {len(gene_clusters_dict)} "
                                f"gene clusters anvi'o identified in your pangenome, and you can find the FASTA files for them in your logs "
                                f"above. SAD DAY FOR EVERYONE.", header="GENE CLUSTERS SEQUENCE ALIGNMENT WARNING")
            else:
               self.run.warning(f"SOME BAD NEWS :/ The alignment of sequences failed for {len(unsuccessful_alignments)} of the total "
                                f"{len(gene_clusters_dict)} gene clusters anvi'o identified in your pangenome :/ This often happens "
                                f"when you have extremely long genes (often due to improper gene calling across scaffolds) or gene "
                                f"clusters with extremely large number of genes (often happens with transposons or other mobile elements). "
                                f"This is not the end of the world, since most things will continue to work. But you will be missing some "
                                f"things for these gene clusters, including gene cluster functional and geometric homogeneity calculations "
                                f"(which will have the value of `-1` to help you parse them out later in your downstream analyses). Even "
                                f"the gene clusters failed at the alignment step will show up in all your downstream analyses. But the fact "
                                f"that some of your gene clusters will not have any alignments may affect your science downstream depending "
                                f"on what you wish to do with them (if you are not sure, come to anvi'o Discord, and tell us about what you "
                                f"want to do with your pangenome so we can discuss more). Additionally, if you would like to see what is going "
                                f"wrong with those gene clusters that are not aligned well, you can re-run the program with the `--debug` "
                                f"flag, upon which anvi'o will anvi'o will store all the sequences in gene clusters that filed to align "
                                f"into temporary FASTA files and will print out their locations. Here is the name of those gene clusters: "
                                f"{', '.join(unsuccessful_alignments)}.", header="GENE CLUSTERS SEQUENCE ALIGNMENT WARNING")

        # populate the pan db with results
        gene_clusters_dict = self.process_gene_clusters(gene_clusters_dict)

        # store gene clusters dict into the db
        self.store_gene_clusters(gene_clusters_dict)

        # generate a hierarchical clustering of gene clusters (or don't)
        self.gen_hierarchical_clustering_of_gene_clusters()

        # generate orderings of gene_clusters based on synteny of genes
        self.gen_synteny_based_ordering_of_gene_clusters(gene_clusters_dict)

        # populate layers additional data and orders
        self.populate_layers_additional_data_and_orders()

        # work with gene cluster homogeneity index
        self.populate_gene_cluster_homogeneity_index(gene_clusters_dict, gene_clusters_failed_to_align=unsuccessful_alignments)

        # let people know if they have too much data for their own comfort
        if len(gene_clusters_dict) > 20000 or len(self.genomes) > 150:
            if len(gene_clusters_dict) > 20000 and len(self.genomes) > 150:
                _ = "gene clusters and genomes"
            elif len(gene_clusters_dict) > 20000:
                _ = "gene clusters"
            else:
                _ = "genomes"

            self.run.warning(f"It seems you have a lot of {_} in this pan database :) It is all good! But please be aware that you may "
                             f"run into performance issues when you try to interactively visaulize these data using `anvi-display-pan`. "
                             f"In some cases it may even be impossible to do it, in fact. This is largely because the part of the "
                             f"anvi'o workflow to offer interactive access to a pangenomes is not designed to accommodate very"
                             f" large number of {_}, but rather enable in-depth exploratory analyses of pangenomes interactively. You still "
                             f"can work with large pangenomes via the command line utilities and do a lot of science with them. If you "
                             f"are unable to work with the interactive interface and it is critical for you, you have multiple options, "
                             f"You can use the `--min-occurrence` flag to reduce the number of gene clusters, or use the program "
                             f"`anvi-dereplicate-genomes` in an attempt to reduce the number of redundant genomes in your analysis. "
                             f"If you are unsure what would be the best game plan for you, you can consider coming to the anvi'o Discord "
                             f"channel and consult the opinion of the anvi'o community. Despite all these, it is still a good idea to run "
                             f"`anvi-display-pan` and see what it says first.", lc="cyan", header="FRIENDLY WARNING")

        # done
        self.run.info_single(f"Your pangenome is ready with a total of {pp(len(gene_clusters_dict))} gene clusters across "
                             f"{len(self.genomes)} genomes 🎉", mc="green", nl_after=1)


        self.run.quit()


class Pangraph():
    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # important Anvi'o artifacts that are necessary for successfull execution
        self.pan_db = A('pan_db')
        self.external_genomes_txt = A('external_genomes')
        self.genomes_storage_db = A('genomes_storage')

        # additional optinal input variables
        self.testing_yaml = A('testing_yaml')
        self.pan_graph_json = A('pan_graph_json')
        self.gene_additional_data = A('gene_additional_data')
        self.gc_additional_data = A('gc_additional_data')

        # additional variables for presetting the UI values
        self.max_edge_length_filter = A('max_edge_length_filter')
        self.gene_cluster_grouping_threshold = A('gene_cluster_grouping_threshold')
        self.groupcompress = A('gene_cluster_grouping_compression')
        self.ungroup = A('ungrouping_area')
        self.pass_filter = 0

        # additional variables for special cases e.g. a user wants to tune the tool
        # in a very specific direction
        # TODO: skip genomes might be a very useful function
        self.priority_genome = A('priority_genome')
        
        if A('genome_names'):
            self.genomes_names = A('genome_names').split(',')
        else:
            self.genomes_names = []

        if A('genome_reverse'):
            self.genomes_reverse = A('genome_reverse').split(',')
        else:
            self.genomes_reverse = []

        # the different data storage related variables e.g. input and output files
        # TODO: DB storage is not yet implemented -> will be GRAPH.db at one point
        self.skip_storing_in_pan_db = True
        self.json_output_file_path = A('output_pan_graph_json')

        self.output_graphml = A('output_graphml')
        self.output_yaml = A('output_testing_yaml')
        self.output_summary = A('output_dir_summary')
        self.output_graphics = A('output_dir_graphics')

        self.output_raw_gc_additional_data = A('output_raw_gc_additional_data')
        self.output_raw_gene_additional_data = A('output_raw_gene_additional_data')

        # learn what gene annotation sources are present across all genomes if we
        # are running things in normal mode
        if not self.testing_yaml and not self.pan_graph_json:
            self.functional_annotation_sources_available = DBInfo(self.genomes_storage_db, expecting='genomestorage').get_functional_annotation_sources()
        else:
            self.functional_annotation_sources_available = []

        # dictionary containing the layer information that will be saved in the
        # json file to be used in the front end
        self.data_table_dict = {}

        # this is the dictionary that wil keep all data that is going to be loaded
        # from anvi'o artifacts
        self.gene_synteny_data_dict = {}
        self.genomes = []

        # these are the main graph data structures for the different steps of the tools execution
        self.initial_graph = nx.DiGraph()
        self.pangenome_graph = nx.DiGraph()
        self.edmonds_graph = nx.DiGraph()
        self.ancest = nx.DiGraph()

        # additional self variables that are used during the tools execution
        self.single_copy_core = set()
        self.grouping = {}
        self.project_name = ''
        self.global_y = 0
        self.global_x = 1
        self.global_x_offset = 0
        self.k = 0
        self.genome_gc_occurence = {}
        self.ghost = 0
        self.debug = False
        self.jsondata = {}
        self.position = {}
        self.offset = {}
        self.x_list = {}
        self.path = {}
        self.edges = {}
        self.layers = []
        self.removed = set()

    def sanity_check(self):

        """Sanity check for incompatible settings, like skip storing and no output path"""

        if self.skip_storing_in_pan_db and not self.json_output_file_path:
            raise ConfigError("You are initializing the Pangraph class with `--skip-storing-in-pan-db` without an `--output-pan-graph-json` "
                            "parameter for the graph results to be stored. Please set an output file path so anvi'o has at least one "
                            "way to store results.")
    
        # if not self.genomes_names:
        #     raise ConfigError("Please be nice enough to use a genome name that is present in your dataset, you know, to let the code take the wheel.")


    def process(self):
        """Primary driver function for the class"""

        self.sanity_check()

        # the pathway over the pan_graph_json serves the purpose of updating the
        # json file
        if self.pan_graph_json:
            self.load_graph_from_json_file()

        if not self.pan_graph_json and self.testing_yaml:
            # skip sanity check EVERYTHING in case of testing mode
            self.prepare_yaml()

        if not self.pan_graph_json and not self.testing_yaml:
            # populate self.gene_synteny_data_dict
            self.get_gene_synteny_data_dict()

        if not self.pan_graph_json and self.output_yaml:
            self.export_gene_synteny_to_yaml()

        if not self.pan_graph_json and not self.output_yaml:
            # contextualize paralogs
            # TODO Incorporate gene direction
            self.run_contextualize_paralogs_algorithm()

            # build graph
            self.build_graph()

            # reconnect open leaves in the graph to generate
            # a flow network from left to right
            self.run_tree_to_flow_network_algorithm()

            self.prepare_synteny_graph()

        if not self.output_yaml:
            self.run_synteny_layout_algorithm()

        if not self.pan_graph_json and not self.output_yaml and not self.testing_yaml:
            # run Alex's layout algorithm
            self.generate_data_table()

        if not self.output_yaml and self.output_raw_gc_additional_data:
            self.get_additional_gc_layer_table()

        if not self.output_yaml and self.output_raw_gene_additional_data:
            self.get_additional_gene_layer_table()

        if not self.pan_graph_json and not self.output_yaml and self.gc_additional_data:
            self.add_additional_gc_layer_values()

        if not self.pan_graph_json and not self.output_yaml and self.gene_additional_data:
            self.add_additional_gene_layer_values()

        if not self.pan_graph_json and not self.output_yaml and self.output_summary:
            self.get_genomic_motifs()

        if not self.pan_graph_json and not self.output_yaml:
            self.get_json_dict_for_graph()

        if self.pan_graph_json:
            self.update_json_dict()

        if not self.output_yaml:
            # store network in the database
            self.store_network()


    def prepare_yaml(self):
        yaml_genomes = utils.get_yaml_as_dict(self.testing_yaml)
        self.project_name = 'YAML TEST'

        for genome in yaml_genomes.keys():
            self.gene_synteny_data_dict[genome] = {}

            if genome in self.genomes_names:
                self.genomes.append(genome)
                self.run.info_single(f"Loaded genome {genome}.")
            else:
                self.run.info_single(f"Skipped genome {genome}.")

            for i, gcs in enumerate(yaml_genomes[genome]):

                contig = genome + '_' + str(i).zfill(5)

                self.gene_synteny_data_dict[genome][contig] = {}

                for j, gc in enumerate(gcs.split(' ')):
                    if gc.endswith('!'):
                        gc = gc[:-1]
                        direction = 'r'
                    else:
                        direction = 'f'

                    self.gene_synteny_data_dict[genome][contig][j] = {
                        'gene_cluster_name': gc,
                        'gene_cluster_id': '',
                        'direction': direction,
                        'rev_compd': 'False',
                        'max_num_paralogs': 0
                    }


    def export_gene_synteny_to_yaml(self):

        self.run.warning(None, header="Export the gene synteny to yaml file", lc="green")
        gene_synteny_yaml_dict = {}
        for genome in self.gene_synteny_data_dict.keys():
            gene_synteny_yaml_dict[genome] = []
            for contig in self.gene_synteny_data_dict[genome].keys():
                contig_data = self.gene_synteny_data_dict[genome][contig]

                contig_string = []

                for gene_call in contig_data.keys():
                    if contig_data[gene_call]['direction'] == 'f':
                        contig_string += [contig_data[gene_call]['gene_cluster_name']]
                    else:
                        contig_string += [contig_data[gene_call]['gene_cluster_name'] + '!']

                gene_synteny_yaml_dict[genome] += [' '.join(contig_string)]

        utils.save_dict_as_yaml(gene_synteny_yaml_dict, self.output_yaml)

        self.run.info_single("Done")


    def get_gene_synteny_data_dict(self):
        """A function to reduce a comprehensive data structure from anvi'o artifacts for
           downstream analyses.
        """
        self.run.warning(None, header="Loading data from database", lc="green")

        filesnpaths.is_file_tab_delimited(self.external_genomes_txt)

        if not utils.is_all_columns_present_in_TAB_delim_file(["name","contigs_db_path"], self.external_genomes_txt):
            raise ConfigError("Your external genomes file does not seem to contain that anvi'o expects to find "
                              "in an external genomes file :/")

        pan_db = dbops.PanSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        self.project_name = pan_db.p_meta['project_name']

        pan_db.init_gene_clusters()
        pan_db.init_gene_clusters_functions_summary_dict()
        pan_db.init_items_additional_data()

        gene_cluster_dict = pan_db.gene_callers_id_to_gene_cluster
        additional_info_cluster = pan_db.items_additional_data_dict

        external_genomes = pd.read_csv(self.external_genomes_txt, header=0, sep="\t", names=["name","contigs_db_path"])
        external_genomes.set_index("name", inplace=True)

        for genome, contigs_db_path in external_genomes.iterrows():

            if genome in self.genomes_names or not self.genomes_names:

                self.genomes.append(genome)

                args = argparse.Namespace(contigs_db=contigs_db_path.item())
                contigs_db = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

                caller_id_cluster = gene_cluster_dict[genome]
                caller_id_cluster_df = pd.DataFrame.from_dict(caller_id_cluster, orient="index", columns=["gene_cluster_name"]).rename_axis("gene_caller_id").reset_index()
                caller_id_cluster_df["gene_cluster_id"] = ""

                contigs_db.init_functions()
                gene_function_calls_df = pd.DataFrame.from_dict(contigs_db.gene_function_calls_dict, orient="index").rename_axis("gene_caller_id").reset_index()

                all_gene_calls = caller_id_cluster_df['gene_caller_id'].values.tolist()
                genes_in_contigs_df = pd.DataFrame.from_dict(contigs_db.get_sequences_for_gene_callers_ids(all_gene_calls, include_aa_sequences=True, simple_headers=True)[1], orient="index").rename_axis("gene_caller_id").reset_index()
                additional_info_df = pd.DataFrame.from_dict(additional_info_cluster, orient="index").rename_axis("gene_cluster_name").reset_index()
                additional_info_df.rename(columns={entry: entry + "_known" for entry in self.functional_annotation_sources_available}, inplace=True)

                joined_contigs_df = caller_id_cluster_df.merge(genes_in_contigs_df, on="gene_caller_id", how="left").merge(gene_function_calls_df, on="gene_caller_id", how="left").merge(additional_info_df, on="gene_cluster_name", how="left")

                if genome in self.genomes_reverse:
                    joined_contigs_df.sort_values(["contig", "start", "stop"], axis=0, ascending=True, inplace=True)
                else:
                    joined_contigs_df.sort_values(["contig", "start", "stop"], axis=0, ascending=False, inplace=True)

                joined_contigs_df.set_index(["contig", "gene_caller_id"], inplace=True)

                self.gene_synteny_data_dict[genome] = joined_contigs_df.fillna("None").groupby(level=0).apply(lambda df: df.xs(df.name).to_dict("index")).to_dict()

                self.run.info_single(f"Loaded genome {genome}.")

            else:
                self.run.info_single(f"Skipped genome {genome}.")

        if not self.genomes:
            raise ConfigError(f"Please keep at least one genome in the dataset. With the current setting you skip over all.")
        else:
            if not self.genomes_names:
                self.genomes_names = self.genomes

        self.run.info_single("Done")


    def context_distance(self, a, b):
        r_val = 0
        f_val = 0
        div = len(a)

        if not a[int(len(a) / 2)] == b[int(len(b) / 2)]: 
            raise ConfigError(f"This is very weird and should not have happend, we are very sorry :(.")

        if set(self.genome_gc_occurence[tuple(a)].keys()).isdisjoint(set(self.genome_gc_occurence[tuple(b)].keys())):
            for i, (n, m) in enumerate(zip(a, b)):

                if not i == int(len(a) / 2):
                    if (n[0] == '-' and m[0] == '-') or (n[0] == '+' and m[0] == '+'):
                        if n[1] != m[1]:
                            r_val += 0.5
                    elif n[0] == '-' or m[0] == '-' or n[0] == '+' or m[0] == '+':
                        r_val += 0
                    elif n == m:
                        r_val += 1

            for i, (n, m) in enumerate(zip(a, b[::-1])):

                if not i == int(len(a) / 2):
                    if (n[0] == '-' and m[0] == '-') or (n[0] == '+' and m[0] == '+'):
                        if n[1] != m[1]:
                            f_val += 0.5
                    elif n[0] == '-' or m[0] == '-' or n[0] == '+' or m[0] == '+':
                        f_val += 0
                    elif n == m:
                        f_val += 1

            if r_val >= f_val:
                return (1 - r_val / div)
            else:
                return (1 - f_val / div)

        elif tuple(a) == tuple(b):
            return (0)

        else:
            return (1)


    def context_split(self, label, gcp):

        if len(label) == 1:
            result = list(zip(label, [1]))
            return (result)

        else:
            X = cdist(np.asarray(label), np.asarray(label), metric=self.context_distance)
            condensed_X = squareform(X)

            Z = linkage(condensed_X, 'ward')

            for t in sorted(set(sum(Z.tolist(), [])), reverse=True):
                clusters = fcluster(Z, t, criterion='distance')

                valid = True
                for c in set(clusters.tolist()):
                    pos = np.where(clusters == c)[0]

                    for i, j in it.combinations(pos, 2):

                        if X[i][j] == 1.:
                            valid = False

                if valid is True:
                    break

            result = list(zip(label, clusters))

            if self.output_graphics:
                fig = plt.figure(figsize=(25, 10))
                ax = plt.axes()
                dn = dendrogram(Z, ax=ax, labels=label, orientation='right')
                ax.axvline(linestyle='--', x=t)
                plt.tight_layout()
                fig.savefig(self.output_graphics + gcp + '.pdf')

            return (result)


    def run_contextualize_paralogs_algorithm(self):
        """A function that resolves the graph context of paralogs based on gene synteny information across genomes"""
        self.run.warning(None, header="Select paralog context", lc="green")

        unresolved = True
        solved = set()

        genome_gc_order_values = {}
        genome_gc_order_left_placeholder = {}
        genome_gc_order_right_placeholder = {}

        while unresolved:

            unresolved = False
            drop = set()

            for genome in self.gene_synteny_data_dict.keys():
                for contig in self.gene_synteny_data_dict[genome].keys():
                    genome_gc_order = [self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_name"] for gene_call in self.gene_synteny_data_dict[genome][contig].keys()]
                
                    if self.k == 0:
                        if genome not in genome_gc_order_values:
                            genome_gc_order_values[genome] = {contig: genome_gc_order}
                            genome_gc_order_left_placeholder[genome] = {contig: '-0'}
                            genome_gc_order_right_placeholder[genome] = {contig: '+0'}
                        else:
                            items = list(genome_gc_order_values[genome].values())
                            times = items.count(genome_gc_order)
                            if times == 0:
                                genome_gc_order_values[genome][contig] = genome_gc_order
                                genome_gc_order_left_placeholder[genome][contig] = '-0'
                                genome_gc_order_right_placeholder[genome][contig] = '+0'
                            else:
                                genome_gc_order_values[genome][contig] = genome_gc_order
                                genome_gc_order_left_placeholder[genome][contig] = '-' + chr(ord('a')+times-1)
                                genome_gc_order_right_placeholder[genome][contig] = '+' + chr(ord('a')+times-1)

                    left_place = genome_gc_order_left_placeholder[genome][contig]
                    right_place = genome_gc_order_right_placeholder[genome][contig]

                    for i in range(0, len(genome_gc_order)):

                        start = (i - self.k) if (i - self.k) >= 0 else 0
                        stop = (i + self.k + 1) if (i + self.k + 1) <= len(genome_gc_order) else len(genome_gc_order)
                        entry = genome_gc_order[start:stop]

                        if len(entry) == 1 + (2 * self.k):
                            gc_k = tuple(entry)
                        elif start == 0 and stop == len(genome_gc_order):
                            gc_k = tuple([left_place] * (self.k - i) + entry + [right_place] * ((i + self.k + 1) - len(genome_gc_order)))
                        elif start == 0:
                            gc_k = tuple([left_place] * (self.k - i) + entry)
                        elif stop == len(genome_gc_order):
                            gc_k = tuple(entry + [right_place] * ((i + self.k + 1) - len(genome_gc_order)))
                        else:
                            raise ConfigError(f"The search frame of entry {entry} is malformed. It is not tuple of size {1 + (2 * self.k)}"
                                              f"neiter is the frame overlapping with the beginning and the end of a very short contig"
                                              f"nor with either of those. This is weird and should not happen, please check your input data"
                                              f"and make sure there is nothing very weird going on there. Sorry :/")

                        if len(gc_k) != 1 + (2 * self.k):
                            raise ConfigError(f"Wow, for some unforseeable reason the length of the search frame is not matching the"
                                              f"current iterration. I thought this is impossible and I'm very curious on how you accomplished"
                                              f"that. Jokes aside errors can happen and I probably forgot some specific exception. Sorry...")

                        gc = gc_k[int(len(gc_k) / 2)]
                        if gc not in solved:

                            if gc_k not in self.genome_gc_occurence.keys():
                                self.genome_gc_occurence[gc_k] = {genome: 1}

                            else:
                                if genome not in self.genome_gc_occurence[gc_k].keys():
                                    self.genome_gc_occurence[gc_k][genome] = 1
                                else:
                                    self.genome_gc_occurence[gc_k][genome] += 1

            for gc, genome_gc_frequency in self.genome_gc_occurence.items():
                if max(genome_gc_frequency.values()) > 1:

                    unresolved = True
                    drop.add(gc[int(len(gc)/2)])

                elif self.k == 0 and len(genome_gc_frequency) == len(self.genomes):
                    self.single_copy_core.add(gc)

            if self.k == 0:
                self.paralog_dict = copy.deepcopy(self.genome_gc_occurence)

            if unresolved:
                keys = list(self.genome_gc_occurence.keys())
                for gc in keys:
                    if gc[int(len(gc)/2)] in drop:
                        self.genome_gc_occurence.pop(gc)

                solved = set([gc[int(len(gc)/2)] for gc in self.genome_gc_occurence.keys()])
                self.k += 1

        g_cleaned = {}
        gcs = set([gc[int(len(gc)/2)] for gc in self.genome_gc_occurence.keys()])

        for gcp in gcs:
            label = [gc for gc in self.genome_gc_occurence.keys() if gc[int(len(gc)/2)] == gcp]
            context = self.context_split(label, gcp)

            for name, cluster in context:
                g_cleaned[name] = gcp + '_' + str(cluster)

        self.run.info_single(f"{pp(self.k+1)} iteration(s) to expand {pp(len(self.paralog_dict.keys()))} GCs to {pp(len(self.genome_gc_occurence.keys()))} GCs without paralogs")
        syn_calls = 0

        for genome in self.gene_synteny_data_dict.keys():

            for contig in self.gene_synteny_data_dict[genome].keys():
                genome_gc_order = [(self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_name"], gene_call) for gene_call in self.gene_synteny_data_dict[genome][contig].keys()]

                left_place = genome_gc_order_left_placeholder[genome][contig]
                right_place = genome_gc_order_right_placeholder[genome][contig]

                for i in range(0, len(genome_gc_order)):
                    start = i-self.k if i-self.k >= 0 else 0
                    stop = i+self.k+1 if i+self.k+1 <= len(genome_gc_order) else len(genome_gc_order)
                    entry = [item[0] for item in genome_gc_order[start:stop]]
                    gene_call = genome_gc_order[i][1]
                    name = genome_gc_order[i][0]

                    if len(entry) == 1 + (2 * self.k):
                        gc_k = tuple(entry)
                    elif start == 0 and stop == len(genome_gc_order):
                        gc_k = tuple([left_place] * (self.k - i) + entry + [right_place] * ((i + self.k + 1) - len(genome_gc_order)))
                    elif start == 0:
                        gc_k = tuple([left_place] * (self.k - i) + entry)
                    elif stop == len(genome_gc_order):
                        gc_k = tuple(entry + [right_place] * ((i + self.k + 1) - len(genome_gc_order)))
                    else:
                        raise ConfigError(f"The search frame of entry {entry} is malformed. It is not tuple of size {1 + (2 * self.k)}"
                                          f"neiter is the frame overlapping with the beginning and the end of a very short contig"
                                          f"nor with either of those. This is weird and should not happen, please check your input data"
                                          f"and make sure there is nothing very weird going on there. Sorry :/")

                    for j in range(0, self.k+1):

                        gc_group = gc_k[int(len(gc_k)/2)-j:int(len(gc_k)/2)+j+1]

                        if gc_group in self.genome_gc_occurence.keys():

                            self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_id"] = g_cleaned[gc_group]
                            syn_calls += 1
                            break

                        else:
                            pass

        num_calls = 0
        for _, value in self.genome_gc_occurence.items():
            num_calls += sum(value.values())

        if num_calls != syn_calls:
            raise ConfigError(f"It looks like {abs(num_calls - syn_calls)} calls were not  captured by the algorithm. Proceeding from"
                              f"here is a very bad idea. We will try to do better next time. For now please check your data to make"
                              f"sure you did not try some very crazy stuff.")

        self.run.info_single("Done")

    def add_node_to_graph(self, gene_cluster, name, info):

        if not self.initial_graph.has_node(gene_cluster):
            self.initial_graph.add_node(
                gene_cluster,
                name=name,
                pos=(0, 0),
                weight=1,
                layer={},
                genome=info
            )

        else:
            self.initial_graph.nodes[gene_cluster]['weight'] += 1
            self.initial_graph.nodes[gene_cluster]['genome'].update(info)


    def add_edge_to_graph(self, gene_cluster_i, gene_cluster_j, info, weight_add=0):

        if self.priority_genome in info.keys():
            weight_add += 100
        # elif gene_cluster_i in self.single_copy_core and gene_cluster_j in self.single_copy_core:
        #     weight_add += 2
        # elif gene_cluster_i in self.single_copy_core or gene_cluster_j in self.single_copy_core:
        #     weight_add += 1
        else:
            weight_add += 0

        draw = {genome: {'gene_call': -1} for genome in info.keys()}

        if not self.initial_graph.has_edge(*(gene_cluster_i, gene_cluster_j)):
            self.initial_graph.add_edge(
                *(gene_cluster_i, gene_cluster_j),
                weight=1 + weight_add,
                genome=draw,
                bended=[],
                direction='R'
            )

        else:
            self.initial_graph[gene_cluster_i][gene_cluster_j]['weight'] += 1 + weight_add
            self.initial_graph[gene_cluster_i][gene_cluster_j]['genome'].update(draw)


    def build_graph(self):

        factor = 1.0 / 2
        decisison_making = {}

        # Unfortunately this part here is very arbitiary but necessary. Find better way later!
        for genome in self.genomes_names:
            decisison_making[genome] = factor
            factor /= 2

        self.run.warning(None, header="Building directed gene cluster graph G", lc="green")

        for genome in self.gene_synteny_data_dict.keys():

            for contig in self.gene_synteny_data_dict[genome].keys():
                gene_cluster_kmer = []
                for gene_call in self.gene_synteny_data_dict[genome][contig].keys():
                    gene_cluster_kmer.append((self.gene_synteny_data_dict[genome][contig][gene_call]['gene_cluster_id'], self.gene_synteny_data_dict[genome][contig][gene_call]['gene_cluster_name'], {genome: {'contig':contig, 'gene_call':gene_call, **self.gene_synteny_data_dict[genome][contig][gene_call]}}))

                if len(gene_cluster_kmer) > 1:
                    gene_cluster_pairs = map(tuple, zip(gene_cluster_kmer, gene_cluster_kmer[1:]))
                    first_pair = next(gene_cluster_pairs)

                    weight_add = 0 #decisison_making[genome]

                    self.add_node_to_graph(first_pair[0][0], first_pair[0][1], first_pair[0][2])
                    self.add_node_to_graph(first_pair[1][0], first_pair[1][1], first_pair[1][2])
                    self.add_edge_to_graph(first_pair[0][0], first_pair[1][0], first_pair[1][2], weight_add)

                    for gene_cluster_pair in gene_cluster_pairs:

                        weight_add = 0 #decisison_making[genome]

                        self.add_node_to_graph(gene_cluster_pair[1][0], gene_cluster_pair[1][1], gene_cluster_pair[1][2])
                        self.add_edge_to_graph(gene_cluster_pair[0][0], gene_cluster_pair[1][0], gene_cluster_pair[1][2], weight_add)

                else:
                    self.add_node_to_graph(gene_cluster_kmer[0][0], gene_cluster_kmer[0][1], gene_cluster_kmer[0][2])

        self.run.info_single(f"Adding {pp(len(self.initial_graph.nodes()))} nodes and {pp(len(self.initial_graph.edges()))} edges to G")

        if self.pass_filter == 0:
            self.run.info_single("Setting algorithm to 'no pass filter'")
        else:
            self.run.info_single(f"Setting algorithm to 'remove edges < {pp(self.pass_filter)}'")

        edges_to_remove = [(i, j) for i,j,k in self.initial_graph.edges(data=True) if k['weight'] < self.pass_filter]

        self.run.info_single(f"Removed {pp(len(edges_to_remove))} edges due to pass filter.")

        self.initial_graph.remove_edges_from(edges_to_remove)

        connectivity = nx.is_connected(self.initial_graph.to_undirected())
        self.run.info_single(f"Connectivity is {connectivity}")

        if connectivity == False:

            weakly_components = list(nx.weakly_connected_components(self.initial_graph))
            self.run.info_single(f"Graph contains {pp(len(weakly_components))} independant components.")

            self.pangenome_graph = nx.DiGraph(self.initial_graph.subgraph(max(weakly_components, key=len)))
            self.run.info_single(f"Keeping longest subgraph with {pp(len(self.pangenome_graph.nodes()))} nodes and {pp(len(self.pangenome_graph.edges()))} edges")

            connectivity = nx.is_connected(self.pangenome_graph.to_undirected())
            self.run.info_single(f"Connectivity is {connectivity}")
        else:
            self.pangenome_graph = nx.DiGraph(self.initial_graph)

        self.run.info_single("Done")
        self.run.warning(None, header="Building maximum branching graph M of G", lc="green")

        selfloops = list(nx.selfloop_edges(self.pangenome_graph))
        self.run.info_single(f"Found and removed {pp(len(selfloops))} selfloop edge(s)")
        self.pangenome_graph.remove_edges_from(selfloops)

        edmonds = nx.algorithms.tree.branchings.Edmonds(self.pangenome_graph)
        self.edmonds_graph = edmonds.find_optimum(
            attr="weight",
            default=1,
            kind="max",
            style="arborescence",
            preserve_attrs=False,
            partition=None,
        )

        if not nx.algorithms.tree.recognition.is_arborescence(self.edmonds_graph):
            self.run.info_single('No maximum aborescence. Entering fallback mode.')
            edmonds_sub_graph = max(nx.weakly_connected_components(self.edmonds_graph), key=len)
            self.edmonds_graph = nx.DiGraph(self.edmonds_graph.subgraph(edmonds_sub_graph))
            self.pangenome_graph = nx.DiGraph(self.pangenome_graph.subgraph(edmonds_sub_graph))
            if nx.algorithms.tree.recognition.is_arborescence(self.edmonds_graph):
                self.run.info_single('Found aborescence.')

                #TODO calculate number of missing edges
                self.run.info_single(f'{len(self.edmonds_graph.nodes())-len(edmonds_sub_graph)} nodes removed to capture synteny.')
                self.run.info_single('Be warned the quality of the pangraph might have suffered.')
                self.run.info_single('Proceed.')
            else:
                raise ConfigError(f"I'm very sorry to inform you that your data is not solvable by the current version of"
                                  f"anvi'o pangraph. The fallback mode tried to solve your dataset by sacrificing some"
                                  f"of the included information, but at this scale it will not lead to a acceptable result :(")

        else:
            self.run.info_single('Found aborescence. Proceed.')

        nx.set_edge_attributes(self.edmonds_graph, {(i, j): d for i, j, d in self.pangenome_graph.edges(data=True) if (i, j) in self.edmonds_graph.edges()})
        nx.set_node_attributes(self.edmonds_graph, {k: d for k, d in self.pangenome_graph.nodes(data=True) if k in self.edmonds_graph.nodes()})

        add_start = []
        for node in self.edmonds_graph.nodes():
            if len(list(self.edmonds_graph.predecessors(node))) == 0:
                add_start.append(node)

        for graph in [self.pangenome_graph, self.edmonds_graph]:
            graph.add_node(
                'start',
                name='start',
                pos=(0, 0),
                weight=len(self.genomes),
                layer={},
                genome={genome: {'gene_call': -1, 'contig': '', 'direction': ''} for genome in self.genomes}
            )

            for u in add_start:

                weight = graph.nodes[u]['weight']
                genomes = graph.nodes[u]['genome'].keys()

                graph.add_edge(
                    *('start', u),
                    genome={genome: {'gene_call': -1} for genome in genomes},
                    weight=weight,
                    bended=[],
                    direction='R'
                )

        self.run.info_single(f"Removing {pp(len(self.pangenome_graph.edges()) - len(self.edmonds_graph.edges()))} edges from G to create M")
        self.run.info_single("Done")

    def find_next_branching_point(self, current):
        if current != 'start':
            pred = list(self.edmonds_graph.predecessors(current))[0]
            pred_successors = list(self.edmonds_graph.successors(pred))
            if len(pred_successors) > 1:
                return(pred)
            else:
                self.find_next_branching_point(pred)


    def mean_edmonds_graph_path_weight(self, source, target):
        path = nx.shortest_path(G=self.edmonds_graph, source=source, target=target, weight='weight')
        path_weight = nx.path_weight(G=self.edmonds_graph, path=path, weight='weight')
        return(path_weight)


    def get_edge(self, node_i, node_j, reverse = False):

        if reverse == False:
            pangenome_graph_edge_data = self.pangenome_graph.get_edge_data(node_i, node_j)
            return(node_i, node_j, pangenome_graph_edge_data)
        else:
            pangenome_graph_edge_data = {y:z if y != 'direction' else 'L' for y,z in self.pangenome_graph.get_edge_data(node_i, node_j).items()}
            return(node_j, node_i, pangenome_graph_edge_data)


    def get_leaves(self, current_branch_successor, edmonds_graph_successors):
        leaves = set()
        while True:
            successors = edmonds_graph_successors[current_branch_successor]

            if successors:
                if len(successors) > 1:
                    for successor in successors:
                        leaves.update(self.get_leaves(successor, edmonds_graph_successors))
                    break
                else:
                    current_branch_successor = successors[0]
            else:
                leaves.update(set([current_branch_successor]))
                break

        return(leaves)


    def edge_check(self, node_i, node_j, data):

        new_data = copy.deepcopy(data)

        if self.edmonds_graph.has_edge(node_i, node_j):
            old_data = self.edmonds_graph[node_i][node_j]
            new_data['genome'].update(old_data['genome'])
            new_data['weight'] += old_data['weight']
            new_data['direction'] = 'B'

        return(new_data)


    def run_tree_to_flow_network_algorithm(self):

        self.run.warning(None, header="Building flow network F from M and G", lc="green")

        edmonds_graph_edges = set(self.edmonds_graph.edges())
        edmonds_graph_nodes = set(self.edmonds_graph.nodes())

        pangenome_graph_edges = set(self.pangenome_graph.edges())
        pangenome_graph_nodes = set(self.pangenome_graph.nodes())

        edmonds_graph_removed_edges = pangenome_graph_edges - edmonds_graph_edges
        edmonds_graph_end = max([(self.mean_edmonds_graph_path_weight('start', node), node) for node in edmonds_graph_nodes if len(list(self.edmonds_graph.successors(node))) == 0])[1]

        for graph in [self.pangenome_graph, self.edmonds_graph]:

            graph.add_node(
                'stop',
                name='stop',
                pos=(0, 0),
                weight=len(self.genomes),
                layer={},
                genome={genome: {'gene_call': -1, 'contig': '', 'direction': ''} for genome in self.genomes}
            )
            graph.add_edge(
                *(edmonds_graph_end, 'stop'),
                genome={genome: {'gene_call': -1} for genome in self.genomes},
                weight=len(self.genomes),
                bended=[],
                direction='R'
            )

        edmonds_graph_nodes.add('stop')
        edmonds_graph_predecessors = {edmonds_graph_node: list(self.edmonds_graph.predecessors(edmonds_graph_node))[0] for edmonds_graph_node in edmonds_graph_nodes if edmonds_graph_node != 'start'}
        edmonds_graph_successors = {edmonds_graph_node: list(self.edmonds_graph.successors(edmonds_graph_node)) for edmonds_graph_node in edmonds_graph_nodes}
        pangenome_graph_successors = {pangenome_graph_node: list(self.pangenome_graph.successors(pangenome_graph_node)) for pangenome_graph_node in pangenome_graph_nodes}
        edmonds_graph_distances = {edmonds_graph_node: self.mean_edmonds_graph_path_weight('start', edmonds_graph_node) for edmonds_graph_node in edmonds_graph_nodes}

        i = 0
        resolved_nodes = set(['stop'])

        pred = ''
        x = 0
        self.progress.new("Solving complex graph")
        while len(resolved_nodes) != len(pangenome_graph_nodes) + 1 or x < 1:
            if len(resolved_nodes) == len(pangenome_graph_nodes) + 1:
                x += 1

            i += 1

            visited_nodes = set(['stop'])
            current_node = 'stop'
            while len(visited_nodes) != len(pangenome_graph_nodes) + 1:

                self.progress.update(f"{str(len(resolved_nodes)).rjust(len(str(len(pangenome_graph_nodes) + 1)), ' ')} / {len(pangenome_graph_nodes) + 1}")

                if pred:
                    current_branch_root = pred
                    pred = ''
                else:
                    current_branch_root = edmonds_graph_predecessors[current_node]

                current_forward_connected = []
                current_backward_connected = []
                successor_branch_leaves = set()

                for current_branch_successor in edmonds_graph_successors[current_branch_root]:
                    if current_branch_successor not in visited_nodes and current_branch_successor != current_node:
                        successor_branch_leaves.update(self.get_leaves(current_branch_successor, edmonds_graph_successors))

                if not successor_branch_leaves:
                    current_node = current_branch_root
                else:
                    current_node = max([(edmonds_graph_distances[successor_branch_leaf], successor_branch_leaf) for successor_branch_leaf in successor_branch_leaves])[1]

                if current_node in resolved_nodes:
                    connected = True
                else:
                    connected = False

                if connected != True or x == 1:
                    for current_node_successor in pangenome_graph_successors[current_node]:

                        if current_node_successor in resolved_nodes:

                            if current_node_successor in nx.ancestors(self.edmonds_graph, current_node) or (current_node_successor not in visited_nodes and current_node_successor not in resolved_nodes):
                                if (current_node, current_node_successor) in edmonds_graph_removed_edges:
                                    current_backward_connected.append(current_node_successor)

                            else:
                                if (current_node, current_node_successor) in edmonds_graph_removed_edges:
                                    current_forward_connected.append(current_node_successor)
                                    connected = True

                                else:
                                    connected = True

                    if connected == False:
                        if len(list(self.pangenome_graph.successors(current_node))) == 0:
                            pangenome_graph_edge_data = {
                                'genome':{genome: {'gene_call': -1} for genome in self.genomes},
                                'weight':len(self.genomes),
                                'bended': [],
                                'direction': 'R'
                            }

                            new_data = self.edge_check(current_node, 'stop', pangenome_graph_edge_data)
                            self.edmonds_graph.add_edge(current_node, 'stop', **new_data)
                            connected = True

                    if connected == True:
                        for current_forward in current_forward_connected:
                            node_i, node_j, data = self.get_edge(current_node, current_forward, reverse = False)

                            new_data = self.edge_check(node_i, node_j, data)
                            self.edmonds_graph.add_edge(node_i, node_j, **new_data)
                            edmonds_graph_removed_edges.remove((current_node, current_forward))

                        for current_backward in current_backward_connected:
                            node_i, node_j, data = self.get_edge(current_node, current_backward, reverse = True)

                            new_data = self.edge_check(node_i, node_j, data)
                            self.edmonds_graph.add_edge(node_i, node_j, **new_data)
                            edmonds_graph_removed_edges.remove((current_node, current_backward))

                        resolved_nodes.add(current_node)

                    else:
                        if current_backward_connected:

                            number = max([(self.pangenome_graph.get_edge_data(current_node, backward)['weight'], i) for (i, backward) in enumerate(current_backward_connected)])[1]

                            node_i, node_j, data = self.get_edge(current_node, current_backward_connected[number], reverse = True)
                            self.edmonds_graph.remove_edge(edmonds_graph_predecessors[current_node], current_node)

                            new_data = self.edge_check(node_i, node_j, data)
                            self.edmonds_graph.add_edge(node_i, node_j, **new_data)

                            edmonds_graph_removed_edges.remove((current_node, current_backward_connected[number]))
                            edmonds_graph_removed_edges.add((edmonds_graph_predecessors[current_node], current_node))

                            edmonds_graph_successors[edmonds_graph_predecessors[current_node]].remove(current_node)
                            edmonds_graph_successors[current_backward_connected[number]] += [current_node]

                            pred = edmonds_graph_predecessors[current_node]

                            edmonds_graph_predecessors.pop(current_node, None)
                            edmonds_graph_predecessors[current_node] = current_backward_connected[number]

                            edmonds_graph_distances[current_node] = self.mean_edmonds_graph_path_weight('start', current_node)

                            resolved_nodes.add(current_node)

                visited_nodes.add(current_node)

                if not nx.is_directed_acyclic_graph(self.edmonds_graph):
                    raise ConfigError(f"Oh no. It looks like your graph is so complex or includes a motif I haven't seen before"
                                      f"therefore the reattachement algorithm itself included a loop to the graph. We had multiple"
                                      f"sanity checks to prevent this but unfortunatly nobody is perfect. We will include more"
                                      f"checks in the next version. Sorry :/")

        self.progress.end()

        remaining_stops = [node for node in self.edmonds_graph.nodes() if self.edmonds_graph.out_degree(node) == 0 and node != 'stop']
        self.run.info_single(f"{pp(i)} iterations to solve the graph")

        for stop in remaining_stops:

            pangenome_graph_edge_data = {
                'genome':{genome: {'gene_call': -1} for genome in self.genomes},
                'weight':len(self.genomes),
                'bended': [],
                'direction': 'R'
            }

            new_data = self.edge_check(stop, 'stop', pangenome_graph_edge_data)
            self.edmonds_graph.add_edge(stop, 'stop', **new_data)

            edmonds_graph_successors[stop] += ['stop']

        if not nx.is_directed_acyclic_graph(self.edmonds_graph):
            raise ConfigError(f"Oh no. It looks like your graph is so complex or includes a motif I haven't seen before"
                              f"therefore the reattachement algorithm itself included a loop to the graph. We had multiple"
                              f"sanity checks to prevent this but unfortunatly nobody is perfect. We will include more"
                              f"checks in the next version. Sorry :/")

        self.run.info_single("Done")


    def load_graph_from_json_file(self):

        self.run.warning(None, header="Import graph F from json file", lc="green")

        filesnpaths.is_file_json_formatted(self.pan_graph_json)

        self.jsondata = json.load(open(self.pan_graph_json))

        self.global_x = int(self.jsondata["infos"]["meta"]["global_x"])
        self.priority_genome = self.jsondata["infos"]["priority_genome"]

        self.global_y = 0
        self.ancest = nx.DiGraph()
        self.x_list = {}
        self.position = {}
        self.edges = {}

        for x in range(0, self.global_x+1):
            self.x_list[x] = []

        for node in self.jsondata["elements"]["nodes"]:
            name = self.jsondata["elements"]["nodes"][node]["name"]
            weight = self.jsondata["elements"]["nodes"][node]["weight"]
            genome = self.jsondata["elements"]["nodes"][node]["genome"]
            x = self.jsondata["elements"]["nodes"][node]["position"]["x"]

            self.x_list[x].append(node)
            self.position[node] = (x, -1)
            self.ancest.add_node(node, name=name, pos=(0,0), weight=weight, layer={}, genome=genome)

        for edge in self.jsondata["elements"]["edges"]:
            source = self.jsondata["elements"]["edges"][edge]["source"]
            target = self.jsondata["elements"]["edges"][edge]["target"]
            genome = self.jsondata["elements"]["edges"][edge]["genome"]
            direction = self.jsondata["elements"]["edges"][edge]["direction"]

            self.ancest.add_edge(source, target, weight=-1, genome=genome, bended=[], direction=direction)

        self.run.info_single("Done")


    def prepare_synteny_graph(self):

        for x, generation in enumerate(nx.topological_generations(self.edmonds_graph)):
            self.x_list[x] = generation
            for node in generation:
                self.position[node] = (x, -1)

        self.global_x = x

        self.ancest = nx.DiGraph(self.edmonds_graph)
        nx.set_edge_attributes(self.ancest, values=-1, name='weight')


    def run_synteny_layout_algorithm(self):

        self.run.warning(None, header="Calculating coordinates on the nodes from F", lc="green")

        layout_graph_nodes = list(self.ancest.nodes())
        layout_graph_successors = {layout_graph_node: list(self.ancest.successors(layout_graph_node)) for layout_graph_node in layout_graph_nodes}

        n_removed = 0
        ghost = 0
        for x in range(self.global_x-1, 0, -1):
            for node in self.x_list[x]:
                node_x_position = self.position[node][0]

                change = []
                for successor in layout_graph_successors[node]:
                    if successor != 'stop':
                        successor_x_position = self.position[successor][0]

                        if successor_x_position <= node_x_position:
                            raise ConfigError(f"The node {node} is succeded by the node {successor} which is in front of the node {node}"
                                              f"that does not make a lot of sense and should not happen. We are sorry for the inconvenience :/")
                        else:
                            change.append((successor_x_position, successor))

                if change:
                    if min(change)[0] > 1:
                        new_x_position = min([(x, n) for (x, n) in change if x > 1])[0] - 1
                        self.position[node] = (new_x_position, -1)

                    for (position, extend_successor) in change:
                        x_difference = position - new_x_position

                        if x_difference <= self.max_edge_length_filter or self.max_edge_length_filter == -1:

                            path_list = [node]

                            for i in range(1, x_difference):
                                path_list += ['Ghost_' + str(ghost)]
                                self.position['Ghost_' + str(ghost)] = (new_x_position + i, -1)
                                ghost += 1

                            path_list += [extend_successor]

                            self.edges[(node, extend_successor)] = (path_list, self.ancest.edges[node, extend_successor])
                            self.ancest.remove_edge(node, extend_successor)

                            if len(path_list) == 2:
                                self.ancest.add_edges_from(map(tuple, zip(path_list, path_list[1:])), weight=-1)
                            else:
                                self.ancest.add_edges_from(map(tuple, zip(path_list, path_list[1:])), weight=-0.5)

                        else:
                            n_removed += 1
                            self.removed.add((node, extend_successor))

        for i, j in self.ancest.edges():

            if self.position[j][0] - self.position[i][0] != 1 and i != 'start' and j != 'stop' and (i,j) not in self.removed:
                raise ConfigError(f"Hmmm. This situation would create a very weird looking connection."
                                  f"The ede {(i, j)} is longer than it should be. I don't know what created"
                                  f"this but we will work on a solution on the next release. Sorry :(")

        if self.max_edge_length_filter == -1:
            self.run.info_single("Setting algorithm to 'keep all edges'")
        else:
            self.run.info_single(f"Setting algorithm to 'max edge length < {pp(self.max_edge_length_filter)}'")

        self.run.info_single(f"Removed {pp(n_removed)} edge(s) due to length cutoff")

        longest_path = nx.bellman_ford_path(G=self.ancest, source='start', target='stop', weight='weight')
        m = set(longest_path)

        starts = [node for node in self.ancest.nodes() if self.ancest.in_degree(node) == 0 if node != 'start']
        for st in starts:
            self.ancest.add_edge('start', st, weight=-1)
        dfs_list = list(nx.dfs_edges(self.ancest, source='start'))

        group = 0
        groups = {}
        groups_rev = {}

        # TODO Currently no seperation between unequal genome context
        if self.gene_cluster_grouping_threshold == -1:
            self.run.info_single("Setting algorithm to 'no grouping'")
        else:
            self.run.info_single(f"Setting algorithm to 'Grouping single connected chains size > {pp(self.gene_cluster_grouping_threshold)}'")

        for node_v, node_w in dfs_list:
            if node_v != 'start' and self.ancest.in_degree(node_v) == 1 and self.ancest.out_degree(node_v) == 1 and self.ancest.in_degree(node_w) == 1 and self.ancest.out_degree(node_w) == 1:
                
                if self.ungroup[0] != -1 and self.ungroup[1] != -1 and self.position[node_v][0] > self.ungroup[0] and self.position[node_v][0] < self.ungroup[1] and self.position[node_w][0] > self.ungroup[0] and self.position[node_w][0] < self.ungroup[1]:
                    continue
                else:
                    if node_v not in groups_rev.keys():
                        group_name = 'GCG_' + str(group).zfill(8)
                        groups[group_name] = [node_v, node_w]
                        groups_rev[node_v] = group_name
                        groups_rev[node_w] = group_name
                        group += 1

                    else:
                        group_name = groups_rev[node_v]
                        groups[group_name] += [node_w]
                        groups_rev[node_w] = group_name

        for label, condense_nodes in groups.items():

            condense_nodes = [node for node in condense_nodes if not node.startswith('Ghost_')]

            if len(condense_nodes) >= self.gene_cluster_grouping_threshold and self.gene_cluster_grouping_threshold != -1:
                self.grouping[label] = condense_nodes

        self.run.info_single(f"Grouped {pp(len(sum(self.grouping.values(), [])))} nodes in {pp(len(self.grouping.keys()))} groups")

        for st in starts:
            self.ancest.remove_edge('start', st)

        self.ancest.remove_edges_from(self.removed)

        branches = {}
        sortable = []
        for g in groups.keys():
            branch = groups[g]

            if not set(branch).isdisjoint(m) and not set(branch).issubset(m):
                raise ConfigError(f"A group is neither disjoint from the main path nor subset of the main path"
                                  f"we should not continue from here as this is not something that should happen.")
                # print('Sanity Error. Code 10.')
                # break

            elif set(branch).isdisjoint(m):

                start = self.position[branch[0]][0]
                length = len(branch)

                if start in branches.keys():
                    if length in branches[start].keys():
                        if branches[start][length].keys():
                            num = max(branches[start][length].keys()) + 1
                        else:
                            num = 1

                        branches[start][length][num] = branch
                    else:
                        num = 1
                        branches[start][length] = {num: branch}
                else:
                    num = 1
                    branches[start] = {length: {num: branch}}

                sortable += [(start, length, num)]

        left_nodes = set(self.ancest.nodes()) - set(groups_rev.keys())
        for n in left_nodes:

            if not set([n]).isdisjoint(m) and not set([n]).issubset(m):
                raise ConfigError(f"A group is neither disjoint from the main path nor subset of the main path"
                                  f"we should not continue from here as this is not something that should happen.")
                # print('Sanity Error. Code 11.')
                # break

            elif set([n]).isdisjoint(m):
                start = self.position[n][0]
                length = 1
                if n != 'start' and n != 'stop':
                    if start in branches.keys():
                        if length in branches[start].keys():
                            if branches[start][length].keys():
                                num = max(branches[start][length].keys()) + 1
                            else:
                                num = 1

                            branches[start][length][num] = [n]
                        else:
                            num = 1
                            branches[start][length] = {num: [n]}
                    else:
                        num = 1
                        branches[start] = {length: {num: [n]}}

                    sortable += [(start, length, num)]

        used = set()
        finished = set()

        y_new = 0
        for node in longest_path:
            x_pos = self.position[node][0]
            self.position[node] = (x_pos, y_new)
            used.add((x_pos, y_new))
            finished.add(node)

        stack = [longest_path]
        while stack:

            current = stack[0]

            remove = True
            for i,j,k in sorted(sortable, key=lambda x: (x[1], x[0]), reverse = False):
                branch = branches[i][j][k]
                branch_pred = set(self.ancest.predecessors(branch[0]))
                branch_succ = set(self.ancest.successors(branch[-1]))
                if (not branch_pred.isdisjoint(set(current))) or (not branch_succ.isdisjoint(set(current))) or (not branch_pred.isdisjoint(set(current)) and not branch_succ.isdisjoint(set(current))):

                    remove = False
                    sortable.remove((i,j,k))
                    y_new = max(sum([[self.position[ypred][1] for ypred in branch_pred], [self.position[ysucc][1] for ysucc in branch_succ]], []))

                    stack = [branch] + stack
                    while True:
                        repeat = False
                        for xnew in range(i, i+j):
                            if (xnew, y_new) in used:
                                repeat = True

                        if repeat == False:
                            break
                        else:
                            y_new += 1

                    for node in branch:
                        x_pos = self.position[node][0]
                        self.position[node] = (x_pos, y_new)
                        used.add((x_pos, y_new))
                        finished.add(node)

                        self.global_y = y_new if y_new > self.global_y else self.global_y

                    break

            if remove == True:
                stack.remove(current)

        if len(set(self.position.values())) != len(self.position.values()):
            # print(len(self.position.values()) - len(set(self.position.values())))
            raise ConfigError(f"No no no no. Something went very wrong here. Some nodes overlap in the UI."
                              f"We don't want this, we definitely don't want this...")

        for edge_i, edge_j in self.edges.keys():
            path_list, infos = self.edges[(edge_i, edge_j)]

            self.ancest.add_edge(edge_i, edge_j, **infos)
            self.ancest[edge_i][edge_j]['weight'] = sum([1 if value != self.priority_genome else 10 for value in self.ancest[edge_i][edge_j]['genome'].keys()])
            self.ancest[edge_i][edge_j]['bended'] = [self.position[p] for p in path_list[1:-1]]
            self.ancest.remove_nodes_from(path_list[1:-1])

        for node in self.ancest.nodes():
            self.ancest.nodes[node]['pos'] = self.position[node]

            self.offset[node] = self.position[node][0]

        x_positions_list = []

        for group_name in self.grouping.keys():

            group = self.grouping[group_name]

            group_size = len(group)
            group_size_compressed = round(group_size * self.groupcompress)
            compressed_factor = group_size_compressed if group_size_compressed == 0 else group_size_compressed - 1

            node_distance_factor = compressed_factor / (group_size - 1)

            start_x = self.ancest.nodes[group[0]]['pos'][0]

            for i, node in enumerate(group):
                self.offset[node] = round(start_x + i * node_distance_factor)

            compressed_length = int(self.offset[group[-1]] - self.offset[group[0]])

            x_positions_list += [i for i in range(start_x, start_x + compressed_length + 1)]

        x_positions_list += list(self.offset.values())
        x_positions = set(x_positions_list)
        empty_spots = []

        for x in range(self.global_x+1):
            if x not in x_positions:
                empty_spots.append(x)

        for node in self.ancest.nodes():
            decrease = 0
            x = self.offset[node]
            for e in empty_spots:
                if x > e:
                    decrease += 1

                    if e == empty_spots[-1]:
                        self.offset[node] = x - decrease
                        break

                else:
                    self.offset[node] = x - decrease
                    break

        for edge_i, edge_j, data in self.ancest.edges(data=True):
            if edge_i != 'start' and edge_j != 'stop':
                if data['bended']:
                    for i, (x, y) in enumerate(data['bended']):
                        decrease = 0
                        for e in empty_spots:
                            if x > e-1:
                                decrease += 1
                            else:
                                data['bended'][i] = (x - decrease, y)
                                break
                        else:
                            data['bended'][i] = (x - decrease, y)

                else:
                    if self.offset[edge_j] - self.offset[edge_i] != 1:

                        y = self.position[edge_i][1]
                        x = self.offset[edge_i] + 1

                        while x < self.offset[edge_j]:
                            data['bended'].append((x, y))
                            x += 1

        self.global_x_offset = self.offset['stop']

        self.run.info_single(f"Final graph {pp(len(self.ancest.nodes()))} nodes and {pp(len(self.ancest.edges()))} edges")
        
        self.run.info_single("Done")


    def generate_data_table(self):

        self.run.warning(None, header="Import layer values inherited from the pangenome", lc="green")

        max_paralogs = 0

        for gene_cluster in self.ancest.nodes():

            if gene_cluster != 'start' and gene_cluster != 'stop':
                node = self.ancest.nodes()[gene_cluster]
                num = len(node['genome'].keys())
                max_paralog_list = []
                num_dir_r = 0
                num_dir_l = 0

                func_homogeneity_list = []
                geo_homogeneity_list = []
                comb_homogeneity_list = []

                for genome in node['genome'].keys():
                    info = node['genome'][genome]

                    max_paralog_list.append(info['max_num_paralogs'])
                    func_homogeneity_list.append(info['functional_homogeneity_index'])
                    geo_homogeneity_list.append(info['geometric_homogeneity_index'])
                    comb_homogeneity_list.append(info['combined_homogeneity_index'])

                    if info['direction'] == 'r':
                        num_dir_r += 1
                    else:
                        num_dir_l += 1

                max_paralogs = max(max_paralog_list) if max(max_paralog_list) > max_paralogs else max_paralogs

                node['layer'] = {
                    'Paralogs': max(max_paralog_list) - 1,
                    'Direction': 1 - max(num_dir_r, num_dir_l)/num,
                    'Functional_Homogeneity': 1 - max(func_homogeneity_list),
                    'Geometric_Homogeneity': 1 - max(geo_homogeneity_list),
                    'Combined_Homogeneity': 1 - max(comb_homogeneity_list)
                }

        self.data_table_dict['Functional_Homogeneity'] = {
            'max': 1,
            'min': 0
        }
        self.data_table_dict['Geometric_Homogeneity'] = {
            'max': 1,
            'min': 0
        }
        self.data_table_dict['Combined_Homogeneity'] = {
            'max': 1,
            'min': 0
        }
        self.data_table_dict['Paralogs'] = {
            'max': max_paralogs - 1 if max_paralogs - 1 != 0 else 1,
            'min': 0
        }
        self.data_table_dict['Direction'] = {
            'max': 0.5,
            'min': 0
        }

        self.run.info_single("Done")


    def get_additional_gene_layer_table(self):

        self.run.warning(None, header="Export empty dataframe for genecall annotations", lc="green")

        gene_layer_dict = {}

        i = 0
        for node, data in self.ancest.nodes(data=True):
            if node != 'start' and node != 'stop':
                for genome in data['genome'].keys():
                    contig = data['genome'][genome]['contig']
                    genecall = data['genome'][genome]['gene_call']

                    gene_layer_dict[i] = {'genome': genome, 'contig': contig, 'genecall': genecall, 'value': 0}
                    i += 1

        gene_layer_df = pd.DataFrame.from_dict(gene_layer_dict, orient='index')

        gene_layer_df.to_csv(self.output_raw_gene_additional_data, index=False)

        self.run.info_single("Done")


    def get_additional_gc_layer_table(self):

        self.run.warning(None, header="Export empty dataframe for gc annotations", lc="green")

        gc_layer_dict = {}
        included = set()

        i = 0
        for node, data in self.ancest.nodes(data=True):
            if node != 'start' and node != 'stop':
                name = data['name']
                if name not in included:
                    gc_layer_dict[i] = {'genecluster': name, 'value': 0}
                    included.add(name)
                    i += 1

        gc_layer_df = pd.DataFrame.from_dict(gc_layer_dict, orient='index')

        gc_layer_df.to_csv(self.output_raw_gc_additional_data, index=False)

        self.run.info_single("Done")


    def add_additional_gene_layer_values(self):

        self.run.warning(None, header="Appending layer values from external gene data", lc="green")

        df = pd.read_csv(self.gene_additional_data)
        df.set_index(['genome', 'contig', 'genecall'], inplace=True)
        layer_names = list(df.columns)
        layer_max = {layer_name: 0 for layer_name in layer_names}

        for node, data in self.ancest.nodes(data=True):

            if node != 'start' and node != 'stop':

                for layer_name in layer_names:
                    value_list = []

                    for genome in data['genome'].keys():
                        contig = data['genome'][genome]['contig']
                        genecall = data['genome'][genome]['gene_call']

                        value = df.loc[(genome, contig, genecall)][layer_name].item()
                        value_list.append(value)

                    value_sum = sum(value_list) / len(value_list)

                    data['layer'][layer_name] = value_sum
                    layer_max[layer_name] = layer_max[layer_name] if layer_max[layer_name] > value_sum else value_sum

        for layer_name in layer_names:
            self.data_table_dict[layer_name] = {
                'max': layer_max[layer_name],
                'min': 0
            }

        self.run.info_single("Done")


    def add_additional_gc_layer_values(self):
        self.run.warning(None, header="Appending layer values from external gc data", lc="green")

        df = pd.read_csv(self.gc_additional_data)
        df.set_index(['genecluster'], inplace=True)
        layer_names = list(df.columns)
        layer_max = {layer_name: 0 for layer_name in layer_names}

        for node, data in self.ancest.nodes(data=True):
            if node != 'start' and node != 'stop':
                name = data['name']
                for layer_name in layer_names:

                    value = df.loc[name][layer_name].item()
                    data['layer'][layer_name] = value
                    layer_max[layer_name] = layer_max[layer_name] if layer_max[layer_name] > value else value

        for layer_name in layer_names:

            self.data_table_dict[layer_name] = {
                'max': layer_max[layer_name],
                'min': 0
            }

        self.run.info_single("Done")


    def walk_one_step(self, G, current, nodes_position_dict, visited):
        successors = [successor for successor in G.successors(current) if successor not in visited]
        predecessors = [predecessor for predecessor in G.predecessors(current) if predecessor not in visited]
        current_position_x, current_position_y = nodes_position_dict[current]
        
        if len(successors) > 0:
            nodes_position = [(nodes_position_dict[successor][0] - current_position_x, successor) for successor in successors]
            successor_position, successor = sorted(nodes_position)[0]
            successor_position_x, successor_position_y = nodes_position_dict[successor] 
       
            return(successor_position_x, successor_position_y, successor)

        elif len(predecessors) > 0:
            nodes_position = [(current_position_x - nodes_position_dict[predecessor][0], predecessor) for predecessor in predecessors]
            predecessor_position, predecessor = sorted(nodes_position, reverse=True)[0]            
            predecessor_position_x, predecessor_position_y = nodes_position_dict[predecessor]

            return(predecessor_position_x, predecessor_position_y, predecessor)

        else:
            return(-1, -1, '')


    def is_node_core(self, node):
        if set(self.ancest.nodes()[node]['genome'].keys()) == set(self.genomes_names):
            return(True)
        else:
            return(False)


    def get_genomic_motifs(self):

        self.run.warning(None, header="Generate pangenome graph summary tables", lc="green")

        graph_min = self.ancest.nodes()['start']['pos'][0] + 1
        graph_max = self.ancest.nodes()['stop']['pos'][0] - 1

        region_id = 0
        regions = {}
        region_sides = {}

        for genome in self.genomes_names:
            # print(genome)

            G_subset_nodes = [n for n,d in self.ancest.nodes(data='genome') if genome not in d.keys()]
            G_subset_edges = [(i,j) for i,j,d in self.ancest.edges(data='genome') if genome not in d.keys()]
            G3 = self.ancest.copy()
            G3.remove_nodes_from(['start', 'stop'] + G_subset_nodes)
            G3.remove_edges_from(G_subset_edges)

            region_type = {}
            regions_reverse = []
            nodes_position_dict = dict(G3.nodes(data="pos"))

            starting_points = []
            for node in list(G3.nodes()):
                if len(list(G3.successors(node))) <= 1 and len(list(G3.predecessors(node))) == 0:
                    starting_points += [node]

            starting_points = sorted(starting_points, key=lambda x: nodes_position_dict[x][0]) 

            # Prepare first starting point in data structure
            node = starting_points[0]
            starting_points.remove(node)

            node_position_x, node_position_y = nodes_position_dict[node]
            former_position = node_position_x

            region_id += 1

            # gene cluster core?
            if self.is_node_core(node):
                # start new core region
                regions[region_id] = [(node_position_x, node_position_y, node)]
                is_core = True
            else:
                # start new accessory region
                regions[region_id] = [(node_position_x, node_position_y, node)]
                is_core = False

            visited = set([node])
            mode = ''

            while True:
                node_position_x, node_position_y, node = self.walk_one_step(G3, node, nodes_position_dict, visited)

                if node:
                    mode = ''
                else:
                    if starting_points:
                        node = starting_points[0]
                        starting_points.remove(node)
                        node_position_x, node_position_y = nodes_position_dict[node]
                        direction = ''
                        mode = 'break'
                    else:
                        break
                
                # extend current region?
                if node_position_x == former_position + 1 or node_position_x == (former_position - graph_max) + 1:
                    # extend core region?
                    if is_core:
                        # gene cluster core?
                        if self.is_node_core(node):
                            # extend core region with core gene cluster
                            regions[region_id] += [(node_position_x, node_position_y, node)]
                            is_core = True
                        else:
                            # end core region
                            region_type[region_id] = ('core', regions[region_id][0][0], regions[region_id][-1][0])
                            region_id += 1

                            # start accessory region
                            regions[region_id] = [(node_position_x, node_position_y, node)]
                            is_core = False
                    else:
                        # gene cluster core?
                        if self.is_node_core(node):
                            # end accessory region
                            region_type[region_id] = ('accessory', regions[region_id][0][0], regions[region_id][-1][0])
                            region_id += 1

                            # start core region
                            regions[region_id] = [(node_position_x, node_position_y, node)]
                            is_core = True
                        else:
                            # extend accessory region
                            regions[region_id] += [(node_position_x, node_position_y, node)]
                            is_core = False

                elif node_position_x == former_position:
                    if is_core:
                        region_type[region_id] = ('core', regions[region_id][0][0], regions[region_id][-1][0])
                    else:
                        region_type[region_id] = ('accessory', regions[region_id][0][0], regions[region_id][-1][0])

                    region_id += 1

                    if self.is_node_core(node):
                        # start new core region
                        regions[region_id] = [(node_position_x, node_position_y, node)]
                        is_core = True
                    else:
                        # start new accessory region
                        regions[region_id] = [(node_position_x, node_position_y, node)]
                        is_core = False

                # new forward bridge region? 
                elif (node_position_x > former_position + 1 and abs(node_position_x - former_position) <= graph_max / 2) or (node_position_x < former_position + 1 and abs(former_position - node_position_x) > graph_max / 2):
                    # is the region to end core?
                    if is_core:  
                        region_type[region_id] = ('core', regions[region_id][0][0], regions[region_id][-1][0])
                    else:
                        region_type[region_id] = ('accessory', regions[region_id][0][0], regions[region_id][-1][0])

                    region_id += 1

                    if node_position_x > former_position + 1:
                        regions[region_id] = [(i, '', '') for i in range(former_position + 1, node_position_x, 1)]
                    else:
                        regions[region_id] = [(i, '', '') for i in range(former_position + 1, graph_max + 1, 1)] + [(i, '', '') for i in range(graph_min, node_position_x, 1)]

                    if mode == 'break':
                        region_type[region_id] = ('break', regions[region_id][0][0], regions[region_id][-1][0])
                    else:
                        region_type[region_id] = ('forward', regions[region_id][0][0], regions[region_id][-1][0])

                    region_id += 1

                    # gene cluster core?
                    if self.is_node_core(node):
                        # start new core region
                        regions[region_id] = [(node_position_x, node_position_y, node)]
                        is_core = True
                    else:
                        # start new accessory region
                        regions[region_id] = [(node_position_x, node_position_y, node)]
                        is_core = False

                # new backward bridge region?
                elif (node_position_x < former_position + 1 and abs(former_position - node_position_x) <= graph_max / 2) or (node_position_x > former_position + 1 and abs(node_position_x - former_position) > graph_max / 2):
                    # is the region to end core?
                    if is_core:  
                        region_type[region_id] = ('core', regions[region_id][0][0], regions[region_id][-1][0])
                    else:
                        region_type[region_id] = ('accessory', regions[region_id][0][0], regions[region_id][-1][0])

                    region_id += 1

                    if node_position_x < former_position + 1:
                        regions[region_id] = [(i, '', '') for i in range(former_position - 1, node_position_x, -1)]
                    else:
                        regions[region_id] = [(i, '', '') for i in range(former_position - 1, graph_min - 1, -1)] + [(i, '', '', '', genome) for i in range(graph_max, node_position_x, -1)]

                    region_type[region_id] = ('reverse', regions[region_id][0][0], regions[region_id][-1][0])

                    # mark a backward bridge region
                    regions_reverse += [(regions[region_id][0][0], regions[region_id][-1][0])]
                    region_id += 1
                    
                    # gene cluster core?
                    if self.is_node_core(node):
                        # start new core region
                        regions[region_id] = [(node_position_x, node_position_y, node)]
                        is_core = True
                    else:
                        # start new accessory region
                        regions[region_id] = [(node_position_x, node_position_y, node)]
                        is_core = False

                else:
                    print('A motif of this type is not yet implemented as we never saw it before.')

                visited.add(node)
                former_position = node_position_x

            if is_core:  
                region_type[region_id] = ('core', regions[region_id][0][0], regions[region_id][-1][0])
            else:
                region_type[region_id] = ('accessory', regions[region_id][0][0], regions[region_id][-1][0])       

            for region_id_inner, region_tuple in region_type.items():
                # overlap = False
                # for region_reverse in regions_reverse:

                #     if region_tuple[1] == region_reverse[0] and region_tuple[2] == region_reverse[1]:
                #         overlap = True
                    
                #     elif region_tuple[1] <= region_tuple[2] and region_reverse[0] >= region_reverse[1]:
                #         # reverse: 15 -> 12
                #         # current: 13 -> 16
                        
                #         #    ------------
                #         # ------------
                #         # 12 13 14 15 16
            
                #         # left_difference: abs(12 - 13) = 1
                #         # right_difference: abs(15 - 16) = 1
                #         # region_compare_size: (16 - 12) + 1 = 5
                #         # region_overlap_size = 5 - (1 +  1) = 3
            
                #         # region_reverse_size: (15 - 12) + 1 = 4
                #         # region_current_size: (16 - 13) + 1 = 4
            
                #         # test: 3 >= (4 * (1/2)) = True
            
                #         left_difference = abs(min(region_tuple[1], region_tuple[2]) - min(region_reverse[0], region_reverse[1]))
                #         right_difference = abs(max(region_tuple[1], region_tuple[2]) - max(region_reverse[0], region_reverse[1]))
                #         region_compare_size = (max(region_tuple[1], region_tuple[2], region_reverse[0], region_reverse[1]) - min(region_tuple[1], region_tuple[2], region_reverse[0], region_reverse[1])) + 1
                #         region_overlap_size = region_compare_size - (left_difference +  right_difference)
            
                #         region_reverse_size = (max(region_reverse[0], region_reverse[1]) - min(region_reverse[0], region_reverse[1])) + 1
                #         region_current_size = (max(region_tuple[1], region_tuple[2]) - min(region_tuple[1], region_tuple[2])) + 1
            
                #         if region_overlap_size >= (region_current_size / 2):
                #             overlap = True

                #     else:
                #         print('A rearrangement event over the start and end border of the graph is not yet implemented')
                
                # if overlap:
                #     region_sides[region_id_inner] = {'start': region_tuple[1], 'end': region_tuple[2], 'types': region_tuple[0], 'consensus': 1.0, 'motif': 'reassortment'}
                # else:
                region_sides[region_id_inner] = {'start': region_tuple[1], 'end': region_tuple[2], 'types': region_tuple[0], 'genome': genome, 'consensus': -1, 'motif': ''}
        
        tab = pd.DataFrame.from_dict(region_sides, orient='index')
        tab.rename_axis("region_id", inplace=True)
        tab.sort_values(['start', 'end'], inplace=True)
        
        drop_rows = []
        groups = tab.groupby(['genome', 'types'])
        # come up with grouping like 0-10, 0-6, 7-10 end up in one group 
        # easy search for all connections from i = 0 -> 5, 8, 10 take the highest
        # i = 0 ... 10 then pick all fragments inside of 0 and 10 e.g. 0 ... 5, 4 ... 8
        # then check every position for at most every genome once before combine
        for name, group in groups:

            i = 0
            j = i + 1
            while i < len(group) - 1:

                row_i = group.iloc[i]
                start_i, end_i, types_i, genome_i, consensus_i, motif_i = row_i
                index_i = row_i.name
        
                row_j = group.iloc[j]
                start_j, end_j, types_j, genome_j, consensus_j, motif_j = row_j
                index_j = row_j.name
                
                if (end_i + 1) == start_j:
                    end_new = max(end_i, end_j)

                    tab.at[index_i, 'end'] = end_new
                    drop_rows += [index_j]
                    
                    regions[index_i] += regions[index_j]
                    regions.pop(index_j)
                    tab.drop(index=index_j, inplace=True)

                    j += 1

                else:
                    i = j + 1
                    j = i + 1

        # return(regions, all_nodes)
        tab.reset_index(drop=False, inplace=True)
        tab['motif_id'] = tab.loc[:, 'region_id']

        # for key in regions.keys():
        #     regions[key] = set(regions[key])
        
        groups = tab.groupby(['start', 'end'])
        for name, group in groups:
        
            types_list = list(group['types'])
        
            if len(set(types_list)) == 1 and len(types_list) <= len(self.genomes_names):
        
                joining_list = sorted([(i, j, True) if regions[group.at[i,'motif_id']] == regions[group.at[i,'motif_id']] else (i, j, False) for i, j in it.combinations(group.index, 2)])
                
                for index_i, index_j, join_boolean in joining_list:
                    if join_boolean:
                        region_id_i, start_i, end_i, types_i, consensus_i, genome_i, motif_i, motif_id_i = tab.loc[index_i]
                        region_id_j, start_j, end_j, types_j, consensus_j, genome_j, motif_j, motif_id_j = tab.loc[index_j]
                    
                        tab.at[index_j, 'motif_id'] = motif_id_i
                        tab.at[index_j, 'motif_id'] = motif_id_i

                        tab.at[index_i, 'consensus'] = 1.0
                        tab.at[index_j, 'consensus'] = 1.0
                        
                        tab.at[index_i, 'motif'] = 'hypervariability' if types_i == 'accessory' and len(types_list) == len(self.genomes_names) else types_i
                        tab.at[index_j, 'motif'] = 'hypervariability' if types_i == 'accessory' and len(types_list) == len(self.genomes_names) else types_i

                        # regions[region_id_i] += regions[region_id_j]
                        
            elif set(types_list) == set(['forward', 'accessory']) and len(types_list) <= len(self.genomes_names):
        
                forward_group = group[group["types"] == 'forward']
                accessory_group = group[group["types"] == 'accessory']
                
                joining_equal = all([True if regions[accessory_group.at[i,'motif_id']] == regions[accessory_group.at[i,'motif_id']] else False for i, j in it.combinations(accessory_group.index, 2)])
                joining_list = sorted([(i, j, True) if regions[group.at[i,'motif_id']] == regions[group.at[i,'motif_id']] else (i, j, False) for i, j in it.combinations(group.index, 2)])
                
                for index_i, index_j, join_boolean in joining_list:
                    if join_boolean:
                        region_id_i, start_i, end_i, types_i, consensus_i, genome_i, motif_i, motif_id_i = tab.loc[index_i]
                        region_id_j, start_j, end_j, types_j, consensus_j, genome_j, motif_j, motif_id_j = tab.loc[index_j]
                        
                        tab.at[index_j, 'motif_id'] = motif_id_i
                        tab.at[index_j, 'motif_id'] = motif_id_i

                        if joining_equal:
                            a = len(forward_group)
                            b = len(accessory_group)
                            if a >= b:
                                consensus = round(a / (a+b), 3)
                                motif = 'insertion'
                            else:
                                consensus = round(b / (a+b), 3)
                                motif = 'deletion'
                        else:
                            consensus = 1.0
                            motif = 'hypervariability'

                        tab.at[index_i, 'consensus'] = consensus
                        tab.at[index_j, 'consensus'] = consensus

                        tab.at[index_i, 'motif'] = motif
                        tab.at[index_j, 'motif'] = motif

                        # regions[region_id_i] += regions[region_id_j]
                
            elif set(types_list) == set(['break', 'accessory']) and len(types_list) <= len(self.genomes_names):
                
                joining_equal = all([True if regions[accessory_group.at[i,'motif_id']] == regions[accessory_group.at[i,'motif_id']] else False for i, j in it.combinations(accessory_group.index, 2)])
                joining_list = sorted([(i, j, True) if regions[group.at[i,'motif_id']] == regions[group.at[i,'motif_id']] else (i, j, False) for i, j in it.combinations(group.index, 2)])
                
                for index_i, index_j, join_boolean in joining_list:
                    if join_boolean:
                        region_id_i, start_i, end_i, types_i, consensus_i, genome_i, motif_i, motif_id_i = tab.loc[index_i]
                        region_id_j, start_j, end_j, types_j, consensus_j, genome_j, motif_j, motif_id_j = tab.loc[index_j]
                    
                        tab.at[index_j, 'motif_id'] = motif_id_i
                        tab.at[index_j, 'motif_id'] = motif_id_i

                        if joining_equal:
                            consensus = 1.0
                            motif = 'broken core'
                        else:
                            consensus = 1.0
                            motif = 'hypervariability'

                        tab.at[index_i, 'consensus'] = consensus
                        tab.at[index_j, 'consensus'] = consensus

                        tab.at[index_i, 'motif'] = motif
                        tab.at[index_j, 'motif'] = motif

                        # regions[region_id_i] += regions[region_id_j]

            else:
                continue

        tab.sort_values(['start', 'end'], inplace=True)
        tab.reset_index(drop=True, inplace=True)

        groups = tab.groupby(['motif_id'])
        nodes_dict_id = 0
        nodes_dict = {}
        for motif_id, group in groups:

            for index, row in group.iterrows():

                region_id = row['region_id']
                genome = row['genome']

                for region in regions[region_id]:
                    node_position_x, node_position_y, node = region

                    if node:
                        nodes_dict[nodes_dict_id] = {'motif_id': motif_id, 'node_position_x': node_position_x, 'node_position_y': node_position_y, 'gene_cluster_id': node, 'genome': genome}
                        nodes_dict_id += 1

        tab2 = pd.DataFrame.from_dict(nodes_dict, orient='index').set_index(['gene_cluster_id', 'genome'])
        tab1 = tab.drop(['region_id'], axis=1).set_index(['motif_id', 'genome'])
        
        nodes_dict = {n:d for n,d in self.ancest.nodes(data='genome')}
        tab3 = pd.DataFrame.from_dict(
            {(i,j): nodes_dict[i][j]
            for i in nodes_dict.keys()
            for j in nodes_dict[i].keys()
            if i != 'start' and i != 'stop'}, orient='index')
        
        tab3.drop('gene_cluster_id', inplace=True, axis=1)
        tab3.rename_axis(['gene_cluster_id', 'genome'], inplace=True)

        for source in self.functional_annotation_sources_available:
            tab3[source] = tab3[source].apply(lambda x: ('None', 'None', 'None') if x == 'None' else x)
            tab3[[source + '_ID', source + 'TEXT', source + '_E_VALUE']] = pd.DataFrame(tab3[source].values.tolist(), index=tab3.index)
            tab3.drop(source, inplace=True, axis=1)

        tab3.fillna('None', inplace=True)
        tab4 = pd.merge(tab2, tab3, how='left', on=['gene_cluster_id', 'genome'])
        
        if len(tab2) != len(tab3):
            self.run.info_single('Genecall entries are missing from the result table!')
            
        tab1.to_csv(self.output_summary + '/motifs.tsv', sep='\t')
        tab4.to_csv(self.output_summary + '/genes.tsv', sep='\t')

        self.run.info_single("Done")


    def update_json_dict(self):
        self.run.warning(None, header="Update JSON dictionary for export", lc="green")

        self.jsondata["infos"]['meta']['global_y'] = self.global_y
        self.jsondata["infos"]['meta']['global_x_offset'] = self.global_x_offset
        self.jsondata["infos"]['max_edge_length_filter'] = self.max_edge_length_filter
        self.jsondata["infos"]['gene_cluster_grouping_threshold'] = self.gene_cluster_grouping_threshold
        self.jsondata["infos"]['groupcompress'] = self.groupcompress
        self.jsondata["infos"]['ungroup'] = self.ungroup
        self.jsondata["infos"]['layout_graph']['edges'] = len(self.ancest.edges())
        self.jsondata["infos"]['data'] = list(self.ancest.graph.items())
        self.jsondata["infos"]['directed'] = self.ancest.is_directed()
        self.jsondata["infos"]['groups'] = self.grouping
        self.jsondata['infos']['grouped'] = {'nodes': len(self.ancest.nodes()), 'edges': len(self.ancest.edges())}

        for i, j in self.ancest.nodes(data=True):
            self.jsondata["elements"]["nodes"][str(i)]["position"]['y'] = j['pos'][1]
            self.jsondata["elements"]["nodes"][str(i)]["position"]['x_offset'] = self.offset[i]

        for edge in self.jsondata["elements"]["edges"]:
            source = self.jsondata["elements"]["edges"][edge]["source"]
            target = self.jsondata["elements"]["edges"][edge]["target"]
            self.jsondata["elements"]["edges"][edge]["shown"] = 1

            if (source, target) in self.removed:
                self.jsondata["elements"]["edges"][edge]["shown"] = 0
                self.jsondata["elements"]["edges"][edge]["bended"] = ""

            elif self.ancest.has_edge(source, target):
                self.jsondata["elements"]["edges"][edge]["bended"] = [{'x': x, 'y': y} for x, y in self.ancest[source][target]["bended"]] if len(self.ancest[source][target]["bended"]) != 0 else ""

        self.run.info_single("Done")


    # TODO rework that section for better debugging and add more features as an example fuse start and top
    # together or remove both so the graph becomes a REAL circle. Aside from that there is a bug in the remove
    # edges section for (k,o) in circular edges and for (k,o) in pangenome edges. Change and test while reworking.
    def get_json_dict_for_graph(self):
        self.run.warning(None, header="Calculating syntenous distances", lc="green")

        # NOTE: Any change in `self.jsondata` will require the pangraph JSON in anvio.tables.__init__
        #       to incrase by one (so that the new code that works with the new structure requires
        #       previously generated JSON to be recomputed).

        nodes_all = len(self.ancest.nodes())
        edges_all = len(self.ancest.edges())

        X = np.zeros([len(self.genomes_names), len(self.genomes_names)])
        for genome_i, genome_j in it.combinations(self.genomes_names, 2):
            nodes_similar = 0
            edges_similar = 0
            nodes_unsimilar = 0
            edges_unsimilar = 0
            for _, data in self.ancest.nodes(data=True):
                if genome_i in data['genome'].keys() and genome_j in data['genome'].keys():
                    nodes_similar += 1
                elif genome_i in data['genome'].keys() or genome_j in data['genome'].keys():
                    nodes_unsimilar += 1
            for _, _, data in self.ancest.edges(data=True):
                if genome_i in data['genome'].keys() and genome_j in data['genome'].keys():
                    edges_similar += 1
                elif genome_i in data['genome'].keys() or genome_j in data['genome'].keys():
                    edges_unsimilar += 1

            i = self.genomes_names.index(genome_i)
            j = self.genomes_names.index(genome_j)

            elements_similar = nodes_similar + edges_similar
            elements_unsimilar = nodes_unsimilar + edges_unsimilar
            elements_all = nodes_all + edges_all

            X[i][j] = elements_unsimilar / (elements_similar + elements_unsimilar)
            X[j][i] = elements_unsimilar / (elements_similar + elements_unsimilar)

            self.run.info_single(f"{genome_i} and {genome_j}: {round(X[i][j], 3)}")

        condensed_X = squareform(X)
        Z = linkage(condensed_X, 'ward')
        tree = to_tree(Z, False)
        newick = clustering.get_newick(tree, tree.dist, self.genomes_names)

        if self.output_graphics:
            fig = plt.figure(figsize=(25, 10))
            ax = plt.axes()
            dn = dendrogram(Z, ax=ax, labels=self.genomes_names, orientation='right')
            plt.tight_layout()
            fig.savefig(self.output_graphics + 'phylogenie' + '.pdf')

        self.run.info_single("Done.")
        self.run.warning(None, header="Creating JSON dictionary for export", lc="green")

        self.jsondata["infos"] = {
            'meta': {
                'title': self.project_name,
                'version': anvio.__pangraph__version__,
                'global_x': self.global_x,
                'global_y': self.global_y,
                'global_x_offset': self.global_x_offset
            },
            'genomes': {genome: 'on' for genome in self.genomes},
            'functional_annotation_sources_available': self.functional_annotation_sources_available,
            'num_genomes': len(self.genomes),
            'max_edge_length_filter': self.max_edge_length_filter,
            'gene_cluster_grouping_threshold': self.gene_cluster_grouping_threshold,
            'groupcompress': self.groupcompress,
            'ungroup': self.ungroup,
            'priority_genome': self.priority_genome,
            'newick': newick,
            'layout_graph': {
                'nodes': len(self.ancest.nodes()),
                'edges': len(self.ancest.edges())
            },
            'data': list(self.ancest.graph.items()),
            'directed': self.ancest.is_directed(),
            'groups': self.grouping,
            'grouped': {
                'nodes': len(self.ancest.nodes()),
                'edges': len(self.ancest.edges())
            },
            'layers_data': self.data_table_dict
        }

        self.jsondata['elements'] = {
            'nodes': {},
            'edges': {}
        }

        for i, j in self.ancest.nodes(data=True):

            self.jsondata["elements"]["nodes"][str(i)] = {
                "name": j["name"],
                "weight": j["weight"],
                "layer": j["layer"],
                "position": {
                    'x': j['pos'][0],
                    'y': j['pos'][1],
                    'x_offset': self.offset[i]
                },
                "genome": {key: {item: j["genome"][key][item] for item in ['gene_call', 'contig', 'direction']} for key in j["genome"].keys()}
            }

        for l, (k, o, m) in enumerate(self.ancest.edges(data=True)):

            self.jsondata["elements"]["edges"]['E_' + str(l).zfill(8)] = {
                "source": k,
                "target": o,
                "shown": 1,
                "genome": m["genome"],
                "weight": m["weight"],
                "direction": m["direction"],
                "bended": [{'x': x, 'y': y} for x, y in m["bended"]] if len(m["bended"]) != 0 else ""
            }

        self.run.info_single("Done.")


    def store_network(self):
        """Function to store final graph structure in a pan-db and/or JSON flat text output file"""

        self.run.warning(None, header="Exporting pangenome graph to JSON", lc="green")

        if self.output_graphml:

            graphml = nx.DiGraph()
            graphml.add_nodes_from([n for n, d in self.ancest.nodes(data=True)])
            graphml.add_edges_from([(i, j) for i, j, d in self.ancest.edges(data=True)])

            nx.write_graphml_lxml(graphml, self.output_graphml)
            self.run.info("GraphML output file", os.path.abspath(self.output_graphml))

        if self.json_output_file_path:

            filesnpaths.is_output_file_writable(self.json_output_file_path)

            with open(self.json_output_file_path, 'w') as output:
                output.write(json.dumps(self.jsondata, indent=2))
            self.run.info("JSON output file", os.path.abspath(self.json_output_file_path))
        else:
            self.run.info("JSON output file", "Skipped (but OK)", mc='red')

        if not self.skip_storing_in_pan_db:
            raise ConfigError("The storage of graph data in pan-db is not yet implemented :/")
        else:
            self.run.info("Into the pan-db", "Skipped (but OK)", mc='red')

        self.run.info_single("Done")