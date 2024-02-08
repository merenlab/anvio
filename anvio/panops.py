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
import time

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

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
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

        miscdata.TableForItemAdditionalData(self.args, r=terminal.Run(verbose=False)).add(d, ['functional_homogeneity_index', 'geometric_homogeneity_index', 'combined_homogeneity_index'], skip_check_names=True)


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


    def process(self):
        # start by processing the additional params user may have passed for the blast step
        self.process_additional_params()

        # load genomes from genomes storage
        self.load_genomes()

        # check sanity
        self.sanity_check()

        # gen pan_db
        self.generate_pan_db()

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
        self.pan_db = A('pan_db')
        self.external_genomes_txt = A('external_genomes')
        self.genomes_storage_db = A('genomes_storage')
        self.max_edge_length_filter = A('max_edge_length_filter')
        self.gene_cluster_grouping_threshold = A('gene_cluster_grouping_threshold')
        self.debug = anvio.DEBUG

        # data storage related variables
        self.skip_storing_in_pan_db = True # FIXME: DB storage is not yet implemented
        self.json_output_file_path = A('output_file')

        # learn what gene annotation sources are present across all genomes
        self.functional_annotation_sources_available = DBInfo(self.genomes_storage_db, expecting='genomestorage').get_functional_annotation_sources()

        self.priority_genome = 'wPip'

        # this is the dictionary that wil keep all data that is going to be loaded
        # from anvi'o artifacts
        self.gene_synteny_data_dict = {}
        self.genome_coloring = {}

        self.initial_graph = nx.DiGraph()
        self.pangenome_graph = nx.DiGraph()
        self.edmonds_graph = nx.DiGraph()
        self.ancest = nx.DiGraph()

        self.fusion_events = 0
        self.removed_edges = 0

        self.leaf_path = []
        self.grouping = {}

        self.global_y = 0
        self.global_x = 1
        self.k = 0
        self.genome_gc_occurence = {}
        self.ghost = 0
        self.debug = False

        self.position = {}
        self.x_list = {}
        self.path = {}
        self.edges = []


    def sanity_check(self):
        if self.skip_storing_in_pan_db and not self.json_output_file_path:
            raise ConfigError("You are initializing the Pangraph class with `--skip-storing-in-pan-db` without an `--output-file` "
                              "parameter for the graph results to be stored. Please set an output file path so anvi'o has at least one "
                              "way to store results.")

        if self.json_output_file_path:
            filesnpaths.is_output_file_writable(self.json_output_file_path)

        if self.external_genomes_txt:
            filesnpaths.is_file_tab_delimited(self.external_genomes_txt)

            if not utils.is_all_columns_present_in_TAB_delim_file(["name","contigs_db_path"], self.external_genomes_txt, including_first_column=True):
                raise ConfigError("Your external genomes file does not seem to contain that anvi'o expects to find "
                                  "in an external genomes file :/")

        # make sure the pan-db and genomes-storage-db are compatible
        utils.is_pan_db_and_genomes_storage_db_compatible(self.pan_db, self.genomes_storage_db)


    def process(self):
        """Primary driver function for the class"""

        # sanity check EVERYTHING
        self.sanity_check()

        # populate self.gene_synteny_data_dict
        self.get_gene_synteny_data_dict()

        # contextualize paralogs
        # TODO Incorporate gene direction
        self.contextualize_paralogs()

        # build graph
        self.build_graph()

        # reconnect open leaves in the graph to generate
        # a flow network from left to right
        # self.run_tree_to_flow_network_algorithm()
        
        ##### START OF ALEX TESTING AREA
        self.new_run_tree_to_flow_network_algorithm()
        ##### END OF ALEX TESTING AREA

        # process edges and nodes to extract unique paths
        # from the nework
        # TODO Edge direction overlay
        # self.calculate_component_paths()

        # run Alex's layout algorithm
        # TODO Rework the algorithm
        # self.run_synteny_layout_algorithm()

        # condense gene clusters into groups
        # self.condense_gene_clusters_into_groups()

        # store network in the database
        self.store_network()


    def get_gene_synteny_data_dict(self):
        """A function to roduce a comprehensive data structure from anvi'o artifacts for
           downstream analyses.
        """
        self.run.warning(None, header="Loading data from database", lc="green")

        pan_db = dbops.PanSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

        pan_db.init_gene_clusters()
        pan_db.init_gene_clusters_functions_summary_dict()
        gene_cluster_dict = pan_db.gene_callers_id_to_gene_cluster

        external_genomes = pd.read_csv(self.external_genomes_txt, header=0, sep="\t", names=["name","contigs_db_path"])
        external_genomes.set_index("name", inplace=True)

        for genome, contigs_db_path in external_genomes.iterrows():

            if genome not in self.genome_coloring.keys():
                self.genome_coloring[genome] = "on"

            if self.genome_coloring[genome] != "off":
                args = argparse.Namespace(contigs_db=contigs_db_path.item())
                contigs_db = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))

                caller_id_cluster = gene_cluster_dict[genome]
                caller_id_cluster_df = pd.DataFrame.from_dict(caller_id_cluster, orient="index", columns=["gene_cluster_name"]).rename_axis("gene_caller_id").reset_index()
                caller_id_cluster_df["gene_cluster_id"] = ""
                caller_id_cluster_df["max_paralog"] = 0
                caller_id_cluster_df["draw"] = self.genome_coloring[genome]

                contigs_db.init_functions()
                gene_function_calls_df = pd.DataFrame.from_dict(contigs_db.gene_function_calls_dict, orient="index", columns=self.functional_annotation_sources_available).rename_axis("gene_caller_id").reset_index()

                all_gene_calls = caller_id_cluster_df['gene_caller_id'].values.tolist()
                genes_in_contigs_df = pd.DataFrame.from_dict(contigs_db.get_sequences_for_gene_callers_ids(all_gene_calls, include_aa_sequences=True, simple_headers=True)[1], orient="index", columns=["contig", "start", "stop", "direction", "partial", "call_type", "source", "version", "sequence", "length", "rev_compd", "aa_sequence", "header"]).rename_axis("gene_caller_id").reset_index()

                joined_contigs_df = caller_id_cluster_df.merge(genes_in_contigs_df, on="gene_caller_id", how="left").merge(gene_function_calls_df, on="gene_caller_id", how="left")
                joined_contigs_df.sort_values(["contig", "start", "stop"], axis=0, ascending=True, inplace=True)
                joined_contigs_df.set_index(["contig", "gene_caller_id"], inplace=True)

                self.gene_synteny_data_dict[genome] = joined_contigs_df.fillna("None").groupby(level=0).apply(lambda df: df.xs(df.name).to_dict("index")).to_dict()

            self.run.info_single(f"{contigs_db_path.item().split('/')[-1]}")

        self.run.info_single("Done")


    def contextualize_paralogs(self):
        """A function that resolves the graph context of paralogs based on gene synteny information across genomes"""
        self.run.warning(None, header="Select paralog context", lc="green")

        unresolved = True
        solved = set()

        genome_gc_order_values = {}
        genome_gc_order_placeholder = {}

        while unresolved:

            unresolved = False
            drop = set()

            for genome in self.gene_synteny_data_dict.keys():
                for contig in self.gene_synteny_data_dict[genome].keys():
                    genome_gc_order = [self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_name"] for gene_call in self.gene_synteny_data_dict[genome][contig].keys()]

                    if self.k == 0:
                        if genome not in genome_gc_order_values:
                            genome_gc_order_values[genome] = {contig: genome_gc_order}
                            genome_gc_order_placeholder[genome] = {contig: ''}
                        else:
                            items = list(genome_gc_order_values[genome].values())
                            times = items.count(genome_gc_order)
                            if times == 0:
                                genome_gc_order_values[genome][contig] = genome_gc_order
                                genome_gc_order_placeholder[genome][contig] = ''
                            else:
                                genome_gc_order_values[genome][contig] = genome_gc_order
                                genome_gc_order_placeholder[genome][contig] = chr(ord('a')+times-1)

                    place = genome_gc_order_placeholder[genome][contig]

                    for i in range(0, len(genome_gc_order)):

                        start = (i - self.k) if (i - self.k) >= 0 else 0
                        stop = (i + self.k + 1) if (i + self.k + 1) <= len(genome_gc_order) else len(genome_gc_order)
                        entry = genome_gc_order[start:stop]
                        # gene_call = genome_gc_order[i][1]
                        # name = genome_gc_order[i][0]

                        if len(entry) == 1 + (2 * self.k):
                            gc_k = tuple(entry)
                        elif start == 0 and stop == len(genome_gc_order):
                            gc_k = tuple([place] * (self.k - i) + entry + [place] * ((i + self.k + 1) - len(genome_gc_order)))
                        elif start == 0:
                            gc_k = tuple([place] * (self.k - i) + entry)
                        elif stop == len(genome_gc_order):
                            gc_k = tuple(entry + [place] * ((i + self.k + 1) - len(genome_gc_order)))
                        else:
                            print("Sanity Error.")
                            exit()

                        if len(gc_k) != 1 + (2 * self.k):
                            print("Sanity Error.")
                            exit()

                        gc = gc_k[int(len(gc_k) / 2)]
                        if gc not in solved:
                            
                            if gc_k not in self.genome_gc_occurence.keys():
                                self.genome_gc_occurence[gc_k] = {genome: 1}

                            elif gc_k in self.genome_gc_occurence.keys():
                                if genome not in self.genome_gc_occurence[gc_k].keys():
                                    self.genome_gc_occurence[gc_k][genome] = 1
                                else:
                                    self.genome_gc_occurence[gc_k][genome] += 1

                            else:
                                # if genome not in self.genome_gc_occurence[gc_k[::-1]].keys():
                                #     self.genome_gc_occurence[gc_k[::-1]][genome] = 1
                                # else:
                                #     self.genome_gc_occurence[gc_k[::-1]][genome] += 1
                                print('Sanity Error.')
                                exit()

            for gc, genome_gc_frequency in self.genome_gc_occurence.items():
                if max(genome_gc_frequency.values()) > 1:

                    unresolved = True
                    drop.add(gc[int(len(gc)/2)])

            self.run.info_single(f"Iteration #{str(self.k)}: {pp(len(self.genome_gc_occurence))} GCs containing {len(drop)} paralogs")
            
            if self.k == 0:
                self.paralog_dict = copy.deepcopy(self.genome_gc_occurence)

            if unresolved:
                keys = list(self.genome_gc_occurence.keys())
                for gc in keys:
                    if gc[int(len(gc)/2)] in drop:
                        self.genome_gc_occurence.pop(gc)

                solved = set([gc[int(len(gc)/2)] for gc in self.genome_gc_occurence.keys()])
                self.k += 1

        syn_calls = 0
        for genome in self.gene_synteny_data_dict.keys():

            for contig in self.gene_synteny_data_dict[genome].keys():
                genome_gc_order = [(self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_name"], gene_call) for gene_call in self.gene_synteny_data_dict[genome][contig].keys()]

                place = genome_gc_order_placeholder[genome][contig]

                for i in range(0, len(genome_gc_order)):
                    start = i-self.k if i-self.k >= 0 else 0
                    stop = i+self.k+1 if i+self.k+1 <= len(genome_gc_order) else len(genome_gc_order)
                    entry = [item[0] for item in genome_gc_order[start:stop]]
                    gene_call = genome_gc_order[i][1]
                    name = genome_gc_order[i][0]

                    if len(entry) == 1 + (2 * self.k):
                        gc_k = tuple(entry)
                    elif start == 0 and stop == len(genome_gc_order):
                        gc_k = tuple([place] * (self.k - i) + entry + [place] * ((i + self.k + 1) - len(genome_gc_order)))
                    elif start == 0:
                        gc_k = tuple([place] * (self.k - i) + entry)
                    elif stop == len(genome_gc_order):
                        gc_k = tuple(entry + [place] * ((i + self.k + 1) - len(genome_gc_order)))
                    # else:
                    #     print("Problem!")

                    for j in range(0, self.k+1):

                        gc_group = gc_k[int(len(gc_k)/2)-j:int(len(gc_k)/2)+j+1]

                        if gc_group in self.genome_gc_occurence.keys():

                            self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_id"] = ','.join(gc_group)
                            self.gene_synteny_data_dict[genome][contig][gene_call]["max_paralog"] = self.paralog_dict[tuple([name])][genome]
                            syn_calls += 1
                            break

                        # elif gc_group[::-1] in self.genome_gc_occurence.keys():

                        #     self.gene_synteny_data_dict[genome][contig][gene_call]["gene_cluster_id"] = ','.join(gc_group[::-1])
                        #     self.gene_synteny_data_dict[genome][contig][gene_call]["max_paralog"] = self.paralog_dict[tuple([name])][genome]
                        #     break

                        else:
                            pass

        num_calls = 0
        for _, value in self.genome_gc_occurence.items():            
            num_calls += sum(value.values()) 

        if num_calls != syn_calls:
            print("Sanity Error.")
            exit()

        self.run.info_single("Done")

    # ANCHOR Node adding
    def add_node_to_graph(self, gene_cluster, name, info):

        if not self.initial_graph.has_node(gene_cluster):
            self.initial_graph.add_node(
                gene_cluster,
                name=name,
                pos=(0, 0),
                weight=1,
                genome=info
            )

        else:
            self.initial_graph.nodes[gene_cluster]['weight'] += 1
            self.initial_graph.nodes[gene_cluster]['genome'].update(info)


    # ANCHOR Edge adding
    def add_edge_to_graph(self, gene_cluster_i, gene_cluster_j, info):

        if self.priority_genome in info.keys():
            weight_add = 10
        else:
            weight_add = 0

        draw = {genome: {y: info[genome][y] for y in info[genome].keys() if y == 'draw'} for genome in info.keys()}

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

    # TODO Should reverse genes also be connected in reverse?
    def build_graph(self):
        """FIXME"""

        self.run.warning(None, header="Building directed gene cluster graph G", lc="green")

        for genome in self.gene_synteny_data_dict.keys():

            for contig in self.gene_synteny_data_dict[genome].keys():
                gene_cluster_kmer = []
                for gene_call in self.gene_synteny_data_dict[genome][contig].keys():
                    gene_cluster_kmer.append((self.gene_synteny_data_dict[genome][contig][gene_call]['gene_cluster_id'], self.gene_synteny_data_dict[genome][contig][gene_call]['gene_cluster_name'], {genome: {'contig':contig, 'gene_call':gene_call, **self.gene_synteny_data_dict[genome][contig][gene_call]}}))

                if len(gene_cluster_kmer) > 1:
                    gene_cluster_pairs = map(tuple, zip(gene_cluster_kmer, gene_cluster_kmer[1:]))
                    first_pair = next(gene_cluster_pairs)

                    self.add_node_to_graph(first_pair[0][0], first_pair[0][1], first_pair[0][2])
                    self.add_node_to_graph(first_pair[1][0], first_pair[1][1], first_pair[1][2])
                    self.add_edge_to_graph(first_pair[0][0], first_pair[1][0], first_pair[1][2])

                    for gene_cluster_pair in gene_cluster_pairs:

                        self.add_node_to_graph(gene_cluster_pair[1][0], gene_cluster_pair[1][1], gene_cluster_pair[1][2])
                        self.add_edge_to_graph(gene_cluster_pair[0][0], gene_cluster_pair[1][0], gene_cluster_pair[1][2])

                else:
                    self.add_node_to_graph(gene_cluster_kmer[0][0], gene_cluster_kmer[0][1], gene_cluster_kmer[0][2])

        self.run.info_single(f"Adding {pp(len(self.initial_graph.nodes()))} nodes and {pp(len(self.initial_graph.edges()))} edges to G")
        connectivity = nx.is_connected(self.initial_graph.to_undirected())
        self.run.info_single(f"Connectivity is {connectivity}")

        if connectivity == False:
            self.pangenome_graph = nx.DiGraph(self.initial_graph.subgraph(max(nx.weakly_connected_components(self.initial_graph), key=len)))
            self.run.info_single(f"Keeping Subgraph with {pp(len(self.pangenome_graph.nodes()))} nodes and {pp(len(self.pangenome_graph.edges()))} edges")

            connectivity = nx.is_connected(self.pangenome_graph.to_undirected())
            self.run.info_single(f"Connectivity is {connectivity}")
        else:
            self.pangenome_graph = nx.DiGraph(self.initial_graph)

        self.run.info_single("Done")

        self.run.warning(None, header="Building maximum branching graph M of G", lc="green")

        add_start = []
        for node in self.pangenome_graph.nodes():
            if len(list(self.pangenome_graph.predecessors(node))) == 0:
                add_start.append(node)

        self.pangenome_graph.add_node(
            'start',
            name='start',
            pos=(0, 0),
            weight=len(self.genome_coloring.keys()),
            genome={genome: {'draw': 'on'} for genome in self.genome_coloring.keys()}
        )

        for u in add_start:

            weight = self.pangenome_graph.nodes[u]['weight']
            genomes = self.pangenome_graph.nodes[u]['genome'].keys()

            self.pangenome_graph.add_edge(
                *('start', u),
                genome={genome: {'draw': 'on'} for genome in genomes},
                weight=weight,
                bended=[],
                direction='R'
            )

        # ANCHOR Edmonds Algorithm
        selfloops = list(nx.selfloop_edges(self.pangenome_graph))
        self.run.info_single(f"Found and removed {pp(len(selfloops))} selfloop edge(s).")
        self.pangenome_graph.remove_edges_from(selfloops)

        self.edmonds_graph = nx.algorithms.tree.branchings.maximum_spanning_arborescence(self.pangenome_graph, attr="weight")
        nx.set_edge_attributes(self.edmonds_graph, {(i, j): d for i, j, d in self.pangenome_graph.edges(data=True) if (i, j) in self.edmonds_graph.edges()})
        nx.set_node_attributes(self.edmonds_graph, {k: d for k, d in self.pangenome_graph.nodes(data=True) if k in self.edmonds_graph.nodes()})

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

            # if edmonds_graph_successors[current_branch_successor]:
            #     for successor in edmonds_graph_successors[current_branch_successor]:
            #         leaves.update(self.get_leaves(successor, edmonds_graph_successors))
            # else:
            #     leaves.update(set([current_branch_successor]))
            
        return(leaves)  

    def new_run_tree_to_flow_network_algorithm(self):

        self.run.warning(None, header="Building flow network F from M and G", lc="green")

        edmonds_graph_edges = set(self.edmonds_graph.edges())
        edmonds_graph_nodes = set(self.edmonds_graph.nodes())

        pangenome_graph_edges = set(self.pangenome_graph.edges())
        pangenome_graph_nodes = set(self.pangenome_graph.nodes())

        edmonds_graph_removed_edges = pangenome_graph_edges - edmonds_graph_edges
        edmonds_graph_end = max([(self.mean_edmonds_graph_path_weight('start', node), node) for node in edmonds_graph_nodes if len(list(self.pangenome_graph.successors(node))) == 0])[1]

        for graph in [self.pangenome_graph, self.edmonds_graph]:

            graph.add_node(
                'stop',
                name='stop',
                pos=(0, 0),
                weight=len(self.genome_coloring.keys()),
                genome={genome: {'draw': 'on'} for genome in self.genome_coloring.keys()}
            )
            graph.add_edge(
                *(edmonds_graph_end, 'stop'),
                genome={genome: {'draw': 'on'} for genome in self.genome_coloring.keys()},
                weight=len(self.genome_coloring.keys()),
                bended=[],
                direction='R'
            )

        # edmonds_graph_clean = nx.DiGraph(self.edmonds_graph)

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

                # print(i, len(resolved_nodes), "/", len(pangenome_graph_nodes) + 1)
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

                # print(successor_branch_leaves)
                
                if not successor_branch_leaves:
                    current_node = current_branch_root
                else:
                    current_node = max([(edmonds_graph_distances[successor_branch_leaf], successor_branch_leaf) for successor_branch_leaf in successor_branch_leaves])[1]
            
                if current_node in resolved_nodes:
                    connected = True
                else:
                    connected = False
                
                if connected != True or x == 1:
                    # print(current_node)
                    for current_node_successor in pangenome_graph_successors[current_node]:
                        # print(current_node, current_node_successor, current_node_successor in nx.ancestors(self.edmonds_graph, current_node))
                        
                        if current_node_successor in resolved_nodes:

                            # self.edmonds_graphAYBE BOTH?
                            # if current_node_successor in nx.ancestors(self.edmonds_graph, current_node):
                            # if current_node_successor not in visited_nodes:
                            if current_node_successor in nx.ancestors(self.edmonds_graph, current_node) or current_node_successor not in visited_nodes:
                                # print('possible reverse')
                                if (current_node, current_node_successor) in edmonds_graph_removed_edges:
                                    current_backward_connected.append(current_node_successor)
                                    # print('reverse')
                
                            else:
                                # print('possible forward')
                                if (current_node, current_node_successor) in edmonds_graph_removed_edges:
                                    current_forward_connected.append(current_node_successor)
                                    # print('forward')
                                    connected = True
                
                                else:
                                    # print('exist')
                                    connected = True
                
                    if connected == False:
                        # print(current_node, list(self.pangenome_graph.successors(current_node)))
                        
                        if len(list(self.pangenome_graph.successors(current_node))) == 0:
                            pangenome_graph_edge_data = {
                                'genome':{genome: {'draw': 'on'} for genome in self.genome_coloring.keys()},
                                'weight':len(self.genome_coloring.keys()),
                                'bended': [],
                                'direction': 'R'
                            }
                        
                            self.edmonds_graph.add_edge(current_node, 'stop', **pangenome_graph_edge_data)
                            connected = True
                
                    # print("connected is ", connected)
                    if connected == True:
                        for current_forward in current_forward_connected:
                            node_i, node_j, data = self.get_edge(current_node, current_forward, reverse = False)
                            self.edmonds_graph.add_edge(node_i, node_j, **data)
                            edmonds_graph_removed_edges.remove((current_node, current_forward))
                            
                        for current_backward in current_backward_connected:
                            node_i, node_j, data = self.get_edge(current_node, current_backward, reverse = True)
                            self.edmonds_graph.add_edge(node_i, node_j, **data)
                            edmonds_graph_removed_edges.remove((current_node, current_backward))
                            
                        resolved_nodes.add(current_node)
                
                    else:
                        if current_backward_connected:
                
                            # print('backward event', current_node, current_backward_connected)

                            # number = random.randint(0,len(current_backward_connected)-1)
                            number = max([(self.pangenome_graph.get_edge_data(current_node, backward)['weight'], i) for (i, backward) in enumerate(current_backward_connected)])[1]

                            node_i, node_j, data = self.get_edge(current_node, current_backward_connected[number], reverse = True)
                            # print('add ', node_i, ' to ', node_j)
                            # print('remove ', edmonds_graph_predecessors[current_node], ' to ', current_node)
                            self.edmonds_graph.remove_edge(edmonds_graph_predecessors[current_node], current_node)
                            self.edmonds_graph.add_edge(node_i, node_j, **data)
                            
                            edmonds_graph_removed_edges.remove((current_node, current_backward_connected[number]))
                            edmonds_graph_removed_edges.add((edmonds_graph_predecessors[current_node], current_node))
                
                            edmonds_graph_successors[edmonds_graph_predecessors[current_node]].remove(current_node)
                            edmonds_graph_successors[current_backward_connected[number]] += [current_node]

                            pred = edmonds_graph_predecessors[current_node]

                            edmonds_graph_predecessors.pop(current_node, None)
                            edmonds_graph_predecessors[current_node] = current_backward_connected[number] 
                
                            edmonds_graph_distances[current_node] = self.mean_edmonds_graph_path_weight('start', current_node)

                            # print("removed", edmonds_graph_removed_edges, "\nsuccessors", edmonds_graph_successors, "\npredecessors", edmonds_graph_predecessors)
                
                            resolved_nodes.add(current_node)
                            
                            # break
            
                visited_nodes.add(current_node)
            
                if not nx.is_directed_acyclic_graph(self.edmonds_graph):
                    print('Sanity Error.')
                    # print(pred, current_node)
                    exit()

        self.progress.end()

        remaining_stops = [node for node in self.edmonds_graph.nodes() if self.edmonds_graph.out_degree(node) == 0 and node != 'stop']
        self.run.info_single(f"{i} iterations to solve the graph.") 

        for stop in remaining_stops:

            pangenome_graph_edge_data = {
                'genome':{genome: {'draw': 'on'} for genome in self.genome_coloring.keys()},
                'weight':len(self.genome_coloring.keys()),
                'bended': [],
                'direction': 'R'
            }
        
            self.edmonds_graph.add_edge(stop, 'stop', **pangenome_graph_edge_data)

            edmonds_graph_successors[stop] += ['stop']

        if not nx.is_directed_acyclic_graph(self.edmonds_graph):
            print('Sanity Error.')
            exit()

        for x, generation in enumerate(nx.topological_generations(self.edmonds_graph)):
            self.x_list[x] = generation
            for node in generation:
                self.position[node] = (x, -1)

        remove = []
        for edge_i, edge_j in self.edmonds_graph.edges():
            if self.position[edge_j][0] - self.position[edge_i][0] > self.max_edge_length_filter:
                remove.append((edge_i, edge_j))
                print('y')
        self.edmonds_graph.remove_edges_from(remove)

        self.global_x = x

        self.ancest = nx.DiGraph(self.edmonds_graph)
        nx.set_edge_attributes(self.ancest, values=-1, name='weight')
        # self.ancest.remove_node('stop')
        layout_graph_nodes = list(self.ancest.nodes())
        layout_graph_successors = {layout_graph_node: list(self.ancest.successors(layout_graph_node)) for layout_graph_node in layout_graph_nodes}

        ghost = 0
        for y in range(self.global_x-1, 0, -1):
            for node in self.x_list[y]:
                node_x_position = self.position[node][0]

                change = []
                for successor in layout_graph_successors[node]:
                    if successor != 'stop':
                        successor_x_position = self.position[successor][0]
                        
                        if successor_x_position <= node_x_position:
                            print('Sanity Error.')
                            exit()
                        else:
                            change.append((successor_x_position, successor))

                if change:
                    if min(change)[0] > 1:
                        new_x_position = min([(x, n) for (x, n) in change if x > 1])[0] - 1
                        self.position[node] = (new_x_position, -1)

                    for (position, extend_successor) in change:
                        x_difference = position - new_x_position

                        path_list = [node]

                        for i in range(1, x_difference):
                            path_list += ['Ghost_' + str(ghost)]
                            self.position['Ghost_' + str(ghost)] = (new_x_position + i, -1)
                            ghost += 1

                        path_list += [extend_successor]

                        self.ancest.remove_edge(node, extend_successor)
                        
                        self.edges.append(path_list)
                        self.ancest.add_edges_from(map(tuple, zip(path_list, path_list[1:])), weight=-1)

        for i, j in self.ancest.edges():

            if self.position[j][0] - self.position[i][0] != 1 and i != 'start' and j != 'stop':
                
                # print(i, j, self.position[i][0], self.position[j][0])
                print('Sanity Error.')
                exit()

        # for pred in self.ancest.predecessors('stop'):
        #     print(pred)

        longest_path = nx.bellman_ford_path(G=self.ancest, source='start', target='stop', weight='weight')
        # print(len(longest_path))
        m = set(longest_path)

        dfs_list = list(nx.dfs_edges(self.ancest, source='start'))

        group = 0
        degree = dict(self.ancest.degree())
        groups = {}
        groups_rev = {}

        for node_v, node_w in dfs_list:
            if node_v != 'start' and degree[node_v] == 2 and degree[node_w] == 2:

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

        branches = {}
        # branches_rev = {}
        sortable = []
        for g in groups.keys():
            # print([(n, self.position[n]) for n in groups[g]],"\n")
            branch = groups[g]

            # for br in branch:
            #     branches_rev[br] = branch

            if not set(branch).isdisjoint(m) and not set(branch).issubset(m):
                print('Sanity Error.')
                break

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

            # branches_rev[node] = node

            if not set([n]).isdisjoint(m) and not set([n]).issubset(m):
                print('Sanity Error.')
                break

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


        # for i,j,k in sorted(sortable, key=lambda x: (x[1], x[0]), reverse = False):
        #     br = branches[i][j][k]
        #     if set(br).issubset(m):
        #         print(br)
        
        used = {}

        y_new = 0
        for node in longest_path:
            x_pos = self.position[node][0]
            self.position[node] = (x_pos, y_new)
            used[x_pos] = y_new

        stack = [longest_path]
        while stack:

            current = stack[0]

            # print(current)

            remove = True
            for i,j,k in sorted(sortable, key=lambda x: (x[1], x[0]), reverse = False):
                branch = branches[i][j][k]
                branch_pred = set(self.ancest.predecessors(branch[0]))
                branch_succ = set(self.ancest.successors(branch[-1]))
                if not branch_pred.isdisjoint(set(current)) or not branch_succ.isdisjoint(set(current)) or (not branch_pred.isdisjoint(set(current)) and not branch_succ.isdisjoint(set(current))):
                    
                    # print(branch, i, j, k)

                    remove = False
                    sortable.remove((i,j,k))

                    stack = [branch] + stack
                    y_new = max([used[s] for s in range(i, i+j)]) + 1

                    for node in branch:
                        x_pos = self.position[node][0]
                        self.position[node] = (x_pos, y_new)

                        # print((x_pos, y_new))

                        used[x_pos] = y_new

                        self.global_y = y_new if y_new > self.global_y else self.global_y

                    break

            if remove == True:
                stack.remove(current)

        #     print('\n')
            
        # print(sortable)
                
        if len(set(self.position.values())) != len(self.position.values()):
            print('Sanity Error.')
            exit()

        nx.set_edge_attributes(self.ancest, {(i, j): d for i, j, d in self.edmonds_graph.edges(data=True)})

        for edge in self.edges:            
            self.ancest.add_edge(edge[0], edge[-1], **self.edmonds_graph[edge[0]][edge[-1]])
            self.ancest[edge[0]][edge[-1]]['bended'] = [self.position[p] for p in edge[1:-1]]
            self.ancest.remove_nodes_from(edge[1:-1])

        nx.set_node_attributes(self.ancest, {k: d for k, d in self.edmonds_graph.nodes(data=True)})

        for node in self.ancest.nodes():
            self.ancest.nodes[node]['pos'] = self.position[node]

        # self.ancest.remove_edge('start', 'stop')

        self.run.info_single(f"Final graph {len(self.ancest.nodes())} nodes and {len(self.ancest.edges())} edges.")
        self.run.info_single(f"Done.")

        # self.position['start'] = (0, 0)
        # for z in range(1, x):
        #     print(z)
        #     if z in branches.keys():
        #         for item in sorted(branches[z].keys(), reverse=False):
        #             for num in branches[z][item]:
        #                 branch = branches[z][item][num]
        #                 predecessors = self.ancest.predecessors(branch[0])
        #                 start_y, predec = min([(self.position[pred][1], pred) for pred in predecessors])

        #                 print(start_y)
        #                 print(branch)

        #                 repeat = True
        #                 while repeat == True:
        #                     repeat = False
        #                     for i, node in enumerate(branch):
        #                         if (i+z, start_y) in self.position.values():
        #                             repeat = True
        #                             start_y += 1
        #                             print('this')
        #                             break
        #                         else:
        #                             self.position[node] = (i+z, start_y)
        #                             print('add ', node, ' ', i+z, start_y)

        # exit()

                # if possible_x_change != -1 and do_change == True:
                #     change.append((node, possible_x_change))

            # for (node, x_change) in change:
            #     self.position[node] = (x_change, 0)

        # print(self.position[',GC_00000004,GC_00000425'])

        # for x, generation in enumerate(nx.topological_generations(self.edmonds_graph)):
        #     nodes = {}
        #     self.x_list[x] = generation
        #     for node in generation:
        #         self.position[node] = (x, 0)
        #         node_list = node.split(',')

        #         if node_list[int(len(node_list)/2)] in nodes.keys():

        #             found = False
        #             for contractor in nodes[node_list[int(len(node_list)/2)]]:

        #                 intersection = set(self.edmonds_graph.nodes()[contractor]['genome'].keys()).intersection(set(self.edmonds_graph.nodes()[node]['genome'].keys()))
        #                 if not intersection:

        #                     self.edmonds_graph.nodes()[contractor]['weight'] += self.edmonds_graph.nodes()[node]['weight']
        #                     self.edmonds_graph.nodes()[contractor]['genome'].update(self.edmonds_graph.nodes()[node]['genome'])

        #                     nx.contracted_nodes(self.edmonds_graph, contractor, node, copy=False)
        #                     self.fusion_events += 1
        #                     found = True
        #                     break

        #             if found == False:
        #                 nodes[node_list[int(len(node_list)/2)]] += [node]

        #         else:
        #             nodes[node_list[int(len(node_list)/2)]] = [node]

        # edmonds_graph_edges = list(self.edmonds_graph.edges())
        # for i, j in edmonds_graph_edges:
        #     if abs(self.position[j][0] - self.position[i][0]) > self.max_edge_length_filter and self.max_edge_length_filter != -1:
        #         self.edmonds_graph.remove_edge(i, j)
        #         self.removed_edges += 1

        # self.run.info_single(f"{self.fusion_events} fusion events.") 
        # self.run.info_single(f"{self.removed_edges} removed edges due to length cutoff.")
        # self.run.info_single(f"Done.")

    # TODO Speed up component path finding (multithreading)
    # def calculate_component_paths(self):

    #     self.run.warning(None, header="Extracting component paths from F", lc="green")

    #     self.progress.new("Solving Path")

    #     edmonds_graph_edges = list(self.edmonds_graph.edges())
    #     number = len(str(len(edmonds_graph_edges)))

    #     j = 0
    #     for i, (node_i, node_j) in enumerate(edmonds_graph_edges):

    #         self.progress.update(f"{str(i).rjust(number, ' ')} / {len(edmonds_graph_edges)}")

    #         if nx.has_path(self.edmonds_graph, 'start', node_i) and nx.has_path(self.edmonds_graph, node_j, 'stop'):

    #             path_leaf = nx.shortest_path(self.edmonds_graph, 'start', node_i, method='bellman-ford')
    #             path_succ = nx.shortest_path(self.edmonds_graph, node_j, 'stop', method='bellman-ford')
    #             full_path = path_leaf + path_succ

    #             value = nx.path_weight(self.edmonds_graph, full_path, 'weight')/len(full_path)

    #             self.leaf_path.append((value, full_path))

    #         else:
    #             j += 1

    #     self.progress.end()

    #     self.run.info_single(f"Removed {j} unsolvable edges.")
    #     self.run.info_single("Done.")

    # # ANCHOR Sub path calculation script
    # def calculate_unknown_edges(self, path, known):

    #     ancest_nodes = list(self.ancest.nodes())

    #     unknown_edges = []
    #     sub_edges = []

    #     for k, o in map(tuple, zip(path, path[1:])):
    #         if not (k, o) in known:

    #             if k in ancest_nodes and o in ancest_nodes:
    #                 if sub_edges:
    #                     unknown_edges.append(sub_edges)
    #                     sub_edges = []

    #                 unknown_edges.append([(k, o)])

    #             elif o in ancest_nodes:
    #                 sub_edges.append((k, o))
    #                 unknown_edges.append(sub_edges)
    #                 sub_edges = []

    #             else:
    #                 sub_edges.append((k, o))

    #         else:
    #             if sub_edges:
    #                 unknown_edges.append(sub_edges)
    #                 sub_edges = []

    #     if sub_edges:
    #         unknown_edges.append(sub_edges)
    #         sub_edges = []

    #     return(unknown_edges)

    # # ANCHOR Main position calculation
    # # TODO Recalculate Topo Coordinated
    # # It is possible that the topological x positions have to be recalculated here as
    # # change of the graph can happen after the first topological sorting by removing / adding
    # # more edges.
    # def run_synteny_layout_algorithm(self):

    #     self.run.warning(None, header="Calculating graph P node positions", lc="green")

    #     self.ancest.add_edge("start", "stop", weight=1)
    #     known = set([('start', 'stop')])

    #     paths = [value for value in sorted(self.leaf_path, key=lambda x: x[0], reverse=True)]
    #     self.progress.new("Running path")

    #     number = len(str(len(paths)))
    #     for i, (_, path) in enumerate(paths):
    #         self.progress.update(f"{str(i).rjust(number, ' ')} / {len(paths)}")

    #         # try:
    #         unknown_edges = self.calculate_unknown_edges(path, known)

    #         # print("1", unknown_edges)

    #         for sub_edges in unknown_edges:

    #             # print("2", sub_edges)

    #             known.update(sub_edges)

    #             node_start = sub_edges[0][0]
    #             node_stop = sub_edges[-1][1]

    #             sub_path = path[path.index(node_start): path.index(node_stop)+1]

    #             node_first_x, node_first_y = self.position[sub_path[1]]
    #             node_last_x, node_last_y = self.position[sub_path[-2]]

    #             node_start_x, node_start_y = self.position[node_start]
    #             node_stop_x, node_stop_y = self.position[node_stop]

    #             # print("3", node_start, node_start_x, node_start_y)
    #             # print("4", node_stop, node_stop_x, node_stop_y)

    #             for z in range(node_start_x, node_stop_x+1):
    #                 if sub_path[z-node_start_x] not in self.x_list[z]:

    #                     sub_path = sub_path[:z-node_start_x] + ["Ghost_" + str(self.ghost)] + sub_path[z-node_start_x:]

    #                     self.x_list[z].append("Ghost_" + str(self.ghost))

    #                     self.ghost += 1

    #             # print("5", self.x_list)

    #             curr_path = []
    #             for s in sub_path:
    #                 if not s.startswith('Ghost_'):
    #                     if len(curr_path) > 1:
    #                         curr_path.append(s)
    #                         self.edges.append(curr_path)
    #                     curr_path = [s]
    #                 else:
    #                     curr_path.append(s)

    #             # print("6", self.x_list)

    #             sub_path = sub_path[1:-1]

    #             # print("7", sub_path)

    #             if node_start == 'start':
    #                 next_y = self.y_shifting(sub_path, node_first_x-1, node_stop_x, node_start_y, node_stop_y)
    #             elif node_stop == 'stop':
    #                 next_y = self.y_shifting(sub_path, node_start_x, node_last_x+1, node_start_y, node_stop_y)
    #             elif node_start == 'start' and node_stop == 'stop':
    #                 next_y = self.y_shifting(sub_path, node_first_x-1, node_last_x+1, node_start_y, node_stop_y)
    #             else:
    #                 next_y = self.y_shifting(sub_path, node_start_x, node_stop_x, node_start_y, node_stop_y)

    #             # print("8", next_y)

    #             self.add_new_edges(sub_path, next_y, node_start, node_stop, node_start_x, node_stop_x)

    #             to_be_removed = set()
    #             if node_start == 'start':
    #                 t = node_start_x + 1
    #                 for node in sub_path:
    #                     if node.startswith('Ghost_'):
    #                         self.x_list[t].remove(node)
    #                         self.ancest.remove_node(node)

    #                         to_be_removed.add(node)
    #                         self.position.pop(node)
    #                         self.path.pop(node)
    #                         t += 1
    #                     else:
    #                         break

    #             if node_stop == 'stop':
    #                 s = node_stop_x - 1
    #                 for node in sub_path[::-1]:
    #                     if node.startswith('Ghost_'):
    #                         self.x_list[s].remove(node)
    #                         self.ancest.remove_node(node)

    #                         to_be_removed.add(node)
    #                         self.position.pop(node)
    #                         self.path.pop(node)
    #                         s -= 1
    #                     else:
    #                         break

    #             # print("9", self.x_list)

    #             cut_sub_path = [node for node in sub_path if node not in to_be_removed]

    #             for node in cut_sub_path:
    #                 self.path[node] = cut_sub_path

    #             # print("10", self.path)

    #             # print("11", self.position)

    #             # print("12", sub_path)

    #             # exit()

    #         # except Exception as error:
    #         #     print('Sanity Error')
    #         #     exit()

    #     self.progress.end()

    #     if len(list(self.position.values())) != len(set(self.position.values())):  
    #         print(len(list(self.position.values())) - len(set(self.position.values())))

    #     nx.set_edge_attributes(self.ancest, {(i, j): d for i, j, d in self.edmonds_graph.edges(data=True)})

    #     for edge in self.edges:
    #         if edge[-1] != 'stop' and edge[0] != 'start':
    #             self.ancest.add_edge(edge[0], edge[-1], **self.edmonds_graph[edge[0]][edge[-1]])
    #             self.ancest[edge[0]][edge[-1]]['bended'] = [self.position[p] for p in edge[1:-1]]
    #             self.ancest.remove_nodes_from(edge[1:-1])

    #     nx.set_node_attributes(self.ancest, {k: d for k, d in self.edmonds_graph.nodes(data=True)})

    #     for node in self.ancest.nodes():
    #         self.ancest.nodes[node]['pos'] = self.position[node]

    #     self.ancest.remove_edge('start', 'stop')

    #     self.run.info_single(f"Final graph {len(self.ancest.nodes())} nodes and {len(self.ancest.edges())} edges.")
    #     self.run.info_single(f"Done.")

        # print(self.position['GC_00000071,GC_00000004,GC_00000125'])

    # ANCHOR Gene Cluster grouping
    # TODO Degree is calculated by pangenome graph not edmonds graph probably not a bad idea due
    # to easier adding of additional edges.
    # def condense_gene_clusters_into_groups(self):

    #     self.run.warning(None, header="Grouping GCs to gene cluster groups (GCGs)", lc="green")

    #     if self.gene_cluster_grouping_threshold == -1:
    #         self.run.info_single("Setting algorithm to 'no grouping'")
    #     else:
    #         self.run.info_single(f"Setting algorithm to 'Grouping single connected chains size > {str(self.gene_cluster_grouping_threshold)}'")

    #     dfs_list = list(nx.dfs_edges(self.ancest, source='start'))

    #     group = 0
    #     degree = dict(self.ancest.degree())
    #     groups = {}
    #     groups_rev = {}

    #     for node_v, node_w in dfs_list:

    #         print(node_v, self.ancest.nodes[node_v].keys())

    #         if node_v != 'start' and node_w != 'stop' and degree[node_v] == 2 and degree[node_w] == 2 and set(self.ancest.nodes[node_v]['genome'].keys()) == set(self.ancest.nodes[node_w]['genome'].keys()):

    #             if node_v not in groups_rev.keys():
    #                 group_name = 'GCG_' + str(group).zfill(8)
    #                 groups[group_name] = [node_v, node_w]
    #                 groups_rev[node_v] = group_name
    #                 groups_rev[node_w] = group_name
    #                 group += 1

    #             else:
    #                 group_name = groups_rev[node_v]
    #                 groups[group_name] += [node_w]
    #                 groups_rev[node_w] = group_name

    #     for label, condense_nodes in groups.items():

    #         if len(condense_nodes) >= self.gene_cluster_grouping_threshold and self.gene_cluster_grouping_threshold != -1:
    #             self.grouping[label] = condense_nodes

    #     self.run.info_single(f"Created {len(groups.keys())} groups.")
    #     self.run.info_single("Done")

    # ANCHOR y-shifting script
    # TODO Wrong hierarchy bug:
    # Sometimes very small branches are included on top of way longer and higher weighted ones I'm currently
    # not completely sure why this occures and have to solve it.
    # def y_shifting(self, sub_path, node_start_x, node_stop_x, node_start_y, node_stop_y):
    #     current_start_x = node_start_x + 1
    #     current_stop_x = node_stop_x - 1
    #     current_path_length = (current_stop_x - current_start_x) - 1
    #     current_y = max(node_start_y, node_stop_y) + 1

    #     next_y = -1

    #     increase_layer = []

    #     while current_y <= self.global_y + 1:

    #         node = ''
    #         current_layer_start_x = self.global_x
    #         current_layer_stop_x = 0
    #         layer_branches = []
    #         sub_branch = []
    #         layer_size = 0

    #         z = current_start_x
    #         while z <= current_stop_x:
    #             for check in self.x_list[z]:
    #                 if check not in sub_path and self.position[check] == (z, current_y):
    #                     node = check

    #                     if not sub_branch or node not in sub_branch:

    #                         sub_branch = self.path[node]

    #                         if (current_y, sub_branch) not in layer_branches:
    #                             layer_branches.append((current_y, sub_branch))

    #                             sub_branch_start_x, _ = self.position[sub_branch[0]]
    #                             sub_branch_stop_x, _ = self.position[sub_branch[-1]]

    #                             layer_size += sub_branch_stop_x - sub_branch_start_x + 1

    #                             current_layer_start_x = sub_branch_start_x if sub_branch_start_x < current_layer_start_x else current_layer_start_x
    #                             current_layer_stop_x = sub_branch_stop_x if sub_branch_stop_x > current_layer_stop_x else current_layer_stop_x

    #                             current_start_x = sub_branch_start_x if sub_branch_start_x < current_start_x else current_start_x
    #                             current_stop_x = sub_branch_stop_x if sub_branch_stop_x > current_stop_x else current_stop_x

    #                             z = current_layer_stop_x

    #             z += 1

    #         if layer_size > current_path_length:

    #             if next_y == -1:
    #                 next_y = current_y

    #         if next_y != -1:
    #             increase_layer.extend(layer_branches)

    #         if not sub_branch:
    #             if current_y != max(node_start_y, node_stop_y):
    #                 break
    #             else:
    #                 current_y += 1
    #         else:
    #             current_y += 1

    #     for _, layer_branch in sorted(increase_layer, reverse=True):
    #         for node_branch in layer_branch:
    #             node_branch_x = self.position[node_branch][0]
    #             node_branch_y = self.position[node_branch][1]

    #             self.position[node_branch] = (node_branch_x, node_branch_y + 1)

    #     if next_y == -1:
    #         next_y = current_y

    #     self.global_y = current_y if current_y > self.global_y else self.global_y
    #     self.global_x = self.position["stop"][0] if self.position["stop"][0] > self.global_x else self.global_x

    #     return(next_y)

    # # ANCHOR adding new nodes and edges
    # # TODO take a closer look at the len(sub_path) == 0 situation for now added if sub_path
    # # keep in mind this can become a error later
    # def add_new_edges(self, sub_path, next_y, node_start, node_stop, node_start_x, node_stop_x):
    #     if sub_path:
    #         if sub_path[0].startswith('Ghost_'):
    #             self.ancest.add_edge(node_start, sub_path[0], weight=0)

    #         else:
    #             self.ancest.add_edge(node_start, sub_path[0], weight=1)

    #         for i, new in enumerate(sub_path, 1):
    #             self.path[new] = sub_path

    #             if new == sub_path[-1]:
    #                 self.position[new] = (node_stop_x - 1, next_y)
    #             else:
    #                 self.position[new] = (node_start_x + i, next_y)

    #             if new != sub_path[-1]:
    #                 if sub_path[i].startswith('Ghost_') or new.startswith('Ghost_') or (sub_path[i].startswith('Ghost_') and new.startswith('Ghost_')):
    #                     self.ancest.add_edge(new, sub_path[i], weight=0)
    #                 else:
    #                     self.ancest.add_edge(new, sub_path[i], weight=1)

    #         if sub_path[-1].startswith('Ghost_'):
    #             self.ancest.add_edge(sub_path[-1], node_stop, weight=0)
    #         else:
    #             self.ancest.add_edge(sub_path[-1], node_stop, weight=1)

    # ANCHOR Converting network to JSON data
    # TODO rework that section for better debugging and add more features as an example fuse start and top
    # together or remove both so the graph becomes a REAL circle. Aside from that there is a bug in the remove
    # edges section for (k,o) in circular edges and for (k,o) in pangenome edges. Change and test while reworking.
    def get_json_dict_for_graph(self):
        self.run.warning(None, header="Exporting network to JSON", lc="green")

        jsondata = {}

        instances_pangenome_graph = 0
        for node, attr in self.pangenome_graph.nodes(data=True):
            instances_pangenome_graph += len(list(attr['genome'].keys()))

        instances_ancest_graph = 0
        for node, attr in self.ancest.nodes(data=True):
            instances_ancest_graph += len(list(attr['genome'].keys()))

        instances_pangenome_graph = 0
        for node, attr in self.pangenome_graph.nodes(data=True):
            if node != 'start' and node != 'stop':
                instances_pangenome_graph += len(list(attr['genome'].keys()))

        instances_ancest_graph = 0
        for node, attr in self.ancest.nodes(data=True):
            if node != 'start' and node != 'stop':
                instances_ancest_graph += len(list(attr['genome'].keys()))

        self.run.info_single(f"Total fraction of recovered genecall information {round((instances_ancest_graph/instances_pangenome_graph)*100, 3)}%")
        self.run.info_single(f"Total fraction of recovered geneclusters {round((len(self.ancest.nodes())+self.fusion_events-2)/(len(self.pangenome_graph.nodes())-2)*100, 3)}%")

        self.ancest.remove_nodes_from(['start', 'stop'])

        # NOTE: Any change in `jsondata` will require the pangraph JSON in anvio.tables.__init__
        #       to incrase by one (so that the new code that works with the new structure requires
        #       previously generated JSON to be recomputed).
        jsondata["infos"] = {'meta': {'version': anvio.__pangraph__version__,
                                      'global_x': self.global_x,
                                      'global_y': self.global_y},
                             'genomes': self.genome_coloring,
                             'functional_annotation_sources_available': self.functional_annotation_sources_available,
                             'num_genomes': len(self.genome_coloring),
                             'max_edge_length_filter': self.max_edge_length_filter,
                             'gene_cluster_grouping_threshold': self.gene_cluster_grouping_threshold,
                             'original': {'nodes': len(self.pangenome_graph.nodes()),
                                          'edges': len(self.pangenome_graph.edges()),
                                          'instances': instances_pangenome_graph},
                             'visualization': {'nodes': len(self.ancest.nodes()),
                                               'edges': len(self.ancest.edges()),
                                               'instances': instances_ancest_graph},
                             'data': list(self.ancest.graph.items()),
                             'directed': self.ancest.is_directed(),
                             'groups': self.grouping}

        jsondata['infos']['grouped'] = {'nodes': len(self.ancest.nodes()), 'edges': len(self.ancest.edges())}

        jsondata["elements"] = {"nodes": {}, "edges": {}}

        for i, j in self.ancest.nodes(data=True):
            # if i != 'start' and i != 'stop':
            jsondata["elements"]["nodes"][str(i)] = {
                "name": j["name"],
                "weight": j["weight"],
                "genome": j["genome"],
                "position": {
                    'x': j['pos'][0],
                    'y': j['pos'][1]
                }
            }

        for l, (k, o, m) in enumerate(self.ancest.edges(data=True)):
            # if k != 'start' and k != 'stop' and o != 'start' and o != 'stop':
            jsondata["elements"]["edges"]['E_' + str(l).zfill(8)] = {
                "source": k,
                "target": o,
                "genome": m["genome"],
                "weight": m["weight"],
                "direction": m["direction"],
                "bended": [{'x': x, 'y': y} for x, y in m["bended"]] if len(m["bended"]) != 0 else ""
            }

        self.run.info_single("Done.")

        return(jsondata)

    def store_network(self):
        """Function to store final graph structure in a pan-db and/or JSON flat text output file"""

        self.run.warning(None, header="Storing the graph structure", lc="green")

        if self.json_output_file_path:
            with open(self.json_output_file_path, 'w') as output:
                output.write(json.dumps(self.get_json_dict_for_graph(), indent=2))
            self.run.info("JSON output file", os.path.abspath(self.json_output_file_path))
        else:
            self.run.info("JSON output file", "Skipped (but OK)", mc='red')

        if not self.skip_storing_in_pan_db:
            raise ConfigError("The storage of graph data in pan-db is not yet implemented :/")
        else:
            self.run.info("Into the pan-db", "Skipped (but OK)", mc='red')


    def get_genome_gc_fused(self):
        return(self.genome_gc_fused)