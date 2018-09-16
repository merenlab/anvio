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
import multiprocessing
import pandas as pd
from itertools import chain

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.tables.miscdata as miscdata

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
        self.skip_homogeneity = self.args.skip_homogeneity
        self.quick_homogeneity = self.args.quick_homogeneity
        self.align_with = A('align_with')
        self.overwrite_output_destinations = A('overwrite_output_destinations')
        self.debug = anvio.DEBUG
        self.min_percent_identity = A('min_percent_identity')
        self.gene_cluster_min_occurrence = A('min_occurrence')
        self.mcl_inflation = A('mcl_inflation')
        self.sensitive = A('sensitive')
        self.minbit = A('minbit')
        self.use_ncbi_blast = A('use_ncbi_blast')
        self.exclude_partial_gene_calls = A('exclude_partial_gene_calls')
        self.description_file_path = A('description')
        self.skip_hierarchical_clustering = A('skip_hierarchical_clustering')
        self.enforce_hierarchical_clustering = A('enforce_hierarchical_clustering')

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
                       'diamond_sensitive': self.sensitive,
                       'minbit': self.minbit,
                       'exclude_partial_gene_calls': self.exclude_partial_gene_calls,
                       'gene_alignments_computed': False if self.skip_alignments else True,
                       'genomes_storage_hash': self.genomes_storage.get_storage_hash(),
                       'project_name': self.project_name,
                       'gene_clusters_clustered': False,
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
            raise ConfigError("Please set a project name using the `--project-name` parameter, and be prepared to see\
                               it around as anvi'o will use it for multiple things, such as setting the output directory\
                               and naming various output files including the database file that will be generated at the\
                               end of the process. If you set your own output directory name, you can have multiple\
                               projects in it and all of those projects can use the same intermediate files whenever\
                               possible.")

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
            raise ConfigError("Minimum percent identity must be between 0%% and 100%%. Although your %.2f%% is\
                               pretty cute, too." % self.min_percent_identity)


        if len([c for c in list(self.genomes.values()) if 'genome_hash' not in c]):
            raise ConfigError("self.genomes does not seem to be a properly formatted dictionary for\
                               the anvi'o class Pangenome.")

        if self.enforce_hierarchical_clustering and self.skip_hierarchical_clustering:
            raise ConfigError("You are confusing anvi'o :/ You can't tell anvi'o to skip hierarchical clustering\
                               while also asking it to enforce it.")

        if self.description_file_path:
            filesnpaths.is_file_plain_text(self.description_file_path)
            self.description = open(os.path.abspath(self.description_file_path), 'rU').read()

        self.pan_db_path = self.get_output_file_path(self.project_name + '-PAN.db')


    def run_diamond(self, unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict):
        diamond = Diamond(unique_AA_sequences_fasta_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads, overwrite_output_destinations=self.overwrite_output_destinations)

        diamond.names_dict = unique_AA_sequences_names_dict
        diamond.search_output_path = self.get_output_file_path('diamond-search-results')
        diamond.tabular_output_path = self.get_output_file_path('diamond-search-results.txt')

        diamond.sensitive = self.sensitive

        return diamond.get_blast_results()


    def run_blast(self, unique_AA_sequences_fasta_path, unique_AA_sequences_names_dict):
        self.run.warning("You elected to use NCBI's blastp for amino acid sequence search. Running blastp will \
                          be significantly slower than DIAMOND (although, anvi'o developers are convinced that \
                          you *are* doing the right thing, so, kudos to you).")
        blast = BLAST(unique_AA_sequences_fasta_path, run=self.run, progress=self.progress,
                          num_threads=self.num_threads, overwrite_output_destinations=self.overwrite_output_destinations)

        blast.names_dict = unique_AA_sequences_names_dict
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
        self.progress.update('(initial pass of the serach results to set the self bit scores ...)')
        for line in open(blastall_results):
            fields = line.strip().split('\t')

            try:
                query_id, subject_id, perc_id, aln_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score = \
                    [mapping[i](fields[i]) for i in range(0, len(mapping))]
            except Exception as e:
                self.progress.end()
                raise ConfigError("Something went wrong while processing the blastall output file in line %d.\
                                    Here is the error from the uppoer management: '''%s'''" % (line_no, e))
            line_no += 1
            all_ids.add(query_id)
            all_ids.add(subject_id)

            if query_id == subject_id:
                self_bit_scores[query_id] = bit_score

        self.progress.end()

        ids_without_self_search = all_ids - set(self_bit_scores.keys())
        if len(ids_without_self_search):
            search_tool = 'BLAST' if self.use_ncbi_blast else 'DIAMOND'
            self.run.warning("%s did not retun search results for %d of %d the amino acid sequences in your input FASTA file.\
                              Anvi'o will do some heuristic magic to complete the missing data in the search output to recover\
                              from this. But since you are a scientist, here are the amino acid sequence IDs for which %s\
                              failed to report self search results: %s." \
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
                raise ConfigError("Something horrible happened. This can only happend if you started a new analysis with\
                                    additional genomes without cleaning the previous work directory. Sounds familiar?")

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

        for gene_cluster in gene_clusters:
            self.view_data[gene_cluster] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.view_data_presence_absence[gene_cluster] = dict([(genome_name, 0) for genome_name in self.genomes])
            self.additional_view_data[gene_cluster] = {'num_genes_in_gene_cluster': 0, 'num_genomes_gene_cluster_has_hits': 0, 'SCG': 0, 'max_num_paralogs': 0}

            for gene_entry in gene_clusters_dict[gene_cluster]:
                genome_name = gene_entry['genome_name']

                self.view_data[gene_cluster][genome_name] += 1
                self.view_data_presence_absence[gene_cluster][genome_name] = 1
                self.additional_view_data[gene_cluster]['num_genes_in_gene_cluster'] += 1

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
                self.run.warning("You have %s gene_clusters, which exceeds the number of gene_clusters anvi'o is comfortable to cluster. But\
                                  since you have used the flag `--enforce-hierarchical-clustering`, anvi'o will attempt\
                                  to create a hierarchical clustering of your gene_clusters anyway. It may take a bit of \
                                  time. Pour yourself a coffee. Or go to a nice vacation. See you in 10 mins, or next year \
                                  or never." % pp(len(gene_clusters_dict)))
            else:
                self.run.warning("It seems you have %s gene clusters in your pangenome. This exceeds the soft limit\
                                  of %s for anvi'o to attempt to create a hierarchical clustering of your gene clusters\
                                  (which becomes the center tree in all anvi'o displays). If you want a hierarchical\
                                  clustering to be done anyway, please see the flag `--enforce-hierarchical-clustering`." \
                                            % (pp(len(gene_clusters_dict)), pp(self.max_num_gene_clusters_for_hierarchical_clustering)))
                self.skip_hierarchical_clustering = True

        ########################################################################################
        #                           STORING FILTERED DATA IN THE DB
        ########################################################################################
        table_structure=['gene_cluster'] + sorted(self.genomes.keys())
        table_types=['text'] + ['numeric'] * len(self.genomes)
        TablesForViews(self.pan_db_path).create_new_view(
                                        data_dict=self.view_data,
                                        table_name='gene_cluster_frequencies',
                                        table_structure=table_structure,
                                        table_types=table_types,
                                        view_name = 'gene_cluster_frequencies')

        TablesForViews(self.pan_db_path).create_new_view(
                                        data_dict=self.view_data_presence_absence,
                                        table_name='gene_cluster_presence_absence',
                                        table_structure=table_structure,
                                        table_types=table_types,
                                        view_name = 'gene_cluster_presence_absence')

        item_additional_data_table = miscdata.TableForItemAdditionalData(self.args)
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
           from the same contig. Of course, the synteny does not mean much for genmes that
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

            dbops.add_items_order_to_db(self.pan_db_path, order_name, ','.join(gene_clusters_order_based_on_genome_synteny), order_data_type_newick=False, run=self.run)

        gene_cluster_gene_cluster_edges = []
        # network description of gene_cluster-gene_cluster relationshops given the gene synteny.
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
                                                      distance=self.distance, linkage=self.linkage, run=self.run, progress=self.progress)


    def populate_gene_cluster_homogeneity_index(self, gene_clusters_dict):
        if self.skip_alignments:
            self.run.warning('Skipping homogeneity calculations because gene clusters are not alligned.')
            return

        if self.skip_homogeneity:
            self.run.warning("Skipping homogeneity calculations per the '--skip-homogeneity' flag.")
            return
        
        pan = dbops.PanSuperclass(args=self.args)
        gene_cluster_names = set(list(gene_clusters_dict.keys()))

        d = pan.compute_homogeneity_indices_for_gene_clusters(gene_cluster_names=gene_cluster_names, num_threads=self.num_threads)

        if d is None:
            self.run.warning("Anvi'o received an empty dictionary for homogeneity indices. Not good :/ Returning empty handed,\
                              without updating anything in the pan database...")
            return

        miscdata.TableForItemAdditionalData(self.args).add(d, ['functional_homogeneity_index', 'geometric_homogeneity_index'], skip_check_names=True)


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

        layers_additional_data_keys.extend(['num_genes', 'avg_gene_length', 'num_genes_per_kb'])

        for genome_name in self.genomes:
            new_dict = {}
            for key in layers_additional_data_keys:
                new_dict[key] = self.genomes[genome_name][key]

            layers_additional_data_dict[genome_name] = new_dict

        # summarize gene cluster stats across genomes
        layers_additional_data_keys.extend(['num_gene_clusters', 'singleton_gene_clusters'])
        for genome_name in self.genomes:
            layers_additional_data_dict[genome_name]['num_gene_clusters'] = 0
            layers_additional_data_dict[genome_name]['singleton_gene_clusters'] = 0
            for gene_cluster in self.view_data_presence_absence:
                genomes_contributing_to_gene_cluster = [t[0] for t in self.view_data_presence_absence[gene_cluster].items() if t[1]]

                # tracking singletons
                if len(genomes_contributing_to_gene_cluster) == 1 and genomes_contributing_to_gene_cluster[0] == genome_name:
                    layers_additional_data_dict[genome_name]['singleton_gene_clusters'] += 1

                # tracking the total number of gene clusters
                if self.view_data_presence_absence[gene_cluster][genome_name]:
                    layers_additional_data_dict[genome_name]['num_gene_clusters'] += 1

        self.progress.end()

        miscdata.TableForLayerOrders(self.args).add(layer_orders_data_dict)
        miscdata.TableForLayerAdditionalData(self.args).add(layers_additional_data_dict, layers_additional_data_keys)


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

        if self.skip_alignments and self.align_with:
            raise ConfigError("You are asking anvi'o to skip aligning sequences within your gene clusters, and then you \
                               are also asking it to use '%s' for aligning sequences within your gene clusters. It is easy \
                               to ignore this and skip the alignment, but anvi'o gets nervous when it realizes her users are \
                               being inconsistent. Please make up your mind, and come back as the explicit person you are" \
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

        self.run.info('gene clusters info', '%d gene_clusters stored in the database' % len(gene_clusters_dict))


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
                    raise ConfigError("Something horrible happened. This can only happen if you started a new analysis with\
                                        additional genomes without cleaning the previous work directory. Sounds familiar?")

                gene_clusters_dict[gene_cluster].append({'gene_caller_id': int(gene_caller_id), 'gene_cluster_id': gene_cluster, 'genome_name': genome_name, 'alignment_summary': ''})

        self.progress.end()

        return gene_clusters_dict


    def compute_alignments_for_gene_clusters(self, gene_clusters_dict):
        if self.skip_alignments:
            self.run.warning('Skipping gene alignments.')
            return gene_clusters_dict

        # we run "select aligner" to print the citation information (the actual selection is
        # done in the `alignment_worker` down below)
        aligners.select(self.align_with)

        self.progress.new('Aligning amino acid sequences for genes in gene clusters')
        self.progress.update('...')
        gene_cluster_names = list(gene_clusters_dict.keys())

        # we only need to align gene clusters with more than one sequence
        non_singleton_gene_cluster_names = [g for g in gene_cluster_names if len(gene_clusters_dict[g]) > 1]
        num_non_singleton_gene_clusters = len(non_singleton_gene_cluster_names)

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
        while received_gene_clusters < num_non_singleton_gene_clusters:
            try:
                gene_clusters_item = output_queue.get()

                if gene_clusters_item:
                    # worker returns None if there is nothing to align
                    # we do not need to owerwrite it to gene_clusters_dict
                    gene_clusters_dict[gene_clusters_item['name']] = gene_clusters_item['entry']

                if self.debug:
                    print(json.dumps(gene_clusters_item, indent=2))

                received_gene_clusters += 1
                self.progress.update("Processed %d of %d non-singleton GCs in %d threads." %
                    (received_gene_clusters, num_non_singleton_gene_clusters, self.num_threads))

            except KeyboardInterrupt:
                print("Anvi'o profiler recieved SIGINT, terminating all processes...")
                break

        for worker in workers:
            worker.terminate()

        self.progress.end()

        return gene_clusters_dict


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

            # sometimes alignments fail, and because pangenomic analyses can take forever,
            # everything goes into the trash bin. to prevent that, here we have a try/except
            # block with lots of warnings if something goes wrong.
            try:
                alignments = aligner(run=r).run_stdin(gene_sequences_in_gene_cluster)
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
                    with open(temp_file_path, 'w') as output:
                        for tpl in gene_sequences_in_gene_cluster:
                            output.write('>%s\n%s\n' % (tpl[0], tpl[1]))
                    debug_info = "The %d sequences in gene cluster %s are stored in the temporary file '%s'" % \
                                        (len(gene_sequences_in_gene_cluster), gene_cluster_name, temp_file_path)
                else:
                    debug_info = "If you re-run your last command with a `--debug` flag, anvi'o will generate more\
                                  information for you about the contenets of this gene cluster (but if you are seeing\
                                  millions of these warnings, it may not be a good idea since with the `--debug` flag\
                                  anvi'o will generate a FASTA file in a temporary directory with the contents of the\
                                  gene cluster, and will not attempt to delete them later)."

                run.warning("VERY BAD NEWS. The alignment of seqeunces with '%s' in the gene cluster '%s' failed\
                             for some reason. Since the real answer to 'why' is too deep in the matrix, there is\
                             no reliable solution for anvi'o to find it for you, BUT THIS WILL AFFECT YOUR SCIENCE\
                             GOING FORWARD, SO YOU SHOULD CONSIDER ADDRESSING THIS ISSUE FIRST. %s" % \
                                                       (aligner.__name__, gene_cluster_name, debug_info), nl_before=1)


            output = {'name': gene_cluster_name, 'entry': copy.deepcopy(gene_clusters_dict[gene_cluster_name])}
            for gene_entry in output['entry']:
                gene_entry['alignment_summary'] = utils.summarize_alignment(alignments['%s_%d' % (gene_entry['genome_name'], gene_entry['gene_caller_id'])])

            output_queue.put(output)


    def process(self):
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
        gene_clusters_dict = self.compute_alignments_for_gene_clusters(gene_clusters_dict)

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
        self.populate_gene_cluster_homogeneity_index(gene_clusters_dict)
        
        # done
        self.run.info('log file', self.run.log_file_path)
        self.run.quit()
