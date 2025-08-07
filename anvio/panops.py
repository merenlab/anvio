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
import yaml
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt

from itertools import chain
from scipy.optimize import curve_fit

# multiprocess is a fork of multiprocessing that uses the dill serializer instead of pickle
# using the multiprocessing module directly results in a pickling error in Python 3.10 which
# goes like this:
#
#   >>> AttributeError: Can't pickle local object 'SOMEFUNCTION.<locals>.<lambda>' multiprocessing
#
import multiprocess as multiprocessing

import anvio
import anvio.tables as t
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.clustering as clustering
import anvio.filesnpaths as filesnpaths
import anvio.tables.miscdata as miscdata

from anvio.dbinfo import DBInfo
from anvio.drivers.mcl import MCL
from anvio.drivers import Aligners
from anvio.drivers.blast import BLAST
from anvio.drivers.diamond import Diamond

from anvio.genomestorage import GenomeStorage
from anvio.tables.views import TablesForViews
from anvio.tables.states import TablesForStates
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.geneclusters import TableForGeneClusters
from anvio.tables.pangraphdata import TableForNodes, TableForEdges

from anvio.directedforce import DirectedForce
from anvio.topologicallayout import TopologicalLayout
from anvio.syntenygenecluster import SyntenyGeneCluster
from anvio.pangenomegraphmaster import PangenomeGraphManager


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


class RarefactionAnalysis:
    """Takes in a pangenome, calculates rarefaction curves and Heaps' Law fit to assess the openness of the pangenome.

        >>> import argparse
        >>> args = argparse.Namespace(pan_db="PATH/TO/PAN.db", iterations=100, output_file='rarefaction_curves.svg')
        >>> rarefaction_curves = RarefactionAnalysis(args)
        >>> k, alpha = rarefaction_curves.process()

    A client of this class is the program `anvi-compute-rarefaction-curves`
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.pan_db_path = A('pan_db')
        self.iterations = A('iterations') or 100
        self.output_file_prefix = A('output_file_prefix')
        self.skip_output_files= A('skip_output_files')

        # to be filled from the pan-db
        self.gene_cluster_data = None
        self.unique_genomes = None
        self.num_genomes = None

        # Load gene cluster data
        self.load_data()


    def load_data(self):
        """Load gene cluster data from pan-db, setup essential variables, sanity check."""

        utils.is_pan_db(self.pan_db_path)

        # let's learn a few things form the pan-db.
        pan_db = dbops.PanDatabase(self.pan_db_path)
        self.gene_cluster_data = pan_db.db.get_some_columns_from_table(t.pan_gene_clusters_table_name, "gene_cluster_id, genome_name", as_data_frame=True)
        self.unique_genomes = self.gene_cluster_data["genome_name"].unique()
        self.num_genomes = len(self.unique_genomes)
        self.pan_project_name = pan_db.meta['project_name']
        pan_db.disconnect()

        # here we will set the output file prefix (so we can define all the output file names already)
        if self.skip_output_files:
            self.run.info("Output files", "None will be generated as per user request", mc='red')
        else:
            if self.output_file_prefix:
                self.run.info("Output file prefix", self.output_file_prefix, mc='green')
            else:
                self.run.info("Output file prefix", f"{self.pan_project_name} (automatically set by anvi'o)", mc="cyan")
                self.output_file_prefix = self.pan_project_name

        # output file paths -- whether they will be used or not
        J = lambda x: None if self.skip_output_files else os.path.join(self.output_file_prefix.rstrip('/') + '-' + x)
        self.rarefaction_curves_figure = J('rarefaction-curves.svg')
        self.rarefaction_pangenome_txt = J('rarefaction-pangenome-averages.txt')
        self.iterations_pangenome_txt = J('rarefaction-pangenome-iterations.txt')
        self.rarefaction_core_txt = J('rarefaction-core-averages.txt')
        self.iterations_core_txt = J('rarefaction-core-iterations.txt')

        # if output files will be produced, let's make sure the user has the write
        # permissions to these destinations
        if not self.skip_output_files:
            filesnpaths.is_output_file_writable(self.rarefaction_pangenome_txt)

        # some insights into what's up on the terminal
        self.run.info("Number of genomes found", self.num_genomes)
        self.run.info("Number of iterations to run", self.iterations)


    def calc_rarefaction_curve(self, target="all"):
        """Calculates rarefaction curves for all gene clusters (target='all') or core gene clusters (target='core')."""

        if target not in ['all', 'core']:
            raise ConfigError("The target variable must either be set to 'all', or 'core'. It is not negotiable!")

        results = []
        iteration_results = []  # we store individual individual iteration values to be able to plot them later

        for n in range(1, self.num_genomes + 1):
            sampled_cluster_counts = []

            for i in range(self.iterations):
                self.progress.increment()
                self.progress.update(f"Processing {target} gene clusters in {n} genomes of {self.num_genomes} total at teration {i + 1} of {self.iterations}")
                sampled_genomes = np.random.choice(self.unique_genomes, n, replace=False)
                sampled_data = self.gene_cluster_data[self.gene_cluster_data["genome_name"].isin(sampled_genomes)]

                cluster_counts = sampled_data.groupby("gene_cluster_id")["genome_name"].nunique()

                if target == "core":
                    count = sum(cluster_counts == n)
                else:
                    count = cluster_counts.shape[0]

                sampled_cluster_counts.append(count)
                iteration_results.append([n, count])  # Store each iteration result

            results.append([n, np.mean(sampled_cluster_counts), np.std(sampled_cluster_counts)])

        df_summary = pd.DataFrame(results, columns=["num_genomes", "avg_num_gene_clusters", "standard_deviation"])
        df_iterations = pd.DataFrame(iteration_results, columns=["num_genomes", "GeneClusters"])

        return df_summary, df_iterations


    def heap_law(self, x, k, alpha):
        """Heaps' Law function: V(N) = K * N^alpha"""

        return k * np.power(x, alpha)


    def fit_heaps_law(self, rarefaction_data):
        """Fits Heaps' Law parameters to rarefaction data."""

        x_data = rarefaction_data["num_genomes"]
        y_data = rarefaction_data["avg_num_gene_clusters"]

        # Fit using non-linear least squares
        popt, _ = curve_fit(self.heap_law, x_data, y_data, p0=[1, 0.5])
        k, alpha = popt
        self.run.info("Heaps' Law parameters estimated", f"K={k:.4f}, alpha={alpha:.4f}")

        return k, alpha


    def store_results_as_txt(self):
        """Generate text files for GC gain averages per num genome and each iteration of sampling"""

        self.rarefaction_pangenome.to_csv(self.rarefaction_pangenome_txt, sep ='\t', index=False)
        self.iterations_pangenome.to_csv(self.iterations_pangenome_txt, sep ='\t', index=False)
        self.rarefaction_core.to_csv(self.rarefaction_core_txt, sep ='\t', index=False)
        self.iterations_core.to_csv(self.iterations_core_txt, sep ='\t', index=False)


    def store_results_as_svg(self):
        """Stores a nice visualization of the rarefaction curves"""

        # Generate fitted values for plotting
        x_fit = np.linspace(1, self.num_genomes, 100)
        y_fit = self.heap_law(x_fit, self.k, self.alpha)

        plt.figure(figsize=(10, 6))

        # Plot individual iteration points with transparency
        sns.scatterplot(x="num_genomes", y="GeneClusters", data=self.iterations_pangenome, color="blue", alpha=0.05)
        sns.lineplot(x="num_genomes", y="avg_num_gene_clusters", data=self.rarefaction_pangenome, color="blue", label="All gene clusters")

        sns.scatterplot(x="num_genomes", y="GeneClusters", data=self.iterations_core, color="red", alpha=0.05)
        sns.lineplot(x="num_genomes", y="avg_num_gene_clusters", data=self.rarefaction_core, color="red", label="Core gene clusters")

        # Overlay Heaps’ Law fit
        plt.plot(x_fit, y_fit, color="green", linestyle="dashed", label=f"Heaps’ Law Fit (K={self.k:.2f}, α={self.alpha:.2f})")

        plt.title(f"Rarefaction Curves with Heaps' Law Fit (with {self.iterations} iterations)")
        plt.xlabel("Number of Genomes")
        plt.ylabel("Number of Gene Clusters")
        plt.legend()
        plt.grid()
        plt.savefig(self.rarefaction_curves_figure)
        plt.close()


    def store_output_files(self):
        if self.skip_output_files:
            # no output files for you
            return

        self.progress.new("Storing results")
        self.progress.update('.. as text')
        self.store_results_as_txt()
        self.progress.update('.. as SVG')
        self.store_results_as_svg()
        self.progress.end()

        self.run.warning(None, header="OUTPUT FILES", lc="cyan")
        self.run.info("Rarefaction curves", self.rarefaction_curves_figure)
        self.run.info("GC gain per genome for core (averages)", self.rarefaction_core_txt)
        self.run.info("GC gain per genome for core (each iteration)", self.iterations_core_txt)
        self.run.info("GC gain per genome for all (averages)", self.rarefaction_pangenome_txt)
        self.run.info("GC gain per genome for all (each iteration)", self.iterations_pangenome_txt)


    def process(self):
        """Calculates rarefaction curves, plots the results into self.output_file, and returns K and alpha for Heaps' Law fit."""

        # get all the data needed to calculate Heaps' Law fit and visualize things
        self.progress.new("Calculating Rarefaction Curves", progress_total_items=(self.num_genomes * self.iterations * 2))
        self.progress.update('...')
        self.rarefaction_pangenome, self.iterations_pangenome = self.calc_rarefaction_curve(target="all")
        self.rarefaction_core, self.iterations_core = self.calc_rarefaction_curve(target="core")
        self.progress.end()

        # Fit Heaps' Law
        self.k, self.alpha = self.fit_heaps_law(self.rarefaction_pangenome)

        # store output files
        self.store_output_files()

        return (self.k, self.alpha)


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


# ANCHOR - PangenomeGraph
class PangenomeGraph():
    """The major backend class to create, solve and layout pangenome graphs in anvi'o.
    Please read through individual steps for more in-depth explanaitions on algorithms
    and data structures.
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # ANVI'O INPUTS
        self.pan_db_path = A('pan_db')
        self.genomes_storage = A('genomes_storage')
        self.external_genomes_txt = A('external_genomes')
        self.pan_graph_yaml = A('pan_graph_yaml')
        self.project_name = A('project_name')

        # learn the project name from the pan-db if the user did not
        # provide another
        if self.pan_db_path:
            pan_db = dbops.PanDatabase(self.pan_db_path)
            self.gene_alignments_computed = pan_db.meta['gene_alignments_computed']
        else:
            self.gene_alignments_computed = False

        if not self.project_name:
            if pan_db:
                self.project_name = pan_db.meta['project_name']
            else:
                raise ConfigError("You need to explicitly define a `--project-name` for this "
                                  "run (anvi'o would have figured it out for you, but you don't "
                                  "even have a pan-db).")

        if self.pan_graph_yaml:
            with open(self.pan_graph_yaml) as file:
                self.yaml_file = yaml.safe_load(file)
        else:
            self.yaml_file = {}

        if A('genome_names'):
            self.genome_names = A('genome_names').split(',')
        elif self.external_genomes_txt:
            self.genome_names = pd.read_csv(self.external_genomes_txt, header=0, sep="\t")['name'].to_list()
        elif self.pan_graph_yaml:
            self.genome_names = list(self.yaml_file.keys())
        else:
            self.genome_names = []

        # ANVI'O OUTPUTS
        self.output_dir = A('output_dir')
        self.pan_graph_db_path = os.path.join(self.output_dir, self.project_name + '-PAN-GRAPH.db')
        self.output_synteny_gene_cluster_dendrogram = A('output_synteny_gene_cluster_dendrogram')
        self.output_hybrid_genome = A('output_hybrid_genome')
        self.circularize = A('circularize')
        self.just_do_it = A('just_do_it')

        # ANVI'O FLAGS
        self.start_node = []
        self.start_gene = A('start_gene')
        self.start_column = A('start_column')
        self.min_contig_chain = A('min_contig_chain')

        self.n = A('n')
        self.alpha = A('alpha')
        self.beta = A('beta')
        self.gamma = A('gamma')
        self.delta = A('delta')
        self.min_k = A('min_k')
        self.inversion_aware = A('inversion_aware')

        self.max_edge_length_filter = A('max_edge_length_filter')
        self.gene_cluster_grouping_threshold = A('gene_cluster_grouping_threshold')
        self.groupcompress = A('grouping_compression')
        self.priority_genome = A('priority_genome')
        self.load_state = A('load_state')
        self.import_values = A('import_values').split(',') if A('import_values') else []

        # STANDARD CLASS VARIABLES
        self.version = anvio.__pangraph__version__
        self.functional_annotation_sources_available = DBInfo(self.genomes_storage, expecting='genomestorage').get_functional_annotation_sources() if self.genomes_storage else []
        self.seed = None
        self.pangenome_graph = PangenomeGraphManager()
        self.pangenome_data_df = pd.DataFrame()

        self.newick = ''
        self.meta = {}
        self.bins = {}
        self.states = {}


    def summarize_pangenome_graph(self):
        self.run.warning(None, header="Generate pangenome graph summary tables", lc="green")

        node_positions, edge_positions, node_groups = TopologicalLayout().run_synteny_layout_algorithm(F=self.pangenome_graph.graph)
        self.pangenome_graph.set_node_positions(node_positions)
        region_sides_df, nodes_df, gene_calls_df = self.pangenome_graph.summarize()
        additional_info = pd.merge(region_sides_df.reset_index(drop=False), nodes_df.reset_index(drop=False), how="left", on="region_id").set_index('syn_cluster')

        for index, line in additional_info.iterrows():
            self.pangenome_graph.graph.nodes()[index]['layer']['backbone'] = True if line["motif"] == "BB" else False

        gene_calls_df.to_csv(os.path.join(self.output_dir, 'gene_calls_df.tsv'), sep='\t')
        region_sides_df.to_csv(os.path.join(self.output_dir, 'region_sides_df.tsv'), sep='\t')
        nodes_df.to_csv(os.path.join(self.output_dir, 'nodes_df.tsv'), sep='\t')

        X = len(set(nodes_df.reset_index().query('region_id == -1')['x'].tolist())) / len(set(nodes_df.reset_index()['x'].tolist()))
        A = len(set(nodes_df.reset_index().query('region_id == -1')['syn_cluster'].tolist())) / len(set(nodes_df.reset_index()['syn_cluster'].tolist()))
        complexity_value = 1 - (X + A) / 2

        with open(os.path.join(self.output_dir, 'complexity_value.txt'), 'w') as file:
            file.write(str(complexity_value))

        self.run.info_single(f"Pangenome graph complexity is {round(complexity_value, 3)}.")
        self.run.info_single(f"Exported gene calls table to {os.path.join(self.output_dir, 'gene_calls_df.tsv')}.")
        self.run.info_single(f"Exported region table to {os.path.join(self.output_dir, 'region_sides_df.tsv')}.")
        self.run.info_single(f"Exported nodes table to {os.path.join(self.output_dir, 'nodes_df.tsv')}.")
        self.run.info_single("Done.")


    def layout_pangenome_graph(self):
        self.run.warning(None, header="Running maximum force layout algorithm", lc="green")

        node_positions, edge_positions, node_groups = TopologicalLayout().run_synteny_layout_algorithm(
            F=self.pangenome_graph.graph,
            gene_cluster_grouping_threshold=self.gene_cluster_grouping_threshold,
            groupcompress=self.groupcompress,
        )

        x_max = max([x for x,y in node_positions.values()])
        y_max = max([y for x,y in node_positions.values()])
        self.run.info_single(f"Pangenome graph length = {x_max}.")
        self.run.info_single(f"Pangenome graph height = {y_max}.")
        if y_max <= len(self.genome_names) * 2 :
            self.run.info_single("Pangenome graph height and length, looks fine :)")
        else:
            self.run.info_single("A high amount of layering might affect the readability.")

        self.pangenome_graph.set_edge_positions(edge_positions)
        self.pangenome_graph.set_node_positions(node_positions)
        self.pangenome_graph.set_node_groups(node_groups)

        length_info = {
            '<10': 0,
            '<50': 0,
            '<100': 0,
            '<500': 0,
            '>500': 0,
        }

        for edge_i, edge_j in self.pangenome_graph.graph.edges():
            length = self.pangenome_graph.graph[edge_i][edge_j]['length']
            for length_filter in length_info.keys():
                if eval(str(length) + length_filter) == True:
                    length_info[length_filter] += 1
                    break

        self.run.info_single("Summary of edge length distribution:")
        for length_filter, number in length_info.items():
            self.run.info_single(f"{length_filter}: {number} edge(s)")

        long_edges = self.pangenome_graph.cut_edges(self.max_edge_length_filter)

        self.run.info_single(f"Removed {len(long_edges)} edges due to user defined length cutoff.")
        self.run.info_single("Done.")


    def print_settings(self):
        """Print settings"""
        self.run.warning(None, header="Loading settings to create anvi'o pangenome graph from scratch", lc="green")
        self.run.info_single(f"Circularize genomes: {self.circularize}.")
        self.run.info_single(f"Graph starting gene: {self.start_gene}.")
        self.run.info_single(f"Functional annotation source for starting gene: {self.start_column}.")
        self.run.info_single(f"Minimum number of synteny clusters in contig: {self.min_contig_chain}.")
        self.run.info_single(f"Global context comparison window size: {self.n}.")
        self.run.info_single(f"Global context comparison treshold value alpha: {self.alpha}.")
        self.run.info_single(f"Local context comparison gap to gene value beta: {self.beta}.")
        self.run.info_single(f"Local context comparison gap to gap value gamma: {self.gamma}.")
        self.run.info_single(f"General context comparison max treshold dela: {self.delta}.")
        self.run.info_single(f"Synteny gene cluster min k: {self.min_k}.")
        self.run.info_single(f"Higher inversion awareness: {self.inversion_aware}.")
        self.run.info_single(f"Priority genome: {self.priority_genome}.")
        self.run.info_single("The remaining settings are not affecting the pangenome graph core creation.")
        self.run.info_single("Done.")


    def get_pangenome_graph_from_scratch(self):
        """Populates `self.pangenome_graph` from a pan-db (or a YAML file)"""

        SynGC = SyntenyGeneCluster(self.args)

        self.pangenome_data_df = SynGC.get_data_from_YAML() if self.pan_graph_yaml else SynGC.get_data_from_pan_db()

        if self.start_gene and self.start_column:
            if self.start_column in self.pangenome_data_df.columns:
                start_syn_cluster = self.pangenome_data_df[self.pangenome_data_df[self.start_column].str.contains(self.start_gene)]['syn_cluster'].to_list()
                start_syn_type = self.pangenome_data_df[self.pangenome_data_df[self.start_column].str.contains(self.start_gene)]['syn_cluster_type'].to_list()
                self.start_node += set(start_syn_cluster)

                if len(set(start_syn_cluster)) > 1:
                    self.run.info_single("There is more than one occurance of your start gene of preference, I really hope you know what you are doing.")

                if len(start_syn_cluster) != len(self.genome_names):
                    self.run.info_single("The number of genomes in the dataset does not equal the number of occurances of the start gene. Weird.")

                if any(node != 'core' for node in start_syn_type):
                    self.run.info_single("At least one occurence of a start gene is not a core synteny cluster.")
            else:
                self.run.info_single("The column were we should search for your start gene does not exist...")

        self.create_pangenome_graph()


    def process(self):
        """Main processing method for pangenome graph analysis and creation"""

        # make sure the output directory is there
        filesnpaths.gen_output_directory(self.output_dir)

        # a round of sanity check
        self.sanity_check()

        # display some settings if applicable
        self.print_settings()

        # figure out if we will get the pangenome graph from
        # a user-provided JSON file, or from sctratch
        self.get_pangenome_graph_from_scratch()

        # generate flat text file summaries for downstream
        # analyses
        self.summarize_pangenome_graph()

        # calculate the display
        self.layout_pangenome_graph()

        # if the user has not provided a tree file,
        # calculate one from the graph properties
        if not self.newick:
            self.newick = self.pangenome_graph.calculate_graph_distance(self.output_dir)

        # FIXME: not currently engaged
        if self.output_hybrid_genome:
            self.pangenome_graph.generate_hybrid_genome(self.output_dir)

        # generate pan-graph-db and populate it with information
        self.generate_pan_graph_db()


    # TODO needs more sanity checks!
    def sanity_check(self):
        pass


    def get_default_state(self):
        """Calculates a default state for the pangenome graph"""

        x_max = max([data['position'][0] for node, data in self.pangenome_graph.graph.nodes(data=True)])
        y_max = max([data['position'][1] for node, data in self.pangenome_graph.graph.nodes(data=True)])

        for i, j, data in self.pangenome_graph.graph.edges(data=True):
            if data['route']:
                for x, y in data['route']:
                    y_max = y if y > y_max else y_max

        full_radius = int(180 * (45 * x_max) / (math.pi * 270))

        tracks_radius = int((2 * full_radius / 3))
        inner = int((1 * full_radius / 3))

        tracks_layer = int(tracks_radius / (3/2 * len(self.genome_names) + (5/2)))

        inner_margin = int(tracks_layer / 2)
        backbone = int(tracks_layer / 2)
        arrow = int(tracks_layer / 2)
        search = int(tracks_layer / 2)

        label = int(arrow * 0.25)
        disty = int(tracks_layer / y_max)

        state = {'rearranged_color': '#8FF0A4',
                 'accessory_color': '#DC8ADD',
                 'paralog_color': '#FFA348',
                 'singleton_color': '#99C1F1',
                 'core_color': '#BCBCBC',
                 'trna_color': '#F66151',
                 'layer_color': '#F5F5F5',
                 'non_back_color': '#F8E45C',
                 'back_color': '#3D70A0',
                 'flexsaturation': True,
                 'arrow': arrow,
                 'flexarrow': True,
                 'backbone': backbone,
                 'flexbackbone': True,
                 **{'flex' + layer: False for layer in self.import_values},
                 **{layer: 0 for layer in self.import_values},
                 **{'flex' + genome + 'layer': True for genome in self.genome_names},
                 **{genome + 'layer': tracks_layer for genome in self.genome_names},
                 'flextree': False,
                 'tree_length': 500,
                 'tree_offset': 100,
                 'tree_thickness': 3,
                 'distx': 45,
                 'disty': disty,
                 'size': 15,
                 'circ': 5,
                 'edge': 5,
                 'flexlinear': False,
                 'line': 5,
                 'label': label,
                 'search_hit': search,
                 'inner_margin': inner_margin,
                 'outer_margin': 0,
                 'inner': inner,
                 'angle' : 270,
                 'flexcondtr': True if self.gene_cluster_grouping_threshold != -1 else False,
                 'condtr': self.gene_cluster_grouping_threshold,
                 'flexmaxlength': True if self.max_edge_length_filter != -1 else False,
                 'maxlength': self.max_edge_length_filter,
                 'flexgroupcompress': True if self.groupcompress != 1.0 else False,
                 'groupcompress': self.groupcompress,
                 **{'flex' + genome: True for genome in self.genome_names},
                 **{genome: '#000000' for genome in self.genome_names}
        }

        return state


    def generate_pan_graph_db(self):
        """Generates an empty pan-graph-db and populates it with essential information"""

        # generate an empty pan-graph-db
        meta_values = {
            'project_name': self.project_name,
            'state': self.load_state,
            'version': self.version,
            'genomes_storage_hash': GenomeStorage(self.genomes_storage, storage_hash=None, genome_names_to_focus=self.genome_names).get_storage_hash(),
            'priority_genome': self.priority_genome,
            'genome_names': ','.join(self.genome_names),
            'gene_alignments_computed': self.gene_alignments_computed,
            'gene_function_sources': ','.join(self.functional_annotation_sources_available),
        }

        dbops.PanGraphDatabase(self.pan_graph_db_path, quiet=False).create(meta_values)

        # add a default state
        TablesForStates(self.pan_graph_db_path).store_state('default', json.dumps(self.get_default_state()))

        # populate nodes in pan-graph-db
        self.store_nodes_in_pan_graph_db()

        # populate edges in pan-graph-db
        self.store_edges_in_pan_graph_db()

        # store items additional data
        self.update_pan_graph_db_with_items_additional_data()

        # store layer orders (the newick tree computed from the graph)
        self.update_pan_graph_db_with_layer_orders()


    def update_pan_graph_db_with_layer_orders(self):
        """Adds the newick tree calculated from the graph into the pan-graph-db"""

        args = argparse.Namespace(pan_or_profile_db=self.pan_graph_db_path, target_data_table="layer_orders")
        miscdata.TableForLayerOrders(args, r=terminal.Run(verbose=False)).add({"default": {'data_type': 'newick', 'data_value': self.newick}}, skip_check_names=True)


    def update_pan_graph_db_with_items_additional_data(self):
        """Updates the pan-graph-db with additional node information"""

        data = {}
        keys = set([])

        # let's start with backbone
        nodes = dict(self.pangenome_graph.graph.nodes(data=True))
        for node in nodes:
            data[node] = nodes[node]['layer']
            keys.update(nodes[node]['layer'].keys())

        # (...)

        args = argparse.Namespace(pan_or_profile_db=self.pan_graph_db_path, target_data_table="items")
        miscdata.TableForItemAdditionalData(args, r=terminal.Run(verbose=False)).add(data, list(keys), skip_check_names=True)


    def store_nodes_in_pan_graph_db(self):
        nodes = dict(self.pangenome_graph.graph.nodes(data=True))

        self.progress.new('Storing syn gene cluster nodes in pan-graph-db')
        self.progress.update('...')

        table_for_nodes = TableForNodes(self.pan_graph_db_path, run=self.run, progress=self.progress)
        for node in nodes:
            node_entry = {
                'node_id': node,
                'node_type': nodes[node]['type'],
                'gene_cluster_id': nodes[node]['gene_cluster'],
                'gene_calls_json': json.dumps(nodes[node]['gene_calls']),
                'alignment_summary': json.dumps(nodes[node]['alignment'])
            }

            table_for_nodes.add(node_entry)

        self.progress.end()

        table_for_nodes.store()

        pan_graph_db = dbops.PanGraphDatabase(self.pan_graph_db_path, quiet=True)
        pan_graph_db.db.set_meta_value('num_nodes', len(nodes))
        pan_graph_db.disconnect()


    def store_edges_in_pan_graph_db(self):
        """"""

        edges = {'E_' + str(edge_id).zfill(8): {'source': edge_i, 'target': edge_j, **data} for edge_id, (edge_i, edge_j, data) in enumerate(self.pangenome_graph.graph.edges(data=True))}

        self.progress.new('Storing syn gene cluster edges in pan-graph-db')
        self.progress.update('...')

        table_for_edges = TableForEdges(self.pan_graph_db_path, run=self.run, progress=self.progress)
        for edge in edges:
            edge_entry = {'edge_id': edge,
                          'source': edges[edge]['source'],
                          'target': edges[edge]['target'],
                          'weight': edges[edge]['weight'],
                          'directions': json.dumps(edges[edge]['directions']),
                          'route': json.dumps(edges[edge]['route'])}
            table_for_edges.add(edge_entry)

        self.progress.end()

        table_for_edges.store()

        pan_graph_db = dbops.PanGraphDatabase(self.pan_graph_db_path, quiet=True)
        pan_graph_db.db.set_meta_value('num_edges', len(edges))
        pan_graph_db.disconnect()


    def create_pangenome_graph(self):
        """
        Here the pangenome graph is created. The initial graph is a cyclic graph holding the
        synteny information of the genes present in the genomes. The pangenome graph is a
        non-cyclic graph featured with (x,y) positions and grouped nodes for a layout
        representation.

        Parameters
        ==========
        None

        Returns
        =======
        self.pangenome_graph: PangenomeGraphManager Object
        """

        pan_db = dbops.PanSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        pan_db.init_gene_clusters()
        
        # 2. step: Fill self.pangenome_graph with nodes and edges based on the synteny data
        self.run.warning(None, header="Initalizing pangenome graph and filling with nodes and edges.", lc="green")
        factor = 0.00000001 / 2
        decisison_making = {}

        # TODO I guess it would be better to only use this on sigle edges not in general, like find nodes with equal edges and then give a smallest numpy number bonus.
        # Unfortunately this part here is very arbitiary but necessary. Find better way later!
        # INCLUDE CORE AND NOT CORE INFORMATION HERE!
        for genome in self.genome_names:
            decisison_making[genome] = factor
            factor /= 2

        add_layers = False
        if self.import_values:
            if set(self.import_values).issubset(self.pangenome_data_df.columns) and set(self.pangenome_data_df[self.import_values].dtypes.astype(str).values.tolist()).issubset(['int64', 'float64']):
                self.run.info_single(f"Entries {', '.join(self.import_values)} will be added as optional layers.")
                add_layers = True

        number_gene_calls = {}
        for genome, genome_group in self.pangenome_data_df.groupby(["genome"]):
            extra_connections = []

            for contig, group in genome_group.groupby(["contig"]):
                group.reset_index(drop=False, inplace=True)
                group.sort_values(["position"], axis=0, ascending=True, inplace=True)

                syn_cluster_tuples = list(map(tuple, group[['index', 'syn_cluster', 'syn_cluster_type', 'gene_cluster', 'gene_caller_id']].values.tolist()))

                group.set_index('index', inplace=True)

                # TODO just multiplying by 100 is a bit arbitrary as well...
                if genome == self.priority_genome:
                    add_weight = 1.0 * 100
                else:
                    add_weight = 0

                if len(syn_cluster_tuples) >= self.min_contig_chain:

                    extra_connections += [syn_cluster_tuples[0], syn_cluster_tuples[-1]]
                    add_weight += decisison_making[genome]

                    if len(syn_cluster_tuples) > 1:
                        syn_cluster_tuple_pairs = map(tuple, zip(syn_cluster_tuples, syn_cluster_tuples[1:]))
                        for syn_cluster_tuple_pair in syn_cluster_tuple_pairs:
                            index_i, syn_cluster_i, syn_cluster_type_i, gene_cluster_i, gene_caller_id_i = syn_cluster_tuple_pair[0]
                            index_j, syn_cluster_j, syn_cluster_type_j, gene_cluster_j, gene_caller_id_j = syn_cluster_tuple_pair[1]

                            node_attributes_i = {
                                'gene_cluster': gene_cluster_i,
                                'gene_calls': {genome: gene_caller_id_i},
                                'type': syn_cluster_type_i,
                                'layer': group[self.import_values].loc[index_i].to_dict() if add_layers else {}
                            }

                            node_attributes_j = {
                                'gene_cluster': gene_cluster_j,
                                'gene_calls': {genome: gene_caller_id_j},
                                'type': syn_cluster_type_j,
                                'layer': group[self.import_values].loc[index_j].to_dict() if add_layers else {}
                            }

                            edge_attributes = {
                                'weight': 1.0 + add_weight,
                                'directions': {genome: 'R'}
                            }

                            self.pangenome_graph.add_node_to_graph(syn_cluster_i, node_attributes_i)
                            self.pangenome_graph.add_node_to_graph(syn_cluster_j, node_attributes_j)
                            self.pangenome_graph.add_edge_to_graph(syn_cluster_i, syn_cluster_j, edge_attributes)

                    else:
                        index_i, syn_cluster_i, syn_cluster_type_i, gene_cluster_i, gene_caller_id_i = syn_cluster_tuples[0]

                        node_attributes_i = {
                            'gene_cluster': gene_cluster_i,
                            'gene_calls':{genome: gene_caller_id_i},
                            'type': syn_cluster_type_i,
                            'layer': group[self.import_values].loc[index_i].to_dict() if add_layers else {}
                        }

                        self.pangenome_graph.add_node_to_graph(syn_cluster_i, node_attributes_i)

                    if genome not in number_gene_calls:
                        number_gene_calls[genome] = len(syn_cluster_tuples)
                    else:
                        number_gene_calls[genome] += len(syn_cluster_tuples)
                else:
                    self.run.info_single(f"Skipped {contig} due to small contig size.")

            if self.circularize:

                extra_connections = extra_connections[1:] + [extra_connections[0]]
                num = 0
                while num < len(extra_connections):

                    index_i, extra_connections_syn_i, extra_connections_type_i, extra_connections_gc_i, extra_connections_id_i = extra_connections[num]
                    index_j, extra_connections_syn_j, extra_connections_type_j, extra_connections_gc_j, extra_connections_id_j = extra_connections[num+1]

                    num += 2

                    edge_attributes = {
                        'weight': 1.0 + add_weight,
                        'directions': {genome: 'R'}
                    }

                    if extra_connections_syn_i != extra_connections_syn_j:
                        self.pangenome_graph.add_edge_to_graph(extra_connections_syn_i, extra_connections_syn_j, edge_attributes)

        for genome in number_gene_calls:
            self.run.info_single(f"Added {number_gene_calls[genome]} gene calls from {genome}.")

        num_syn_cluster = len(self.pangenome_data_df['syn_cluster'].unique())
        num_graph_nodes = len(self.pangenome_graph.graph.nodes())

        self.run.info_single(f"Added {num_graph_nodes} nodes and {len(self.pangenome_graph.graph.edges())} edges to pangenome graph.")
        self.run.info_single("Done.")

        if num_syn_cluster != num_graph_nodes:
            self.run.info_single(f"It looks like {abs(num_syn_cluster - num_graph_nodes)} nodes were not added to the graph. Proceeding from "
                              f"here might cause some trouble. Maybe you also set a minimum contig size, in this case great job user, go ahead :)")

        # 3. step: Check connectivity of the graph
        self.pangenome_graph.run_connectivity_check()

        # 4. step: Find edges to reverse to create maxmimum directed force
        self.run.warning(None, header="Running maximum force calculation on pangenome graph", lc="green")

        selfloops = list(nx.selfloop_edges(self.pangenome_graph.graph))
        self.run.info_single(f"Found and removed {len(selfloops)} selfloop edge(s)")
        self.pangenome_graph.graph.remove_edges_from(selfloops)

        # changed_edges, removed_nodes, removed_edges = DirectedForce().return_optimum_complexity()
        changed_edges = DirectedForce().return_optimum_complexity(
            self.pangenome_graph.graph,
            max_iterations=1,
            start_node=self.start_node
        )

        self.pangenome_graph.reverse_edges(changed_edges)
        # self.pangenome_graph.graph.remove_edges_from(removed_edges)
        # self.pangenome_graph.graph.remove_nodes_from(removed_nodes)

        # self.pangenome_data_df.drop(self.pangenome_data_df.loc[self.pangenome_data_df['syn_cluster'].isin(removed_nodes)].index, inplace=True)
        # self.run.info_single(f"The pangenome graph is now a connected non-cyclic graph.")
        if len(changed_edges) == 0:
            self.run.info_single("This does look weird good but maybe you have a perfect dataset without the need of any edge reversal.")
        self.run.info_single("Done.")

        for node, data in self.pangenome_graph.graph.nodes(data=True):
            node_alignment_summaries = {}
            node_alignment_lengths = []
            node_alignments = {}
            for genome_name, gene_caller_id in data['gene_calls'].items():
                if self.gene_alignments_computed:
                    genome_alignments = pan_db.gene_clusters_gene_alignments[genome_name]
                    if gene_caller_id in genome_alignments:
                        alignment_summary = genome_alignments[gene_caller_id]
                        node_alignment_summaries[genome_name] = alignment_summary

                        alignment_summary_list = alignment_summary.split('|')
                        start = alignment_summary_list[0]
                        summary_code = list(map(int, alignment_summary_list[1:]))

                        if start == '-':
                            sequence = sum(summary_code[1::2]) * 'N'
                        else:
                            sequence = sum(summary_code[0::2]) * 'N'

                        alignment = utils.restore_alignment(sequence, alignment_summary)
                        node_alignment_lengths += [len(alignment)]
                    else:
                        alignment_summary = ''
                        node_alignment_summaries[genome_name] = alignment_summary

                        sequence = ''
                        alignment = ''
                        node_alignment_lengths += [0]

                    node_alignments[genome_name] = alignment

            if len(set(node_alignment_lengths)) != 1:
                raise ConfigError("Your alignments have a different length? Oh boy that's not something we like.")

            cleaned_alignments = {genome_name: '' for genome_name in node_alignments.keys()}

            for i in range(0, node_alignment_lengths[0]):
                summary_code_list = [alignment[i] for genome_name, alignment in node_alignments.items()]
                if not all(a == '-' for a in summary_code_list):
                    for genome_name, alignment in node_alignments.items():
                        cleaned_alignments[genome_name] += alignment[i]

            data['alignment'] = {genome_name: utils.summarize_alignment(alignment) if alignment else '' for genome_name, alignment in cleaned_alignments.items()}