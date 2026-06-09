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

from itertools import chain, combinations
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

from anvio.splitter import LocusSplitter
from anvio.genomestorage import GenomeStorage
from anvio.tables.views import TablesForViews
from anvio.tables.states import TablesForStates
from anvio.errors import ConfigError, FilesNPathsError
from anvio.genomedescriptions import GenomeDescriptions
from anvio.tables.geneclusters import TableForGeneClusters
from anvio.tables.genefunctions import TableForGeneFunctions
from anvio.tables.pangraphdata import TableForNodes, TableForEdges, TableForRegions, TableForGenomeDistances

from anvio import panaai
from anvio.topologicallayout import TopologicalLayout
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
P = terminal.pluralize

additional_param_sets_for_sequence_search = {'diamond'   : '--masking 0',
                                             'ncbi_blast': ''}


class PangenomeGraphSubGraph:
    """Takes in a pangenome graph, and exports the genomic loci between two nodes in it as contigs databases.

        >>> import argparse
        >>> args = argparse.Namespace(pan_graph_db="PATH/TO/PAN-GRAPH.db", graph_nodes="NODE_X,NODE_Y", output_dir="OUTPUT_DIR")
        >>> subgraph = PangenomeGraphSubGraph(args)
        >>> subgraph.export()

    A client of this class is the program `anvi-export-pan-subgraph`
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.pan_graph_db_path = A('pan_graph_db')
        self.graph_nodes = A('graph_nodes').split(',') if A('graph_nodes') else None
        self.output_dir = A('output_dir')
        self.external_genomes_file_path = A('external_genomes')

        if not self.graph_nodes:
            raise ConfigError("This program is useless without the `--graph-nodes` parameter :/")

        if not self.pan_graph_db_path:
            raise ConfigError("Please send a pangenome graph database")

        if len(self.graph_nodes) != 2:
            raise ConfigError(f"The `--graph-nodes` parameter must be set to two node names that are separated by a comma :/ "
                              f"Your parameter, '{A('graph_nodes')}', does not really comply with that.")

        utils.is_pan_graph_db(self.pan_graph_db_path)

        filesnpaths.check_output_directory(self.output_dir)


    def export(self):
        """Export the genomic loci between self.node_names from every genome involved in pangenome graph"""

        # get an instance of PanGraphSuperclass
        pangraph = dbops.PanGraphSuperclass(self.args)
        pangraph.init_synteny_gene_clusters()

        missing_nodes = [node for node in self.graph_nodes if node not in pangraph.synteny_gene_cluster_names]
        if len(missing_nodes) == 2:
            raise ConfigError(f"Neither of the nodes you requested, '{self.graph_nodes[0]}' and '{self.graph_nodes[1]}', are "
                              f"found in the pangenome graph database (congratulations) :(")
        elif len(missing_nodes) == 1:
            raise ConfigError(f"One of the nodes you requested, '{missing_nodes[0]}', is not found in the pangenome graph database :(")
        else:
            pass

        # learn the genome names from the external genomes file and make sure the genome names in
        # the pangenome graph db are consistent with those.
        g = GenomeDescriptions(self.args, run=terminal.Run(verbose=False), progress=self.progress)
        g.load_genomes_descriptions(skip_functions=True, init=False)

        missing_genomes = [genome_name for genome_name in pangraph.genome_names if genome_name not in g.genomes]
        if len(missing_genomes):
            raise ConfigError(f"The following genomes are found in the pangenome graph database, but not in the "
                              f"external genomes file: {', '.join(missing_genomes)}. So anvi'o is confuse "
                              f"and not sure how to continue :(")

        self.run.info('Pangenome graph database', pangraph.p_meta['project_name'])
        self.run.info("Pan graph database", self.pan_graph_db_path)
        self.run.info("Nodes to export", ', '.join(self.graph_nodes))
        self.run.info("Loci", '')

        d = {}
        for genome_name in pangraph.genome_names:
            d[genome_name] = []

            for graph_node in self.graph_nodes:
                if not len(pangraph.synteny_gene_clusters[graph_node][genome_name]):
                    raise ConfigError(f"The curent implementation of this tool requires the graph nodes of interest to "
                                      f"correspond to SynGCs that are present in all genomes (so we can select what is "
                                      f"between them in each genome easily). Unfortunately, the graph node '{graph_node}' "
                                      f"does not have any genes from the genome '{genome_name}'.")

                d[genome_name].append(pangraph.synteny_gene_clusters[graph_node][genome_name][0])

            d[genome_name] = sorted(d[genome_name])

            # we know which genes we are interested in for the genome, let's report it to the user
            # before moving on
            self.run.info_single(f"{d[genome_name][0]} to {d[genome_name][1]} ({P('gene', d[genome_name][1] - d[genome_name][0])}) for {genome_name}", level=2)

        # at this stage we have everything we need stored in `d` and `g` to start exporting loci
        # from each contigs database. let's start by generating the output directory
        filesnpaths.gen_output_directory(self.output_dir, delete_if_exists=True)

        self.progress.new("Exporting", progress_total_items=len(pangraph.genome_names))
        for genome_name in pangraph.genome_names:
            progress.update(f"Working on {genome_name} ...", increment=1)

            contigs_db_path = g.genomes[genome_name]['contigs_db_path']
            first_gene_call = d[genome_name][0]
            second_gene_call = d[genome_name][1]

            # build the args for LocusSplitter
            locus_args = argparse.Namespace(contigs_db=contigs_db_path,
                                            gene_caller_ids=f"{first_gene_call},{second_gene_call}",
                                            flank_mode=True,
                                            output_dir=self.output_dir,
                                            output_file_prefix=genome_name,
                                            delimiter=',',
                                            never_reverse_complement=True,
                                            include_fasta_output=False)

            # let's go
            locus_splitter = LocusSplitter(locus_args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
            locus_splitter.process()

        self.progress.end()


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
        self.project_name = A('project_name')
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

        # next, we figure out where to keep the intermediate data files
        user_pan_db_path = A('output_file') or A('pan_db')
        if user_pan_db_path:
            if not user_pan_db_path.endswith('-PAN.db'):
                raise ConfigError("Sorry. The output file names for anvi'o pan-db artifats must end with '-PAN.db'. No exceptions, no exclusions, "
                                  "no creative interpretations. Anvi'o: freedom in data analyses, tyranny in output file suffixes ✊")
            self.pan_db_path = user_pan_db_path
        else:
            # if the user did not specify an output directory for the pan-db, just put it next to the genomes storage
            # file with a name that includes the project name.
            self.pan_db_path = os.path.join(os.path.dirname(self.genomes_storage_path), self.project_name + '-PAN.db')

        # only to keep intermediate data files
        self.intermediate_data_dir = None
        self.remove_intermediate_data_dir_at_the_end = False
        if A('output_dir') or A('intermediate_data_dir'):
            self.intermediate_data_dir = A('output_dir') or A('intermediate_data_dir')
        else:
            self.intermediate_data_dir = filesnpaths.get_temp_directory_path(just_the_path=True)
            self.remove_intermediate_data_dir_at_the_end = True
        self.intermediate_data_dir = os.path.abspath(self.intermediate_data_dir)

        # when it is time to organize gene_clusters
        self.linkage = A('linkage') or constants.linkage_method_default
        self.distance = A('distance') or constants.distance_metric_default

        self.genomes = None
        self.log_file_path = None

        # to be filled during init:
        self.amino_acid_sequences_dict = {}
        self.view_data = {}
        self.view_data_presence_absence = {}
        self.additional_view_data = {}
        self.aligner = None

        # we don't know what we are about
        self.description = None


    def cleanup(self):
        self.run.quit()

        self.run.warning(None, header="CLEANUP (OR LACKTHEREOF)", lc="cyan")

        if self.remove_intermediate_data_dir_at_the_end:
            if anvio.DEBUG:
                self.run.info_single(f"The intermediate data at {self.intermediate_data_dir} is kept because of the `--debug` flag. "
                                     f"You can inspect these data, or re-use them with the `--intermediate-data-dir` parameter with "
                                     f"`anvi-pan-genome`.", level=0, nl_after=2, mc='cyan')
            elif self.intermediate_data_dir and os.path.exists(self.intermediate_data_dir):
                self.run.info_single(f"'{self.intermediate_data_dir}', the intermediate data directory anvi'o used for this analysis, "
                                     f"is now being cleaned up. But you *could* keep it for any reason, such as wanting to debug anvi'o, "
                                     f"or to re-run `anvi-pan-genome` with different parameters using the same search results by simply "
                                     f"defining an explicit output directory path for intermediate data files using the very aptly named "
                                     f"parameter `--intermediate-data-dir`. All is good, but just FYI.", level=0, nl_after=2, mc='cyan')
                import shutil
                shutil.rmtree(self.intermediate_data_dir)
        else:
            self.run.info_single(f"No cleanup! The intermediate data at {self.intermediate_data_dir} is available for you to re-run the "
                                 f"`anvi-pan-genome` with the `--intermediate-data-dir` parameter.", level=0, nl_after=2, mc='cyan')


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
        output_file_path = os.path.join(self.intermediate_data_dir, file_name)

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
        # deal with the intermediate data directory
        try:
            filesnpaths.is_file_exists(self.intermediate_data_dir)
        except FilesNPathsError:
            filesnpaths.gen_output_directory(self.intermediate_data_dir, delete_if_exists=self.overwrite_output_destinations)

        filesnpaths.is_output_dir_writable(self.intermediate_data_dir)

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
                             "failed to report self search results: %s."
                                                    % (search_tool, len(ids_without_self_search), len(all_ids),
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
                                 "clustering to be done anyway, please see the flag `--enforce-hierarchical-clustering`."
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

        # write gene cluster stats and SCG data as items additional data in named groups.
        # we say skip_check_names=True, simply because there is no gene_clusters table has not been
        # generated yet, but the check names functionality in dbops looks for the gene clsuters table to
        # be certain. it is not a big deal here, since we absoluely know what gene cluster names we are
        # working with.
        stats_args = argparse.Namespace(**{**vars(self.args), 'target_data_group': 'gene_cluster_stats'})
        stats_keys = ['num_genomes_gene_cluster_has_hits', 'num_genes_in_gene_cluster', 'max_num_paralogs']
        miscdata.TableForItemAdditionalData(stats_args, r=terminal.Run(verbose=False)).add(self.additional_view_data, stats_keys, skip_check_names=True)

        scg_args = argparse.Namespace(**{**vars(self.args), 'target_data_group': 'SCG'})
        miscdata.TableForItemAdditionalData(scg_args, r=terminal.Run(verbose=False)).add(self.additional_view_data, ['SCG'], skip_check_names=True)

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

            dbops.do_hierarchical_clustering_of_items(self.pan_db_path,
                                                      updated_clustering_configs,
                                                      database_paths={'PAN.db': self.pan_db_path},
                                                      default_clustering_config=constants.pan_default,
                                                      distance=self.distance,
                                                      linkage=self.linkage,
                                                      run=terminal.Run(verbose=False),
                                                      progress=self.progress)


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

        homogeneity_args = argparse.Namespace(**{**vars(self.args), 'target_data_group': 'homogeneity'})
        homogeneity_keys = ['functional_homogeneity_index', 'geometric_homogeneity_index', 'combined_homogeneity_index']
        miscdata.TableForItemAdditionalData(homogeneity_args, r=terminal.Run(verbose=False)).add(d, homogeneity_keys, skip_check_names=True)

        aai_args = argparse.Namespace(**{**vars(self.args), 'target_data_group': 'AAI'})
        aai_keys = ['AAI_min', 'AAI_max', 'AAI_avg']
        miscdata.TableForItemAdditionalData(aai_args, r=terminal.Run(verbose=False)).add(d, aai_keys, skip_check_names=True)


    def populate_layers_additional_data_and_orders(self):
        self.progress.new('Layers additional data and orders')
        self.progress.update('Copmputing the hierarchical clustering of the (transposed) view data')

        layer_orders_data_dict = {}
        for clustering_tuple in [('gene_cluster presence absence', self.view_data_presence_absence), ('gene_cluster frequencies', self.view_data)]:
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
            raise ConfigError("There must be at least two genomes for this workflow to work. You have like '%d' of them :/"
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
                              "being inconsistent. Please make up your mind, and come back as the explicit person you are"
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

        # deal with the intermediate data output directory
        self.cleanup()

        # done
        self.run.log_file_path = None
        self.run.info("The new pan-db", self.pan_db_path, nl_after=1)

        self.run.info_single(f"Your pangenome is ready with a total of {pp(len(gene_clusters_dict))} gene clusters across "
                             f"{len(self.genomes)} genomes 🎉", mc="green", nl_after=1)


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

        # we seem to like longer messages in this class.
        self.run.width = 60

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        # ANVI'O INPUTS
        self.pan_db_path = A('pan_db')
        self.genomes_storage = A('genomes_storage')
        self.external_genomes_txt = A('external_genomes')
        self.project_name = A('project_name')

        # learn the project name from the pan-db if the user did not
        # provide another
        self.genome_names_user_focused = bool(A('genome_names'))
        if A('genome_names'):
            if filesnpaths.is_file_exists(A('genome_names'), dont_raise=True):
                self.genome_names = utils.get_column_data_from_TAB_delim_file(A('genome_names'), column_indices=[0], expected_number_of_fields=1)[0]
            else:
                self.genome_names = [g.strip() for g in A('genome_names').split(',')]
        elif self.external_genomes_txt:
            filesnpaths.is_file_tab_delimited(self.external_genomes_txt, expected_number_of_fields=2)
            self.genome_names = pd.read_csv(self.external_genomes_txt, header=0, sep="\t")['name'].to_list()
        else:
            self.genome_names = []

        if self.pan_db_path:
            if filesnpaths.is_file_exists(self.pan_db_path, dont_raise=False):
                self.pan_db = dbops.PanDatabase(self.pan_db_path)
                self.gene_alignments_computed = self.pan_db.meta['gene_alignments_computed']

                self.pan_super = dbops.PanSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
                self.pan_super.init_gene_clusters()

                # --genome-names must be a subset of the pan-db's genome list:
                # the pan-graph is computed downstream of the pan, so the focus
                # set can only narrow what the pan already knows about.
                if A('genome_names'):
                    pan_genomes = set(self.pan_super.genome_names)
                    absent = sorted(set(self.genome_names) - pan_genomes)
                    if absent:
                        head = ', '.join(absent[:3]) + ('...' if len(absent) > 3 else '')
                        raise ConfigError(f"{len(absent)} name(s) in `--genome-names` are not "
                                          f"in the pan-db's genome list (the pan-graph must be "
                                          f"a subset of the pan's genomes): {head} :/")
        else:
            self.pan_db = None
            self.pan_super = None
            self.gene_alignments_computed = False

        if self.genomes_storage:
            if filesnpaths.is_file_exists(self.genomes_storage, dont_raise=False):
                self.genomes_storage_hash = GenomeStorage(self.genomes_storage, storage_hash=None, genome_names_to_focus=self.genome_names).get_storage_hash()
        else:
            self.genomes_storage_hash = None

        if not self.project_name:
            if self.pan_db:
                self.project_name = self.pan_db.meta['project_name']
            else:
                raise ConfigError("You need to explicitly define a `--project-name` for this "
                                  "run (anvi'o would have figured it out for you, but you don't "
                                  "even have a pan-db).")

        # ANVI'O OUTPUTS
        user_pan_graph_db_path = A('output_file') or A('pan_graph_db')
        if user_pan_graph_db_path:
            self.pan_graph_db_path = user_pan_graph_db_path
        else:
            self.pan_graph_db_path = os.path.join('.', self.project_name + '-PAN-GRAPH.db')

        self.output_dir = os.path.dirname(self.pan_graph_db_path) or '.'

        self.just_do_it = A('just_do_it')

        # ANVI'O FLAGS
        self.min_contig_chain = A('min_contig_chain')
        self.no_include_non_coding_genes = bool(A('no_include_non_coding_genes'))
        self.no_remerge = bool(A('no_remerge'))
        self.remerge_max_length = A('remerge_max_length') if A('remerge_max_length') is not None else -1
        self.max_edge_length_filter = A('max_edge_length_filter')
        self.gene_cluster_grouping_threshold = A('gene_cluster_grouping_threshold')
        self.groupcompress = A('grouping_compression')
        self.component = A('component') if A('component') is not None else 0
        self.region_scope = A('region_scope') or 'global'
        self.load_state = A('load_state')
        self.import_values = A('import_values').split(',') if A('import_values') else []

        # NEMESIS / AAI ENGINE PARAMETERS
        self.locality_window = A('locality_window')
        self.min_window_completeness = A('min_window_completeness')
        self.min_line_pair_hits = A('min_line_pair_hits')
        self.orientation_tie_threshold = A('orientation_tie_threshold')
        self.min_orientation_score = A('min_orientation_score')
        self.orientation_demotion_strategy = A('orientation_demotion_strategy')
        self.ranking_components = A('ranking_components')
        self.ranking_mean = A('ranking_mean')
        self.minbit_floor = A('minbit_floor')
        self.decision_floor = A('decision_floor')
        self.support_floor = A('support_floor')
        self.decision_tie_score = A('decision_tie_score')
        self.decision_boundary_score = A('decision_boundary_score')
        self.min_ranking_score = A('min_ranking_score')
        self.fusion_top_bucket_k = A('fusion_top_bucket_k')
        self.fusion_seed = A('fusion_seed')

        # OPTIONAL DESCRIPTION (from master): stored as pan-graph-db metadata.
        description_file_path = A('description')
        if description_file_path:
            filesnpaths.is_file_plain_text(description_file_path)
            self.description = open(os.path.abspath(description_file_path), 'r').read()
        else:
            self.description = ''

        # STANDARD CLASS VARIABLES
        self.version = anvio.__pangraph__version__
        self.functional_annotation_sources_available = (DBInfo(self.genomes_storage, expecting='genomestorage').get_functional_annotation_sources() or []) if self.genomes_storage else []
        self.seed = None
        self.pangenome_graph = PangenomeGraphManager()

        self.newick = ''
        self.meta = {}
        self.bins = {}
        self.states = {}

        # Populated by create_pangenome_graph (engine bookkeeping for downstream consumers).
        self.lines = None
        self.line_names = None
        self.line_to_genome = None
        self.in_g_flip = None
        self.layers_data = {}


    def summarize_pangenome_graph(self):
        self.run.warning(None, header="GENERATING SUMMARY TABLES", lc="green")

        # bounds-check --component against the actual number of weakly connected components
        n_components = sum(1 for _ in nx.weakly_connected_components(self.pangenome_graph.graph))
        if self.component >= n_components:
            raise ConfigError(f"You asked for component {self.component}, but the graph only has "
                              f"{n_components} component(s) (valid indices are 0 to {n_components - 1}). "
                              f"Pass a smaller value to --component.")

        # Lay out and summarize EVERY component. self.component is the
        # default *display* component; it no longer gates analysis. Every
        # node gets real positions and a real backbone classification.
        self.pangenome_graph.layout_all_components(
            self.gene_cluster_grouping_threshold, self.groupcompress)

        region_sides_df, backbone_by_node = self.pangenome_graph.summarize_all_components(scope=self.region_scope)

        # Default backbone to None as a safety net (the per-component
        # summarize pass should cover every node; None only persists for
        # the rare case where it doesn't, which is more honest than 0).
        for _node, data in self.pangenome_graph.graph.nodes(data=True):
            data['layer'] = data['layer'] | {'backbone': backbone_by_node.get(_node, None)}

        # stash the region-level summary for `generate_pan_graph_db` to persist into
        # the pan_graph_regions table (translating BR/VR -> backbone/variable on the way in)
        self.region_sides_df = region_sides_df

        self.run.info_single(f"{len(region_sides_df)} region(s) summarized across {n_components} component(s); "
                             f"backbone/variable labels and region IDs attached to nodes.")

    def layout_pangenome_graph(self):
        self.run.warning(None, header="Running maximum force layout algorithm", lc="green")

        node_positions, edge_positions, node_groups = TopologicalLayout().run_synteny_layout_algorithm(
            F=self.pangenome_graph.graph,
            gene_cluster_grouping_threshold=self.gene_cluster_grouping_threshold,
            groupcompress=self.groupcompress,
            component=self.component,
        )

        x_max = max([x for x,y in node_positions.values()])
        y_max = max([y for x,y in node_positions.values()])
        self.run.info_single(f"Pangenome graph length = {x_max}.")
        self.run.info_single(f"Pangenome graph height = {y_max}.")
        if y_max <= len(self.genome_names) * 2:
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
        self.run.warning(None, header="SETTINGS", lc="green")
        if self.genome_names_user_focused:
            n = len(self.genome_names)
            preview = ', '.join(self.genome_names[:3]) + ('...' if n > 3 else '')
            self.run.info("Focused genomes (--genome-names)", f"{n} ({preview})")
        else:
            self.run.info("Focused genomes (--genome-names)", 'all (no focus)')
        self.run.info("Minimum number of genes per contig", self.min_contig_chain)
        self.run.info("Remove non-coding genes", self.no_include_non_coding_genes)
        self.run.info("Skip remerge step", self.no_remerge)
        self.run.info("Remerge LCA/LCD-asymmetry cap",
                      'disabled' if self.remerge_max_length < 0 else self.remerge_max_length)
        self.run.info("Component to layout", self.component)
        self.run.info("Region scope (BR/VR denominator)", self.region_scope)
        self.run.info("Locality window", self.locality_window)
        self.run.info("Min window completeness", self.min_window_completeness)
        self.run.info("Min line-pair hits", self.min_line_pair_hits)
        self.run.info("Orientation tie threshold", self.orientation_tie_threshold)
        self.run.info("Min orientation score", self.min_orientation_score)
        self.run.info("Orientation demotion strategy", self.orientation_demotion_strategy)
        self.run.info("Ranking components", self.ranking_components)
        self.run.info("Ranking mean", self.ranking_mean)
        self.run.info("Minbit floor", self.minbit_floor)
        self.run.info("Decision floor", self.decision_floor)
        self.run.info("Support floor", self.support_floor)
        self.run.info("Decision tie/unlabeled score", self.decision_tie_score)
        self.run.info("Decision boundary score", self.decision_boundary_score)
        self.run.info("Min ranking score", self.min_ranking_score)
        self.run.info("Fusion top-bucket k", self.fusion_top_bucket_k)
        self.run.info("Fusion seed", self.fusion_seed)
        self.run.info("Max edge length filter", self.max_edge_length_filter)
        self.run.info("Gene-cluster grouping threshold", self.gene_cluster_grouping_threshold)
        self.run.info("Grouping compression", self.groupcompress)
        self.run.info("Load state", self.load_state)
        self.run.info("Import values", ','.join(self.import_values) if self.import_values else '')
        self.run.info("Just do it", self.just_do_it)


    def process(self):
        """Main processing method for pangenome graph analysis and creation"""

        # make sure the output directory is there
        filesnpaths.gen_output_directory(self.output_dir)

        # fail fast if the user has pointed --pan-graph-db at a non-`.db` filename
        # or at a path that already exists (so we don't burn CPU before crashing
        # inside `PanGraphDatabase.touch()`)
        dbops.is_db_ok_to_create(self.pan_graph_db_path, 'pan-graph')

        # a round of sanity check
        self.sanity_check()

        # display some settings if applicable
        self.print_settings()

        # delegate the graph build to the AAI engine and mirror the result
        # into self.pangenome_graph (a PangenomeGraphManager).
        self.create_pangenome_graph()

        # collapse engine over-splits before per-node alignments are
        # built, so add_layers walks the unified gene_calls on survivors.
        self.remerge_nodes()

        self.add_layers()

        # compute region-level summaries; results are stashed on `self` so they can be
        # written to the pan_graph_regions table during db generation
        self.summarize_pangenome_graph()

        # calculate the display
        # self.layout_pangenome_graph()

        # compute the graph-based pairwise genome distance matrix; the newick
        # tree may end up empty if all distances are zero (identical genomes).
        self.distance_matrix = pd.DataFrame()
        if not self.newick:
            self.newick, self.distance_matrix, _ = self.pangenome_graph.calculate_graph_distance()

        # generate pan-graph-db and populate it with information
        self.generate_pan_graph_db()

        # Let the user know if we ended up with a zero-variation graph
        if not self.newick:
            self.run.warning("Your pangenome graph has no structure, which means the genomes you are working with "
                             "have no gene-level variation. All the files are still generated (because that's how "
                             "anvi'o rolls), but it will be a wasted effort to visualize these data since you will "
                             "not see anything worth noting :/ If you were expecting variation, please double-check "
                             "your input genomes, or take a look at the conventional pangenome using the program "
                             "`anvi-display-pan` for good measure.", header="⚠️ SILLY GENOMES WARNING ⚠️", lc="red")


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

        distx = 45
        full_radius = int(180 * (distx * x_max) / (math.pi * 270))
        tracks_radius = int((2 * full_radius / 3))
        inner = int((1 * full_radius / 3))
        tracks_layer = int(tracks_radius / (3/2 * len(self.genome_names) + (5/2)))

        inner_margin = int(tracks_layer / 2)
        backbone = int(tracks_layer / 2)
        arrow = tracks_layer
        search = int(tracks_layer / 2)

        label = int(arrow * 0.25)

        state = {
            'drawing': {
                'type': 'circular',
                'inner_radius': inner,
                'start_angle': 0,
                'end_angle': 270,
                'node_x_spacing': distx,
                'node_y_spacing': 120
            },
            'nodes': {
                'radius': 15,
                'outline_width': 5,
                'fade_by_prevalence': True,
                'type_colors': {
                    'core': '#BCBCBC',
                    'rearrangement': '#8FF0A4',
                    'accessory': '#DC8ADD',
                    'duplication': '#FFA348',
                    'singleton': '#99C1F1',
                    'rna': '#ECFA28'
                }
            },
            'edges': {
                'width': 5
            },
            'graph_layout': {
                'grouping_enabled': self.gene_cluster_grouping_threshold != -1,
                'grouping_threshold': self.gene_cluster_grouping_threshold,
                'max_edge_length_enabled': True,
                'max_edge_length': self.max_edge_length_filter if self.max_edge_length_filter != -1 else 1000,
                'group_compression_enabled': self.groupcompress != 1.0,
                'group_compression': self.groupcompress,
                'component': self.component
            },
            'layers': {
                'backbone': {
                    'visible': True,
                    'height': backbone,
                    'backbone_color': '#C4C4C4',
                    'variable_region_color': '#CC5252'
                },
                'orientation_arrow': {
                    'visible': True,
                    'height': arrow
                },
                'search': {
                    'hit_height': search
                }
            },
            'layers_tree': {
                'visible': True,
                'height': tracks_layer,
                'offset': int(inner_margin / 2),
                'line_width': 10
            },
            'genome_tracks': {
                'line_width': 5,
                'background_color': '#F5F5F5',
                'genomes': {
                    genome: {
                        'color': '#000000',
                        'show': True,
                        'track_height': tracks_layer,
                        'show_track': True
                    } for genome in self.genome_names
                }
            },
            'imported_layers': {
                layer: {
                    'visible': False,
                    'height': 0
                } for layer in self.import_values
            },
            'labels': {
                'font_size': label,
                'offset': int(inner_margin / 2),
                'position_tick_count': 20
            },
            'margins': {
                'inner': inner_margin,
                'outer': 0
            },
            'region_labels': {
                'font_size': 13,
                'min_width_px': 80,
                'distance': 2
            },
            'bins': {
                'show_labels': True,
                'label_orientation': 'natural',
                'label_font_size': 19.5,
                'ring_height': 4,
                'ring_opacity': 0.8,
                'show_edges': True,
                'edge_thickness': 4,
                'edge_color': '#FFFFFF',
                'edge_opacity': 1.0
            }
        }

        return state


    def generate_pan_graph_db(self):
        self.run.warning(None, header="STORING PAN-GRAPH-DB", lc="green")

        """Generates an empty pan-graph-db and populates it with essential information"""

        # generate an empty pan-graph-db with meta values that capture the settings and
        # and the data that went into the graph construction so we have them for downstream
        # analyeses and for the user to be able to revisit them if/when they need to.
        meta_values = {
            # identity & versioning
            'anvio_version': anvio.__version__,
            'version': self.version,
            'project_name': self.project_name,
            'description': self.description,
            # genome provenance
            'genomes_storage_hash': self.genomes_storage_hash,
            'genome_names': ','.join(self.genome_names),
            'num_genomes': len(self.genome_names),
            'gene_alignments_computed': self.gene_alignments_computed,
            'gene_function_sources': ','.join(self.functional_annotation_sources_available),
            # AAI engine parameters
            'min_contig_chain': self.min_contig_chain,
            'no_include_non_coding_genes': self.no_include_non_coding_genes,
            'no_remerge': self.no_remerge,
            'remerge_max_length': self.remerge_max_length,
            'locality_window': self.locality_window,
            'min_window_completeness': self.min_window_completeness,
            'min_line_pair_hits': self.min_line_pair_hits,
            'orientation_tie_threshold': self.orientation_tie_threshold,
            'min_orientation_score': self.min_orientation_score,
            'orientation_demotion_strategy': self.orientation_demotion_strategy,
            'ranking_components': self.ranking_components,
            'ranking_mean': self.ranking_mean,
            'minbit_floor': self.minbit_floor,
            'min_ranking_score': self.min_ranking_score,
            'fusion_top_bucket_k': self.fusion_top_bucket_k,
            'fusion_seed': self.fusion_seed,
            # layout & simplification parameters
            'max_edge_length_filter': self.max_edge_length_filter,
            'gene_cluster_grouping_threshold': self.gene_cluster_grouping_threshold,
            'grouping_compression': self.groupcompress,
            'component': self.component,
        }

        dbops.PanGraphDatabase(self.pan_graph_db_path, run=self.run, progress=self.progress, quiet=False).create(meta_values)

        # add a default state
        TablesForStates(self.pan_graph_db_path).store_state('default', json.dumps(self.get_default_state(), indent=4))

        # populate nodes in pan-graph-db
        self.store_nodes_in_pan_graph_db()

        # populate edges in pan-graph-db
        self.store_edges_in_pan_graph_db()

        # populate regions (backbone / variable, with CVS) in pan-graph-db
        self.store_regions_in_pan_graph_db()

        # populate pairwise graph-based genome distances in pan-graph-db
        self.store_genome_distances_in_pan_graph_db()

        # store items additional data
        self.update_pan_graph_db_with_items_additional_data()

        # store layer orders (the newick tree computed from the graph)
        self.update_pan_graph_db_with_layer_orders()


    def store_regions_in_pan_graph_db(self):
        """Persists `region_sides_df` (computed by `summarize_pangenome_graph`) into the
           pan_graph_regions table, translating the BR/VR shorthand into 'backbone' /
           'variable' for human-friendly downstream consumption."""

        if not hasattr(self, 'region_sides_df') or self.region_sides_df is None or self.region_sides_df.empty:
            self.run.info_single("No regions to store in the pan-graph-db (skipping).")
            return

        self.progress.new('Storing regions in pan-graph-db')
        self.progress.update('...')

        region_label = {'BR': 'backbone', 'VR': 'variable'}

        df = self.region_sides_df.reset_index(drop=False).rename(columns={'region': 'region_type'})
        df['region_type'] = df['region_type'].map(lambda x: region_label.get(x, x))

        table_for_regions = TableForRegions(self.pan_graph_db_path, run=self.run, progress=self.progress)
        for _, row in df.iterrows():
            entry = {col: row[col] for col in t.pan_graph_regions_table_structure}
            table_for_regions.add(entry)

        self.progress.end()

        table_for_regions.store()


    def store_genome_distances_in_pan_graph_db(self):
        """Persists the pairwise genome distance matrix into the pan_graph_genome_distances
           table. Both directions (A,B) and (B,A) are stored so consumers can index either
           way without reconstruction."""

        if not hasattr(self, 'distance_matrix') or self.distance_matrix is None or self.distance_matrix.empty:
            self.run.info_single("No graph-based genome distances to store in the pan-graph-db (skipping).")
            return

        self.progress.new('Storing genome distances in pan-graph-db')
        self.progress.update('...')

        table_for_distances = TableForGenomeDistances(self.pan_graph_db_path, run=self.run, progress=self.progress)
        for genome_a in self.distance_matrix.index:
            for genome_b in self.distance_matrix.columns:
                if genome_a == genome_b:
                    continue
                table_for_distances.add({'genome_a': genome_a,
                                         'genome_b': genome_b,
                                         'distance': float(self.distance_matrix.loc[genome_a, genome_b])})

        self.progress.end()

        table_for_distances.store()


    def update_pan_graph_db_with_layer_orders(self):
        """Adds the newick tree calculated from the graph into the pan-graph-db"""

        # Only add newick tree if one was generated (not empty for identical genomes)
        if self.newick:
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

        self.progress.new('Storing syn gene cluster nodes in pan-graph-db')
        self.progress.update('...')

        table_for_nodes = TableForNodes(self.pan_graph_db_path, run=self.run, progress=self.progress)
        for node, data in self.pangenome_graph.graph.nodes(data=True):
            node_entry = {
                'node_id': node,
                'node_type': data['type'],
                'region_id': data.get('region_id') if data.get('region_id') is not None else '',
                'gene_cluster_id': data['gene_cluster'],
                'synteny_position_json': json.dumps(data['synteny']),
                'gene_calls_json': json.dumps(data['gene_calls']),
                'alignment_summary': json.dumps(data['alignment']),
                'node_x': data['position'][0],
                'node_y': data['position'][1],
                'component_id': int(data.get('component_id', 0)),
            }

            table_for_nodes.add(node_entry)

        self.progress.end()

        table_for_nodes.store()

        pan_graph_db = dbops.PanGraphDatabase(self.pan_graph_db_path, run=self.run, progress=self.progress, quiet=True,)
        pan_graph_db.db.set_meta_value('num_nodes', len(self.pangenome_graph.graph.nodes()))
        pan_graph_db.disconnect()


    def store_edges_in_pan_graph_db(self):
        """"""

        self.progress.new('Storing syn gene cluster edges in pan-graph-db')
        self.progress.update('...')

        table_for_edges = TableForEdges(self.pan_graph_db_path, run=self.run, progress=self.progress)
        for edge_i, edge_j, data in self.pangenome_graph.graph.edges(data=True):
            edge_entry = {
                'edge_id': data['name'],
                'source': edge_i,
                'target': edge_j,
                'weight': data['weight'],
                'genomes_json': json.dumps(data['genomes']),
            }
            table_for_edges.add(edge_entry)

        self.progress.end()

        table_for_edges.store()

        pan_graph_db = dbops.PanGraphDatabase(self.pan_graph_db_path, run=self.run, progress=self.progress, quiet=True)
        pan_graph_db.db.set_meta_value('num_edges', len(self.pangenome_graph.graph.edges(data=True)))
        pan_graph_db.disconnect()


    def create_pangenome_graph(self):
        """Build the pangenome graph by delegating to ``PangenomeAAIEngine``
        (see :mod:`anvio.panaai`) and mirror the result into
        ``self.pangenome_graph`` (a :class:`PangenomeGraphManager`).

        In panmode the parent gene-cluster map is pulled from the pan-db and
        passed to the engine; in non-panmode the engine still runs but only
        ``core`` / ``accessory`` / ``singleton`` node types are produced and
        alignments are unavailable downstream.
        """

        # In panmode, get the (genome, gene_caller_id) -> gene_cluster map.
        gene_clusters = None
        if self.pan_super is not None:
            gene_clusters = {}
            for genome, gid_map in self.pan_super.gene_callers_id_to_gene_cluster.items():
                for gid, gc_name in gid_map.items():
                    gene_clusters[(genome, gid)] = gc_name
            self.run.info('Gene clusters loaded from pan-db', len(gene_clusters))
        else:
            self.run.info_single('No pan-db given — alignments will be unavailable in the output '
                                 'and node types collapse to core/accessory/singleton.',
                                 mc='yellow')

        # Drive the AAI engine.
        engine = panaai.PangenomeAAIEngine(self.args, r=self.run, p=self.progress)
        G, lines, line_names, line_to_genome, in_g_flip = engine.process(
            gene_clusters=gene_clusters)

        # Stash bookkeeping for downstream consumers (add_layers, summarize,
        # calculate_graph_distance).
        self.lines = lines
        self.line_names = line_names
        self.line_to_genome = line_to_genome
        self.in_g_flip = in_g_flip
        self.genome_calls = engine.genome_calls

        # Per-edge genome sets, derived by walking committed lines.
        edge_genomes = self._compute_edge_genome_sets(
            G, lines, line_names, line_to_genome, in_g_flip)

        # Per-gene synteny position: a global per-genome counter walked in
        # graph-traversal order (a line is iterated reversed when in_g_flip
        # is True). Within a contig, consecutive positions follow forward
        # graph edges so the JS edge_synteny[source][target] lookup hits;
        # between contigs the lookup naturally misses and the JS flushes
        # the genome track, which is the correct rendering.
        gene_to_synteny_pos = {}
        lines_by_genome = {}
        for li, name in enumerate(line_names):
            if li not in in_g_flip:
                continue
            lines_by_genome.setdefault(line_to_genome[name], []).append(li)
        for genome, line_indices in lines_by_genome.items():
            line_indices.sort(key=lambda li: line_names[li])
            pos = 0
            for li in line_indices:
                tokens = lines[li]
                seq = reversed(tokens) if in_g_flip[li] else tokens
                name = line_names[li]
                for tok in seq:
                    gene_to_synteny_pos[f"{name}:{tok}"] = pos
                    pos += 1

        # Initialize layers_data so add_layers can populate it later.
        self.layers_data = {}

        # Mirror G into self.pangenome_graph.
        self.pangenome_graph = PangenomeGraphManager(run=self.run, progress=self.progress)

        for node, attrs in G.nodes(data=True):
            gene_calls = {}
            synteny = {}
            for endpoint in attrs.get('genes', ()):
                line_name, _, gid_str = endpoint.rpartition(':')
                if not gid_str.isdigit():
                    continue
                genome = line_to_genome.get(line_name)
                if genome is None:
                    continue
                gene_calls[genome] = int(gid_str)
                pos = gene_to_synteny_pos.get(endpoint)
                if pos is not None:
                    synteny[genome] = pos

            parent_gc = node.rsplit('_', 1)[0] if gene_clusters is not None else ''

            self.pangenome_graph.add_node_to_graph(node, {
                'gene_cluster': parent_gc,
                'gene_calls': gene_calls,
                'synteny': synteny,
                'type': attrs.get('type', ''),
                'component_id': int(attrs.get('component_id', 0)),
            })

        for u, v in G.edges():
            genomes = sorted(edge_genomes.get((u, v), set()))
            self.pangenome_graph.add_edge_to_graph(u, v, {
                'weight': float(len(genomes)),
                'genomes': genomes,
            })

        # Assign edge names in iteration order.
        for edge_id, (u, v) in enumerate(self.pangenome_graph.graph.edges()):
            self.pangenome_graph.graph[u][v]['name'] = 'E_' + str(edge_id).zfill(8)

        # We keep ALL components per design decision §4.8; the layout/summarize
        # passes filter by component_id when needed. No connectivity check.
        self.run.info('Graph nodes', self.pangenome_graph.graph.number_of_nodes(), mc='green')
        self.run.info('Graph edges', self.pangenome_graph.graph.number_of_edges())


    def _compute_edge_genome_sets(self, G, lines, line_names, line_to_genome, in_g_flip):
        """Walk each committed line in its placement orientation and
        accumulate ``(super_node_u, super_node_v) -> {genomes}`` for every
        consecutive pair on the line.

        Lines that did not get committed are skipped (they contributed no
        edges to G).
        """
        gene_to_super = {}
        for n, attrs in G.nodes(data=True):
            for endpoint in attrs.get('genes', ()):
                gene_to_super[endpoint] = n

        edge_genomes = {}
        for line_idx, flip in in_g_flip.items():
            name = line_names[line_idx]
            genome = line_to_genome.get(name)
            if genome is None:
                continue
            tokens = lines[line_idx]
            seq = list(reversed(tokens)) if flip else list(tokens)
            prev = None
            for tok in seq:
                super_node = gene_to_super.get(f"{name}:{tok}")
                if super_node is None:
                    continue
                if prev is not None and prev != super_node:
                    edge_genomes.setdefault((prev, super_node), set()).add(genome)
                prev = super_node

        return edge_genomes


    def add_layers(self):
        """Populate the ``alignment`` attribute on every super-node and copy
        user-requested numeric layers from CONTIGS.dbs into ``node['layer']``.

        In panmode, per-genome alignments are pulled from
        ``self.pan_super.gene_clusters_gene_alignments`` (the existing pan-db
        alignment machinery), padded to a common length, and re-summarized.
        In non-panmode the pan-db is unavailable, so we emit empty alignment
        strings and move on. After alignments, the columns named in
        ``self.import_values`` are aggregated per super-node by mean across
        the super-node's genes; non-numeric values are skipped.
        """
        self.run.warning("Adding per-node alignments.",
                         header="ADDING LAYERS", lc="green")

        if not self.gene_alignments_computed or self.pan_super is None:
            self.run.info_single("No pan-db (or no alignments computed) — alignment summaries "
                                 "will be empty.")
            for _node, data in self.pangenome_graph.graph.nodes(data=True):
                data['alignment'] = {g: '' for g in data['gene_calls']}
            self._import_layer_values()
            return

        # Nodes with alignments of differing lengths (we pad and warn).
        diff_len_nodes = []
        # Nodes where no alignments could be recovered for >1 genome.
        no_alignment_nodes = []

        for node, data in self.pangenome_graph.graph.nodes(data=True):
            node_alignments = {}
            node_alignment_lengths = []

            for genome_name, gene_caller_id in data['gene_calls'].items():
                genome_alignments = self.pan_super.gene_clusters_gene_alignments.get(genome_name, {})
                if gene_caller_id in genome_alignments:
                    alignment_summary = genome_alignments[gene_caller_id]
                    parts = alignment_summary.split('|')
                    start = parts[0]
                    summary_code = list(map(int, parts[1:]))

                    if start == '-':
                        sequence = sum(summary_code[1::2]) * 'N'
                    else:
                        sequence = sum(summary_code[0::2]) * 'N'

                    alignment = utils.restore_alignment(sequence, alignment_summary)
                    node_alignment_lengths.append(len(alignment))
                else:
                    alignment = ''
                    node_alignment_lengths.append(0)

                node_alignments[genome_name] = alignment

            valid_alignment_lengths = {l for l in node_alignment_lengths if l > 0}
            if not valid_alignment_lengths:
                data['alignment'] = {g: '' for g in node_alignments}
                if len(node_alignments) > 1:
                    no_alignment_nodes.append(node)
                continue

            alignment_length = max(valid_alignment_lengths)
            if len(valid_alignment_lengths) != 1:
                diff_len_nodes.append(node)

            # poor man's fix-up: pad shorter alignments with gaps so every
            # string is the same length, then drop columns that are all-gap.
            padded = {g: a.ljust(alignment_length, '-') for g, a in node_alignments.items()}
            cleaned = {g: '' for g in padded}
            for i in range(alignment_length):
                column = [padded[g][i] for g in padded]
                if not all(c == '-' for c in column):
                    for g in padded:
                        cleaned[g] += padded[g][i]

            data['alignment'] = {g: utils.summarize_alignment(a) if a else ''
                                 for g, a in cleaned.items()}

        if diff_len_nodes:
            examples = ", ".join(diff_len_nodes[:10])
            more = "" if len(diff_len_nodes) <= 10 else ", ..."
            self.run.warning(f"Alignments have differing lengths for {len(diff_len_nodes)} node(s) "
                             f"(e.g., {examples}{more}). This can happen when the upstream alignment "
                             f"step failed for individual gene clusters. The shorter sequences were "
                             f"padded with gap characters so downstream processing can continue, but "
                             f"the padding may lead to undesirable (or misleading) visualization of "
                             f"the node's sequence content.")

        if no_alignment_nodes:
            examples = ", ".join(no_alignment_nodes[:10])
            more = "" if len(no_alignment_nodes) <= 10 else ", ..."
            self.run.warning(f"No alignments found for {len(no_alignment_nodes)} multi-gene node(s) "
                             f"(e.g., {examples}{more}). Storing empty alignments and continuing.")
        else:
            self.run.info_single("Alignments were found for all nodes and successfully added.")

        self._import_layer_values()


    def remerge_nodes(self):
        """Collapse same-parent-GC super-nodes that look like an engine over-split.

        For every parent gene cluster that produced >= 2 super-nodes, walk
        pairs and merge them when **all four** guards pass:

        1. Same ``component_id`` -- this pass NEVER bridges weakly connected
           components; cross-component homologs are left alone (use the
           engine's tuning knobs if you want them pulled in).
        2. ``nx.lowest_common_ancestor`` is not one of the two nodes -- if
           it is, the pair is in-series along the same path (real
           duplication-like topology) and must not be collapsed.  ``None``
           (no common ancestor) is acceptable: the pair is parallel.
        3. Disjoint ``gene_calls`` genomes -- merging requires the engine's
           same-genome-conflict invariant (<=1 gene per genome per node).
        4. Disjoint ``synteny`` genomes -- redundant with (3) in the
           current flow but kept as defense in depth.

        When a pair merges, the lower-numbered survivor X absorbs Y:
        ``gene_calls`` / ``synteny`` / ``layers_data`` are unioned, edges
        are rewired *manually* (predecessor/successor edges of Y are added
        onto X, with the ``genomes`` list unioned and ``weight`` refreshed
        to ``float(len(genomes))``), and Y is removed.  This avoids
        ``nx.contracted_nodes`` which would silently truncate parallel
        edge attributes.

        After merging, survivors that started as ``rearrangement`` get
        retyped to ``core`` (if every genome is now represented) or
        ``accessory`` (if the merged genome set matches the parent GC's
        full genome coverage in the pan-db).  ``rna`` / ``duplication`` /
        ``singleton`` survivors keep their type.

        Component IDs are NOT recomputed -- guard (1) keeps every merge
        intra-component, so the engine's ``component_id`` attribute stays
        valid.
        """
        if self.no_remerge:
            self.run.info_single('Remerge step skipped (`--no-remerge`).')
            return

        self.run.warning("Remerging nodes. This is extremely useful in case of highly sensitive graph creation "
                         "settings. The algorithm will attempt to find e.g. false rearrangement nodes and join "
                         "them together as a single synteny gene cluster. Use `--no-remerge` to disable this "
                         "step if you experience unexpected results. Merges are blocked across weakly "
                         "connected components and across in-series (ancestor/descendant) pairs.",
                         header="REMERGING SENSITIVE NODES", lc="green")

        graph = self.pangenome_graph.graph

        # Bucket nodes by parent gene cluster.
        gc_to_syns = {}
        for syn, data in graph.nodes(data=True):
            gc_to_syns.setdefault(data.get('gene_cluster', ''), []).append(syn)

        # Precompute parent_gc -> {unique genomes covered} (pan-mode only;
        # in non-panmode every node is its own parent_gc so the outer
        # bucketing loop never finds a pair anyway).
        parent_gc_genomes = {}
        if self.pan_super is not None:
            for genome, gid_map in self.pan_super.gene_callers_id_to_gene_cluster.items():
                for _gid, gc_name in gid_map.items():
                    parent_gc_genomes.setdefault(gc_name, set()).add(genome)

        n_genomes_total = len(self.genome_names)
        original_num_nodes = graph.number_of_nodes()
        merged_pairs = 0
        new_core_num = 0
        new_accessory_num = 0
        n_rejected_asymmetry = 0

        # O(1) reverse view used to compute the lowest common descendant:
        # LCD(a, b) in `graph` == LCA(a, b) in the reverse direction.
        reverse_view = graph.reverse(copy=False) if self.remerge_max_length >= 0 else None

        # Count pairs (the real unit of work) rather than GCs -- most GCs are
        # singletons and skip instantly, while a handful of fat GCs spend all
        # the time in nx.lowest_common_ancestor inside `combinations`.
        total_pairs = sum(len(syns) * (len(syns) - 1) // 2
                          for syns in gc_to_syns.values() if len(syns) >= 2)
        self.progress.new('Remerging nodes', progress_total_items=total_pairs)
        self.progress.update('...')

        for gc, syns in gc_to_syns.items():
            if len(syns) < 2:
                continue
            for syn_a, syn_b in combinations(syns, 2):
                self.progress.increment()
                self.progress.update(f"{merged_pairs} pair(s) merged so far")
                # Either side may have been absorbed by a prior merge.
                if syn_a not in graph or syn_b not in graph:
                    continue

                node_a = graph.nodes[syn_a]
                node_b = graph.nodes[syn_b]

                # Guard 1: components must match (cheapest check first).
                if node_a.get('component_id', 0) != node_b.get('component_id', 0):
                    continue

                # Guard 2: LCA must not be either of the two (in-series).
                lca = nx.lowest_common_ancestor(graph, syn_a, syn_b)
                if lca == syn_a or lca == syn_b:
                    continue

                # Guard 2b: LCA/LCD-asymmetry cap. Reject pairs whose
                # shortest-path distances from the LCA (or to the LCD) differ
                # by more than `remerge_max_length` edges. Each check is
                # applied independently when its anchor exists; whichever
                # exist must ALL pass. The cap is skipped only when it is
                # disabled (< 0) or when NEITHER anchor exists.
                if self.remerge_max_length >= 0:
                    lcd = nx.lowest_common_ancestor(reverse_view, syn_a, syn_b)
                    reject = False

                    if lca is not None:
                        try:
                            dist_a = nx.shortest_path_length(graph, lca, syn_a)
                            dist_b = nx.shortest_path_length(graph, lca, syn_b)
                            if abs(dist_a - dist_b) > self.remerge_max_length:
                                reject = True
                        except nx.NetworkXNoPath:
                            # Real LCA in a DAG always has a path; defensive.
                            reject = True

                    if not reject and lcd is not None:
                        try:
                            dist_a = nx.shortest_path_length(graph, syn_a, lcd)
                            dist_b = nx.shortest_path_length(graph, syn_b, lcd)
                            if abs(dist_a - dist_b) > self.remerge_max_length:
                                reject = True
                        except nx.NetworkXNoPath:
                            reject = True

                    if reject:
                        n_rejected_asymmetry += 1
                        continue

                # Guard 3 / 4: disjoint genome membership.
                if not set(node_a['gene_calls']).isdisjoint(node_b['gene_calls']):
                    continue
                if not set(node_a['synteny']).isdisjoint(node_b['synteny']):
                    continue

                # Pick survivor X = lower numeric suffix.
                a_idx = int(syn_a.rsplit('_', 1)[1])
                b_idx = int(syn_b.rsplit('_', 1)[1])
                if a_idx <= b_idx:
                    syn_x, syn_y = syn_a, syn_b
                else:
                    syn_x, syn_y = syn_b, syn_a
                node_x = graph.nodes[syn_x]
                node_y = graph.nodes[syn_y]

                # Union node-level dicts.
                node_x['gene_calls'] = {**node_x['gene_calls'], **node_y['gene_calls']}
                node_x['synteny'] = {**node_x['synteny'], **node_y['synteny']}

                # Layers data (only present when `--import-values` was used).
                if syn_y in self.layers_data:
                    sink = self.layers_data.setdefault(syn_x, {})
                    for k, v in self.layers_data[syn_y].items():
                        bucket = sink.setdefault(k, [])
                        bucket.extend(v if isinstance(v, list) else [v])
                    del self.layers_data[syn_y]

                # Edge rewire: predecessors of Y -> X.
                for pred in list(graph.predecessors(syn_y)):
                    if pred == syn_x:
                        continue
                    y_edge = graph[pred][syn_y]
                    y_genomes = set(y_edge.get('genomes', []))
                    if graph.has_edge(pred, syn_x):
                        x_edge = graph[pred][syn_x]
                        merged_g = sorted(set(x_edge.get('genomes', [])) | y_genomes)
                        x_edge['genomes'] = merged_g
                        x_edge['weight'] = float(len(merged_g))
                    else:
                        graph.add_edge(pred, syn_x, **dict(y_edge))

                # Edge rewire: successors of Y -> X.
                for succ in list(graph.successors(syn_y)):
                    if succ == syn_x:
                        continue
                    y_edge = graph[syn_y][succ]
                    y_genomes = set(y_edge.get('genomes', []))
                    if graph.has_edge(syn_x, succ):
                        x_edge = graph[syn_x][succ]
                        merged_g = sorted(set(x_edge.get('genomes', [])) | y_genomes)
                        x_edge['genomes'] = merged_g
                        x_edge['weight'] = float(len(merged_g))
                    else:
                        graph.add_edge(syn_x, succ, **dict(y_edge))

                graph.remove_node(syn_y)
                merged_pairs += 1

                # Retype only if survivor began as `rearrangement`.
                if node_x.get('type') == 'rearrangement':
                    n_calls = len(node_x['gene_calls'])
                    if n_calls == n_genomes_total:
                        node_x['type'] = 'core'
                        new_core_num += 1
                    elif gc in parent_gc_genomes and n_calls == len(parent_gc_genomes[gc]):
                        node_x['type'] = 'accessory'
                        new_accessory_num += 1

        self.progress.end()

        self.run.info_single(f"{merged_pairs} pair(s) merged; "
                             f"{original_num_nodes - graph.number_of_nodes()} node(s) removed.")
        self.run.info_single(f"{new_core_num} node(s) retyped 'rearrangement' -> 'core'.")
        self.run.info_single(f"{new_accessory_num} node(s) retyped 'rearrangement' -> 'accessory'.")
        if self.remerge_max_length >= 0:
            self.run.info_single(f"{n_rejected_asymmetry} pair(s) rejected "
                                 f"(LCA/LCD-asymmetry > {self.remerge_max_length}).")

        if not nx.is_directed_acyclic_graph(graph):
            raise ConfigError("Cyclic graphs are not implemented. Remerge produced a cycle, which means "
                              "one of the guards above failed -- please report this with the dataset.")


    def _import_layer_values(self):
        """Copy the numeric per-gene columns named in ``self.import_values``
        (from each genome's ``genes_in_contigs_dict``) onto each super-node's
        ``layer`` dict, aggregated by mean across the super-node's genes.
        Values that cannot be coerced to float are skipped; a column with no
        usable values for a node is omitted from that node's layer.
        """
        if not self.import_values:
            return

        gene_lookup = {}
        for genome, calls in (self.genome_calls or {}).items():
            gene_lookup[genome] = {c['gene_callers_id']: c for c in calls}

        missing_columns = set()
        per_column_seen = {col: False for col in self.import_values}

        for _node, data in self.pangenome_graph.graph.nodes(data=True):
            if 'layer' not in data or data['layer'] is None:
                data['layer'] = {}
            for col in self.import_values:
                numeric_values = []
                for genome, gid in data['gene_calls'].items():
                    gene_info = gene_lookup.get(genome, {}).get(gid)
                    if gene_info is None or col not in gene_info:
                        if gene_info is not None and col not in gene_info:
                            missing_columns.add(col)
                        continue
                    raw = gene_info[col]
                    if raw is None:
                        continue
                    try:
                        numeric_values.append(float(raw))
                    except (TypeError, ValueError):
                        continue
                if numeric_values:
                    data['layer'][col] = sum(numeric_values) / len(numeric_values)
                    per_column_seen[col] = True

        skipped = [c for c, seen in per_column_seen.items() if not seen]
        if skipped:
            self.run.warning(f"No numeric values found for imported column(s): "
                             f"{', '.join(sorted(skipped))}. Those layers will be empty.")

        imported = [c for c, seen in per_column_seen.items() if seen]
        if imported:
            self.run.info("Imported per-node layers", ', '.join(imported))


class FragmentedGeneAnnotator():
    """Identifies fragmented genes (pseudogenes) across a pangenome.

    In pangenomes, a gene that is intact in some genomes may be split into two or more
    adjacent open reading frames in others -- typically due to a premature stop codon
    introduced by a point mutation or transposon insertion. The MCL algorithm correctly
    groups the fragments with the full-length gene into one gene cluster, but their
    presence creates spurious singleton nodes in the pangenome graph.

    This class detects such cases by looking for multiple genes from the same genome
    within a single gene cluster that are adjacent on the same contig. It then compares
    the lengths of these fragments against the full-length representative (the longest
    gene in the cluster from a genome where the gene is not split) and annotates them
    under a 'PSEUDO_GENES' function source in each relevant contigs database.

    Labels assigned:
        - `fragmented_gene`: the longest fragment in a genome, likely still functional
          (only if it is >= min_full_length_ratio of the full-length representative)
        - `gene_fragment`: shorter fragments, or all fragments if none meets the
          length threshold
    """

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        self.pan_db_path = A('pan_db')
        self.genomes_storage_path = A('genomes_storage')
        self.external_genomes_path = A('external_genomes')
        self.min_full_length_ratio = A('min_full_length_ratio') or 0.50
        self.max_combined_length_ratio = A('max_combined_length_ratio') or 1.20
        self.skip_reporting = A('skip_reporting') or False
        self.report_only = A('report_only') or False
        self.find_stray_fragments = A('find_stray_fragments') or False
        self.annotation_source = A('annotation_source')

        if not self.pan_db_path:
            raise ConfigError("You must provide a pan database path.")

        if not self.genomes_storage_path:
            raise ConfigError("You must provide a genomes storage path.")

        if not self.external_genomes_path:
            raise ConfigError("You must provide an external genomes file.")

        utils.is_pan_db_and_genomes_storage_db_compatible(self.pan_db_path, self.genomes_storage_path)


    def process(self):
        """Main entry point for fragmented gene annotation."""

        # we import here to avoid circular imports since genomedescriptions imports dbops
        # which imports panops
        import anvio.genomedescriptions as genomedescriptions
        genome_desc_args = argparse.Namespace(external_genomes=self.external_genomes_path, internal_genomes=None,
                                              skip_checking_genome_hashes=False, just_do_it=False, gene_caller=None,
                                              list_hmm_sources=False, list_available_gene_names=False)
        self.genome_descriptions = genomedescriptions.GenomeDescriptions(genome_desc_args, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))
        self.genome_descriptions.load_genomes_descriptions(skip_functions=True, init=False)

        # initialize pan superclass and gene clusters
        pan_args = argparse.Namespace(pan_db=self.pan_db_path, genomes_storage=self.genomes_storage_path)
        self.pan_super = dbops.PanSuperclass(pan_args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        self.pan_super.init_gene_clusters()

        # if the user wants gene cluster functions in the report, initialize them now
        if self.annotation_source:
            self.pan_super.init_gene_clusters_functions()

            if self.annotation_source not in self.pan_super.gene_clusters_function_sources:
                available_sources = ', '.join(sorted(self.pan_super.gene_clusters_function_sources)) if self.pan_super.gene_clusters_function_sources else 'None'
                raise ConfigError(f"The annotation source '{self.annotation_source}' is not available in the genomes storage. "
                                  f"Here are the sources that are available: {available_sources}.")

        # initialize genomes storage for sequence access
        self.genomes_storage = GenomeStorage(self.genomes_storage_path, run=terminal.Run(verbose=False), progress=terminal.Progress(verbose=False))

        # load genes_in_contigs_dict for each genome so we know contig membership and positions
        self.genes_in_contigs = {}
        self.contig_gene_order = {}  # {genome_name: {contig: [gene_ids sorted by start]}}
        for genome_name in self.genome_descriptions.genomes:
            contigs_db_path = self.genome_descriptions.genomes[genome_name]['contigs_db_path']
            contigs_db_args = argparse.Namespace(contigs_db=contigs_db_path)
            contigs_super = dbops.ContigsSuperclass(contigs_db_args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
            self.genes_in_contigs[genome_name] = contigs_super.genes_in_contigs_dict

            # pre-compute sorted gene order per contig for fast adjacency lookups
            genes_by_contig = {}
            for gene_id, gene_info in self.genes_in_contigs[genome_name].items():
                contig = gene_info['contig']
                if contig not in genes_by_contig:
                    genes_by_contig[contig] = []
                genes_by_contig[contig].append(gene_id)

            for contig in genes_by_contig:
                genes_by_contig[contig].sort(key=lambda g: self.genes_in_contigs[genome_name][g]['start'])

            self.contig_gene_order[genome_name] = genes_by_contig

        gene_clusters = self.pan_super.gene_clusters

        # find fragmentation events across all gene clusters
        self.run.warning(None, header="IDENTIFYING FRAGMENTED GENES", lc="green")
        self.run.info_single("Please read the documentation of this program to familiarize yourelf with its "
                             "assumptions and how to make sense of the results displayed below. You can find "
                             "the documentation at https://anvio.org/m/anvi-annotate-fragmented-genes",
                             level=0, mc='green')
        self.run.info('Num genomes', len(self.genome_descriptions.genomes), nl_before=1)
        self.run.info('Num gene clusters', len(gene_clusters))
        self.run.info('Min full-length ratio', self.min_full_length_ratio)
        self.run.info('Max combined length ratio', self.max_combined_length_ratio)
        self.run.info('Search for stray fragments', self.find_stray_fragments)
        self.run.info('Annotation source for report', self.annotation_source or 'None')
        self.run.info('Report only', self.report_only)

        # annotations_per_genome will be {genome_name: {entry_counter: {gene_callers_id, source, accession, function, e_value}}}
        annotations_per_genome = {g: {} for g in self.genome_descriptions.genomes}
        entry_counter_per_genome = {g: 0 for g in self.genome_descriptions.genomes}

        # collect all fragmentation events via in-cluster adjacency
        in_cluster_events = self.scan_in_cluster_fragmentation(gene_clusters)

        # if requested, also look for stray out-of-frame fragments adjacent to truncated genes
        stray_events = []
        if self.find_stray_fragments:
            # build a set of gene_callers_ids already flagged by the in-cluster scan so we
            # don't double-annotate them
            already_flagged = set()
            for _, fragmentation_events, _, _, _ in in_cluster_events:
                for genome_name, adjacent_group in fragmentation_events:
                    for gene_id in adjacent_group:
                        already_flagged.add((genome_name, gene_id))

            stray_events = self.scan_stray_fragment_events(gene_clusters, already_flagged)

        all_events = in_cluster_events + stray_events
        stray_gene_cluster_ids = set(gc_id for gc_id, _, _, _, _ in stray_events)

        total_fragmented_genes = 0
        total_gene_fragments = 0
        total_stray_fragmented_genes = 0
        total_stray_gene_fragments = 0
        gene_clusters_with_fragmentation = 0

        # report and annotate
        for gene_cluster_id, fragmentation_events, reference_length, reference_genome, reference_gene_id in all_events:
            gene_clusters_with_fragmentation += 1

            if not self.skip_reporting:
                # get consensus function for this gene cluster if an annotation source was provided
                gc_function = None
                if self.annotation_source:
                    _acc, _func = self.pan_super.get_gene_cluster_function_summary(gene_cluster_id, self.annotation_source)
                    gc_function=(_acc.split('!!!')[0] if _acc else _acc, _func.split('!!!')[0] if _func else _func)

                self.report_gene_cluster(gene_cluster_id, fragmentation_events, reference_length, reference_genome, reference_gene_id, gc_function=gc_function)

            for genome_name, adjacent_group in fragmentation_events:
                gene_lengths = []
                for gene_callers_id in adjacent_group:
                    gene_info = self.genes_in_contigs[genome_name][gene_callers_id]
                    gene_length = gene_info['stop'] - gene_info['start']
                    gene_lengths.append((gene_callers_id, gene_length))

                gene_lengths.sort(key=lambda x: x[1], reverse=True)

                longest_gene_id, longest_length = gene_lengths[0]
                ratio = longest_length / reference_length

                for gene_callers_id, gene_length in gene_lengths:
                    frag_ratio = gene_length / reference_length

                    is_stray = gene_cluster_id in stray_gene_cluster_ids

                    if gene_callers_id == longest_gene_id and ratio >= self.min_full_length_ratio:
                        label = 'fragmented_gene'
                        function_text = (f"Putative fragmented gene ({ratio * 100:.1f}% of full-length), "
                                         f"based on a homologous gene in {reference_genome} with gene caller id {reference_gene_id}")
                        total_fragmented_genes += 1
                        if is_stray:
                            total_stray_fragmented_genes += 1
                    else:
                        label = 'gene_fragment'
                        function_text = (f"Putative gene fragment ({frag_ratio * 100:.1f}% of full-length), "
                                         f"based on a homologous gene in {reference_genome} with gene caller id {reference_gene_id}")
                        total_gene_fragments += 1
                        if is_stray:
                            total_stray_gene_fragments += 1

                    entry_id = entry_counter_per_genome[genome_name]
                    annotations_per_genome[genome_name][entry_id] = {
                        'gene_callers_id': gene_callers_id,
                        'source': 'PSEUDO_GENES',
                        'accession': label,
                        'function': function_text,
                        'e_value': 0,
                    }
                    entry_counter_per_genome[genome_name] += 1

        num_non_singleton_gene_clusters = sum(1 for gc_id in gene_clusters if sum(1 for g in gene_clusters[gc_id] if gene_clusters[gc_id][g]) > 1)
        pct_with_fragmentation = gene_clusters_with_fragmentation / num_non_singleton_gene_clusters * 100 if num_non_singleton_gene_clusters else 0

        self.run.info('Non-singleton gene clusters', num_non_singleton_gene_clusters, nl_before=1)
        self.run.info('Gene clusters with fragmentation', f"{gene_clusters_with_fragmentation} ({pct_with_fragmentation:.1f}% of non-singleton GCs)")
        self.run.info('Total fragmented genes', total_fragmented_genes)
        self.run.info('Total gene fragments', total_gene_fragments)

        if self.find_stray_fragments:
            self.run.info('Stray fragmented genes', total_stray_fragmented_genes, nl_before=1)
            self.run.info('Stray gene fragments', total_stray_gene_fragments)

        if self.report_only:
            self.run.warning("The --report-only flag is set, so no annotations have been written to any contigs database.",
                             header="REPORT ONLY MODE", lc="yellow")
            return

        # write annotations to each contigs-db
        self.progress.new("Annotating contigs-dbs", progress_total_items=len(self.genome_descriptions.genomes))
        genomes_annotated = 0
        for genome_name in self.genome_descriptions.genomes:
            contigs_db_path = self.genome_descriptions.genomes[genome_name]['contigs_db_path']
            self.progress.update(f"Working on {genome_name} ...", increment=True)
            functions_dict = annotations_per_genome[genome_name]

            if not len(functions_dict):
                continue

            gene_functions_table = TableForGeneFunctions(contigs_db_path, terminal.Run(verbose=False), terminal.Progress(verbose=False))
            gene_functions_table.create(functions_dict)

            genomes_annotated += 1

        self.progress.end()

        self.run.info('Contigs databases annotated', genomes_annotated, nl_before=1)


    def find_fragmentation_events(self, gene_cluster_id):
        """For a given gene cluster, find groups of adjacent genes from the same genome on
        the same contig that likely represent a fragmented gene.

        Returns a list of (genome_name, [gene_callers_id, ...]) tuples, where each inner
        list is a group of 2+ adjacent genes forming a fragmentation event.
        """

        gene_clusters = self.pan_super.gene_clusters
        fragmentation_events = []

        for genome_name in gene_clusters[gene_cluster_id]:
            gene_ids = gene_clusters[gene_cluster_id][genome_name]

            # only genomes with 2+ genes in this cluster can have fragmentation
            if len(gene_ids) < 2:
                continue

            genes_in_contigs = self.genes_in_contigs[genome_name]

            # group genes by contig
            genes_by_contig = {}
            for gene_id in gene_ids:
                if gene_id not in genes_in_contigs:
                    continue
                contig = genes_in_contigs[gene_id]['contig']
                if contig not in genes_by_contig:
                    genes_by_contig[contig] = []
                genes_by_contig[contig].append(gene_id)

            for contig, contig_gene_ids in genes_by_contig.items():
                if len(contig_gene_ids) < 2:
                    continue

                # use precomputed contig gene order for fast adjacency lookups
                all_genes_on_contig = self.contig_gene_order[genome_name].get(contig, [])
                gene_to_position = {g: i for i, g in enumerate(all_genes_on_contig)}

                # find groups of genes that are adjacent (consecutive positions on the contig)
                positions = [(gene_id, gene_to_position[gene_id]) for gene_id in contig_gene_ids]
                positions.sort(key=lambda x: x[1])

                # walk through and build groups of consecutive positions
                current_group = [positions[0][0]]
                for i in range(1, len(positions)):
                    if positions[i][1] == positions[i - 1][1] + 1:
                        current_group.append(positions[i][0])
                    else:
                        if len(current_group) >= 2:
                            fragmentation_events.append((genome_name, current_group))
                        current_group = [positions[i][0]]

                if len(current_group) >= 2:
                    fragmentation_events.append((genome_name, current_group))

        return fragmentation_events


    def get_full_length_reference(self, gene_cluster_id, fragmentation_events):
        """Find the full-length representative for a gene cluster.

        The reference is the longest gene from a genome where the gene cluster has exactly
        one gene (i.e., an unfragmented genome). Returns (length, genome_name, gene_callers_id)
        or (None, None, None) if no suitable reference exists.
        """

        gene_clusters = self.pan_super.gene_clusters

        # collect genomes with fragmentation for this cluster
        fragmented_genomes = set(genome_name for genome_name, _ in fragmentation_events)

        best_length = 0
        best_genome = None
        best_gene_id = None

        for genome_name in gene_clusters[gene_cluster_id]:
            gene_ids = gene_clusters[gene_cluster_id][genome_name]

            if len(gene_ids) != 1:
                continue

            if genome_name in fragmented_genomes:
                continue

            gene_id = gene_ids[0]
            gene_info = self.genes_in_contigs[genome_name][gene_id]
            gene_length = gene_info['stop'] - gene_info['start']

            if gene_length > best_length:
                best_length = gene_length
                best_genome = genome_name
                best_gene_id = gene_id

        if best_length == 0:
            return None, None, None

        return best_length, best_genome, best_gene_id


    def report_gene_cluster(self, gene_cluster_id, fragmentation_events, reference_length, reference_genome, reference_gene_id, gc_function=None):
        """Print a terminal visualization for a gene cluster with fragmentation events.

        Shows each genome's gene(s) as colored bars proportional to length, with fragment
        positioning based on actual start/stop offsets within each genome.
        """

        # local import to avoid pulling ttycolors into every panops consumer
        from anvio.ttycolors import color_text

        gene_clusters = self.pan_super.gene_clusters
        genes_in_contigs = self.genes_in_contigs

        # collect all genomes and their genes for this cluster
        fragmented_genomes = {}
        for genome_name, adjacent_group in fragmentation_events:
            if genome_name not in fragmented_genomes:
                fragmented_genomes[genome_name] = []
            fragmented_genomes[genome_name].append(adjacent_group)

        # identify stray fragment genes that belong to a different gene cluster so
        # the report can show them alongside the truncated gene they were paired with
        stray_genes_info = {}
        for genome_name, adjacent_group in fragmentation_events:
            cluster_gene_ids = set(gene_clusters[gene_cluster_id].get(genome_name, []))
            for g in adjacent_group:
                if g not in cluster_gene_ids:
                    home_cluster = self.pan_super.gene_callers_id_to_gene_cluster.get(genome_name, {}).get(g, '?')
                    stray_genes_info[(genome_name, g)] = home_cluster

        # determine the longest fragment per genome for labeling
        longest_fragment_per_genome = {}
        for genome_name, groups in fragmented_genomes.items():
            for group in groups:
                gene_lengths = [(g, genes_in_contigs[genome_name][g]['stop'] - genes_in_contigs[genome_name][g]['start']) for g in group]
                gene_lengths.sort(key=lambda x: x[1], reverse=True)
                longest_id = gene_lengths[0][0]
                longest_length = gene_lengths[0][1]
                ratio = longest_length / reference_length
                longest_fragment_per_genome.setdefault(genome_name, set())
                if ratio >= self.min_full_length_ratio:
                    longest_fragment_per_genome[genome_name].add(longest_id)

        # determine bar width (terminal characters for the reference gene)
        bar_width = 80

        # collect rows for display
        rows = []

        # determine column widths
        max_genome_len = 0
        max_gene_id_len = 0

        for genome_name in gene_clusters[gene_cluster_id]:
            gene_ids = gene_clusters[gene_cluster_id][genome_name]
            if not gene_ids:
                continue
            if len(genome_name) > max_genome_len:
                max_genome_len = len(genome_name)
            for gene_id in gene_ids:
                gene_id_str = str(gene_id)
                if len(gene_id_str) > max_gene_id_len:
                    max_gene_id_len = len(gene_id_str)

        for (_, gene_id) in stray_genes_info:
            if len(str(gene_id)) > max_gene_id_len:
                max_gene_id_len = len(str(gene_id))

        for genome_name in sorted(gene_clusters[gene_cluster_id].keys()):
            gene_ids = gene_clusters[gene_cluster_id][genome_name]
            if not gene_ids:
                continue

            is_fragmented_genome = genome_name in fragmented_genomes

            if is_fragmented_genome:
                # for fragmented genomes, get all gene positions relative to the group span
                all_fragment_genes = set()
                for group in fragmented_genomes[genome_name]:
                    all_fragment_genes.update(group)

                # include stray fragment neighbors (from other clusters) in the display
                stray_in_genome = {g for (gn, g) in stray_genes_info if gn == genome_name}
                genes_to_show = list(gene_ids) + sorted(stray_in_genome - set(gene_ids))

                for gene_id in sorted(genes_to_show, key=lambda g: genes_in_contigs[genome_name][g]['start']):
                    gene_info = genes_in_contigs[genome_name][gene_id]
                    gene_length = gene_info['stop'] - gene_info['start']

                    if gene_id in all_fragment_genes:
                        # find the group this gene belongs to
                        group_for_gene = None
                        for group in fragmented_genomes[genome_name]:
                            if gene_id in group:
                                group_for_gene = group
                                break

                        # compute span of this fragment group
                        group_starts = [genes_in_contigs[genome_name][g]['start'] for g in group_for_gene]
                        group_stops = [genes_in_contigs[genome_name][g]['stop'] for g in group_for_gene]
                        span_start = min(group_starts)
                        span_length = max(group_stops) - span_start

                        # position within the span, scaled to bar_width
                        if span_length > 0:
                            rel_start = (gene_info['start'] - span_start) / span_length
                            rel_end = (gene_info['stop'] - span_start) / span_length
                        else:
                            rel_start = 0
                            rel_end = 1

                        # scale to reference length proportion
                        span_ratio = span_length / reference_length if reference_length > 0 else 1
                        scaled_bar_width = int(bar_width * min(span_ratio, 1.0))
                        if scaled_bar_width < 1:
                            scaled_bar_width = bar_width

                        bar_start = int(rel_start * scaled_bar_width)
                        bar_end = int(rel_end * scaled_bar_width)
                        if bar_end == bar_start:
                            bar_end = bar_start + 1

                        bar = list('░' * scaled_bar_width + '░' * (bar_width - scaled_bar_width))
                        for i in range(bar_start, min(bar_end, len(bar))):
                            bar[i] = '█'

                        bar_str = ''.join(bar)

                        is_stray = (genome_name, gene_id) in stray_genes_info

                        if gene_id == reference_gene_id and genome_name == reference_genome:
                            color = 'green'
                            label = 'reference'
                        elif gene_id in longest_fragment_per_genome.get(genome_name, set()):
                            color = 'blue'
                            label = 'fragmented_gene'
                        else:
                            color = 'red'
                            label = 'gene_fragment'

                        if is_stray:
                            label += f" (stray, from {stray_genes_info[(genome_name, gene_id)]})"

                        length_pct = f"{gene_length / reference_length * 100:.1f}%"

                        rows.append((genome_name, str(gene_id), color_text(bar_str, color), f"{gene_length:>6} nt  {length_pct:>6}  {label}"))
                    else:
                        # gene in this genome but not part of a fragment group
                        gene_ratio = gene_length / reference_length if reference_length > 0 else 1
                        filled = max(1, int(bar_width * min(gene_ratio, 1.0)))
                        bar_str = '█' * filled + '░' * (bar_width - filled)

                        rows.append((genome_name, str(gene_id), color_text(bar_str, 'gray'), f"{gene_length:>6} nt"))
            else:
                # non-fragmented genome
                for gene_id in gene_ids:
                    gene_info = genes_in_contigs[genome_name][gene_id]
                    gene_length = gene_info['stop'] - gene_info['start']
                    gene_ratio = gene_length / reference_length if reference_length > 0 else 1
                    filled = max(1, int(bar_width * min(gene_ratio, 1.0)))
                    bar_str = '█' * filled + '░' * (bar_width - filled)

                    length_pct = f"{gene_length / reference_length * 100:.1f}%"

                    if gene_id == reference_gene_id and genome_name == reference_genome:
                        color = 'green'
                        label = 'reference'
                    else:
                        color = 'gray'
                        label = ''

                    rows.append((genome_name, str(gene_id), color_text(bar_str, color), f"{gene_length:>6} nt  {length_pct:>6}  {label}"))

        # compute the total content width from known column sizes
        max_info_len = max(len(info) for _, _, _, info in rows) if rows else 0
        content_width = 2 + max_genome_len + 2 + max_gene_id_len + 2 + bar_width + 2 + max_info_len

        # print the report anvi'o way
        run_width = self.run.width
        self.run.width = content_width
        header = f"{gene_cluster_id}"
        self.run.warning(None, header=header)
        self.run.width = run_width

        if self.annotation_source:
            if gc_function:
                func_str = f"{self.annotation_source} Consensus: {gc_function[0]} | {gc_function[1]}"
            else:
                func_str = f"{self.annotation_source} Consensus: Unknown"
            self.run.info_single(func_str, level=0, mc='red', nl_after=1, cut_after=None)

        prev_genome = None
        for genome_name, gene_id_str, bar_str, info in rows:
            if prev_genome is not None and genome_name != prev_genome:
                nl_before=1
            else:
                nl_before=0

            prev_genome = genome_name
            self.run.info_single(f"   {genome_name:<{max_genome_len}}  {gene_id_str:>{max_gene_id_len}}  {bar_str}  {info}",
                                 cut_after=None, nl_before=nl_before, pretty_indentation=False, level=0)


    def scan_in_cluster_fragmentation(self, gene_clusters):
        """Scan all gene clusters for in-cluster fragmentation events.

        This is the standard algorithm: for each gene cluster, look for genomes that
        contribute 2+ adjacent genes on the same contig, which indicates a gene that has
        been split by a premature stop codon.

        Returns a list of tuples:
            (gene_cluster_id, fragmentation_events, reference_length, reference_genome, reference_gene_id)
        """

        all_events = []

        self.progress.new("Scanning gene clusters", progress_total_items=len(gene_clusters))
        for gene_cluster_id in gene_clusters:
            self.progress.update(f"Processing {gene_cluster_id} ...", increment=True)

            fragmentation_events = self.find_fragmentation_events(gene_cluster_id)

            if not fragmentation_events:
                continue

            reference_length, reference_genome, reference_gene_id = self.get_full_length_reference(gene_cluster_id, fragmentation_events)

            if reference_length is None:
                continue

            # filter out groups whose combined gene length far exceeds the reference,
            # which indicates tandem paralogs (gene duplications) rather than a gene
            # split by a premature stop codon. in a true fragmentation event the
            # fragments should sum to roughly the reference length, not 2x or more.
            filtered_events = []
            for genome_name, adjacent_group in fragmentation_events:
                combined_length = sum(self.genes_in_contigs[genome_name][g]['stop'] - self.genes_in_contigs[genome_name][g]['start'] for g in adjacent_group)
                if combined_length <= reference_length * self.max_combined_length_ratio:
                    filtered_events.append((genome_name, adjacent_group))

            if not filtered_events:
                continue

            all_events.append((gene_cluster_id, filtered_events, reference_length, reference_genome, reference_gene_id))

        self.progress.end()

        return all_events


    def scan_stray_fragment_events(self, gene_clusters, already_flagged):
        """Scan for out-of-frame gene fragments that ended up in different gene clusters.

        When a premature stop codon splits a gene and the downstream fragment is in a
        different reading frame, the fragment will not cluster with the original gene. Instead,
        it appears as a short gene in a different gene cluster. This method detects such cases
        by looking for genomes where a gene cluster contains a single gene that is significantly
        shorter than the full-length reference, and then checking whether an adjacent gene on
        the same contig (belonging to a different gene cluster) fills in the missing length.

        Parameters
        ----------
        gene_clusters : dict
            The gene_clusters dict from PanSuperclass.
        already_flagged : set
            Set of (genome_name, gene_callers_id) tuples already identified by the in-cluster
            scan, to avoid double-annotation.

        Returns a list of tuples with the same structure as scan_in_cluster_fragmentation.
        """

        stray_events = []

        # reverse lookup: gene_callers_id -> gene_cluster_id per genome
        gene_to_cluster = self.pan_super.gene_callers_id_to_gene_cluster

        self.progress.new("Scanning for stray fragments", progress_total_items=len(gene_clusters))
        for gene_cluster_id in gene_clusters:
            self.progress.update(f"Processing {gene_cluster_id} ...", increment=True)

            # first, determine the full-length reference for this cluster using genomes
            # that contribute exactly one gene and are not fragmented
            best_length = 0
            best_genome = None
            best_gene_id = None

            for genome_name in gene_clusters[gene_cluster_id]:
                gene_ids = gene_clusters[gene_cluster_id][genome_name]
                if len(gene_ids) != 1:
                    continue

                gene_id = gene_ids[0]
                if (genome_name, gene_id) in already_flagged:
                    continue

                gene_info = self.genes_in_contigs[genome_name][gene_id]
                gene_length = gene_info['stop'] - gene_info['start']

                if gene_length > best_length:
                    best_length = gene_length
                    best_genome = genome_name
                    best_gene_id = gene_id

            if best_length == 0:
                continue

            reference_length = best_length
            reference_genome = best_genome
            reference_gene_id = best_gene_id

            fragmentation_events = []

            for genome_name in gene_clusters[gene_cluster_id]:
                gene_ids = gene_clusters[gene_cluster_id][genome_name]

                # we are looking for genomes with a single, truncated gene in this cluster
                if len(gene_ids) != 1:
                    continue

                gene_id = gene_ids[0]

                if (genome_name, gene_id) in already_flagged:
                    continue

                gene_info = self.genes_in_contigs[genome_name][gene_id]
                gene_length = gene_info['stop'] - gene_info['start']

                # skip if this gene is already close to full length
                if gene_length >= reference_length * self.min_full_length_ratio:
                    continue

                # check adjacent genes on the same contig
                contig = gene_info['contig']
                contig_genes = self.contig_gene_order[genome_name].get(contig, [])
                if not contig_genes:
                    continue

                gene_position = {g: i for i, g in enumerate(contig_genes)}
                if gene_id not in gene_position:
                    continue

                pos = gene_position[gene_id]

                # look at immediate neighbors (upstream and downstream)
                neighbor_ids = []
                if pos > 0:
                    neighbor_ids.append(contig_genes[pos - 1])
                if pos < len(contig_genes) - 1:
                    neighbor_ids.append(contig_genes[pos + 1])

                for neighbor_id in neighbor_ids:
                    if (genome_name, neighbor_id) in already_flagged:
                        continue

                    # the neighbor must be in a *different* gene cluster
                    neighbor_cluster = gene_to_cluster.get(genome_name, {}).get(neighbor_id, None)
                    if neighbor_cluster is None or neighbor_cluster == gene_cluster_id:
                        continue

                    # if the reference genome also has a gene in the neighbor's cluster, the
                    # neighbor is a real independent gene, not a stray fragment (because, and
                    # bear with me here, the genome with the intact full-length gene also has
                    # separate gene in that family .. assumptions assumptions.. but this logic
                    # really fixed the issue of over-identifying bona fide genes that are
                    # distinct asfragments)
                    if reference_genome in gene_clusters.get(neighbor_cluster, {}):
                        continue

                    neighbor_info = self.genes_in_contigs[genome_name][neighbor_id]
                    neighbor_length = neighbor_info['stop'] - neighbor_info['start']

                    # check if the combined span of the truncated gene and its neighbor
                    # approximates the full-length reference (but does not far exceed it,
                    # which would indicate paralogs rather than fragments)
                    combined_length = gene_length + neighbor_length
                    if combined_length >= reference_length * self.min_full_length_ratio and combined_length <= reference_length * self.max_combined_length_ratio:
                        fragmentation_events.append((genome_name, [gene_id, neighbor_id]))
                        # mark these so we don't flag them again from the neighbor's cluster
                        already_flagged.add((genome_name, gene_id))
                        already_flagged.add((genome_name, neighbor_id))
                        break

            if fragmentation_events:
                stray_events.append((gene_cluster_id, fragmentation_events, reference_length, reference_genome, reference_gene_id))

        self.progress.end()

        return stray_events
