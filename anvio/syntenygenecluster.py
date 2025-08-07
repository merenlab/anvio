# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Classes for pan operations.

    anvi-pan-genome is the default client using this module
"""

import os
import re
import yaml
import random
import argparse
import numpy as np
import pandas as pd
import itertools as it
import matplotlib.pyplot as plt

from itertools import chain
from scipy.optimize import curve_fit
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, to_tree

from warnings import simplefilter
simplefilter("ignore", category=pd.errors.PerformanceWarning)

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.dbinfo import DBInfo

from anvio.genomestorage import GenomeStorage

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alexander Henoch"
__email__ = "ahenoch@outlook.de"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print


# TODO implement multithreading here. This will speed up the SynCluster Algorithm
# ANCHOR - SyntenyGeneCluster
class SyntenyGeneCluster():
    """A class to mine a dataframe containing all the information related to gene calls
    e.g. functions, contig, direction, etc. Furthermore the run_contextualize_paralogs_algorithm
    function can add another column to the mined dataframe, containing the syncluster ID.
    This ID is unique per genome and can be used before creating pangenome graphs.
    """

    def __init__(self, args, r=run, p=progress):
        self.args = args
        self.run = r
        self.progress = p

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.pan_db_path = A('pan_db')
        self.external_genomes = A('external_genomes')
        self.genomes_storage = A('genomes_storage')
        self.pan_graph_yaml = A('pan_graph_yaml')
        self.output_dir = A('output_dir')

        self.n = A('n')
        self.alpha = A('alpha')
        self.beta = A('beta')
        self.gamma = A('gamma')
        self.delta = A('delta')
        self.inversion_aware = A('inversion_aware')
        self.min_k = A('min_k')
        self.output_synteny_gene_cluster_dendrogram = A('output_synteny_gene_cluster_dendrogram')

        if self.pan_graph_yaml:
            with open(self.pan_graph_yaml) as file:
                self.yaml_file = yaml.safe_load(file)
        else:
            self.yaml_file = {}

        if A('genome_names'):
            self.genome_names = A('genome_names').split(',')
        elif self.external_genomes:
            self.genome_names = pd.read_csv(self.external_genomes, header=0, sep="\t")['name'].to_list()
        elif self.pan_graph_yaml:
            self.genome_names = list(self.yaml_file.keys())
        else:
            raise ConfigError("Unfortunately we couldn't find an external genomes files, please add one :)")

        if self.pan_graph_yaml:
            self.functional_annotation_sources_available = []
        else:
            self.functional_annotation_sources_available = DBInfo(self.genomes_storage, expecting='genomestorage').get_functional_annotation_sources()
        # self.pangenome_data_df = pd.DataFrame()
        self.contig_identifiers = {}

        self.just_do_it = A('just_do_it')


    def get_num_contigs_and_genome_length(self, file_path):
        fasta = f.SequenceSource(file_path)

        num_contigs = 0
        length = 0

        while next(fasta):
            num_contigs += 1
            length += len(fasta.seq)

        return(num_contigs, length)


    def get_data_from_YAML(self, contextualize_paralogs=True):
        """Create a data tale form the YAML file"""
        i = 0
        pangenome_data_dict = {}
        for genome in self.genome_names:
            current_pos = 0
            for contig_num, contig in enumerate(self.yaml_file[genome]):
                for gene_call, gene_cluster in enumerate(contig):
                    direction = 'l' if gene_cluster.endswith('!') else 'r'
                    gene_cluster = gene_cluster.replace('!', '')

                    pangenome_data_dict[i] = {
                        'position': gene_call,
                        'genome': genome,
                        'gene_cluster': gene_cluster,
                        'gene_caller_id': gene_call,
                        'contig': genome + '_' + str(contig_num),
                        'direction': direction,
                        'start': current_pos,
                        'stop': current_pos + 400
                    }
                    current_pos += 500
                    i += 1

        pangenome_data_df = pd.DataFrame.from_dict(pangenome_data_dict, orient='index').set_index(["genome", "gene_caller_id"])
        self.run.info_single("Done.")

        if contextualize_paralogs:
            return self.run_contextualize_paralogs_algorithm(pangenome_data_df)
        else:
            pangenome_data_df


    def get_data_from_pan_db(self, contextualize_paralogs=True):
        """Major mining function. Can be used outside of the purpose of creating anvi'o
        pangenome graphs.

        Parameters
        ==========

        Returns
        =======
        pangenome_data_df: ps.DataFrame
            Dataframe containing gene call informations present in the anvi'o dbs.
        """

        self.run.warning(None, header="Loading data from database", lc="green")

        filesnpaths.is_file_tab_delimited(self.external_genomes)

        if not utils.is_all_columns_present_in_TAB_delim_file(["name","contigs_db_path"], self.external_genomes, including_first_column=True):
            raise ConfigError("Your external genomes file does not seem to contain what anvi'o expects to find "
                              "in an external genomes file :/")

        pan_db = dbops.PanSuperclass(self.args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
        pan_db.init_gene_clusters()
        pan_db.init_gene_clusters_functions_summary_dict()
        pan_db.init_items_additional_data()

        gene_cluster_dict = pan_db.gene_callers_id_to_gene_cluster
        additional_info_cluster = pan_db.items_additional_data_dict

        external_genomes = pd.read_csv(self.external_genomes, header=0, sep="\t", names=["name","contigs_db_path"])
        external_genomes.set_index("name", inplace=True)

        # TODO trna and rrna gene cluster update also include more types and split trna and rrna in individual sets. Instead of GC_00000000 add the two at the end.
        # TODO Should the reversed once also have reversed gene call orientation?
        pangenome_data_list = []
        for genome, contigs_db_path in external_genomes.iterrows():
            if genome in self.genome_names:
                args = argparse.Namespace(contigs_db=contigs_db_path.item())
                contigs_db = dbops.ContigsSuperclass(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
                contigs_db.init_functions()

                gene_function_calls_df = pd.DataFrame.from_dict(contigs_db.gene_function_calls_dict, orient="index").rename_axis("gene_caller_id").reset_index()

                # all_gene_calls = caller_id_cluster_df['gene_caller_id'].values.tolist()
                # genes_in_contigs_df = pd.DataFrame.from_dict(contigs_db.get_sequences_for_gene_callers_ids(all_gene_calls, include_aa_sequences=True, simple_headers=True)[1], orient="index").rename_axis("gene_caller_id").reset_index()
                genes_in_contigs_df = pd.DataFrame.from_dict(contigs_db.genes_in_contigs_dict, orient="index").rename_axis("gene_caller_id").reset_index()

                trnas = genes_in_contigs_df.query("source.str.contains('RNA')", engine='python')['gene_caller_id'].tolist()
                caller_id_cluster = {**gene_cluster_dict[genome], **{trna:"GC_00000000" for trna in trnas}}

                caller_id_cluster_df = pd.DataFrame.from_dict(caller_id_cluster, orient="index", columns=["gene_cluster"]).rename_axis("gene_caller_id").reset_index()
                # caller_id_cluster_df["syn_cluster"] = ""

                additional_info_df = pd.DataFrame.from_dict(additional_info_cluster, orient="index").rename_axis("gene_cluster").reset_index()
                additional_info_df.drop(self.functional_annotation_sources_available, axis=1, errors='ignore', inplace=True)

                joined_contigs_df = caller_id_cluster_df.merge(genes_in_contigs_df, on="gene_caller_id", how="left").merge(gene_function_calls_df, on="gene_caller_id", how="left").merge(additional_info_df, on="gene_cluster", how="left")

                joined_contigs_df.fillna("None", inplace=True)
                joined_contigs_df['genome'] = genome
                joined_contigs_df.set_index(["genome", "gene_caller_id"], inplace=True)

                for source in self.functional_annotation_sources_available:
                    joined_contigs_df[source] = joined_contigs_df[source].apply(lambda x: ('None', 'None', 'None') if x == 'None' else x)
                    joined_contigs_df[[source + '_ID', source + '_TEXT', source + '_E_VALUE']] = pd.DataFrame(joined_contigs_df[source].values.tolist(), index=joined_contigs_df.index)
                    joined_contigs_df.drop(source, inplace=True, axis=1)

                joined_contigs_df.reset_index(drop=False, inplace=True)

                for contig, group in joined_contigs_df.groupby(["contig"]):
                    group.sort_values(["contig", "start", "stop"], axis=0, ascending=[True, True, True], ignore_index=True, inplace=True)

                    group.rename_axis("position", inplace=True)
                    group.reset_index(drop=False, inplace=True)

                    pangenome_data_list += [group]
                self.run.info_single(f"Successfully mined data from genome {genome}.")
            else:
                self.run.info_single(f"Skipped genome {genome} on users request.")

        pangenome_data_df = pd.concat(pangenome_data_list)
        self.run.info_single("Done.")

        if contextualize_paralogs:
            return self.run_contextualize_paralogs_algorithm(pangenome_data_df)
        else:
            return pangenome_data_df


    # TODO complete distance works like a charm here but I should consider implementing ward distance for cases where both distances are NOT 1.0
    def k_mer_split(self, gene_cluster, gene_cluster_k_mer_contig_positions, gene_cluster_contig_order, single_copy_core):

        synteny_gene_cluster_type = ''

        if len(gene_cluster_k_mer_contig_positions) == 1:
            clusters = np.array([1])
            synteny_gene_cluster_type = 'singleton'
        else:
            X = np.empty([len(gene_cluster_k_mer_contig_positions), len(gene_cluster_k_mer_contig_positions)])
            gene_cluster_k_mer_contig_dict_i = {}
            gene_cluster_k_mer_contig_dict_j = {}

            for i, gene_cluster_k_mer_contig_position_a in enumerate(gene_cluster_k_mer_contig_positions):
                for j, gene_cluster_k_mer_contig_position_b in enumerate(gene_cluster_k_mer_contig_positions):
                    gene_cluster_k_mer_contig_dict_i[i] = gene_cluster_k_mer_contig_position_a[0]
                    gene_cluster_k_mer_contig_dict_j[j] = gene_cluster_k_mer_contig_position_b[0]
                    X[i][j] = self.k_mer_distance(gene_cluster_k_mer_contig_position_a, gene_cluster_k_mer_contig_position_b, gene_cluster_contig_order, single_copy_core)

            np.fill_diagonal(X, 0.0)
            condensed_X = squareform(X)
            Z = linkage(condensed_X, 'complete')

            # Maye the Z.tolist is not the best way to estimate the steps
            for threshold in sorted(set(sum(Z.tolist(), [])), reverse=True):
                clusters = fcluster(Z, threshold, criterion='distance')
                valid = True
                for c in set(clusters.tolist()):
                    pos = np.where(clusters == c)[0]
                    for i, j in it.combinations(pos, 2):
                        if X[i][j] == 1.0:
                            valid = False
                            if gene_cluster_k_mer_contig_dict_i[i] == gene_cluster_k_mer_contig_dict_j[j]:
                                synteny_gene_cluster_type = 'duplication'
                            else:
                                if not synteny_gene_cluster_type:
                                    synteny_gene_cluster_type = 'rearrangement'

                if valid is True:
                    if not synteny_gene_cluster_type:
                        if len(gene_cluster_k_mer_contig_positions) == len(self.genome_names):
                            synteny_gene_cluster_type = 'core'
                        else:
                            synteny_gene_cluster_type = 'accessory'
                    break

        if gene_cluster == "GC_00000000":
            synteny_gene_cluster_type = 'rna'

        labels = []
        synteny_gene_cluster_id_contig_positions = []
        for cluster, (genome, contig, position, gene_caller_id, gene_cluster_kmer) in zip(clusters, gene_cluster_k_mer_contig_positions):
            k = int(len(gene_cluster_kmer)/2)
            gene_cluster_id = gene_cluster_kmer[k] + '_' + str(cluster)
            synteny_gene_cluster_id_contig_positions += [(genome, contig, position, gene_caller_id, gene_cluster_id)]
            labels += [gene_cluster_id + ' ' + str(gene_cluster_kmer)]

        num_cluster = len(set(clusters.tolist()))
        # if len(gene_cluster_k_mer_contig_positions) != 1:
        if self.output_synteny_gene_cluster_dendrogram and num_cluster > 1:
        # if output_synteny_gene_cluster_dendrogram and len(gene_cluster_k_mer_contig_positions) != 1:
            # cmap = plt.cm.tab20(np.linspace(0, 1, num_cluster))
            # colors = [mcolors.rgb2hex(rgb) for rgb in cmap]
            colors = ["#%06x" % random.randint(0, 0xFFFFFF) for _ in range(num_cluster)]
            label_colors = {label: colors[cluster-1] for label, cluster in zip(labels, clusters)}

            fig = plt.figure(figsize=(15+1+k*8, len(label_colors)))
            ax = plt.gca()
            dendrogram(Z, ax=ax, labels=labels, orientation='right')

            #new_labels = []
            y_tick_labels = ax.get_ymajorticklabels()
            for label in y_tick_labels:
                label_text = label.get_text()
                label.set_color(label_colors[label_text])

                #label_match = re.search(r'\((.+)\)', label_text)
                #new_label = label_match.group(0)
                #new_labels += [new_label]

            # ax.set_yticklabels(new_labels)
            plt.tight_layout()
            fig.savefig(os.path.join(self.output_dir, gene_cluster + '.svg'))
            plt.close(fig)

            labels_matrix = [re.search(r'\((.+)\)', label_matrix).group(0) for label_matrix in labels]
            synteny_gene_cluster_matrix = pd.DataFrame(X, index=labels_matrix, columns=labels_matrix)
            synteny_gene_cluster_matrix.to_csv(os.path.join(self.output_dir, gene_cluster + '.tsv'), sep='\t')

        return(synteny_gene_cluster_id_contig_positions, synteny_gene_cluster_type)

    def single_copy_core_context(self, gene_cluster_order, single_copy_core, position):

        left_side_of_kmer = []
        l = 0
        while True:
            if position - l < 0:
                break
            else:
                gene_cluster = gene_cluster_order[position - l]

            if len(left_side_of_kmer) == self.n:
                break
            elif gene_cluster in single_copy_core:

                left_side_of_kmer += [gene_cluster]

            l += 1

        right_side_of_kmer = []
        r = 0
        while True:

            if position + r >= len(gene_cluster_order):
                break
            else:
                gene_cluster = gene_cluster_order[position + r]

            if len(right_side_of_kmer) == self.n:
                break
            elif gene_cluster in single_copy_core:

                right_side_of_kmer += [gene_cluster]

            r += 1

        return(set(left_side_of_kmer), set(right_side_of_kmer))


    # TODO was it a wise choice to remove reverse kmer comparison? I guess yes better to create more nodes then removing necessary ones
    def k_mer_distance(self, gene_cluster_k_mer_contig_position_a, gene_cluster_k_mer_contig_position_b, gene_cluster_contig_order, single_copy_core):
        genome_a, contig_a, position_a, gene_caller_id_a, gene_cluster_kmer_a = gene_cluster_k_mer_contig_position_a
        gene_cluster_order_a = gene_cluster_contig_order[(genome_a, contig_a)]

        gene_cluster_k_mer_left_context_a, gene_cluster_k_mer_right_context_a = self.single_copy_core_context(gene_cluster_order_a, single_copy_core, position_a)

        genome_b, contig_b, position_b, gene_caller_id_b, gene_cluster_kmer_b = gene_cluster_k_mer_contig_position_b
        gene_cluster_order_b = gene_cluster_contig_order[(genome_b, contig_b)]

        gene_cluster_k_mer_left_context_b, gene_cluster_k_mer_right_context_b = self.single_copy_core_context(gene_cluster_order_b, single_copy_core, position_b)

        gene_cluster_k_mer_left_context_min_length = min(len(gene_cluster_k_mer_left_context_a), len(gene_cluster_k_mer_left_context_b))
        gene_cluster_k_mer_right_context_min_length = min(len(gene_cluster_k_mer_right_context_a), len(gene_cluster_k_mer_right_context_b))

        gene_cluster_k_mer_left_context_score = len(gene_cluster_k_mer_left_context_a.intersection(gene_cluster_k_mer_left_context_b)) / gene_cluster_k_mer_left_context_min_length if gene_cluster_k_mer_left_context_min_length != 0 else 0.0
        gene_cluster_k_mer_right_context_score = len(gene_cluster_k_mer_right_context_a.intersection(gene_cluster_k_mer_right_context_b)) / gene_cluster_k_mer_right_context_min_length if gene_cluster_k_mer_right_context_min_length != 0 else 0.0

        if genome_a == genome_b:
            return(1.0)
        elif self.n != 0 and gene_cluster_k_mer_left_context_score < self.alpha and gene_cluster_k_mer_right_context_score < self.alpha:
            return(1.0)
        elif gene_cluster_kmer_a == gene_cluster_kmer_b:
            return(0.0)
        else:
            r_val = 0
            f_val = 0
            div = len(gene_cluster_kmer_a)-1
            # div = len(gene_cluster_kmer_a)

            for i, (n, m) in enumerate(zip(gene_cluster_kmer_a, gene_cluster_kmer_b)):
                if not i == int(len(gene_cluster_kmer_a) / 2):
                    if (n[0] == '-' and m[0] == '-') or (n[0] == '+' and m[0] == '+'):
                        # if n[1:] != m[1:]:
                        r_val += self.beta
                    elif n[0] == '-' or m[0] == '-' or n[0] == '+' or m[0] == '+':
                        r_val += self.gamma
                    elif n == m:
                        r_val += 1.0
                    else:
                        r_val += 0.0

            if self.inversion_aware:
                for i, (n, m) in enumerate(zip(gene_cluster_kmer_a, gene_cluster_kmer_b[::-1])):
                    if not i == int(len(gene_cluster_kmer_a) / 2):
                        if (n[0] == '-' and m[0] == '-') or (n[0] == '+' and m[0] == '+'):
                            # if n[1:] != m[1:]:
                            f_val += self.beta
                        elif n[0] == '-' or m[0] == '-' or n[0] == '+' or m[0] == '+':
                            f_val += self.gamma
                        elif n == m:
                            f_val += 1.0
                        else:
                            f_val += 0.0

                if r_val >= f_val:
                    sim_value = 1.0 - r_val / div
                else:
                    sim_value = 1.0 - f_val / div

            else:
                sim_value = 1.0 - r_val / div

            return(sim_value if sim_value <= self.delta else 1.0)


    def run_contextualize_paralogs_algorithm(self, pangenome_data_df):
        """A function that resolves the graph context of paralogs based on gene synteny
        information across genomes and adds this information to the pangenome_data_df dataframe
        as a new column called syn_cluster. A syn cluster is a subset of a gene cluster
        containing at most one gene call per genome.

        G1: GC1 ----- GC2 ----- GC2 ----- GC3 ----- GC4 ----- GC4 ----- GC5 ----- GC5
        G2: GC1 ----- GC2 ----- GC2 ----- GC3 ----- GC4 ----- GC4 ----- GC4 ----- GC5

        This example shows two genomes that contain multiple paralogous genes. The following
        algorithm reads the context of these paralogous genes to split the related gene cluster
        into the smaller unit of syn clusters based on the surrounding similariy.

        G1: GC1_1 --- GC2_1 --- GC2_2 --- GC3_1 --- GC4_1 ------------- GC4_3 --- GC5_1 --- GC5_2
        G2: GC1_1 --- GC2_1 --- GC2_2 --- GC3_1 --- GC4_1 --- GC4_2 --- GC4_3 --- GC5_1 ---------

        Every genome contains a set of syn clusters without duplications now. This process to split
        these gene clusters is as follows:

        Every gene cluster will be saved as a tuple in a multi layer dictionary containing the
        genomes and number of gene calls related to the gene cluster. If the number for at
        least one genome in the gene cluster is over one, the context of the gene cluster will
        be expanded iteratively.

        Iteration 1:

        (GC1, ): {G1: 1, G2: 1} --> OK
        (GC2, ): {G1: 2, G2: 2} --> not OK
        (GC3, ): {G1: 1, G2: 1} --> OK
        (GC4, ): {G1: 2, G2: 3} --> not OK
        (GC5, ): {G1: 2, G2: 1} --> not OK

        Iteratiob n:

        (GC1, ): {G1: 1, G2: 1}         --> OK
        (GC1, GC2, GC2): {G1: 1, G2: 1} --> OK
        (GC2, GC2, GC3): {G1: 1, G2: 1} --> OK
        (GC2, GC3, GC4): {G1: 1, G2: 1} --> OK
        (GC3, GC4, GC4): {G1: 1, G2: 1} --> OK
        (GC4, GC4, GC4): {G1: 0, G2: 1} --> OK
        (GC4, GC4, GC5): {G1: 1, G2: 1} --> OK
        (GC4, GC5, GC5): {G1: 1, G2: 0} --> OK
        (GC4, GC5,  +0): {G1: 0, G2: 1} --> OK
        (GC5, GC5,  +0): {G1: 1, G2: 0} --> OK

        As soon as all entries in the dictionary have a maximum number of genecalls per genome of one,
        some entries can be combined again to account for slightly different context e.g. early contig end.
        To combine a small tree calculation is used based on the similarity of the tuples and the condition
        of having one gene call per genome at most. The distance calculation is performed in context distance
        and the combination in context_split.

        (GC1, ): {G1: 1, G2: 1}         --> GC1_1
        (GC1, GC2, GC2): {G1: 1, G2: 1} --> GC2_1
        (GC2, GC2, GC3): {G1: 1, G2: 1} --> GC2_2
        (GC2, GC3, GC4): {G1: 1, G2: 1} --> GC3_1
        (GC3, GC4, GC4): {G1: 1, G2: 1} --> GC4_1
        (GC4, GC4, GC4): {G1: 0, G2: 1} --> GC4_2
        (GC4, GC4, GC5): {G1: 1, G2: 1} --> GC4_3
        (GC4, GC5, GC5): {G1: 1, G2: 0} |-> GC5_1
        (GC4, GC5,  +0): {G1: 0, G2: 1} |
        (GC5, GC5,  +0): {G1: 1, G2: 0} --> GC5_2

        Parameters
        ==========

        Returns
        =======
        pangenome_data_df: ps.DataFrame
            Dataframe containing gene call informations present in the anvi'o dbs
            plus the additional column of unipque syn clusters.
        """

        self.run.warning(None, header="Select paralog context", lc="green")

        gene_cluster_positions = {}
        gene_cluster_contig_order = {}
        gene_cluster_id_contig_positions = {}

        groups = pangenome_data_df.groupby(["genome", "contig"])
        ident = 0
        for name, group in groups:
            genome, contig = name

            self.contig_identifiers[contig] = str(ident)
            ident += 1

            group.reset_index(drop=False, inplace=True)
            group.sort_values('position', axis=0, ascending=True, inplace=True)
            gene_cluster_info = group[['gene_cluster', 'gene_caller_id']].values.tolist()
            gene_cluster_contig_order[(genome, contig)] = [gene_cluster for gene_cluster, gene_caller_id in gene_cluster_info]

            for i, (gene_cluster, gene_caller_id) in enumerate(gene_cluster_info):
                if not gene_cluster in gene_cluster_positions:
                    gene_cluster_positions[gene_cluster] = {(genome, contig): [(i, gene_caller_id)]}
                elif not (genome, contig) in gene_cluster_positions[gene_cluster]:
                    gene_cluster_positions[gene_cluster][(genome, contig)] = [(i, gene_caller_id)]
                else:
                    gene_cluster_positions[gene_cluster][(genome, contig)] += [(i, gene_caller_id)]


        gene_cluster_count = {}
        for (genome, contig), gene_cluster_sequence in gene_cluster_contig_order.items():
            for gene_cluster in gene_cluster_sequence:
                if gene_cluster not in gene_cluster_count:
                    gene_cluster_count[gene_cluster] = [genome]
                else:
                    gene_cluster_count[gene_cluster] += [genome]

        single_copy_core = set()
        for gene_cluster, count in gene_cluster_count.items():
            if len(count) == len(self.genome_names) and len(set(count)) == len(self.genome_names):
                single_copy_core.add(gene_cluster)


        j = 0
        self.progress.new("Contextualize gene clusters")
        for z, (gene_cluster, genome_contig_gene_cluster_positions) in enumerate(gene_cluster_positions.items()):
            k = self.min_k
            while True:
                gene_cluster_k_mer_genome_frequency = {}
                gene_cluster_k_mer_contig_positions = []

                for (genome, contig), values in genome_contig_gene_cluster_positions.items():
                    #contig_identifier = str(int(contig.split('_')[-1]))
                    contig_identifier = self.contig_identifiers[contig]
                    gene_cluster_order = gene_cluster_contig_order[(genome, contig)]

                    for position, gene_caller_id in values:

                        left_k_mer = [gene_cluster_order[l] if l >= 0 and l < len(gene_cluster_order) else '-' + contig_identifier for l in range(position - k, position)]
                        right_k_mer = [gene_cluster_order[l] if l >= 0 and l < len(gene_cluster_order) else '+' + contig_identifier for l in range(position, position + k + 1)]
                        gene_cluster_k_mer = tuple(left_k_mer + right_k_mer)

                        if not gene_cluster_k_mer in gene_cluster_k_mer_genome_frequency:
                            gene_cluster_k_mer_genome_frequency[gene_cluster_k_mer] = {genome: 1}
                        elif not genome in gene_cluster_k_mer_genome_frequency[gene_cluster_k_mer]:
                            gene_cluster_k_mer_genome_frequency[gene_cluster_k_mer][genome] = 1
                        else:
                            gene_cluster_k_mer_genome_frequency[gene_cluster_k_mer][genome] += 1

                        gene_cluster_k_mer_contig_positions += [(genome, contig, position, gene_caller_id, gene_cluster_k_mer)]

                for gene_cluster_k_mer, genome_frequency in gene_cluster_k_mer_genome_frequency.items():
                    if max(genome_frequency.values()) > 1:
                        k += 1
                        break

                else:
                    synteny_gene_cluster_id_contig_positions, synteny_gene_cluster_type = self.k_mer_split(gene_cluster, gene_cluster_k_mer_contig_positions, gene_cluster_contig_order, single_copy_core)

                    for genome, contig, position, gene_caller_id, gene_cluster_id in synteny_gene_cluster_id_contig_positions:
                        gene_cluster_id_contig_positions[j] = {'genome': genome, 'contig': contig, 'gene_caller_id': gene_caller_id, 'syn_cluster': gene_cluster_id, 'syn_cluster_type': synteny_gene_cluster_type}
                        j += 1
                    break
            self.progress.update(f"{str(z+1).rjust(len(str(len(gene_cluster_positions))), ' ')} / {len(gene_cluster_positions)}")

        self.progress.end()
        gene_cluster_id_contig_positions_df = pd.DataFrame.from_dict(gene_cluster_id_contig_positions, orient='index')
        pangenome_data_df = pangenome_data_df.merge(gene_cluster_id_contig_positions_df, on=['genome', 'contig', 'gene_caller_id'], how='inner')

        self.run.info_single(f'{len(pangenome_data_df)} gene caller entries.')
        self.run.info_single(f'{len(pangenome_data_df["gene_cluster"].unique())} gene cluster entries.')

        value_counts = pangenome_data_df["syn_cluster_type"].value_counts()
        self.run.info_single(f'{value_counts.get("core", 0)} core synteny gene caller entries.')
        self.run.info_single(f'{value_counts.get("duplication", 0)} paralog synteny gene caller entries.')
        self.run.info_single(f'{value_counts.get("rna", 0)} rna synteny gene caller entries.')
        self.run.info_single(f'{value_counts.get("rearrangement", 0)} rearranged synteny gene caller entries.')
        self.run.info_single(f'{value_counts.get("accessory", 0)} remaining accessory synteny gene caller entries.')
        self.run.info_single(f'{value_counts.get("singleton", 0)} singleton synteny gene caller entries.')
        self.run.info_single(f'{len(pangenome_data_df["syn_cluster"].unique())} synteny gene cluster entries in total.')
        self.run.info_single("Done.")

        if len(pangenome_data_df["syn_cluster"].unique()) > 2 * len(pangenome_data_df["gene_cluster"].unique()):
            if self.just_do_it:
                self.run.info_single("We wanted to inform you, that the number of gene to syn clusters doesn't really line up but you choosed just-do-it so either you hate us, yourself, or you are very sure about what you are doint. GOOD LUCK")
            else:
                raise ConfigError("We are sorry to inform you, that the number of gene to syn clusters doesn't really line up something might have gone wrong.")

        pangenome_data_df.set_index('position').to_csv(os.path.join(self.output_dir, 'synteny_cluster.tsv'), sep='\t')
        self.run.info_single(f"Exported mining table to {os.path.join(self.output_dir, 'synteny_cluster.tsv')}.")
        self.run.info_single("Done.")

        return(pangenome_data_df)
