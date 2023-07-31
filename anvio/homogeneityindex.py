# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    A class to compute homogeneity indices for a dictionary of gene clusters using a novel algorithm
"""

import anvio

import anvio.utils as utils
import anvio.terminal as terminal


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2019, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Mahmoud Yousef"
__email__ = "mahmoudyousef@uchicago.edu"


class HomogeneityCalculator(object):
    def __init__(self, quick_homogeneity=False):
        self.run = terminal.Run()
        self.progress = terminal.Progress()
        self.quick_homogeneity = quick_homogeneity
        self.functional = {}
        self.geometric = {}
        self.overall = {}


    def compute_functional_index(self, gene_cluster_sequences):
        """Given an array of aligned gene sequences of a gene cluster, computes the functional homogeneity index.
           Every amino acid residue of the same residual position in all genes is checked and assigned a similarity score
           based on how biochemically close the two residues are. Greater similarity scores indicate greater biochemical similarities
           between the two residues (as an extension, a greater functional index indicates greater biochemical similarities between the
           functional outcome of these genes)."""
        num_sequences = len(gene_cluster_sequences)
        if num_sequences == 1:
            return 1.0
        elif num_sequences == 0: #this is an error condition
            return 0.0

        num_residues = len(gene_cluster_sequences[0])
        similarity_score = 0
        max_score = 0
        for residue_number in range(0, num_residues):
            residues = []
            for gene_sequence in range(0, (num_sequences)-1):
                residues.append(gene_cluster_sequences[gene_sequence][residue_number])
            for spot_a in range(0,len(residues)):
                for spot_b in range(spot_a+1,len(residues)):
                    amino_acid_residue_1 = residues[spot_a]
                    amino_acid_residue_2 = residues[spot_b]

                    if amino_acid_residue_1 != "-" and amino_acid_residue_2 != "-":
                        max_score += 3
                        if amino_acid_residue_1 == amino_acid_residue_2 and (amino_acid_residue_1 != 'X' and amino_acid_residue_1 != 'J' \
                                                                            and amino_acid_residue_1 != 'B' and amino_acid_residue_1 != 'Z'):
                            similarity_score += 3
                        elif utils.is_amino_acid_functionally_conserved(amino_acid_residue_1,amino_acid_residue_2):
                            similarity_score += 2
                        #elif amino_acid_residue_2 != "-":
                            #similarity_score += 1
                            #max_score -= 3
                        else:
                            similarity_score += 0

        functional_index = similarity_score / max_score

        return functional_index


    def convert_sequences_to_binary_array(self, gene_sequences, bygene = False):
        """This function takes an array of aligned gene sequences of a gene cluster and converts it to a binary array.
           Gaps are represented by 1s. For efficiency purposes, arrays are stored as single-dimensional arrays of binary numbers
           and can be ordered either by residue (every array entry represents a residue position in all genes) or by gene (every
           array entry represents a gene)."""
        num_genes = len(gene_sequences)
        num_residues = len(gene_sequences[0])
        if not bygene: #array will be ordered by residue column
            array = [0] * num_residues
            for gene in range(num_genes):
                for residue in range(num_residues):
                    if gene_sequences[gene][residue] == "-":
                        array[residue] = ((array[residue]) << 1) + 1
                    else:
                        array[residue] = ((array[residue]) << 1) + 0
        else: #array will be ordered by gene, rather than by residue column
            array = [0] * num_genes
            for residue in range(num_residues):
                for gene in range(num_genes):
                    if gene_sequences[gene][residue] == "-":
                        array[gene] = ((array[gene]) << 1) + 1
                    else:
                        array[gene] = ((array[gene]) << 1) + 0

        return array
        #now we have an array of binary numbers that represent the gap-residue pattern of the entire gene cluster


    def compute_geometric_index(self, gene_cluster_sequences, quick_homogeneity=False):
        """This function calculates the geometric homogeneity index for a gene cluster. gene_cluster_sequences is an array of
           aligned gene sequences of that gene cluster.
           This function will compute vertical (residue-level) homogeneity as well as horizontal (gene-level) homogeneity by default.
           Adding to the arguments list 'quick_homogeneity=True' will skip horizontal homogeneity calculations. This will be faster, but
           the resulting geometric index will not factor in the spread of alignment gaps across genes"""
        num_genes = len(gene_cluster_sequences)
        if num_genes == 1:
            return 1.0

        binary_matrix_by_residue = self.convert_sequences_to_binary_array(gene_cluster_sequences)
        num_residues = len(binary_matrix_by_residue)
        residue_uniformity = []

        for col in range(num_residues):
            differences = []
            for counter in range(num_residues):
                if col == counter:
                    continue
                diff = bin(binary_matrix_by_residue[col] ^ binary_matrix_by_residue[counter])[2:].zfill(num_genes) #this has been converted to a string
                number_of_similarities = diff[-num_genes:].count('0')
                #Let's explain what happened - we converted diff into a binary number. .zfill ensured that we have the sufficient number of digits.
                #This line of code creates an array - similariies - containing all of the digits that are 0 (equal to each other in the ^ expresson)
                #number_of_similarities is the length of that array - effectively, the number of 0s in the binary number
                differences.append(number_of_similarities / num_genes)
            residue_uniformity.append(sum(differences) / len(differences))

        by_residue = sum(residue_uniformity) / len(residue_uniformity)

        if quick_homogeneity:
            return by_residue

        binary_matrix_by_gene = self.convert_sequences_to_binary_array(gene_cluster_sequences, bygene=True)
        gene_uniformity = []

        for gene in range(num_genes):
            differences = []
            for counter in range(num_genes):
                if gene == counter:
                    continue
                diff = bin(binary_matrix_by_gene[gene] ^ binary_matrix_by_gene[counter])[2:].zfill(num_residues)
                number_of_similarities = diff[-num_residues:].count('0')
                differences.append(number_of_similarities / num_residues)
            gene_uniformity.append(sum(differences) / len(differences))

        by_gene = sum(gene_uniformity) / len(gene_uniformity)

        geometric_index = (by_residue + by_gene) / 2

        return geometric_index


    def get_homogeneity_dicts(self, gene_clusters_dict):
        """ The main function called by dbops.PanSuperClass. It retrieves the gene clusters dictionary passed to
            the HomogeneityCalculator intiatior and calculates functional and geometric indices for each.
            This function returns three dictionaries - functional, geometric, and a representation of the average -
            with the following structure:
                {gene_cluster_id_1: index_1,
                 gene_cluster_id_2: index_2,
                 etc...
                 }
            Note that this function assumes that all gene sequences have been aligned properly"""

        for gene_cluster in gene_clusters_dict:
            cluster_sequences = []
            genes_in_cluster = gene_clusters_dict[gene_cluster]

            for genome_name in genes_in_cluster:
                gene_caller_ids = genes_in_cluster[genome_name]
                for gene_caller_id in gene_caller_ids:
                    cluster_sequences.append(gene_caller_ids[gene_caller_id])

            fun = self.compute_functional_index(cluster_sequences)
            geo = self.compute_geometric_index(cluster_sequences, self.quick_homogeneity)
            self.functional[gene_cluster] = fun
            self.geometric[gene_cluster] = geo
            if fun == 0 and geo == 0:
                self.overall[gene_cluster] = 0
            else:
                self.overall[gene_cluster] = 2 * fun * geo / (fun + geo)

        return self.functional, self.geometric, self.overall
