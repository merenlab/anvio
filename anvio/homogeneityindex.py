"""
homogeneityindex.py: code for determining homogeneity indices

The Homogeneity Indices calculate the relative homogeneity of all of the genomes in each gene cluster.
It is divided into functional and geometric homogeneity.

Author: Mahmoud Yousef
"""

conserved_amino_acid_groups = {
    'Nonpolar': ['L','V','I','M','C','H','A'],
    'Aromatic': ['F','W','Y'],
    'Bases': ['K','R','H'],
    'Neutral Amines': ['Q, N'],
    'Acids': ['D','E'],
    'Polar and Nonpolar': ['H','Y'],
    'Mostly nonpolar': ['S','T'],
    'B': ['B','N','D'],
    'Z': ['Z','Q','E'],
    'J': ['J','L','I'],
    'None': []
}

amino_acid_property_group = {}
for key in ['A','I','L','V','M','C']:
    amino_acid_property_group[key] = 'Nonpolar'
for key in ['F','W']:
    amino_acid_property_group[key] = 'Aromatic'
for key in ['K','R']:
    amino_acid_property_group[key] = 'Bases'
for key in ['Q', 'N']:
    amino_acid_property_group[key] = 'Neutral Amines'
for key in ['D','E']:
    amino_acid_property_group[key] = 'Acids'
for key in ['H','Y']:
    amino_acid_property_group[key] = 'Polar and Nonpolar'
for key in ['S','T']:
    amino_acid_property_group[key] = 'Mostly nonpolar'
for key in ['G','P','X']:
    amino_acid_property_group[key] = 'None'
amino_acid_property_group['B'] = 'B'
amino_acid_property_group['Z'] = 'Z'
amino_acid_property_group['J'] = 'J'


def is_amino_acid_functionally_conserved(amino_acid_residue_1, amino_acid_residue_2):
    group = amino_acid_property_group[amino_acid_residue_1]
    conserved_group = conserved_amino_acid_groups[group]

    if amino_acid_residue_2 in conserved_group:
        return True
    if group == 'Polar and Nonpolar': #they fall in more than one group, multiple tests needed
        if amino_acid_residue_1 == 'H' and (amino_acid_residue_2 in conserved_amino_acid_groups['Nonpolar'] \
                                            or amino_acid_residue_2 in conserved_amino_acid_groups['Bases']):
            return True
        if amino_acid_residue_1 == 'Y' and (amino_acid_residue_2 in conserved_amino_acid_groups['Aromatic']):
            return True
    return False


class HomogeneityCalculator(object):
    def __init__(self, gene_clusters_dict, quick_homogeneity=False):
        self.gene_clusters_dict = gene_clusters_dict
        self.quick_homogeneity = quick_homogeneity
        self.functional = {}
        self.geometric = {}

        self.cluster_sizes = {}
        self.total_functional = 0
        self.total_geometric = 0
        #the above are foundations for potential future capabilities


    def compute_functional_index(self, gene_cluster_sequences):
        #gene_cluster_sequences is a list of sequences from a single gene cluster
        num_sequences = len(gene_cluster_sequences)
        if num_sequences == 1: 
            return 100
        elif num_sequences == 0: #this is an error condition
            return 0
        
        num_residues = len(gene_cluster_sequences[0])
        similarity_score = 0
        max_score = 0
        for residue_number in range(0, num_residues):
            residues = []
            for gene_sequence in range(0, num_sequences):
                residues.append(gene_cluster_sequences[gene_sequence][residue_number])
            for spot_a in range(0,len(residues)):
                for spot_b in range(spot_a+1,len(residues)):
                    amino_acid_residue_1 = residues[spot_a]
                    amino_acid_residue_2 = residues[spot_b]
                    if amino_acid_residue_1 == "-":
                        amino_acid_residue_1 = residues[spot_b]
                        amino_acid_residue_2 = "-"
                    max_score += 3
                    if amino_acid_residue_1 == amino_acid_residue_2 and (amino_acid_residue_1 != 'X' and amino_acid_residue_1 != 'J' \
                                                                        and amino_acid_residue_1 != 'B' and amino_acid_residue_1 != 'Z'):
                        similarity_score += 3
                    elif is_amino_acid_functionally_conserved(amino_acid_residue_1,amino_acid_residue_2):
                        similarity_score += 2
                    elif amino_acid_residue_2 != "-":
                        similarity_score += 1
                    else:
                        similarity_score += 0
        functional_index = (similarity_score / max_score) * 100
        return functional_index


    def convert_sequences_to_binary_array(self, gene_sequences, bygene = False): #1 indicates gaps
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
        num_genes = len(gene_cluster_sequences)
        if num_genes == 1:
            return 100
        
        binary_matrix_by_residue = self.convert_sequences_to_binary_array(gene_cluster_sequences)
        num_residues = len(binary_matrix_by_residue)
        residue_uniformity = []

        for col in range(num_residues):
            differences = []
            for counter in range(num_residues):
                if col == counter:
                    continue
                diff = bin(binary_matrix_by_residue[col] ^ binary_matrix_by_residue[counter])[2:].zfill(num_genes)
                number_of_similarities = diff[-num_genes:].count('0')
                #Let's explain what happened - we converted diff into a binary number. .zfill ensured that we have the sufficient number of digits.
                #This line of code creates an array - similariies - containing all of the digits that are 0 (equal to each other in the ^ expresson)
                #number_of_similarities is the length of that array - effectively, the number of 0s in the binary number
                differences.append(number_of_similarities / num_genes)
            residue_uniformity.append(sum(differences) / len(differences))
        
        by_residue = 100 * (sum(residue_uniformity) / len(residue_uniformity))

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

        by_gene = 100 * (sum(gene_uniformity) / len(gene_uniformity))

        geometric_index = (by_residue + by_gene) / 2
        return geometric_index


    def get_homogeneity_dicts(self):
        """desc"""

        for gene_cluster in self.gene_clusters_dict:
            cluster_sequences = []
            genes_in_cluster = self.gene_clusters_dict[gene_cluster]

            for genome_name in genes_in_cluster:
                gene_caller_ids = genes_in_cluster[genome_name]
                for gene_caller_id in gene_caller_ids:
                    cluster_sequences.append(gene_caller_ids[gene_caller_id])

            self.functional[gene_cluster] = self.compute_functional_index(cluster_sequences)
            self.geometric[gene_cluster] = self.compute_geometric_index(cluster_sequences, self.quick_homogeneity)
        
        return self.functional, self.geometric
