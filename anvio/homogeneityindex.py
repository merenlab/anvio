"""
homogeneityindex.py: code for determining homogeneity indices

The Homogeneity Indices calculate the relative homogeneity of all of the genomes in each gene cluster.
It is divided into functional and geometric homogeneity.

Author: Mahmoud Yousef
"""

conserved_groups = {
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

residue_id = {}
for key in ['A','I','L','V','M','C']:
    residue_id[key] = 'Nonpolar'
for key in ['F','W']:
    residue_id[key] = 'Aromatic'
for key in ['K','R']:
    residue_id[key] = 'Bases'
for key in ['Q', 'N']:
    residue_id[key] = 'Neutral Amines'
for key in ['D','E']:
    residue_id[key] = 'Acids'
for key in ['H','Y']:
    residue_id[key] = 'Polar and Nonpolar'
for key in ['S','T']:
    residue_id[key] = 'Mostly nonpolar'
for key in ['G','P']:
    residue_id[key] = 'None'
for key in ['X']:
    residue_id[key] = 'None'
residue_id['B'] = 'B'
residue_id['Z'] = 'Z'
residue_id['J'] = 'J'


def is_conserved(r1, r2):
    group = residue_id[r1]
    conserved_group = conserved_groups[group]

    if r2 in conserved_group:
        return True
    if group == 'Polar and Nonpolar': #they fall in more than one group, multiple tests needed
        if r1 == 'H' and (r2 in conserved_groups['Nonpolar'] or r2 in conserved_groups['Bases']):
            return True
        if r1 == 'Y' and (r2 in conserved_groups['Aromatic']):
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
        
        length = len(gene_cluster_sequences[0])
        similarity = 0
        max = 0
        for spot in range(0, length):
            residues = []
            for sequence in range(0, num_sequences):
                residues.append(gene_cluster_sequences[sequence][spot])
            for a in range(0,len(residues)):
                for b in range(0,len(residues)):
                    if a >= b:
                        continue
                    r1 = residues[a]
                    r2 = residues[b]
                    if r1 == "-":
                        r1 = residues[b]
                        r2 = "-"
                    max += 3
                    if r1 == r2 and (r1 != 'X' and r1 != 'J' and r1 != 'B' and r1 != 'Z'):
                        similarity += 3
                    elif is_conserved(r1,r2):
                        similarity += 2
                    elif r2 != "-":
                        similarity += 1
                    else:
                        similarity += 0
        index = (similarity / max) * 100
        return index


    def label_gaps(self, gene_sequences, bygene = False): #1 indicates gaps
        if bygene == False:
            array = [0] * len(gene_sequences[0])
            for i in range(len(gene_sequences)):
                for j in range(len(gene_sequences[0])):
                    if gene_sequences[i][j] == "-":
                        array[j] = ((array[j]) << 1) + 1
                    else:
                        array[j] = ((array[j]) << 1) + 0
        else: #array will be ordered by gene, rather than by column
            array = [0] * len(gene_sequences)
            for i in range(len(gene_sequences[0])):
                for j in range(len(gene_sequences)):
                    if gene_sequences[j][i] == "-":
                        array[j] = ((array[j]) << 1) + 1
                    else:
                        array[j] = ((array[j]) << 1) + 0

        return array
        #now we have an array of binary numbers that represent the gap-residue pattern of each column


    @profile
    def compute_geometric_index(self, gene_cluster_sequences, quick_homogeneity=False): 
        num_genes = len(gene_cluster_sequences)
        if num_genes == 1:
            return 100
        
        grid = self.label_gaps(gene_cluster_sequences)
        length = len(grid)
        residue_uniformity = []

        for col in range(length):
            differences = []
            for counter in range(length):
                if col == counter:
                    continue
                diff = grid[col] ^ grid[counter]
                runsum = len([ones for ones in bin(diff).zfill(num_genes)[-num_genes:] if ones=='0'])
                #runsum = 0
                #for i in range(num_genes):
                #    if abs(diff) % 2 == 0:
                #        runsum += 1
                #    diff = diff >> 1
                differences.append(runsum / num_genes)
            residue_uniformity.append(sum(differences) / len(differences))
        
        by_residue = 100 * (sum(residue_uniformity) / len(residue_uniformity))

        if quick_homogeneity:
            return by_residue

        grid2 = self.label_gaps(gene_cluster_sequences, bygene=True)
        length = len(grid2)
        gene_uniformity = []

        for gene in range(length):
            differences = []
            for counter in range(length):
                if gene == counter:
                    continue
                diff = grid2[gene] ^ grid2[counter]
                runsum = len([ones for ones in bin(diff).zfill(len(grid))[-len(grid):] if ones=='0'])
                #runsum = 0
                #for i in range(len(gene_cluster_sequences[0])):
                #    if abs(diff) % 2 == 0:
                #        runsum += 1
                #    diff = diff >> 1
                differences.append(runsum / len(grid))
            gene_uniformity.append(sum(differences) / len(differences))

        by_gene = 100 * (sum(gene_uniformity) / len(gene_uniformity))

        index = (by_residue + by_gene) / 2
        return index


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
            #self.cluster_sizes[gene_cluster] = len(cluster_sequences)
        
        return self.functional, self.geometric
