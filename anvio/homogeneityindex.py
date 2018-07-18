"""
homogeneityindex.py: code for determining homogeneity index

Homogeneity index is a numerical score between 0 and 100 that indicates the relative homogeneity of genes within a gene cluster

Author: Mahmoud Yousef
"""

conserved_groups = {
    'Nonpolar': ['L','V','I','M','C','H','P'],
    'Aromatic': ['F','W'],
    'Positive': ['K','R'],
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
    residue_id[key] = 'Positive'
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
    return False

def compute_functional_index(gene_cluster_sequences):
    #gene_cluster_sequences is a list of sequences from a single gene cluster
    num_sequences = len(gene_cluster_sequences)
    if num_sequences == 1: 
        return 100
    elif num_sequences == 0: #this is an error condition
        return 0
    
    length = len(gene_cluster_sequences[0])
    similarity = 0
    max = 0
    for j in range(0, length):
        residues = []
        for k in range(0, num_sequences):
            residues.append(gene_cluster_sequences[k][j])
        for a in range(0,len(residues)):
            for b in range(0,len(residues)):
                if a >= b:
                    continue
                r1 = residues[a]
                r2 = residues[b]
                if r1 == "-":
                    max += 1
                    if r2 == "-":
                        similarity += 1
                    else:
                        similarity += 0
                else:
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

def label_gaps(gene_sequences): #1 indicates gaps
    matrix = []
    for row in gene_sequences:
        array = []
        for residue in row:
            if residue == "-":
                array.append("1")
            else:
                array.append("0")
        matrix.append(array)

    return matrix

def find_state_change_probabilities(column):
    gap_to_gap = 0
    gap_to_residue = 0
    residue_to_gap = 0
    residue_to_residue = 0
    total = len(column)
    if total == 1:
        return 0 #singletons will have 0 functional homogeneity
    count = 0

    for index1 in range(total - 1): #This entropy check is faulty
        r1 = column[index1]
        for index2 in range(index1, (total - 1)):
            r2 = column[index2 + 1]
            count += 1
            if r1 == "1":
                if r2 == "1":
                    gap_to_gap += 1
                else:
                    gap_to_residue += 1
            else:
                if r2 == "1":
                    residue_to_gap += 1
                else:
                    residue_to_residue += 1

    prob_gap_to_gap = gap_to_gap / (total - 1)
    prob_gap_to_residue = gap_to_residue / (total - 1)
    prob_residue_to_gap = residue_to_gap / (total - 1)
    prob_residue_to_residue = residue_to_residue / (total - 1)

    return prob_gap_to_residue + prob_residue_to_gap

def compute_structural_index(gene_cluster_sequences):
    grid = label_gaps(gene_cluster_sequences)

    length = len(grid[0])
    entropy = []

    for col in range(length):
        column = []
        for row in range(len(gene_cluster_sequences)):
            column.append(grid[row][col])
        entropy.append(find_state_change_probabilities(column))
    
    return 100 -((sum(entropy) / len(entropy)) * 100)


def compute_homogeneity_index(gene_cluster_sequences, num_genomes, num_unique): #will expand on this soon
    functional = compute_functional_index(gene_cluster_sequences)

    diff = num_genomes - num_unique #Can I use this to modify the functional index? Should we ignore paralogs 

    structural = compute_structural_index(gene_cluster_sequences)
    print(structural)
    return functional 

