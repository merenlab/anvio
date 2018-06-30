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


def compute_homogeneity_index(gene_cluster_sequences):
    #gene_cluster_sequences is a list of sequences from a single gene cluster
    #This only considers genomes that have genes in this cluster.
    num_sequences = len(gene_cluster_sequences)
    if num_sequences == 1:
        return 100
    elif num_sequences == 0:
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
                

