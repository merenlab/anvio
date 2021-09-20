from Bio import SeqIO

import pandas as pd
import numpy as np
import glob
import os.path

# This script checks if there is a genomic SCG in a cluster
# and takes note of it in the misc data file

# remove excessive header lines in concatenated file first
lines = open(snakemake.input.misc_data).readlines()
new_content = [l for l in lines[1:] if not l.startswith('scg_name')]
open(snakemake.input.misc_data, "w").writelines([lines[0]] + new_content)

# Import
#-------
misc_data = pd.read_csv(snakemake.input.misc_data, \
                  sep="\t", \
                  index_col=False)

print(misc_data)
sys.exit()

cluster_rep_index = pd.read_csv(snakemake.params.cluster_rep_index, \
                  sep="\t", \
                  index_col=False, \
                  names=["representative", "cluster_members"])

fasta_df = pd.DataFrame({'header': [], 'sequence': []})

for seq_record in SeqIO.parse(snakemake.input.final_list_of_sequences_for_tree_calculation, "fasta"):
    fasta_df = fasta_df.append({'header': str(seq_record.description), 'sequence': str(seq_record.seq)}, ignore_index=True)


# Clean
# -----
cluster_rep_index_dict = cluster_rep_index.groupby('representative')['cluster_members'].apply(list).to_dict()
seq_in_tree_list = fasta_df.header.to_list()
misc_data['has_genomic_SCG_in_cluster'] = 'no'

# Make split name column for mapping
misc_data['split_name'] = misc_data['scg_name'].astype(str) + '_split_00001'
first_column = misc_data.pop('split_name')
misc_data.insert(0, 'split_name', first_column)

# Export
#-------
misc_data.to_csv(snakemake.output.misc_data_final, \
           sep="\t", \
           index=None, \
           na_rep="NA")
