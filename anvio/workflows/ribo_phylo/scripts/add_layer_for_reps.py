from Bio import SeqIO

import pandas as pd
import numpy as np
import glob
import os.path

# This script determines if an SCG derived from an isolate genomes is
# clustered with other metagenome derived SCGs then takes note of it in
# the misc data file.

# Import
#-------
misc_data = pd.read_csv(snakemake.input.misc_data, \
                  sep="\t", \
                  index_col=False)

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

sup_list = []
for seq in seq_in_tree_list:
  cluster_members_list = cluster_rep_index_dict[seq]
  for external_genome in snakemake.params.external_genomes:
    check = any(external_genome in s for s in cluster_members_list)
    if check is True:
      sup_list.append(seq)

misc_data['has_genomic_SCG_in_cluster'] = np.where(misc_data['new_header'].isin(sup_list), 'yes', 'no')

# Make split name column for mapping
misc_data['split_name'] = misc_data['new_header'].astype(str) + '_split_00001'
first_column = misc_data.pop('split_name')
misc_data.insert(0, 'split_name', first_column)

# Export
#-------
misc_data.to_csv(snakemake.output.misc_data_final, \
           sep="\t", \
           index=None, \
           na_rep="NA")
