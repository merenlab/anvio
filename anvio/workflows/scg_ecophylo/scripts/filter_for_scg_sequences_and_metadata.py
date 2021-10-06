from Bio import SeqIO

import pandas as pd
import numpy as np
import glob
import os.path

# # Import SCG taxonomy results
# #----------------------------
# SCG_taxonomy_results = pd.read_csv(snakemake.input.SCG_taxonomy, \
#                   sep="\t", \
#                   index_col=None)

# # Make gene-callers-id column
# SCG_taxonomy_results['gene_callers_id'] = SCG_taxonomy_results.scg_name.str.split("_",expand=True).iloc[:, -1].astype(int)

# # Export gene-callers-ids to extract SCG NT sequences
# SCG_taxonomy_results['gene_callers_id'].to_csv(snakemake.output.gene_callers_ids, \
#                                        sep="\t", \
#                                        index=None, \
#                                        header=False, \
#                                        na_rep="NA")

# Import import fasta from anvi-get-sequences-for-hmm-hits
#-------------------------------------------------------------------
fasta_df = pd.DataFrame({'header': [], 'sequence': []})

for seq_record in SeqIO.parse(snakemake.input.fasta, "fasta"):
    fasta_df = fasta_df.append({'header': str(seq_record.description), 'sequence': str(seq_record.seq)}, ignore_index=True)


headers = pd.read_csv(snakemake.input.headers, \
                  sep="\t", \
                  index_col=None)

print(headers)
exit()
# Export fasta headers of ribosomal proteins with SCG taxonomy hits
#------------------------------------------------------------------

# make gene_callers_id column for tab fasta
fasta_df['gene_callers_id'] = fasta_df['header'].str.split("gene_callers_id:|\|start:", expand=True)[1].astype(int)

# filter for sequences that have SCG taxonomy results
fasta_df_filtered = fasta_df[fasta_df['gene_callers_id'].isin(SCG_taxonomy_results['gene_callers_id'])]

# Export filtered fasta
fasta = open(snakemake.output.fasta_filtered, 'w')

for index, line in fasta_df_filtered.iterrows():
    header = ">" + line['header'] + "\n"
    sequence = line['sequence'] + "\n"

    fasta.write(header)
    fasta.write(sequence)

fasta.close()

# Create misc tree data
#----------------------

# Add sample column
fasta_df_filtered["sample_name"] = fasta_df_filtered['header'].str.split(" bin_id:|\|source:", expand=True)[1].astype(str)
fasta_df_filtered["sample_name"] = fasta_df_filtered["sample_name"].str.replace("-contigs", "")

# Add project column
fasta_df_filtered['Project'] = pd.np.where(fasta_df_filtered.sample_name.str.contains("ACE"), "ACE",
                               pd.np.where(fasta_df_filtered.sample_name.str.contains("TARA"), "TARA", "no_project"))

# Join with SCG taxonomy data
scg_sequences_misc_data = fasta_df_filtered.merge(SCG_taxonomy_results, on="gene_callers_id", how="inner")
scg_sequences_misc_data = scg_sequences_misc_data[["header", "sample_name", "Project", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]]

# Export misc data
scg_sequences_misc_data.to_csv(snakemake.output.misc_data, \
                                       sep="\t", \
                                       index=None, \
                                       header=False, \
                                       na_rep="NA")