from Bio import SeqIO

import pandas as pd
import numpy as np
import glob
import os.path




# Import fasta header index
#----------------------------
reformat_report = pd.read_csv(snakemake.input.reformat_file_filtered, \
                  sep="\t", \
                  index_col=None)

# Make ribosomal protein column
reformat_report['ribosomal_protein'] = reformat_report['old_header'].str.split("___Bacteria_71", expand=True)[0]

# filter for rows with wildcards from ribosomal_protein_name and sample_name
gene_callers_ids_filtered = reformat_report.loc[(reformat_report['ribosomal_protein'] == snakemake.wildcards.ribosomal_protein_name) & (reformat_report['sample_name'] == snakemake.wildcards.sample_name)]['gene_callers_id']


# Export filter gene-caller-ids
gene_callers_ids_filtered.to_csv(snakemake.output.gene_callers_ids_reps, \
					 sep="\t", \
					 index=None, \
					 na_rep="NA")
