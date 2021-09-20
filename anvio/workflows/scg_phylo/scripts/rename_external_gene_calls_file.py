from Bio import SeqIO

import pandas as pd
import numpy as np
import glob
import os.path

# This script replaces the contig name of the external_gene_calls
# file with the reformated name.

# This script reformats the names from the external_gene_calls.txt file (which comes from anvi-get-sequences-for-gene-calls) with 
# the names of sequences from anvi-estimate-scg-taxonomy

# Import tables
#--------------
external_gene_calls = pd.read_csv(snakemake.input.external_gene_calls_all, \
                                  delim_whitespace=True, \
                                  index_col=False)


headers = pd.read_csv(snakemake.input.headers, \
                      sep="\t", \
                      index_col=False,
                      names=["contig"])

# Make new columns
#-----------------
external_gene_calls[['sample_name', 'contig_number', 'gene_callers_id']] = external_gene_calls['contig'].str.split('_',expand=True)
external_gene_calls['reference_protein_name'] = snakemake.wildcards.reference_protein_name
external_gene_calls['contig'] = external_gene_calls['sample_name'] + '_' + external_gene_calls['reference_protein_name'] + '_' + external_gene_calls['gene_callers_id']

# # Join reformat_report with external_gene_calls
#------------------------------------------------
external_gene_calls_filtered = external_gene_calls.merge(headers, on="contig", how="inner").drop(columns=['reference_protein_name', 'sample_name', 'contig_number', 'reference_protein_name'])


# # Write file
#-------------
external_gene_calls_filtered.to_csv(snakemake.output.external_gene_calls_renamed, \
                                    sep="\t", \
                                    index=False, \
                                    header=True, \
                                    na_rep="NA")