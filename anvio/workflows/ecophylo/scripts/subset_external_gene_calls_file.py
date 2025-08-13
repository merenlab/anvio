import pandas as pd

# This script reformats the names from the external_gene_calls.txt file (which comes from anvi-get-sequences-for-gene-calls) with
# the names of sequences from anvi-estimate-scg-taxonomy

# Import tables
#--------------
external_gene_calls_all = pd.read_csv(snakemake.params.external_gene_calls_all, \
                                      delim_whitespace=True, \
                                      index_col=False)

headers = pd.read_csv(snakemake.input.headers, \
                      sep="\t", \
                      index_col=False,
                      names=["contig"])

# Make new columns
#-----------------
external_gene_calls_subset = external_gene_calls_all[external_gene_calls_all['contig'].isin(headers['contig'].tolist())]

# Make new gene-callers-ids starting from 1:x
#--------------------------------------------
external_gene_calls_subset['gene_callers_id'] = external_gene_calls_subset.index

# Write file
#-------------
external_gene_calls_subset.to_csv(snakemake.output.external_gene_calls_subset, \
                           sep="\t", \
                           index=False, \
                           header=True)
