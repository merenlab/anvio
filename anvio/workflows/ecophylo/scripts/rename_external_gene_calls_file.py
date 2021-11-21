import pandas as pd

# This script reformats the names from the external_gene_calls.txt file (which comes from anvi-get-sequences-for-gene-calls) with
# the names of sequences from anvi-estimate-scg-taxonomy

# Import tables
#--------------
external_gene_calls = pd.read_csv(snakemake.input.external_gene_calls, \
                                  delim_whitespace=True, \
                                  index_col=False)

# Make new columns
#-----------------
external_gene_calls[['name', 'contig_number', 'gene_callers_id']] = external_gene_calls['contig'].str.rsplit('_',2, expand=True)
external_gene_calls['HMM'] = snakemake.wildcards.HMM
external_gene_calls['sample_name'] = snakemake.wildcards.sample_name
external_gene_calls['gene_callers_id_leading_zeros'] = external_gene_calls['gene_callers_id'].apply(lambda x: x.zfill(12))
external_gene_calls['contig'] = external_gene_calls['sample_name'] + '_' + external_gene_calls['HMM'] + '_' + external_gene_calls['gene_callers_id_leading_zeros']

external_gene_calls = external_gene_calls[["gene_callers_id", "contig", "start", "stop", "direction", "partial", "call_type", "source", "version", "aa_sequence"]]

# Write file
#-----------
external_gene_calls.to_csv(snakemake.output.external_gene_calls_renamed, \
                           sep="\t", \
                           index=False, \
                           header=True)
