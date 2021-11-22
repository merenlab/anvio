import pandas as pd

# This script will replace the external-gene-calls names with the simplified names from anvi-script-reformat-fasta
# sample-HMM-number

# Import tables
#--------------
external_gene_calls = pd.read_csv(snakemake.input.external_gene_calls, \
                                  delim_whitespace=True, \
                                  index_col=False)

reformat_report = pd.read_csv(snakemake.input.reformat_file, \
                              sep="\t", \
                              index_col=False, \
                              names=["new_header", "header"])

# Parse input files
#-----------------
# Parse reformat_report to get gene-callers-ids
reformat_report["gene_callers_id"] = reformat_report['header'].str.split("gene_callers_id:|\|start:", expand=True)[1].astype(str)

# Parse external-gene-calls contig column to get gene-callers-ids
external_gene_calls[['name', 'contig_number', 'gene_callers_id']] = external_gene_calls['contig'].str.rsplit('_',2, expand=True)

# Join external-gene-calls with reformat-report on gene-callers-id
#-----------------
external_gene_calls = external_gene_calls.merge(reformat_report, on="gene_callers_id", how="inner")

# Replace contig name from anvi-export-gene-calls with the simplified, new header from anvi-script-reformat-fasta
external_gene_calls = external_gene_calls[["gene_callers_id", "new_header", "start", "stop", "direction", "partial", "call_type", "source", "version", "aa_sequence"]]
external_gene_calls = external_gene_calls.rename(columns={'new_header': 'contig'})

# Write file
#-----------
external_gene_calls.to_csv(snakemake.output.external_gene_calls_renamed, \
                           sep="\t", \
                           index=False, \
                           header=True)
