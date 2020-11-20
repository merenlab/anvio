from Bio import SeqIO

import pandas as pd
import numpy as np
import glob
import os.path


# Import fasta header index
#----------------------------
reformat_report = pd.read_csv(snakemake.input.report_file, \
                  sep="\t", \
                  index_col=None,
                  names=['new_header', 'old_header'])

reformat_report['gene_callers_id'] = reformat_report['old_header'].str.split("gene_callers_id:|\|start:", expand=True)[1].astype(int)
reformat_report['sample_name'] = reformat_report['old_header'].str.split("bin_id:|\-contigs", expand=True)[1]

reformat_report['header'] = reformat_report['sample_name'].map(str) + '_' + reformat_report['gene_callers_id'].map(str)


# Add new names to external_gene_calls files
external_gene_calls = pd.read_csv(snakemake.input.external_gene_calls, \
                  sep="\t", \
                  index_col=None)


external_gene_calls['sample'] = external_gene_calls['contig'].str.rsplit('_', n=2).str.get(0)
external_gene_calls['header'] = external_gene_calls['sample'] + "_" + external_gene_calls['gene_callers_id'].astype(str)
external_gene_calls_filtered = reformat_report[['new_header', 'header']].merge(external_gene_calls, on="header", how="inner").drop(columns=['sample', 'header', 'contig']).rename(columns={'new_header': 'contig'})
external_gene_calls_filtered = external_gene_calls_filtered[['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version']]

external_gene_calls_filtered.to_csv(snakemake.output.external_gene_calls, \
           sep="\t", \
           index=None, \
           na_rep="NA")

# Import import fasta 
#-------------------------------------------------------------------
fasta_df = pd.DataFrame({'header': [], 'sequence': []})

for seq_record in SeqIO.parse(snakemake.input.fna, "fasta"):
    char_list = seq_record.description.split("_")
    header = "_".join(char_list[:6]) + "_" + char_list[7]
    fasta_df = fasta_df.append({'header': header, 'sequence': str(seq_record.seq)}, ignore_index=True)

# # Export
# #------------------------------------------------------------------
fasta_new_headers_df = reformat_report.merge(fasta_df, on="header", how="inner")[['new_header', 'sequence']]

# Export filtered fasta
fasta = open(snakemake.output.fasta, 'w')


for index, line in fasta_new_headers_df.iterrows():
    header = ">" + line['new_header'] + "\n"
    sequence = line['sequence'] + "\n"
    # print(sequence)
    fasta.write(header)
    fasta.write(sequence)

fasta.close()