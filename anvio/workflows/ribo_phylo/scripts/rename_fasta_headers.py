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

# print(reformat_report)

# Import import fasta 
#-------------------------------------------------------------------
fasta_df = pd.DataFrame({'header': [], 'sequence': []})

for seq_record in SeqIO.parse(snakemake.input.fna, "fasta"):
    fasta_df = fasta_df.append({'header': snakemake.wildcards.sample_name + "_" + str(seq_record.description), 'sequence': str(seq_record.seq)}, ignore_index=True)

# print(fasta_df)

# # Export
# #------------------------------------------------------------------
fasta_new_headers_df = reformat_report.merge(fasta_df, on="header", how="inner")[['new_header', 'sequence']]

# print(fasta_new_headers_df)

# Export filtered fasta
fasta = open(snakemake.output.fasta, 'w')
# fasta = open(snakemake.output.fasta_filtered, 'w')


for index, line in fasta_new_headers_df.iterrows():
    header = ">" + line['new_header'] + "\n"
    sequence = line['sequence'] + "\n"

    fasta.write(header)
    fasta.write(sequence)

fasta.close()