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


# Import import fasta as dataframe
#-------------------------------------------------------------------
fasta_df = pd.DataFrame({'new_header': [], 'sequence': []})

for seq_record in SeqIO.parse(snakemake.input.reps, "fasta"):
    fasta_df = fasta_df.append({'new_header': str(seq_record.description), 'sequence': str(seq_record.seq)}, ignore_index=True)


# Join represenative header with reformat report to get original gene-callers-ids for each contigsDB
fasta_reps_df = fasta_df.merge(reformat_report, on="new_header", how="inner")

# Export filtered reformat file
fasta_reps_df.to_csv(snakemake.output.report_file_filtered, \
					 sep="\t", \
					 index=None, \
					 na_rep="NA")
