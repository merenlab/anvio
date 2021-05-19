import pandas as pd
import numpy as np
import glob
import os.path

from Bio import SeqIO
from snakemake.shell import shell

# This script extract sequence headers from fasta files

# Import import fasta as dataframe
fasta_df = pd.DataFrame({'header': [], 'sequence': []})

for seq_record in SeqIO.parse(str(snakemake.input.reps), "fasta"):
    fasta_df = fasta_df.append({'header': str(seq_record.description), 'sequence': str(seq_record.seq)}, ignore_index=True)

# Grab headers
headers = fasta_df.header

# Export filtered reformat file
headers.to_csv(snakemake.output.headers, \
					 sep="\t", \
					 index=None, \
					 na_rep="NA")
