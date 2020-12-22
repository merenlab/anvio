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
print(reformat_report)

taxonomy = pd.read_csv(snakemake.input.taxonomy, \
                  sep="\t", \
                  index_col=None)

taxonomy = taxonomy.rename(columns = {'scg_name': 'sequence_name'})
print(taxonomy)

reformat_report['gene_callers_id'] = reformat_report['old_header'].str.split("gene_callers_id:|\|start:", expand=True)[1].astype(int)

reformat_report['sequence_name'] = snakemake.wildcards.ribosomal_protein_name + '_' + reformat_report['gene_callers_id'].map(str)
reformat_report['sample'] = snakemake.wildcards.sample_name

# # Import import fasta 
# #-------------------------------------------------------------------
# fasta_df = pd.DataFrame({'old_header': [], 'sequence': []})

# for seq_record in SeqIO.parse(str(snakemake.input.faa), "fasta"):
#     # char_list = seq_record.description.split("_")
#     # header = "_".join(char_list[:-2]) + "_" + char_list[-1]
#     fasta_df = fasta_df.append({'old_header': str(seq_record.description), 'sequence': str(seq_record.seq)}, ignore_index=True)

# print("does tis exist")
# print(fasta_df)
# # # # Export
# # # #------------------------------------------------------------------
misc_new_headers_df = reformat_report.merge(taxonomy, on="sequence_name", how="left")[['new_header', 'sample', 't_domain', 't_phylum', 't_order', 't_family', 't_genus', 't_species']]
print(misc_new_headers_df)


misc_new_headers_df.to_csv(snakemake.output.misc_data_final, \
           sep="\t", \
           index=None, \
           na_rep="NA")
# # Export filtered fasta
# fasta = open(snakemake.output.fasta, 'w')


# for index, line in fasta_new_headers_df.iterrows():
#     header = ">" + line['new_header'] + "\n"
#     sequence = line['sequence'] + "\n"
#     # print(sequence)
#     fasta.write(header)
#     fasta.write(sequence)

# fasta.close()