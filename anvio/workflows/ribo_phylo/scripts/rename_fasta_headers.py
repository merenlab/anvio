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

taxonomy = pd.read_csv(snakemake.input.taxonomy, \
                  sep="\t", \
                  index_col=None)

print(snakemake.params.external_genomes)
# Clean
#------
taxonomy = taxonomy.rename(columns = {'bin_name': 'sequence_name'})

reformat_report['gene_callers_id'] = reformat_report['old_header'].str.split("gene_callers_id:|\|start:", expand=True)[1].astype(int)
reformat_report['sequence_name'] = snakemake.wildcards.ribosomal_protein_name + '_' + reformat_report['gene_callers_id'].map(str)
reformat_report['sample'] = snakemake.wildcards.sample_name
reformat_report['contig_db_type'] = np.where(reformat_report['sample'].isin(snakemake.params.external_genomes), 'external_genomes', 'metagenome')

misc_new_headers_df = reformat_report.merge(taxonomy, on="sequence_name", how="left")[['new_header', 'sample', 'contig_db_type', 't_domain', 't_phylum', 't_order', 't_family', 't_genus', 't_species']]

# Export
#------------------------------------------------------------------
misc_new_headers_df.to_csv(snakemake.output.misc_data_final, \
           sep="\t", \
           index=None, \
           na_rep="NA")
