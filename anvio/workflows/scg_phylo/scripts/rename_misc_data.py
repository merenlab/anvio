import pandas as pd
import glob
import os.path

from snakemake.shell import shell


# Import misc data with original headers
# ----------------------------
misc_data_with_original_headers = pd.read_csv(snakemake.input.taxonomy, \
                  sep="\t", \
                  index_col=None)


# Import reformat file with new names for tree calculation (output of anvi-script-reformat-fasta)
# ----------------------------
reformat_report = pd.read_csv(snakemake.input.reformat_report, \
									sep="\t", \
									engine='python', \
									names=["name_new", "header"] \
									)

reformat_report['gene_callers_id'] = reformat_report['header'].str.split("gene_callers_id:|\|start:", expand=True)[1].astype(int)
reformat_report['sample_name'] = reformat_report['header'].str.split("bin_id:|\-contigs", expand=True)[1]

reformat_report['header'] = reformat_report['sample_name'].map(str) + '_' + reformat_report['gene_callers_id'].map(str)

# Join tables to give the misc data the new shortened names for the tree calculation
# ----------------------------
misc_data = misc_data_with_original_headers.merge(reformat_report, on="header", how="inner")

# clean data
misc_data = misc_data[["name_new", "sample_name", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]]
misc_data.columns = ["name", "sample_name", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]


# print(misc_data_with_original_headers)
# print(reformat_report)
# print(misc_data)
# Export final table
# ----------------------------
misc_data.to_csv(snakemake.output.misc_data_final, \
					 sep="\t", \
					 index=None, \
					 na_rep="NA")
