import pandas as pd
import glob
import os.path

# Import SCG taxonomy results
# misc_data_with_old_names  = pd.read_csv(ribosomal_protein + "_protein/" + "03_" + ribosomal_protein +  "_scg_estimate_taxonomy_all.txt", \
# 				  sep="\t", \
# 				  index_col=None)

# print(snakemake.input.misc_data_all)

# Import misc data with original headers
# ----------------------------
misc_data_with_original_headers = pd.read_csv(snakemake.input.misc_data_all, \
                  sep="\t", \
                  index_col=None)


# Import reformat file with new names for tree calculation (output of anvi-script-reformat-fasta)
# ----------------------------
reformat_report = pd.read_csv(snakemake.input.reformat_report_all, \
									sep="\t", \
									engine='python', \
									names=["name_new", "header"] \
									)

# print(misc_data_with_original_headers)
# print(reformat_report)

# Join tables to give the misc data the new shortened names for the tree calculation
# ----------------------------
misc_data = misc_data_with_original_headers.merge(reformat_report, on="header", how="inner")

# clean data
misc_data = misc_data[["name_new", "sample_name", "Project", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]]
misc_data.columns = ["name", "sample_name", "Project", "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]

# Export final table
# ----------------------------
# print("output file should be: ", snakemake.output.misc_data_final)
misc_data.to_csv(snakemake.output.misc_data_final, \
					 sep="\t", \
					 index=None, \
					 na_rep="NA")
