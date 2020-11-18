import pandas as pd
import numpy as np
import glob
import os.path

from snakemake.shell import shell

# Import data
gap_counts_distribution_df = pd.read_csv(snakemake.input.tsv, \
									          sep="\t", \
									          index_col=None)

# FIXME: Here the workflow either uses a user defined num gaps threshold or the default
# value of 1.5*IQR to filter out sequences from the MSA. The problem is that the user
# defined value is currently fixed for every single MSA while the 1.5*IQR is modular for
# each MSA. The user should be able to provide a dict with the name of the MSA and the gap
# gap threshold specific to it i.e. {"Ribosome_L16": "3", "Ribosome_L16": "5"}
if snakemake.params.gap_threshold:

	outlier_threshold = int(snakemake.params.gap_threshold)

else:

	# Calculate outlier threshold
	Q1 = np.percentile(gap_counts_distribution_df.num_gaps, 25)
	Q3 = np.percentile(gap_counts_distribution_df.num_gaps, 75)
	IQR = Q3 - Q1

	outlier_threshold = Q3 + (IQR * 1.5)

if outlier_threshold == 0:
	print("There are no gaps in this alignment")

	shell("anvi-script-reformat-fasta {snakemake.input.fasta} \
                                  -o {snakemake.output}")
else:

	print("The outlier_threshold rule filter_out_outlier_sequences is: " + str(outlier_threshold))
	print("Here is the five num sum: ")
	print(gap_counts_distribution_df.num_gaps.describe())
	# Filter out sequences that have gaps > outlier_threshold
	shell("anvi-script-reformat-fasta {snakemake.input.fasta} \
	                                  -M %d \
	                                  -o {snakemake.output}" % (outlier_threshold))