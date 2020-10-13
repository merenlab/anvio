import pandas as pd
import numpy as np
import glob
import os.path

from snakemake.shell import shell


# Import data
gap_counts_distribution_df = pd.read_csv(snakemake.input.tsv, \
									          sep="\t", \
									          index_col=None)

# Here we need to add some logic where the user can decide to use the 1.5*IQR or any
# other value they choose for their gap filter threshold. It would be great to 
# get this value from the config file... for now 1.5*IQR is hard coded into the workflow

# if threshold == IQR:
# 	Calculate_IQR
# if threshold == user_defined_threshold:
# 	use_that

# Calculate outlier threshold
Q1 = np.percentile(gap_counts_distribution_df.num_gaps, 25)
Q3 = np.percentile(gap_counts_distribution_df.num_gaps, 75)
IQR = Q3 - Q1

outlier_threshold = IQR * 1.5

print("the outlier_threshold is: " + str(outlier_threshold))

# Filter out sequences that have gaps > outlier_threshold
shell("anvi-script-reformat-fasta {snakemake.input.fasta} \
                                  -M %d \
                                  -o {snakemake.output}" % (outlier_threshold))
