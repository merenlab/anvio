This file contains information about your genes. 

It is a tab-delimited text file where each row represents a specific gene and each column provides different information. 

As of now, the only program that returns data in this format is %(anvi-script-gen_stats_for_single_copy_genes.py)s, which returns this information for the single copy core genes in your %(contigs-db)s. 

From left to right, these tell you 
* The source for this gene (ex `Protista_83`)
* The name of the contig that this gene is a part of
* The gene name 
* The e-value (of the HMM hit that was used to find this gene)
