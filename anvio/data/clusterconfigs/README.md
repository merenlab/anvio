The configuration file provides a recipe to anvio to mix multiple sources of information. For instance, when there are multiple metagenomic samples, reliance on coverage to order contigs is very efficient to determine the correct genome bins based on distribution patterns of contigs across samples. However when there are small number of samples from similar environments, the confidence of coverage values decrease quickly. In these cases tetra-nucleotide frequency information could be more useful than coverage. This creates a necessity to combine multiple sources of information. Also, the influence of the pieces of this information should be somehow maintained. For instance, if there are 3 samples, the contribution of TNF should be more compared to how much it would have contributed to the ordering of samples if there were 12 samples. anvio uses these configuration files to mix multiple sources of information by relying on dimension reduction through multidimensional scaling with a distance metric (euclidean, by default), and combining resulting normalized coordinats prior to clustering.

If it there is only one matrix, the scaling step will be ommitted, and the clustering would be done on the single full matrix with all or a subset of the columns specified by `columns_to_use` variable. When there is only one matrix to order contigs, declaring a `ratio` for the only matrix, or declaring the `num_components` variable under the general section would be irrelevant, and for the sake of clarity, the config handler class would raise an exception. This setup is appropriate for ordering contigs based on tetranucleotide-frequency or coverage *alone*. 

More information to come.

# Basic structure of the config file #

Following is a sample config file:

     [general]
     num_components = 16
     output_file = OUTPUT_TREE.txt
     seed = 42
     
     [TETRANUCLEOTIDE-FREQ-MATRIX.txt]
     alias = tnf
     ratio=2
     
     [METADATA-mean_coverage.txt]
     ratio=3
     alias = coverage
     columns_to_use = 204-6M,204-7M,204-9M

There will be more information here.


# Parameters #

* `alias` (string, single word): a one word descriptior of the matrix. This alias will be used for reporting, therefore it is important for it to be brief and accurate. 

* `ratio` (integer): this is an optional parameter that can be set under matrices to  specify how many of the `num_components` (under the 'general' section) should be assigned for a given matrix. if there are more than one matrices, leaving them blank will let the software handling the config file to determine how ratios should be arranged. The default behavior of anvio for merging will be to use TFN and Coverage information. It will increase the influence of Coverage with increasing number of samples. For instance, if there are 16 samples, the influence of TNF will be minimal, if there are two samples, the influence of coverage will be minimal, etc.

* `columns_to_use` (`column_name_1,...,column_name_n`)  this variable defines which columns should be taken into account for scaling. some files may contain multiple sources of information, some of which may be irrelevant for scaling (such as a column of taxanomical strings in the METADATA file for coverage of each contig. the column names of interest (such as sample names in this case) can be listed with this variable. If this does not exist, then all columns that return True for `anvio.constants.IS_ESSENTIAL_FIELD` is considered.  

* `seed` (int)  Seed for reproducable results for testing purposes. If none defined, no seed will be passed to scaling functions. You don't need it if you don't know what it is.   

