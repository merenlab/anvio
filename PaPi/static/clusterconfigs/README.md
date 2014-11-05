* ratio (integer)
this is an optional parameter that can be set under matrices to 
specify how many of the num_components (under the 'general' section)
should be assigned for a given matrix. if there are more than one
matrices, leaving them blank will let the software handling the config
file to determine how ratios should be arranged. The default behavior
of PaPi for merging will be to use TFN and Coverage information. It will
increase the influence of Coverage with increasing number of samples.
For instance, if there are 16 samples, the influence of TNF will be
minimal, if there are two samples, the influence of coverage will be
minimal, etc.


* columns (column_name_1,...,column_name_n)

this variable defines which columns should be taken into account
for scaling. some files may contain multiple sources of information,
some of which may be irrelevant for scaling (such as a column of
taxanomical strings in the METADATA file for coverage of each contig.
the column names of interest (such as sample names in this case) can
be listed with this variable. If this does not exist, then all columns
that return True for PaPi.constants.IS_ESSENTIAL_FIELD is considered.

* skip_scaling (True|False)

If it is True, vectors would be clustered without any scaling step.
If there is no matrix combination, there is no need to scale data,
and clustering can be done on the full matrix. This is the default
behavior for ordering contigs based on tetranucleotide-frequency
or coverage *alone*. When skip_scaling is set to True, it is expected
to be only one matrix, and no num_components and ratio declerations
in the config file.

* seed (int)

Seed for reproducable results for testing purposes. If none defined,
no seed will be passed to scaling functions. You don't need it if you
don't know what it is.
