* influence_ratio (integer)
this number is used to decide how many of the num_components (under
the 'general' section) are  going to be assigned for this file. if
it is not declared, the influence is left to be decided by the
software that handles the config object. 

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
