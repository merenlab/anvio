This program, as one might think, allows you to import a %(collection)s. This allows you to easily import any binning that you've already done into a %(profile-db)s, since the %(bin)ss within that collection will be carried over. 

This information (in the form of a %(collection-txt)s) can either come from another Anvi'o project (using %(anvi-export-collection)s) or you can get the coverage and sequence composion of your data using %(anvi-export-splits-and-coverages)s to bin your contigs with software outside of Anvi'o, then import that data into your database with this program. 

You can run this program like so: 

{{ codestart }}
anvi-import-collection -C my_bins.txt \
                        -p %(profile-db)s \
                        -c %(contigs-db)s 
{{ codestop }}

This will import the collection indicated in `my_bins.txt` into your %(profile-db)s. 

`my_bins.txt` should be a tab-delimited file where the first column lists a split name and the second lists the bin that it is placed in. You can see an example of this [here](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_splits.txt). 

You can also provide this information by listing your contigs instead of your splits (like [this](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_contigs.txt)). Just add the `--contigs-mode` tag. 

You can also provide an information file to describe the source and/or colors of your bins. [This file](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/example_bins_info_file.txt) is an example of such an information file. 


