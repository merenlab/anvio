This file describes the contents of a %(collection)s in a readable format. This is used for taking collections in and out of Anvi'o, or transferring them between Anvi'o projects. 

This is a tab-delimited file that contains two columns. The first lists all of the contigs contained within the bins in this collection. The second lists the associated bin that that contig is placed in. For examples, check out [this page](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_contigs.txt) or below.

{{ codestart }}
contig_1    bin_1
contig_2    bin_1
contig_3    bin_1
contig_4    bin_2
contig_5    bin_3
contig_6    bin_3
{{ codestop }}

You can also list splits instead of contigs in the left column, as seen [here](https://github.com/merenlab/anvio/blob/master/anvio/tests/sandbox/example_files_for_external_binning_results/external_binning_of_splits.txt).
