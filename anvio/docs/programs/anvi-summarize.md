Anvi-summarize lets you look at a **comprehensive overview of your %(collection)s** and its many statistics that anvi'o has calculated. 

It will create a folder (by default called `SUMMARY`) that contains many different summary files, including an HTML output that conviently displays them all for you. The exact contents of this folder will depend on whether you run the program on a %(profile-db)s (i.e., to summarize a %(collection)s of binned contigs, such as metagenome-assembled genomes) or on a %(pan-db)s (i.e., to summarize a %(collection)s of binned gene clusters, such as when you want to compare accessory vs core genome). Due to the extensive set of output files it produces, this program can be useful for sharing information with collaborators, generating supplementary files for manuscripts, and exporting data for use as input to downstream programs/scripts.

## Output files

Regardless of input type, this program always produces an `index.html` file, which you can open in a web browser to view all the summary information in a nicely-formatted interactive webpage.

When run on a %(profile-db)s, this program will:
* produce an overall table of bin statistics (`bins_summary.txt`) like length, GC content, completion and redundancy
* generate a per-bin folder of bin-specific information (in a directory called `bin_by_bin`), including:
    * %(fasta)s files of their contigs
    * information about their %(hmm-hits)s
    * coverage, detection, and other read-recruitment statistics of the bin in each sample stored in the %(profile-db)s
    * domain and taxonomy predictions from single-copy core genes (see %(anvi-run-scg-taxonomy)s)
    * the bin-specific value of each statistic described in `bins_summary.txt`, like length, percent completeness, and redundancy
* generate various tab-delimited matrix files compiling information about all bins across all samples (in a directory called `bins_across_samples`), including read-recruitment statistics and number of Ribosomal RNA annotations per bin (the rRNA info is not described across samples, but happens to live with the other matrix files regardless)
* if you have imported any miscellaneous data into the %(profile-db)s with %(anvi-import-misc-data)s, this information will also be exported (in the directories `misc_data_items` and `misc_data_layers`)

**Confused about the read-recruitment statistics?**

In case you want to learn about the definitions of statistics like coverage, detection, abundance, variability, and so on, you should first read [Mike Lee's explanation of these statistics](https://merenlab.org/2017/05/08/anvio-views/). Our [vocabulary page](https://anvio.org/vocabulary/) might also be helpful. Then, keep in mind that anvi'o computes these values on a per-contig (and per-split) basis. When you run %(anvi-summarize)s, the program will summarize this information for a given bin by taking the average of a statistic's value across all splits in the bin, weighting that average by split length.

When run on a %(pan-db)s, this program will:
* generate a _huge_ table (`[NAME]_gene_clusters_summary.txt`) describing every gene in every gene cluster of your pangenome (even those not in the specified %(collection)s), including:
    * gene-cluster-specific information like the number of genomes contributing to that cluster, maximum number of paralogs in any partipating genome, and cluster homogeneity metrics. 
    * gene-specific information like functional annotations and amino acid sequence
* export any tables of miscellaneous data that were imported into the %(pan-db)s with %(anvi-import-misc-data)s (in the directories `misc_data_items` and `misc_data_layers`)

## Running anvi-summarize 

### Running on a profile database

A standard run of anvi-summarize on a %(profile-db)s will look something like this:

{{ codestart }}
anvi-summarize -c %(contigs-db)s \
               -p %(profile-db)s \
               -o MY_SUMMARY \
               -C %(collection)s
{{ codestop }}

This will name the output directory `MY_SUMMARY` instead of the standard `SUMMARY`. 

When running on a profile database, you also have options to 
* output very accurate (but intensely processed) coverage and detection data for each gene (using `--init-gene-coverages`)
* edit your contig names so that they contain the name of the bin that the contig is in (using `--reformat-contig-names`)
* also display the amino acid sequeunces for your gene calls.  (using `--report-aa-seqs-for-gene-calls`)

### Running on a pan database

When running on a %(pan-db)s, you'll want to instead provide the associated genomes storage database. 

{{ codestart }}
anvi-summarize -g %(genomes-storage-db)s \
               -p %(pan-db)s \
               -C %(collection)s 
{{ codestop }}

You can also choose to display DNA sequences for your gene clusters instead of amino acid sequences with the flag `--report-DNA-sequences`

### Other notes

If you're unsure what collections are in your database, you can run this program with the flag `--list-collections` or by running %(anvi-show-collections-and-bins)s. Don't have a collection at all? If you want to summarize everything in the database, you can generate a default collection of everything by running %(anvi-script-add-default-collection)s.

You can also use the flag `--quick-summary` to get a less comprehensive summary with a much shorter processing time. 
