The output directory produced by %(anvi-summarize)s when run on a %(profile-db)s and %(contigs-db)s pair.

By default the directory is named `[PROJECT]-SUMMARY`. It requires a %(collection)s and provides a comprehensive statistical and sequence-level breakdown of every bin in that collection.

## Output files

### bins_summary.txt

One row per %(bin)s. Columns include bin name, taxon ID (if calculated), total nucleotides, total contigs, N50, GC content, completion, and redundancy.

### bin_by_bin/

A subdirectory per bin, each containing:

- A %(fasta)s file of the bin's contigs
- %(hmm-hits)s information
- Coverage, detection, and other read-recruitment statistics across each sample in the %(profile-db)s
- Domain and taxonomy predictions from single-copy core genes (see %(anvi-run-scg-taxonomy)s)

{:.notice}
In case you want to learn about the definitions of statistics like coverage, detection, abundance, variability, and so on, you should first read [Mike Lee's explanation of these statistics](https://merenlab.org/2017/05/08/anvio-views/). Our [vocabulary page](https://anvio.org/vocabulary/) might also be helpful. Then, keep in mind that anvi'o computes these values on a per-contig (and per-split) basis. When you run %(anvi-summarize)s, the program will summarize this information for a given bin by taking the average of a statistic's value across all splits in the bin, weighting that average by split length.

### bins_across_samples/

Tab-delimited matrix files compiling per-bin statistics across all samples — mean coverage, abundance, variability, and more. See [this post](https://merenlab.org/2017/05/08/anvio-views/) for definitions of these statistics.

### misc_data_layers/ and misc_data_items/

Any miscellaneous data imported into the database pair with %(anvi-import-misc-data)s, exported as %(misc-data-items-txt)s and %(misc-data-layers-txt)s files.

### index.html

An HTML document that formats all summary information for convenient browsing without an anvi'o installation.
