Anvi-summarize lets you export a **comprehensive overview of your data** from an anvi'o database. Depending on the input, it can summarize a %(collection)s of binned contigs (from a %(profile-db)s), a %(collection)s of binned gene clusters (from a %(pan-db)s), or the full contents of a %(pan-graph-db)s. The output is a directory of flat files and an HTML index that conveniently displays them for you. This makes the program useful for sharing information with collaborators, generating supplementary files for manuscripts, and exporting data for use in downstream analyses.

## Output files

### When run on a profile database

Running on a %(profile-db)s produces:

* An overall table of bin statistics (`bins_summary.txt`) like length, GC content, completion and redundancy
* A per-bin folder of bin-specific information (`bin_by_bin/`), including:
    * %(fasta)s files of their contigs
    * information about their %(hmm-hits)s
    * coverage, detection, and other read-recruitment statistics in each sample stored in the %(profile-db)s
    * domain and taxonomy predictions from single-copy core genes (see %(anvi-run-scg-taxonomy)s)
    * the bin-specific value of each statistic in `bins_summary.txt`
* Matrix files compiling information about all bins across all samples (`bins_across_samples/`)
* Miscellaneous data exported from the %(profile-db)s (`misc_data_items/` and `misc_data_layers/`)

**Confused about the read-recruitment statistics?**

In case you want to learn about the definitions of statistics like coverage, detection, abundance, variability, and so on, you should first read [Mike Lee's explanation of these statistics](https://merenlab.org/2017/05/08/anvio-views/). Our [vocabulary page](https://anvio.org/vocabulary/) might also be helpful. Then, keep in mind that anvi'o computes these values on a per-contig (and per-split) basis. When you run %(anvi-summarize)s, the program will summarize this information for a given bin by taking the average of a statistic's value across all splits in the bin, weighting that average by split length.

### When run on a pan database

Running on a %(pan-db)s produces a large table (`[NAME]_gene_clusters_summary.txt`) describing every gene in every gene cluster of your pangenome, including:

* gene-cluster-specific information like the number of genomes contributing to that cluster, maximum number of paralogs in any participating genome, and cluster homogeneity metrics
* gene-specific information like functional annotations and amino acid sequence
* a `bin_name` column indicating which collection bin each gene cluster belongs to (empty when no collection is provided — see below)

It also exports any miscellaneous data imported into the %(pan-db)s (`misc_data_items/` and `misc_data_layers/`).

### When run on a pan-graph database

Running on a %(pan-graph-db)s produces a directory (by default named `[PROJECT]-PAN-GRAPH-SUMMARY`) containing the following files:

* `SYNGCs.txt` — one row per synteny gene cluster (SynGC), with columns for its identity, region membership, node type, number of genomes in which it is present, items additional data, and per-source function consensus annotations
* `GENESxSYNGCs.txt` — long-format table with one row per (genome × gene call), linking each gene call to its SynGC, region, per-gene functional annotations, and amino acid sequence (or DNA sequence if `--report-DNA-sequences` is used; omitted when `--quick-summary` is set)
* `REGIONS.txt` — one row per genomic region (backbone or variable), with all region-level statistics including the composite variability score
* `GENOMES_DIST_MAT.txt` — square genome × genome synteny distance matrix
* `GENOMES_DIST.newick` — newick tree over genomes derived from the distance matrix

It also exports miscellaneous data from the %(pan-graph-db)s (`misc_data_items/` and `misc_data_layers/`).

## Running anvi-summarize

### Running on a profile database

{{ codestart }}
anvi-summarize -c %(contigs-db)s \
               -p %(profile-db)s \
               -o MY_SUMMARY \
               -C %(collection)s
{{ codestop }}

When running on a profile database, you also have options to:

* output very accurate (but intensely processed) coverage and detection data for each gene (using `--init-gene-coverages`)
* edit your contig names so that they contain the name of the bin that the contig is in (using `--reformat-contig-names`)
* also display the amino acid sequences for your gene calls (using `--report-aa-seqs-for-gene-calls`)

### Running on a pan database

A %(collection)s is **optional** when summarizing a %(pan-db)s. Without one, anvi'o will still export the full gene clusters table with all functional annotations — the `bin_name` column will simply be empty. If you have organized your gene clusters into bins using %(anvi-interactive)s or %(anvi-import-collection)s, passing the collection name will populate `bin_name` with the bin each gene cluster belongs to, which makes downstream filtering by bin straightforward.

Run without a collection (exports everything):

{{ codestart }}
anvi-summarize -g %(genomes-storage-db)s \
               -p %(pan-db)s
{{ codestop }}

Run with a collection (adds `bin_name` to the output):

{{ codestart }}
anvi-summarize -g %(genomes-storage-db)s \
               -p %(pan-db)s \
               -C %(collection)s
{{ codestop }}

You can display DNA sequences instead of amino acid sequences with `--report-DNA-sequences`.

### Running on a pan-graph database

A %(collection)s is **optional** when summarizing a %(pan-graph-db)s. Without one, all output files are still produced in full — the `bin_name` column in `SYNGCs.txt` and `GENESxSYNGCs.txt` will simply be empty. If you have organized your SynGCs into bins, passing the collection name will populate `bin_name` in both files, making it easy to filter the output to any bin of interest with a single column filter.

Run without a collection (exports everything):

{{ codestart }}
anvi-summarize -g %(genomes-storage-db)s \
               --pan-graph-db %(pan-graph-db)s
{{ codestop }}

Run with a collection (adds `bin_name` to `SYNGCs.txt` and `GENESxSYNGCs.txt`):

{{ codestart }}
anvi-summarize -g %(genomes-storage-db)s \
               --pan-graph-db %(pan-graph-db)s \
               -C %(collection)s
{{ codestop }}

### Other notes

If you are unsure what collections are in your database, you can run this program with the flag `--list-collections` or by running %(anvi-show-collections-and-bins)s.

You can also use the flag `--quick-summary` to get a less comprehensive summary with a much shorter processing time. For profile-db summaries it skips several heavier computations; for pan-db summaries it omits sequences and annotation text from the gene clusters file; for pan-graph-db summaries it omits sequences from `GENESxSYNGCs.txt`.
