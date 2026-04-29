Anvi-summarize lets you export a **comprehensive overview of your data** from an anvi'o database. Depending on the input, it can summarize a %(collection)s of binned contigs (from a %(profile-db)s), a %(collection)s of binned gene clusters (from a %(pan-db)s), or the full contents of a %(pan-graph-db)s. The output is a directory of flat files and an HTML index that conveniently displays them for you. This makes the program useful for sharing information with collaborators, generating supplementary files for manuscripts, and exporting data for use in downstream analyses.

See also %(anvi-summarize-blitz)s.

## Output files

What this program produces as output depends on its inputs:

* Running it on a %(profile-db)s will produce a %(profile-summary)s output.
* Running it on a %(pan-db)s and %(genomes-storage-db)s will produce a %(pan-summary)s output.
* Running it on a %(pan-graph-db)s and %(genomes-storage-db)s will produce a %(pan-graph-summary)s output.

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
