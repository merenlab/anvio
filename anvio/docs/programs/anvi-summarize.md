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

### Computing discov statistics

When summarizing a profile database, you can optionally compute the distribution of coverage (discov) score for each bin and each contig within each bin. See %(discov-stats)s for a full description of the metric and its parameters.

To enable discov computation, use the `--report-discov` flag:

{{ codestart }}
anvi-summarize -c %(contigs-db)s \
               -p %(profile-db)s \
               -o MY_SUMMARY \
               -C %(collection)s \
               --report-discov
{{ codestop }}

This requires access to the auxiliary data file (`AUXILIARY-DATA.db`) in the same directory as the profile database, and a profile database that includes the number of reads assigned to each split. It produces two additional output files under `bins_across_samples/`:

* `discov_bins.txt` — one row per bin × sample
* `discov_contigs.txt` — one row per contig × sample. If the splits of a contig are in different bins, each of these bins gets a row for the part of the contig that is in it, and the `partial_contig_in_bin` column is set to 1 for these rows.

Both files report the length of each sequence, the number of reads mapping to it (`num_reads`), the number of windows, the observed and expected proportion of covered windows, the normalized spread (_S_), the window depth evenness (_E_), the spread-evenness score (_C_), and the discov score (`dis_cov`).

You can adjust how discov is computed using the following parameters. When neither `--window-length` nor `--window-length-as-percentage` is specified, anvi'o applies context-sensitive defaults: a fixed 10,000 bp window for bins, and a percentage-based window (1%% of contig length, minimum 300 bp) for individual contigs. If you specify window length parameters, they apply to both bins and contigs.

**Window sizing**

* `--window-length INT` — use a fixed window size in bp for all sequences (bins and contigs)
* `--window-length-as-percentage INT` — set window length as a percentage of each sequence's length
* `--min-window-length INT` — minimum window length floor for percentage mode (default for contigs: 300)

**Combining S and E**

* `--alpha FLOAT` — weight of S relative to E, in [0, 1] (default: 0.5)
* `--discov-formula STRING` — `geometric` (C = S^α × E^(1-α)) or `linear` (C = αS + (1-α)E) (default: `geometric`)

The final score is discov = min(1, C + D^2), where D is detection.

{:.notice}
Discov values computed by %(anvi-summarize)s and %(anvi-profile-blitz)s on the same data should be identical, with one exception: the auxiliary data file caps per-nucleotide coverage at 65,535X, so the values can differ slightly for sequences with extremely high coverage. See %(discov-stats)s for details.

### Other notes

If you are unsure what collections are in your database, you can run this program with the flag `--list-collections` or by running %(anvi-show-collections-and-bins)s.

You can also use the flag `--quick-summary` to get a less comprehensive summary with a much shorter processing time. For profile-db summaries it skips several heavier computations; for pan-db summaries it omits sequences and annotation text from the gene clusters file; for pan-graph-db summaries it omits sequences from `GENESxSYNGCs.txt`.

For cases where you want the aggregate bin statistics and per-sample coverage matrices from your profile-db but not per-bin sequence files, you can use the flag `--light-summary` instead.

Just for your reference, in the context of %(profile-db)s summaries, and for a dataset of ~2,000 MAGs described in a 42Gb %(profile-db)s and a 9.5G %(contigs-db)s files, %(anvi-summarize)s took over 240 minutes to run during a test in 2026. With `--light-summary`, the program took about 40 minutes, and with `--quick-summary`, it took about 10 minutes to run on the same dataset.
