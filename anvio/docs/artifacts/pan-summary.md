The output directory produced by %(anvi-summarize)s when run on a %(pan-db)s and a %(genomes-storage-db)s.

By default the directory is named `[PROJECT]-PAN-SUMMARY`. Its central file is a large tab-delimited table that describes every gene in every gene cluster of the pangenome, regardless of whether a %(collection)s was supplied.

## Output files

### [NAME]_gene_clusters_summary.txt

One row per (gene cluster × genome × gene call). Columns include:

- `gene_cluster_id` — the gene cluster the row belongs to
- `bin_name` — the bin the gene cluster was assigned to in the provided %(collection)s, or empty if no collection was given
- `genome_name`, `gene_callers_id` — the genome and gene call for this row
- items additional data keys carried over from the %(pan-db)s
- per-source function accessions and annotations
- `aa_sequence` or `dna_sequence` (omitted when `--quick-summary` is used)

### misc_data_layers/ and misc_data_items/

Any miscellaneous data imported into the %(pan-db)s with %(anvi-import-misc-data)s, exported as tab-delimited files.

### index.html

An HTML document that formats all summary information for convenient browsing without an anvi'o installation.
