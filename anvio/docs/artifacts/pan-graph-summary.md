The output directory produced by %(anvi-summarize)s when run on a %(pan-graph-db)s and a %(genomes-storage-db)s.

By default the directory is named `[PROJECT]-PAN-GRAPH-SUMMARY`. It contains a set of flat files that together describe every synteny gene cluster (SynGC), every gene call, every genomic region, and the genome-level distance structure of the pan-graph.

## Output files

### SYNGCs.txt

One row per SynGC. Columns include:

- `node_id` — the SynGC identifier
- `bin_name` — the bin the SynGC was assigned to in the provided %(collection)s, or empty if no collection was given
- `source_gene_cluster_id` — the gene cluster in the source pangenome from which this SynGC was derived
- `node_type`, `region_id`, `region_type` — structural context in the pan-graph
- `num_genomes_present`, `genomes_present` — which genomes carry this SynGC
- items additional data keys carried over from the %(pan-graph-db)s
- per-source function consensus accessions and annotations

### GENESxSYNGCs.txt

Long-format table with one row per (genome × gene call). Columns include:

- `node_id`, `bin_name`, `source_gene_cluster_id` — links back to `SYNGCs.txt`
- `genome_name`, `gene_caller_id` — the genome and gene call for this row
- `region_id`, `region_type` — region membership
- per-source per-gene-call function accessions and annotations
- `aa_sequence` or `dna_sequence` (omitted when `--quick-summary` is used)

### REGIONS.txt

One row per genomic region (backbone or variable), dumped directly from the pan-graph regions table. Includes all region-level statistics such as the composite variability score, diversity, complexity, and expansion metrics.

### GENOMES_DIST_MAT.txt

Square genome × genome synteny distance matrix.

### GENOMES_DIST.newick

Newick tree over genomes derived from the distance matrix above.

### misc_data_layers/ and misc_data_items/

Any miscellaneous data imported into the %(pan-graph-db)s with %(anvi-import-misc-data)s, exported as tab-delimited files.
