This program builds matrices that compare antiSMASH results across many genomes or metagenomes, reading the antiSMASH annotations straight from contigs databases that have already been processed with %(anvi-run-antismash)s.

It works like `anvi-estimate-metabolism --matrix-format`: point it at your genomes through an %(external-genomes)s file (`-e`), an %(internal-genomes)s file (`-i`), or both, and it reads the antiSMASH function sources (`antiSMASH`, `antiSMASH_ROLE`, `antiSMASH_DOMAIN`, `antiSMASH_REGION`) from each database and pivots them into matrices. An external genome is a whole %(contigs-db)s; an internal genome is a bin, so only that bin's genes are counted, and two bins from the same assembly get their own separate profiles. Because a matrix is a comparison across genomes, a single database (`-c`) is not accepted. There is nothing to run per-genome first beyond the annotation itself, which you can do for a whole set in one command with `anvi-run-antismash -e %(external-genomes)s`.

You choose one of two matrix types:

### Per-cluster matrix (`--regions`)

Rows are BGC **product types** (terpene, NRPS, hserlactone, …), columns are your genomes, and each cell is the **number of regions (clusters)** of that type in that genome, a compact profile of biosynthetic potential across your dataset.

{{ codestart }}
anvi-script-gen-antismash-matrix -e %(external-genomes)s \
                                 --regions \
                                 -O BGC_TYPES
{{ codestop }}

### Per-gene matrix (`--genes`)

Rows are distinct gene **signatures**, described by five columns (`product`, `role`, `function`, `smcog`, and `domains`) that stay constant down each row; columns are your genomes, and each cell is how many genes match that signature. The leading `product` column names the BGC type the gene belongs to, so it lines up directly with the rows of the `--regions` matrix (a gene signature that occurs in more than one cluster type honestly appears as more than one row). Genes that sit inside a cluster but carry no functional annotation of their own are not shown here; they are in the per-gene file described below.

{{ codestart }}
anvi-script-gen-antismash-matrix -e %(external-genomes)s \
                                 --genes \
                                 -O BGC_GENES
{{ codestop }}

### Output

Either mode writes two files of type %(functions-across-genomes-txt)s: a `-FREQUENCY.txt` matrix (counts) and a `-PRESENCE-ABSENCE.txt` matrix (0/1). In both, the descriptor columns come first, followed by one column per genome.

If any of the databases you provide were never annotated with %(anvi-run-antismash)s, the program still runs, but fills that genome's column with `NA` (rather than `0`) and warns you about it. The distinction matters: `NA` means "not annotated, unknown", whereas `0` would wrongly imply "annotated, but no cluster of that kind". Annotate all your databases first for a clean, all-numeric comparison.

### Relating the matrices back to individual genes (`--genes`)

A matrix is a cross-genome summary: each cell is a *count*, so it cannot carry the identifiers that only make sense inside one database: the gene caller id and the region id. To bridge that gap, `--genes` additionally writes a comprehensive per-gene file, `<prefix>-GENE-DETAILS.txt`, with one row per antiSMASH-annotated gene across all your genomes and these columns:

`genome`, `gene_callers_id`, `contig`, `region_id`, `product`, `role`, `function`, `smcog`, `domains`

This one file is the join key for everything downstream:

* `genome` + `gene_callers_id` points back to that genome's %(contigs-db)s, so you can relate an antiSMASH gene to its other annotations (KOfam, COGs, …), since every source shares the same gene caller id (see %(anvi-export-functions)s).
* `region_id` groups genes into their clusters, so you can list all the genes that make up a given BGC region (including accessory genes with no function of their own, which the matrix omits).
* `product` and the signature columns line up with the `--genes` and `--regions` matrices, so you can go from any matrix cell down to the individual genes behind it.
