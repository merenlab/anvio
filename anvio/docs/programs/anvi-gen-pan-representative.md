%(anvi-gen-pan-representative)s generates a **pangenome-supplemented representative genome** as a %(contigs-db)s. The output keeps all contigs from a single representative genome intact, and appends a single supplementary contig that contains one representative gene for every gene cluster in the pangenome that was absent from the representative genome. The result is a single %(contigs-db)s that captures the full gene repertoire of the pangenome.

![](../../images/anvi-gen-pan-representative.png)

By doing so, %(anvi-gen-pan-representative)s offers a tractable solution to a fundamental problem in comparative genomics and genome-resolved metagenomics: selecting a single representative for a set of closely related genomes that does not represent the entirety of the gene pool. Selecting representatives is a necessary step for a broad range of analyses in microbiology, typically performed by (1) clustering all genomes using an arbitrary similarity threshold (such as 95%%), or by identifying a clade of very closely related organisms after a phylogenomic analysis, and then (2) choosing a representative genome to stand in for the entire group in downstream analyses.

The problem %(anvi-gen-pan-representative)s solves arises from the fact that pangenomes almost always contain genes distributed unevenly across member genomes, making it extremely unlikely for any single genome to carry all the genes present in the group. The program addresses this by supplementing the chosen representative with a single additional contig (so called 'the supplementary contig') that carries one representative gene sequence per missing gene cluster. The result is a %(contigs-db)s that retains the full genomic integrity of the representative while extending its gene content to cover the entire pangenome, enabling downstream analyses (read recruitment, functional annotation, phylogenomics, and more) that would otherwise miss genes present in the group but absent from the representative.

## A real example

Here is an example using a *Mycobacterium tuberculosis* pangenome (which is an extremely closed pangenome, as expected from the lifestyle of this particular organism):

![TB pangenome](../../images/anvi-gen-pan-representative-tb-pangenome.png)

The green layer in this anvi'o visualization marks the representative genome, and red selections mark the gene clusters (which contain one or more genes from one or more genomes in the collectoin of all *M. tuberculosis* genomes) that are *absent* in the representative. Running %(anvi-gen-pan-representative)s on this pangenome results in a single %(contigs-db)s, visualization of which shows two contigs:

![TB representative genome](../../images/anvi-gen-pan-representative-tb-pan-contigs-db.png)

While the genome is intact, representative genes from each gene cluster marked red in the previous figure now represented in a single contig where genes can be traced back to the original genome they were encoded. Using %(anvi-export-contigs)s, one can generate a %(fasta)s file for this representative genome with pangenome supplemented contig of missing genes for a comprehensive, but not inflated representatiion of the gene pool of the *Mycobacterium tuberculosis* pangenome considered here.

## Prerequisites

Before running this program you will need:

* A %(pan-db)s and %(genomes-storage-db)s produced by %(anvi-pan-genome)s.

* A %(contigs-db)s for each genome, created with %(anvi-gen-contigs-database)s from a %(contigs-fasta)s with simple deflines. If your FASTA files do not have simple deflines, use %(anvi-script-reformat-fasta)s to fix them first.

* An %(external-genomes)s file listing all genomes and their contigs database paths, which you can generate with %(anvi-script-gen-genomes-file)s.

Optionally, each %(contigs-db)s may carry functional annotations from programs such as %(anvi-run-kegg-kofams)s, %(anvi-run-ncbi-cogs)s, %(anvi-run-hmms)s, or %(anvi-run-scg-taxonomy)s. When present, these annotations are carried over into the output database.

## Basic usage

{{ codestart }}
anvi-gen-pan-representative -p %(pan-db)s \
                            -g %(genomes-storage-db)s \
                            -e %(external-genomes)s \
                            -o PATH/TO/OUTPUT.db
{{ codestop }}

## Selection of the representative genome

By default, the program scores each genome using a weighted combination of genome quality and assembly contiguity:

```
score = alpha × normalized(completeness − redundancy) + (1 − alpha) × (1 − normalized(number of contigs))
```

Both terms are normalized independently across all candidate genomes before combining. The default `--alpha` of `0.8` strongly favors completeness and low redundancy while still penalizing fragmented assemblies. An alpha of `1.0` considers only completeness and redundancy; a value of `0.0` selects only on the basis of the fewest contigs. In the case of a tie, the genome with the most genes wins.

* `--alpha FLOAT` adjusts the weight of the completeness/redundancy term relative to assembly contiguity. Default is `0.8`.

* `--max-num-contigs INT` restricts the candidate pool before scoring. Any genome with more contigs than this threshold is excluded from consideration entirely.

To bypass automatic selection, and exclusively name the representative yourself, you can use the parameter `--representative`:

{{ codestart }}
anvi-gen-pan-representative -p %(pan-db)s \
                            -g %(genomes-storage-db)s \
                            -e %(external-genomes)s \
                            --representative GENOME_NAME \
                            -o PATH/TO/OUTPUT.db
{{ codestop }}

`GENOME_NAME` should match the name of the genome in the %(external-genomes)s file.

## Creation of the supplementary contig

Once the representative genome is chosen, the program builds the supplementary contig using a greedy covering strategy: it repeatedly selects the donor genome that covers the most gene clusters not yet represented, extracts the representative gene sequences for those uncovered gene clusters from that donor genome, then moves on to the next best donor until every gene cluster in the pangenome is accounted for. All extracted gene sequences are concatenated into a single supplementary contig, separated by stretches of `N` characters (20 by default, controlled by `--gap-size`).

By default, each gene is extracted individually from its start codon to its stop codon. Two additional modes offer finer control over what flanking sequence is included:

* `--keep-synteny` Consecutive unrepresented genes from the same source contig are extracted as a single block (i.e., from the start codon of the first gene to the stop codon of the last) and thereby their relative arrangement instead of isolating each gene.


* `--keep-promoter` Extends each extracted block outward to the intergenic boundaries: from the stop codon of the preceding gene (or the contig start if there is none) to the start codon of the following gene (or the contig end if there is none). This captures flanking regulatory regions and implies `--keep-synteny` automatically.

* `--gap-size INT` controls the number of `N` characters inserted between sequences in the supplementary contig. Default is `20`.

Please see the program help menu to see the most up-to-date parameters and output options.
