%(anvi-gen-pan-representative)s generates a **pangenome-supplemented representative genome** as a %(contigs-db)s. The output keeps all contigs from a single representative genome and appends a supplementary contig containing one representative gene for every gene cluster in the pangenome that the representative genome is missing. The result is a single %(contigs-db)s that captures the full gene repertoire of the pangenome.

![](../../images/anvi-gen-pan-representative.png)

#### Prerequisites

Before running this program you will need:
- A %(pan-db)s and %(genomes-storage-db)s produced by %(anvi-pan-genome)s.
- A %(contigs-db)s for each genome, created with %(anvi-gen-contigs-database)s from a %(contigs-fasta)s with simple deflines. If your FASTA files do not have simple deflines, use %(anvi-script-reformat-fasta)s to fix them first.
- An %(external-genomes)s file listing all genomes and their contigs database paths, which you can generate with %(anvi-script-gen-genomes-file)s.

Optionally, each %(contigs-db)s may carry functional annotations from programs such as %(anvi-run-kegg-kofams)s, %(anvi-run-ncbi-cogs)s, %(anvi-run-hmms)s, or %(anvi-run-scg-taxonomy)s. When present, these annotations are carried over into the output database.

#### Basic usage

{{ codestart }}
anvi-gen-pan-representative -p %(pan-db)s \
                            -g %(genomes-storage-db)s \
                            -e %(external-genomes)s \
                            -o PATH/TO/OUTPUT.db
{{ codestop }}

#### How the representative genome is selected

By default, the program scores each genome using a weighted combination of genome quality and assembly contiguity:

```
score = alpha × normalized(completeness − redundancy) + (1 − alpha) × (1 − normalized(number of contigs))
```

Both terms are normalized independently across all candidate genomes before combining. The default `--alpha` of `0.8` strongly favors completeness and low redundancy while still penalizing fragmented assemblies. An alpha of `1.0` considers only completeness and redundancy; a value of `0.0` will select only based on the least number of contigs. In the case of a tie, the genome with the most genes wins 🏆

`--max-num-contigs INT` restricts the candidate pool before scoring. Any genome with more contigs than this threshold is excluded from consideration entirely.

To bypass automatic selection and name the representative yourself, use `--representative`:

{{ codestart }}
anvi-gen-pan-representative -p %(pan-db)s \
                            -g %(genomes-storage-db)s \
                            -e %(external-genomes)s \
                            --representative GENOME_NAME \
                            -o PATH/TO/OUTPUT.db
{{ codestop }}

#### Selection options

`--alpha FLOAT` Weight for the completeness/redundancy term. Default `0.8`.

`--max-num-contigs INT` Exclude genomes with more contigs than this threshold before scoring.

`--representative GENOME_NAME` Skip automatic selection and specify the representative directly.

#### The supplementary contig

Once the representative genome is chosen, the program builds the supplementary contig using a greedy covering strategy: it repeatedly selects the genome that covers the most gene clusters not yet represented, extracts the representative genes per uncovered gene clusters from that genome, then continues with the next genome until every gene cluster in the pangenome is accounted for. All extracted gene sequences are concatenated into a single supplementary contig, separated by stretches of `N` characters (20 by default, controlled by `--gap-size`).

By default each gene is extracted individually from its start codon to its stop codon. Two additional modes offer finer control over what flanking sequence is included:

`--keep-synteny` Consecutive unrepresented genes from the same source contig are extracted as a single block — from the start codon of the first gene to the stop codon of the last — preserving their relative arrangement instead of isolating each gene.
`--keep-promoter` Extends each extracted block outward to the intergenic boundaries: from the stop codon of the preceding gene (or the contig start if there is none) to the start codon of the following gene (or the contig end if there is none). This captures the flanking regulatory regions and implies `--keep-synteny` automatically.

#### Extraction options

`--keep-synteny` Extract consecutive unrepresented genes from the same contig as a single block.

`--keep-promoter` Extend each extracted block to intergenic boundaries; implies `--keep-synteny`.

`--gap-size INT` Number of `N` characters inserted between genes in the supplementary contig. Default `20`.

#### Output options

`--output-file` (required) Path for the output %(contigs-db)s. Must end in `.db`.

`--project-name` A name for the project stored in the output database. If omitted, the name is taken from the input %(pan-db)s.

`--description` Path to a plain text file (Markdown supported) containing a description of the project. The description is stored in the database and rendered in the anvi'o interactive interface and summary outputs.

`--split-length` The length at which long contigs are virtually split for display and analysis purposes. Default is 20,000 bp. Set to `0` or any negative integer to disable splitting entirely.

`--kmer-size` K-mer size used for tetranucleotide frequency calculations in the output contigs database. Default is `4`. Accepted range is 1–5.

`--skip-mindful-splitting` By default anvi'o avoids cutting across gene boundaries when splitting contigs. Use this flag to split strictly at the requested length regardless of gene coordinates. Note that some downstream analyses may not be available for databases created this way.
