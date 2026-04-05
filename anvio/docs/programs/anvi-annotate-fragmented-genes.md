This program **identifies and annotates fragmented genes (pseudogenes) across a pangenome** by scanning gene clusters stored in a %(pan-db)s for evidence of gene fragmentation. The results are stored as %(functions)s in individual %(contigs-db)s files under the function annotation source `PSEUDO_GENES`.

Essentially, this program back propagates insights into gene fragmentation events in individual genomes by considering a pangenome. This allows investigation of the function and context of fragmented genes, their inclusion in formal reporting in interactive interfaces and summary outputs, and enables functional enrichment analyses.

### The problem

In comparative genomics, a gene that is intact in some genomes may be split into two or more adjacent open reading frames in others. Assuming that genomes are high-quality and do not suffer from extensive sequencing errors, fragmented genes occur when a point mutation or a transposon insertion introduces a premature stop codon into a gene, after which the gene caller identifies a new start codon downstream and reports two (or more) shorter genes where there used to be one.

If such fragmentation events are not common across all genomes, such fragmented genes appear as the following alignment in pangenomes:

```
Genome A gene x: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Genome B gene n: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx--------------
Genome B gene m: --------------------------------------xxxxxxxxxxxx

Genome C gene z: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

In this example, gene `n` and gene `m` in Genome B are adjacent on the same contig and together correspond to the full-length gene represented by `x` in Genome A and `z` in Genome C. Due to the adjacency of `n` and `m`, and their overlapping nature with other genes, one can assume that the gene in Genome B was likely split by a mutation that disrupted the original reading frame.

Gene fragmentation can occur in any group, but it is more common in some clades than others. For instance, *Brucella* are known for having many pseudogenes (and in fact this program is coming to life because [Sean Crosson](https://directory.natsci.msu.edu/directory/Profiles/Person/101773), who studies *Brucella* asked for it -- so thank you, Sean!). The anvi'o pangenomics workflow correctly groups the fragments with the full-length gene into a single gene cluster thanks to the all-vs-all BLAST search. But excessive gene fragments create a problem for downstream analyses when the issue is not accounted for: they clutter the analysis results with spurious singletons and obscure truly unique genes.

### The solution

%(anvi-annotate-fragmented-genes)s offers a solution for this issue by scanning every gene cluster in a %(pan-db)s and looking for cases where a single genome contributes two or more genes to the same gene cluster that are **adjacent on the same contig** (i.e., no other gene sits between them). When such a case is found, it compares the lengths of these adjacent fragments against a **full-length reference**, i.e., the longest gene in that cluster from a genome where the gene is not split.

Each fragment is then classified as one of the following:

- **`fragmented_gene`**: The longest fragment in a genome, which likely retains gene function. Or not. It is impossible to say, of course, but the assumption here is that if the gene that has undergone fragmentation, the longest fragment of it is the most likely to continue retaining its function. This label is only assigned if the fragment is at least a certain fraction of the full-length reference, controlled by the `--min-full-length-ratio` flag. The deafult value here is 0.7, which means this label is assigned to a fragment only if it is longer than 70%% of the reference gene.

- **`gene_fragment`**: Shorter fragments that are unlikely to be functional. If even the longest fragment in a genome falls below the length threshold, **all** fragments in that genome are labeled `gene_fragment`.

The results are written as functional annotations under the source name `PSEUDO_GENES` into each relevant %(contigs-db)s. The annotation text for each gene includes the percent coverage relative to the full-length reference and identifies the genome and gene caller ID of that reference, so it remains meaningful even outside the pangenome context and becomes a part of the information in the genome.

### Basic usage

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s
{{ codestop }}

This will scan all gene clusters, report fragmentation events to the terminal, and write `PSEUDO_GENES` annotations into the contigs databases listed in the %(external-genomes)s file.

### Terminal output

When the program runs, it prints a color-coded visualization for each gene cluster where fragmentation was detected. Here is an example output:

![terminal_output](../../images/anvi-annotate-fragmented-genes.png)

In this particular example, four of the gene clusters in the pangenome had fragmented genes. The bars in the report show the relative length of each gene compared to the full-length reference, and their position reflects the actual layout of fragments on the contig. Bars colored in **green** represent the full-length reference gene determined by anvi'o, **blue** represent the longest fragment in a genome (labeled `fragmented_gene`), **red** represent shorter fragments, and **gray** bars represent genes from genomes where the gene is not fragmented. Genome names and anvi'o gene caller ids are also shown to double check things.

The purpose of this report is for you to go back to the pangenome with %(anvi-display-pan)s, search for some of the gene clusters, and inspect them to confirm that you are happy with the result.

If you are satisfied and would like your pangenome to include this information, you will need to restart the pangenomics workflow with these newly annotated %(contigs-db)s files so %(anvi-summarize)s output can include the necessary data for you to be able to do functional enrichment analyses of genes that have `fragmented_gene` annotations.

### Adjusting the length threshold

By default, the longest fragment in a genome must be at least 70%% of the full-length reference to receive the `fragmented_gene` label. You can adjust this threshold:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --min-full-length-ratio 0.50
{{ codestop }}

Setting a lower value is more permissive (more fragments will be labeled `fragmented_gene` rather than `gene_fragment`). Setting a higher value is more conservative.

### Report-only and skip-reporting modes

If you want to see what the program would annotate **without actually writing to any database**, use the `--report-only` flag:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --report-only
{{ codestop }}

Conversely, if you want to annotate silently without the per-gene-cluster terminal visualizations, use `--skip-reporting`:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --skip-reporting
{{ codestop }}

### Re-running the program

If the `PSEUDO_GENES` source already exists in a contigs database (from a previous run), the program will overwrite the existing annotations. This means you can safely re-run the program with different parameters without needing to manually remove old annotations first.
