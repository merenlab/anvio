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

In this example, gene `n` and gene `m` in Genome B are adjacent on the same contig and together correspond to the full-length gene represented by `x` in Genome A and `z` in Genome C. Due to the adjacency of `n` and `m`, and their overlap with other genes in the cluster, one can assume that the gene in Genome B was likely split by a mutation that disrupted the original reading frame.

Gene fragmentation can occur in any group, but it is more common in some clades than others. For instance, *Brucella* are known for having many pseudogenes (and in fact this program is coming to life because [Sean Crosson](https://directory.natsci.msu.edu/directory/Profiles/Person/101773), who studies *Brucella* asked for it -- so thank you, Sean!). The anvi'o pangenomics workflow correctly groups the fragments with the full-length gene into a single gene cluster thanks to the all-vs-all BLAST search. But excessive gene fragments create a problem for downstream analyses when the issue is not accounted for: they clutter the analysis results with spurious singletons and obscure truly unique genes.

### The solution

%(anvi-annotate-fragmented-genes)s offers a solution for this issue by scanning every gene cluster in a %(pan-db)s and looking for cases where a single genome contributes two or more genes to the same gene cluster that are **adjacent on the same contig** (i.e., no other gene sits between them). When such a case is found, it compares the lengths of these adjacent fragments against a **full-length reference**, i.e., the longest gene in that cluster from a genome where the gene is not split.

Each fragment is then classified as one of the following:

- **`fragmented_gene`**: The longest fragment in a genome, which likely retains gene function. Or not. It is impossible to say, of course, but the assumption here is that if a gene has undergone fragmentation, the longest fragment of it is the most likely to continue retaining its function. This label is only assigned if the fragment is at least a certain fraction of the full-length reference, controlled by the `--min-full-length-ratio` flag (default: 0.5, which means the fragment must be longer than 50%% of the reference gene).

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

In this particular example, four of the gene clusters in the pangenome had fragmented genes. The bars in the report show the relative length of each gene compared to the full-length reference, and their position reflects the actual layout of fragments on the contig. **Green** bars represent the full-length reference gene determined by anvi'o, **blue** bars represent the longest fragment in a genome (labeled `fragmented_gene`), **red** bars represent shorter fragments, and **gray** bars represent genes from genomes where the gene is not fragmented. Genome names and anvi'o gene caller ids are also shown to double check things.

You can also include the consensus function annotation for each gene cluster in the report header by providing a functional annotation source with `--annotation-source`:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --annotation-source COG24_FUNCTION
{{ codestop }}

This will display the consensus functional annotation for each gene cluster (i.e., the most common function across all genes in it) for you to have a quick idea about their potential role:

![terminal_output](../../images/anvi-annotate-fragmented-genes-w-functions.png)

The purpose of the terminal output is to give you a sense of the decisions made and how they look in the context of the pangenome in general. It would be useful to go back to the pangenome with %(anvi-display-pan)s, search for some of the gene clusters, and inspect them to confirm that you are happy with the results.

If you are satisfied and would like your pangenome to include this information, you will need to restart the pangenomics workflow with these newly annotated %(contigs-db)s files so %(anvi-summarize)s output can include the necessary data for you to be able to do functional enrichment analyses of genes that have `fragmented_gene` annotations.

### Adjusting the length threshold

By default, the longest fragment in a genome must be at least 50%% of the full-length reference to receive the `fragmented_gene` label. You can adjust this threshold to be more stringent or more permissive with the `--min-full-length-ratio` flag:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --min-full-length-ratio 0.75
{{ codestop }}

Setting a lower value is more permissive (more fragments will be labeled `fragmented_gene` rather than `gene_fragment`). Setting a higher value is more conservative. The latter will risk losing the representation of more fragmented genes in downstream analyses and that's why the default is set to 0.5, but the final call may depend on your survey of the terminal report (so please take time to study your terminal output).

### Distinguishing fragmentation from duplication

Not every pair of adjacent genes in the same gene cluster is a fragmentation event. When a gene has been **duplicated** in tandem (producing two near-full-length paralogs side by side on the contig), both copies end up in the same gene cluster because they share high sequence similarity. Since they are adjacent to one another, they will look like a potential fragmentatin even to our algorithm. The critical insight that will distinguish gene duplication from gene fragmentation will come from the difference between the combined length of adjacent genes compared to the reference: while fragmentatinon will roughly sum to the full lenght of the reference (e.g., a 60%% fragment + a 40%% fragment ≈ 100%% of the reference), duplicated genes will sum to a much larger length than the referenc esince each copy of the gene will be near-full-length, so together they will be closer to ~200%% of the reference (or more, for higher-copy tandem repeats).

By default, %(anvi-annotate-fragmented-genes)s skips any group of adjacent genes whose combined length exceeds 1.2× the full-length reference, treating them as probable paralogs. You can adjust this threshold with the `--max-combined-length-ratio` flag:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --max-combined-length-ratio 1.30
{{ codestop }}

A higher value is more permissive (fewer groups will be excluded as paralogs). A lower value is more conservative. The default of 1.20 allows for some overlap at the fragment boundary while still catching obvious duplications, which appear to be a good idea based on our tests.

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

### Finding stray out-of-frame fragments

By default, this program only detects fragmentation events where both fragments end up in the **same gene cluster**, which happens when the downstream fragment remains in the same reading frame as the original gene. However, if the premature stop codon shifts the reading frame (e.g., a single-nucleotide insertion or deletion rather than a substitution), the downstream fragment will be called in a different frame, and MCL will place it in a **different gene cluster**, since it no longer shares sequence similarity with the original gene. This is what that situation would look like compared to the example before; one gene cluster would look like this

```
Genome A gene x: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Genome B gene n: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx--------------

Genome C gene z: xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

And another one would look like this:

```
Genome B gene m: xxxxxxxxxxxx
```

The program includes an optional flag, `--find-stray-fragments`, to search for these out-of-frame fragments:

{{ codestart }}
anvi-annotate-fragmented-genes -p %(pan-db)s \
                               -g %(genomes-storage-db)s \
                               -e %(external-genomes)s \
                               --find-stray-fragments
{{ codestop }}

When this flag is set, %(anvi-annotate-fragmented-genes)s performs a second scan after the standard in-cluster analysis. For each gene cluster, it looks for genomes where the cluster contains a single gene that is significantly shorter than the full-length reference. It then checks whether an adjacent gene on the same contig (one that belongs to a **different** gene cluster) together with the truncated gene approximates the expected full-length gene. If so, both are annotated as fragments.

This is a more aggressive search, and it may occasionally flag genes that are genuinely short rather than fragmented, so it is off by default. But while the algorithm worked well in our mock datasets, Meren's test with a large *B. fragilis* pangenome in which the program found over 100 gene clusters with fragmented genes, it found 0 stray fragments, so it is safe to assume that its false positive rate will be rather small if any.

{:.notice}
**A note from** {%% include person/display_mini_single.html github="meren" %%}: *The zero strays in a pangenome that contained over 100 regular fragmentation events likely indicates that the process is likely a result of biology rather than bioinformatics. As in, most premature stops likely come from substitutions, not frameshifts (i.e., C-to-T turning CAG (Gln) into TAG (stop) preserves the reading frame, and both fragments stay in-frame, and then BLAST clusters them together, and then the in-cluster scan catches them. It is also possible that most frameshifted downstream sequences often aren't called as genes by Prodigal. Even if a frameshift creates a new reading frame downstream, Prodigal needs to find a valid start codon, a [Shine-Dalgarno-like signal](https://en.wikipedia.org/wiki/Shine–Dalgarno_sequence), and a reasonable ORF length before it calls it a gene. Difficult to know which one is playing a more significant role, but if you are reading these lines, and if you feel that you have an interesting observation from your own pangenome or ideas about why out-of-frame / stray fragments occur in much less frequency compared to in-frame fragments, please let us know and so we can update the code if we are making a mistake, or this section with a better explanation*.

### Re-running the program

If the `PSEUDO_GENES` source already exists in a contigs database (from a previous run), the program will overwrite the existing annotations. This means you can safely re-run the program with different parameters without needing to manually remove old annotations first.
