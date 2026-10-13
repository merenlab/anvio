This program relates changes in tRNA-seq seed abundances to the codon usage of genes or gene functions.

An mRNA transcript is more efficiently translated when its codons have an abundant supply of matching tRNAs. We define a metric called *affinity* to measure the relationship between the pool of tRNA transcripts — tRNA "supply" — and the codon usage of a protein-coding gene or a functional group of genes — tRNA "demand." A positive affinity indicates that the sample's tRNA pool is better matched to the gene's codons than the reference pool, favoring faster translation of the gene, whereas a negative affinity suggests slower translation relative to the reference.

The affinity of a gene averages the log₂-fold change in the abundance of each tRNA isoacceptor relative to a reference, weighted by how much the gene's codons use that isoacceptor. (An isoacceptor is the set of tRNAs that share an anticodon.) Concisely, affinity is the **demand-weighted average change in tRNA supply**.

Affinity is closely related to the tRNA Adaptation Index (tAI). Classic tAI ([dos Reis et al., 2004](https://doi.org/10.1093/nar/gkh834)) weights each codon by the number of tRNA *genes* that read it — a fixed property of the genome. Affinity instead uses *measured* tRNA abundances, which change from sample to sample; in this it follows the expression-based tAI of [Wei et al. (2019)](https://doi.org/10.1038/s41598-019-39369-x) (their "tRNA tpm"). Affinity is also *differential*: it is the log₂ ratio of this expression-based tAI between the sample and the reference. This ratio cancels any bias that scales a tRNA's measured abundance consistently across samples: tRNA-seq undercounts some tRNAs — for instance, where a nucleotide modification or stable structure blocks reverse transcription — but such a bias multiplies that tRNA's abundance by the same factor in the sample and the reference, so it divides out of the log₂ ratio.

## Quick start

Here are examples of commands highlighting some of the program capabilities.

This minimal run on a single genome compares each tRNA-seq sample to a designated reference sample and writes a table of affinities:

{{ codestart }}
anvi-compute-trnaseq-functional-affinity -t %(trnaseq-contigs-db)s \
                                         -s %(seeds-specific-txt)s \
                                         -c %(contigs-db)s \
                                         -r reference_sample \
                                         -o affinity.txt
{{ codestop }}

When no single sample is a natural reference (in a time series, for example), one approach is to compare each sample's tRNA supply to the geometric mean of all samples by using `--reference-mean`:

{{ codestart }}
anvi-compute-trnaseq-functional-affinity -t %(trnaseq-contigs-db)s \
                                         -s %(seeds-specific-txt)s \
                                         -c %(contigs-db)s \
                                         --reference-mean \
                                         -o affinity.txt
{{ codestop }}

tRNA-seq data can be analyzed against multiple genomes, which can be given by an external genomes file:

{{ codestart }}
anvi-compute-trnaseq-functional-affinity -t %(trnaseq-contigs-db)s \
                                         -s %(seeds-specific-txt)s \
                                         -e %(external-genomes)s \
                                         --reference-mean \
                                         -o affinity.txt
{{ codestop }}

## Input data

Seeds, or predicted tRNA transcripts, are the end product of the anvi'o trnaseq %(workflow)s. The program %(anvi-integrate-trnaseq)s matches seeds to the underlying tRNA genes in your (meta)genomes, storing the matches as %(trna-gene-hits)s in the %(trnaseq-contigs-db)s. Seed coverages across the tRNA-seq samples are collected into a single, easily parsed %(seeds-specific-txt)s table by %(anvi-tabulate-trnaseq)s. ("Specific" as opposed to "nonspecific" coverage is explained in the %(trnaseq-profile-db)s artifact.) `anvi-compute-trnaseq-functional-affinity` needs three inputs: the %(trnaseq-contigs-db)s, its %(seeds-specific-txt)s, and a (meta)genomic %(contigs-db)s. It computes the gene and function codon frequencies it needs internally, by the same calculation as %(anvi-get-codon-frequencies)s.

Affinity can be calculated for genes or functional groups of genes (see %(functions)s), using the codon frequencies of the concatenated gene sequences annotated by each function. Anvi'o can annotate contigs databases with protein functions from various sources including COG and KEGG. %(anvi-run-kegg-kofams)s not only annotates KEGG Orthologs (KOs), but also KEGG BRITE hierarchical categories of KOs across biological levels. Affinity values for groups of proteins at different levels of a hierarchy may be useful in addition to affinity values for individual proteins.

The functional groups to analyze are chosen with `--function-sources` (default: `KEGG_BRITE`): give it a list of annotation sources to use specific ones, or use it as a bare flag to include every source annotated in the %(contigs-db)s. For KEGG BRITE, each gene is assigned to its most specific category by default; `--all-brite-categories` instead emits a functional group for every level of the hierarchy. To analyze **individual genes** rather than functions, use `--gene-affinity`: the output is then indexed by `gene_caller_id`, and function selection options no longer apply. Restrict gene mode to specific genes with `--gene-caller-ids`, which requires a single `--contigs-db` (gene-caller IDs are only unique within one database).

Affinity can be calculated using a %(collection)s of bins (e.g., metagenome-assembled genomes), which would have been provided in running %(anvi-integrate-trnaseq)s, and must also be provided using collection options in this program.

When a tRNA-seq seed matches genes in more than one of the analyzed genomes, `--seed-assignment` decides how its coverage is counted:

  - `unambiguous_genome` (the default) keeps only seeds that match a single genome being analyzed by the program.
  - `unambiguous_db` keeps only seeds that are unique within the entire %(trnaseq-contigs-db)s, which can include more matched genomes than are being analyzed by the program.
  - `ambiguous_all` counts a shared seed toward every genome it matches. The affinities of those genomes are then not independent, because they share the same read coverage.
  - `ambiguous_choose` tries to assign each shared seed to one genome — the genome whose *unambiguous* seeds are best covered, provided they beat the runner-up genome's unambiguous seeds by at least a factor of `--min-coverage-ratio` (default 5) in every sample that has unambiguous coverage. Seeds that fail this test are dropped.

This choice only matters when seed-to-genome assignment was not already resolved in running %(anvi-integrate-trnaseq)s.

## Translational affinity metric

### Wobble weights

Like the tRNA Adaptation Index ([dos Reis et al., 2004](https://doi.org/10.1093/nar/gkh834)), affinity weights each codon-anticodon pairing by a *decoding weight* rather than by counting every pairing equally. A decoding weight estimates the selective constraint on the efficiency of codon-anticodon coupling — how efficiently a tRNA anticodon reads a codon it pairs with. Watson-Crick base pairs (A-U, C-G) couple with 100%% efficiency. Wobble pairs are less efficient; for example, G-U has an estimated efficiency of ~37%% ([Sabi and Tuller, 2014](https://doi.org/10.1093/dnares/dsu017)).

Decoding weights can be changed from the defaults to test the sensitivity of affinity measurements to this parameterization: run the program with `--get-default-decoding-weights` to write out the default table, edit the weights, and supply the new table back to the program with `--decoding-weights-txt`. There is uncertainty in wobble weights as they were not measured directly but inferred by optimizing the correlation between tAI and codon usage bias (CUB), assuming that highly expressed genes like ribosomal proteins with high CUB are best adapted to the genomic tRNA repertoire ([Sabi and Tuller, 2014](https://doi.org/10.1093/dnares/dsu017)).

### Reference

The nature of tRNA-seq data requires comparison between samples to measure translational affinity. The user must provide **two or more** tRNA-seq samples in the %(trnaseq-contigs-db)s/%(seeds-specific-txt)s and choose a **reference** against which samples are compared. The reference is the geometric mean of isoacceptor relative abundances over a chosen set of samples, defined by providing exactly one of three mutually exclusive options:

  - `--reference-sample` (`-r`): a single sample is the reference. By default every *other* sample is then analyzed — a sample compared to itself would be trivially zero. This can be useful when one sample is a natural baseline, e.g., an unstressed control in a stress-response experiment, or the first sample of a time series.
  - `--reference-samples`: an explicit set of samples whose geometric mean is the reference. This can be useful when the meaningful baseline is a *group* of samples, e.g., all samples taken at one phase of a diel cycle.
  - `--reference-mean`: the geometric mean of **all** samples is the reference. This is the appropriate choice when no single sample is naturally the reference, e.g., samples were collected evenly over a diel cycle.

The set of samples that *receive* an affinity (the analyzed set) is chosen independently with `--analyzed-samples` (default: all samples); the reference set and analyzed set may overlap. By default the output's sample columns are ordered alphabetically; use `--sample-order` to set an explicit order. (List the sample names available for these options by running the program with just `--seeds-specific-txt` and `--list-samples`.)

Why compare samples at all? Because measurement biases — chiefly in the reverse transcription of tRNA into cDNA — can distort the apparent abundances of different tRNAs within any single sample. Measured relative abundances of, say, tRNA-Arg-ACG, tRNA-Ala-UGC, and tRNA-Tyr-GUA of 20%%, 0.10%%, and 3.0%% in one sample may bear little relation to reality. Two sources of bias include:

  - **Nucleotide modifications.** Some modifications make reverse transcriptase stall or fall off the tRNA template, while others are read through, sometimes introducing a sequencing error. (Those errors are informative — they point to modification sites — and the trnaseq %(workflow)s exploits them.)
  - **tRNA length.** tRNAs have a bimodal length distribution, and the longer ones carry ~15 extra nucleotides in the variable loop. This gives reverse transcriptase more chances to fall off before it reaches the anticodon arm from the 3' end.

We assume these biases are similar across the samples of an experiment, even if the levels of some modifications change. We therefore measure *changes* in each isoacceptor's abundance relative to the reference, rather than trusting absolute abundances. Affinity must be read in this light: a gene is favored or disfavored for translation by how the tRNA pool shifts **relative to the reference**.

Example: Sample 1 (used here as the reference) has tRNA-Arg-ACG, tRNA-Ala-UGC, and tRNA-Tyr-GUA relative abundances of 20%%, 0.10%%, and 3.0%%, respectively; Sample 2 has relative abundances of 30%%, 0.08%%, and 2.5%%; and Sample 3 has relative abundances of 15%%, 0.16%%, and 3.9%%. The fold changes in abundance for Sample 2 versus the reference are 1.5, 0.80, and 0.83, and the changes for Sample 3 are 0.75, 1.6, and 1.3. (With `--reference-mean` the denominator would instead be the geometric mean of all three samples rather than Sample 1.)

### Calculation

The affinity combines these tRNA dynamics with the (static) codon composition of a gene, through two vectors.

The first, **T**, holds the fold change in each isoacceptor's measured relative abundance (sample over reference). In our example, **T** = [1.5, 0.80, 0.83] for Sample 2 and **T** = [0.75, 1.6, 1.3] for Sample 3.

The second vector, **C**, records how heavily the gene's codons use each isoacceptor. An isoacceptor's entry in **C** is the sum of its codons' frequencies in the gene, each weighted by how efficiently the isoacceptor reads that codon. Take the first isoacceptor, Arg-ACG. (Its wobble adenosine is deaminated to inosine, giving Arg-ICG — typically the only tRNA with a wobble inosine in bacteria.) It reads CGU with an efficiency of 1, CGC with ~0.58, and CGA with ~0.12. Suppose Gene A is 1,000 nucleotides long and contains 34 CGU, 21 CGC, and 25 CGA codons; then its Arg-ICG entry is (34/1000)·1 + (21/1000)·0.58 + (25/1000)·0.12 = 0.049. Suppose the entries for Ala-UGC and Tyr-GUA work out to 0.010 and 0.035, so **C** = [0.049, 0.010, 0.035].

Dividing **C** by its total gives the **demand** vector **y** = **C** / Σᵢ Cᵢ = [0.521, 0.106, 0.372] — the fraction of Gene A's (efficiency-weighted) codon usage that each isoacceptor serves.

Affinity **A** is the **demand-weighted average of the log₂ supply fold changes**. Write isoacceptor i's log₂ fold change as **x**ᵢ = log₂(Tᵢ); these form the per-sample **supply** vector **x**, and **A** = Σᵢ yᵢ · xᵢ. For Sample 2, **x** ≈ [0.58, −0.32, −0.26] and **A** ≈ 0.521·0.58 + 0.106·(−0.32) + 0.372·(−0.26) ≈ **+0.17**; for Sample 3, **x** ≈ [−0.42, 0.68, 0.38] and **A** ≈ **0.00**. Sample 2's tRNA pool has the higher (positive) affinity to Gene A — the isoacceptors that increased relative to the reference are the ones Gene A's codons tend to use, favoring its translation; Sample 3 is roughly neutral. Equivalently, 2^**A** is the ratio of Gene A's expression-based tAI in the sample to that in the reference — an exact identity here because both use the same detected isoacceptors and the same demand weights: 2^0.17 ≈ 1.13, i.e. ~13%% higher effective tRNA availability for Gene A in Sample 2 versus reference Sample 1.

Every affinity comes with a **standard error**, written to a sister table (with `-STDERR` inserted before the file extension). The idea is simple: an affinity that rests on poorly covered isoacceptors is less certain, so it gets a larger error bar. The error is propagated from the tRNA-seq read counts themselves — the Poisson counting noise in each isoacceptor's log-supply, including the counts in the reference samples. Isoacceptors not detected in a sample — those whose coverage there is below `--min-coverage` (see Robustness) — are left out of that sample's affinity, i.e. treated as unchanged from the reference.

### Robustness

tRNA-seq experiments, especially from complex samples containing many populations, rarely yield a complete profile of the tRNA species in any organism. The more tRNAs measured, the more robustly affinity estimates the preference of the tRNA pool for the gene. The user can change the minimum number of detected tRNA isoacceptors required to compute an affinity (`--min-isoacceptors`, default 4). An isoacceptor contributes to a sample's affinity only if it clears the coverage threshold in that sample *and* has a defined reference value (it clears the threshold in the reference set), so the number of usable isoacceptors can differ between samples. The count-based standard error (above) measures the resulting uncertainty: an affinity resting on few or poorly covered isoacceptors carries a larger standard error, flagging it as less reliable.

A coverage threshold for detection of a tRNA isoacceptor excludes very low abundance isoacceptors, with the default being a specific coverage of 10 reads at the 3' end (discriminator nucleotide) of the isoacceptor seeds ("specific" as opposed to "nonspecific" coverage is explained in the %(trnaseq-profile-db)s artifact). The isoacceptor's mean coverage across the reference set must clear this threshold (with a single `--reference-sample` this reduces to "that sample's coverage must clear the threshold"). Per-sample isoacceptor measurements below the threshold are dropped. The flag `--shared-isoacceptors` tightens this so that an isoacceptor must clear the threshold in *every* analyzed and reference sample, fixing one common isoacceptor set across all samples — this can be useful when directly comparing raw affinity magnitudes between samples. Note that over many samples this can be very strict and may leave no isoacceptors, in which case the program stops with an explanatory error.

### Standardization

It is useful to compare affinities across genes (or functional groups) to see which are most favored or disfavored by the change in the tRNA pool. The ranking of genes by affinity can differ from sample to sample. It can also help to standardize affinities *within* a sample before comparing one sample to another: the set of isoacceptors measured in a sample affects the overall magnitude of its affinities, so one sample's values may be uniformly larger or smaller than another's even when the gene ranking is similar. `--normalize-affinity` rescales the affinities to make these comparisons cleaner, with three options:

  - `min_max`: rescale each sample's affinities to a 0-to-1 range, so the lowest-affinity gene becomes 0 and the highest becomes 1.
  - `min_max_mean`: apply `min_max`, then subtract the mean so the values are centered on zero (potential range −1 to 1).
  - `magnitude_min_max`: min-max rescale the *absolute* affinities, emphasizing the strength of change regardless of its direction.

## Output

By default the program writes a single tab-delimited **affinity table** to `--output-file`. In the default function mode its leading columns are `genome_name`, `function_source`, `function_accession`, and `function_name`; in gene mode (`--gene-affinity`) they are `genome_name` and `gene_caller_id`. Every remaining column is one analyzed sample (headed by the sample name) and holds that gene's or function's affinity in that sample. A companion **standard-error table** with the same rows and columns is written alongside it, with `-STDERR` inserted before the extension (e.g. `affinity-STDERR.txt`).

A few options reshape these files:

  - `--separate-genomes` writes one file per genome, inserting `-<genome_name>` before the extension; `--separate-function-sources` likewise inserts `-<source>`, giving `<output>-<genome>-<source>.<ext>`. When either is set, nothing is written to the bare `--output-file` path.
  - `--normalize-affinity` (described under Standardization, above) writes an *additional* table per method, inserting the method name (e.g. `-min_max_mean`) before the extension; the raw table is still written under its base name.
  - `--no-raw-affinity` suppresses the raw affinity table and its `-STDERR` companion — useful when only normalized tables or plots are wanted.

Three further options expose the quantities *behind* each affinity, so the calculation can be inspected or reconstructed (all tab-delimited):

  - `--save-isoacceptor-abundance-ratios` writes the per-sample supply ratios **T** (each isoacceptor's relative abundance over the reference) to a `-ABUNDANCE_RATIOS` table, indexed by genome, amino acid, and (modified) anticodon. This table can only be split by genome (`--separate-genomes`), never by function source.
  - `--save-isoacceptor-codon-weights` writes the per-isoacceptor demand **C** (the efficiency-weighted codon usage each isoacceptor serves) to a `-CODON_WEIGHTS` table, with the same rows as the affinity table and one column per anticodon.
  - `--save-codon-frequencies` writes the underlying per-gene or per-function codon counts to a `-CODONS` table.

The `-CODON_WEIGHTS` and `-CODONS` tables are split the same way as the affinity table, by `--separate-genomes` and/or `--separate-function-sources`.

All output files are tab-delimited regardless of the extension given to `--output-file`.

## Isoacceptor contribution analysis

### Supply and demand framing

Each affinity — for one gene in one genome and one sample — combines two lists with one entry per isoacceptor:

  - **Supply** (**x**) is the log₂ fold change in each isoacceptor's relative abundance versus the reference: xᵢ = log₂(Tᵢ), the supply vector from the Calculation section. Supply is *per-sample* — every gene in the same sample shares it. A positive value means that isoacceptor's tRNA rose relative to the reference; a negative value means it fell.
  - **Demand** (**y**) is the gene's relative codon weight (the vector **C** above, scaled so its entries sum to 1 across all of the gene's detected isoacceptors in the union of samples; the subset that enters any one sample's affinity may sum to less than 1, since some isoacceptors can be undetected there). Demand is *per-gene* — it does not depend on the sample. A larger value means the gene leans more heavily on the codons that isoacceptor reads.

The affinity is the demand-weighted sum of the supply changes, A = Σᵢ yᵢ · xᵢ: a positive affinity means the gene's codons favor the tRNAs whose supply rose relative to the reference, and a negative affinity means they favor tRNAs whose supply fell.

### Per-isoacceptor contribution decomposition

Because the affinity is just a demand-weighted sum, it splits *exactly* into one piece per isoacceptor:

**contribution = demand × supply change** (for a given gene, sample, and isoacceptor)

and the pieces add back up to exactly the affinity. Each piece is how much one isoacceptor moved the affinity: large and positive when the gene leans on an isoacceptor whose supply rose, large and negative when it leans on one whose supply fell, and near zero for isoacceptors the gene barely uses. Inspecting the pieces answers, "*which isoacceptors* drove this gene's affinity in this sample?"

An SE-normalized version, **contribution_norm = contribution ÷ SE(affinity)**, is also available. Its purpose is comparison *across* (gene, sample) pairs. A raw contribution of, say, 0.3 means something very different in a precisely measured affinity than in a noisy one; dividing by the affinity's standard error re-expresses each contribution in units of that uncertainty — roughly, how many error bars it is worth — so contributions from well- and poorly-measured pairs can be ranked or pooled on the same footing, without a few noisy pairs dominating simply because their raw numbers are large. (Within a single gene and sample the divisor is the same for every isoacceptor, so it only rescales — it does not change which isoacceptors look important *there*.) Equivalently, a (gene, sample)'s normalized contributions sum to affinity ÷ SE, the affinity's signal-to-noise ratio: the raw contributions partition how *large* an affinity is, while the normalized ones partition how *strong the evidence* for it is. It is reported as NaN wherever the affinity — and therefore its standard error — is undefined, i.e. for a gene with no detected isoacceptors in that sample.

### Outputs

Turn on this analysis with `--save-isoacceptor-contributions`. It produces four kinds of table: one detailed **long** table and three **summary views** (per sample, per gene, and global). Each kind is written for two value types — the raw `contribution` and the SE-normalized `contribution_norm` — and each summary view is also written once per chosen statistic. With the two default statistics (`mean` and `abs_mean`), a default run therefore writes **14 tables per (genome, source)**: for each value type, 1 long table + 3 summary views × 2 statistics = 7 tables, and 7 × 2 = 14.

The **long** tables have one row per (gene, sample, isoacceptor), with the gene/function columns, `trnaseq_sample_name`, `anticodon`, and a `contribution` or `contribution_norm` value. Undetected entries are kept as NaN, so the table has the same shape for every (gene, sample) pair that shares an isoacceptor set.

The three **summary views** collapse the gene dimension, the sample dimension, or both, into wide tables with one column per anticodon:

  - **per sample** summarizes over genes — one row per (genome, source, sample);
  - **per gene** summarizes over samples — one row per gene or function;
  - **global** summarizes over both — one row per (genome, source).

Each cell of a summary view is a statistic over the contributions: `mean` of the signed values, `abs_mean` of the absolute values (useful for ranking isoacceptors by overall influence, regardless of sign), or `std` (their standard deviation). The defaults are `mean` and `abs_mean`; add `std` with `--contribution-statistics`. By default both value types and all four views (`long`, `per_sample`, `per_gene`, `global`) are written; trim any of them with `--contribution-variants`, `--contribution-aggregations`, or `--contribution-statistics`.

In function mode (the default), the summaries are computed separately for each function source. A run annotated with both KEGG and Pfam, for example, gets separate KEGG and Pfam rows for the same (genome, sample); they are not merged into one cross-source value. To merge all sources into a single summary, use `--compare-all-function-sources`.

The files live in a directory `<output_root>-CONTRIBUTIONS/`, where `<output_root>` is `--output-file` with its extension removed. Inside is one subdirectory per view (`long`, `per_sample`, `per_gene`, `global`). Without `--separate-genomes`, each file holds all genomes concatenated; with `--separate-genomes`, each genome gets its own file under a `<genome>/` subdirectory. The file name records the value type (`raw` or `norm`) and, for the summary views, the statistic (e.g. `norm-mean.txt`); with `--separate-function-sources`, the source name is prepended.

So a default run writes files like `<output_root>-CONTRIBUTIONS/long/raw.txt` and `<output_root>-CONTRIBUTIONS/per_sample/norm-mean.txt`. With both `--separate-genomes` and `--separate-function-sources`, the paths become `<output_root>-CONTRIBUTIONS/long/<genome>/<source>-raw.txt` and `<output_root>-CONTRIBUTIONS/per_sample/<genome>/<source>-norm-mean.txt`. A single glob then fetches related tables: every per-sample table for every genome is `CONTRIBUTIONS/per_sample/*/`, and every output for one genome is `CONTRIBUTIONS/*/<genome>/`.

## Plots

With `--plot`, the program renders a heatmap of affinities (one PDF per genome, derived from the `--output-file` path), with functions or genes as rows and tRNA-seq samples as columns. Rows are ordered by a dendrogram that clusters functions/genes by codon composition; add `--plot-sample-dendrogram` to cluster the samples too. By default the heatmap shows only the 25 highest- and 25 lowest-affinity functions/genes, ranked by mean affinity across samples (`--n-highest-affinity` / `--n-lowest-affinity`, each default 25; set both to `-1` to show all functions/genes). For exploring all functions interactively, %(anvi-interactive)s is the better tool: it can load the affinity table together with the Newick trees written by `--save-codon-trees` and `--save-sample-trees`. An affinity table from a previous run can be re-plotted without recomputation by passing it to `--plot-affinity-file`.
