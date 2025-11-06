This program **finds tRNA seed sequences from a set of tRNA-seq samples**.

This program follows %(anvi-trnaseq)s in the [trnaseq-workflow](../../workflows/trnaseq/). %(anvi-trnaseq)s is run on each tRNA-seq sample, producing sample %(trnaseq-db)ss. A tRNA-seq database contains predictions of tRNA sequences, structures, and modification sites in the sample. anvi-merge-trnaseq takes as input the tRNA-seq databases from a set of samples. It compares tRNAs predicted from the samples, finding those in common and calculating their sample coverages. The final tRNA sequences predicted from all samples are called **tRNA seeds** and function like contigs in metagenomic experiments. Seeds are stored in a %(trnaseq-contigs-db)s and sample coverages are stored in a %(trnaseq-profile-db)s. These databases are **variants** of normal %(contigs-db)ss and %(profile-db)ss, performing similar functions in the anvi'o ecosystem but containing somewhat different information.

Most of the heavy computational work in the [trnaseq-workflow](../../workflows/trnaseq/) is performed by %(anvi-trnaseq)s. anvi-merge-trnaseq is meant run relatively quickly, allowing its parameters to be tuned to fit the dataset.

The `anvi-merge-trnaseq --help` menu provides detailed explanations of the parameters controlling the multifacted analyses performed by the program.

## Key parameters

### Number of reported seeds

One key parameter is the number of reported tRNA seed sequences (`--max-reported-trna-seeds`). The default value of 10,000 seeds is more appropriate for a complex microbial community than a pure culture of a bacterial isolate, which should yield a number of tRNA seeds equal to the number of expressed tRNAs, say ~30. Sequence artifacts may be reported in addition to the 30 actual tRNAs with a higher value like 10,000. Artifacts are relatively common despite intensive screening by %(anvi-trnaseq)s and anvi-merge-trnaseq due to nontemplated nucleotides and modification-induced mutations introduced into tRNA-seq reads by reverse transcription. In practice, artifacts are easy to distinguish from true tRNA seeds by analyzing seed coverage in %(anvi-interactive)s and checking seed homology to reference databases, among other measures.

### Modification filters

Other key parameters, `--min-variation` and `--min-third-fourth-nt`, determine the coverage cutoffs that distinguish predicted positions of modified nucleotides from single nucleotide variants. Compared to SNVs, modifications typically produce higher nucleotide variability to three or four different nucleotides. However, modification-induced mutations are often highly skewed to one other nucleotide rather than all three mutant nucleotides. Furthermore, the high coverage of seeds in many tRNA-seq libraries can uncover SNVs with a low-frequency third nucleotide rather than the expected two. Some SNVs that are wrongly called modifications can be easily spotted in %(anvi-interactive)s and the output of %(anvi-plot-trnaseq)s due to covariation at two positions in the seed as a result of base pairing. In other words, SNV frequencies are equivalent at the two base paired positions in every sample, where modification artifacts have no effect on nucleotide variability at another position across the molecule.

## Examples

*Merge two samples.*

{{ codestart }}
anvi-merge-trnaseq trnaseq_database_1 trnaseq_database_2 (...) \
                   -o OUTPUT_DIRECTORY \
                   -n PROJECT_NAME \
{{ codestop }}

*Merge two samples with and without demethylase treatment, giving priority to the demethylase split in calling the underlying nucleotide at modified positions.*

{{ codestart }}
anvi-merge-trnaseq untreated_trnaseq_database demethylase_trnaseq_database (...) \
                   -o OUTPUT_DIRECTORY \
                   -n PROJECT_NAME \
                   --preferred-treatment demethylase
{{ codestop }}
