This program **calculates the pN/pS ratio** for each gene in a %(contigs-db)s and outputs it as a %(pn-ps-data)s artifact.

### What is the pN/pS ratio?

The pN/pS ratio (first described in [Schloissnig et al. 2012](https://doi.org/10.1038/nature11711)) is the ratio of 2 rates: the rates of non-synonymous (pN) and synonymous (pS) **polymorphism**. It is analogous to dN/dS, which is the ratio of rates between non-synonymous (dN) and synonymous **substitutions** between two strains. We calculate pN/pS from allele frequency obtained through SCVs and SAAVs. See the study by [Kiefl et al. 2023](https://www.science.org/doi/10.1126/sciadv.abq4632) for additional information, and [this reproducible workflow](https://merenlab.org/data/anvio-structure/chapter-III/) associated with that study to see use cases.

###  How do I use this program?

First, you will need to run %(anvi-gen-variability-profile)s using the flag `--engine CDN` to get a %(variability-profile-txt)s for SCVs (single codon variants), which we'll name `SCVs.txt` in this example.

Then you can run this program like so:

{{ codestart }}
anvi-get-pn-ps-ratio -V SCVs.txt \
                     -c %(contigs-db)s \
                     -o output_dir
{{ codestop }}

A pN/pS value is calculated for each gene x sample combo. This will result in a directory called `output_dir` that contains several tables that describe each of your genes. See %(pn-ps-data)s for more information.

### Interpreting the output

The output directory will contain five tab-delimited files. By default these are in long format (one row per gene-sample combination), but if you use the `--pivot` flag, they will be in matrix format (genes as rows, samples as columns).

Here is what each file contains:

- **`pNpS.txt`**: The pN/pS ratio for each gene in each sample. This is the primary output of the program. The general right-hand rule for a broad interpretation is that the values greater than 1 suggest diversifying selection, values less than 1 suggest purifying selection, and values around 1 suggest neutral evolution. Some entries may be empty, which means one of two things: (1) the gene had fewer SCVs than the `--minimum-num-variants` threshold set by default (default: 4), so the ratio was not considered reliable enough to report, or (2) both pN and pS were zero, making the ratio undefined. Values of `inf` indicate that pN was greater than zero but pS was zero (i.e., all observed polymorphisms were non-synonymous).

- **`pN.txt`**: The rate of *non-synonymous polymorphism per non-synonymous site* for each gene in each sample. This is the sum of non-synonymous fractions across all SCVs in a gene, normalized by the number of non-synonymous sites in that gene.

- **`pS.txt`**: The same as above but opposite as this one is the rate of *synonymous polymorphism per synonymous site* for each gene in each sample. This is the sum of synonymous fractions across all SCVs in a gene, normalized by the number of synonymous sites in that gene.

- **`num_SCVs.txt`**: The *number of single codon variants observed* for each gene in each sample. This is useful for understanding the statistical support behind each pN/pS estimate. When the `--pivot` flag is used, genes with no SCVs in a given sample will show 0.

- **`potentials.txt`**: The *number of synonymous and non-synonymous sites* for each gene (determined by the gene's codon composition). These values represent the denominators used to normalize pN and pS. A gene's potential depends on which codons it contains: for example, a methionine codon (ATG) contributes zero synonymous sites because any single-nucleotide change to it is non-synonymous, whereas a leucine codon (CTN) contributes more synonymous sites. The potentials also depend on the `--comparison` parameter, since the comparison codon determines what is considered synonymous or non-synonymous.

### Setting the right parameters

Please keep in mind that %(anvi-get-pn-ps-ratio)s has a few filtering choices that are set by default and may become critical for your investigation, including,

- The minimum departure from consensus for a variable position (`--min-departure-from-consensus`).
- The minimum departure from reference for a variable position (`--min-departure-from-reference`).
- The minimum number of SCVs in a grouping (`--minimum-num-variants`).
- The minimum coverage at a variable position (`--min-coverage`).

Please consider finetuning these parameters according to your research question.