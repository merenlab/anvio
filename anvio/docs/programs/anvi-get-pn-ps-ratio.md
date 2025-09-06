This program **calculates the pN/pS ratio** for each gene in a %(contigs-db)s and outputs it as a %(pn-ps-data)s artifact.

### What is the pN/pS ratio?

The pN/pS ratio (first described in [Schloissnig et al. 2012](https://doi.org/10.1038/nature11711)) is the ratio of two rates: the rates of non-synonymous (pN) and synonymous (pS) **polymorphism**. It is analogous to dN/dS, which is the ratio of rates between non-synonymous (dN) and synonymous **substitutions** between two strains. We calculate pN/pS from allele frequency data obtained through SCVs (Single Codon Variants) and SAAVs (Single Amino Acid Variants). 

In molecular evolution, **non-synonymous changes** are nucleotide substitutions that alter the amino acid sequence of a protein, while **synonymous changes** are substitutions that do not change the amino acid due to the degeneracy of the genetic code. The pN/pS ratio provides insights into the selective pressures acting on genes: values significantly greater than 1 may indicate positive selection, while values significantly less than 1 suggest purifying selection. See the study by [Kiefl et al. 2023](https://www.science.org/doi/10.1126/sciadv.abq4632) for additional information, and [this reproducible workflow](https://merenlab.org/data/anvio-structure/chapter-III/) associated with that study to see use cases.

### How do I use this program?

First, you will need to run %(anvi-gen-variability-profile)s using the flag `--engine CDN` to generate a %(variability-profile-txt)s for SCVs (single codon variants), which we'll name `SCVs.txt` in this example.

Then you can run this program as follows:

{{ codestart }}
anvi-get-pn-ps-ratio -V SCVs.txt \
                     -c %(contigs-db)s \
                     -o output_dir
{{ codestop }}

A pN/pS value is calculated for each gene × sample combination. This will result in a directory called `output_dir` that contains several tables describing each of your genes. See %(pn-ps-data)s for more information.

### Other parameters

This program has several default filtering choices that you should pay attention to. You can tune these filter options with the following variables:

- The minimum departure from consensus for a variable position (`--min-departure-from-consensus`).
- The minimum departure from reference for a variable position (`--min-departure-from-reference`).
- The minimum number of SCVs in a grouping (`--minimum-num-variants`).
- The minimum coverage at a variable position (`--min-coverage`).
