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

### Other parameters

This program has some default filtering choices that you should pay mind to. You can tune these filter options with the following variables:

- The minimum departure from consensus for a variable position (`--min-departure-from-consensus`).
- The minimum departure from reference for a variable position (`--min-departure-from-reference`).
- The minimum number of SCVs in a grouping (`--minimum-num-variants`).
- The minimum coverage at a variable position (`--min-coverage`).
