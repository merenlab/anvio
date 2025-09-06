This program generates a matrix of pairwise fixation indices (F<sub>ST</sub>) between your samples.

### What is a fixation index?

As described [in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst), the fixation index is a measure of genetic distance between two populations based on their sequence variants (usually single nucleotide variants, or SNVs). Specifically, the fixation index represents the ratio between the variance in allele frequency between subpopulations and the variance in the total population. 


The fixation index has its own [Wikipedia page](https://en.wikipedia.org/wiki/Fixation_index) and is a special case of [F-statistics](https://en.wikipedia.org/wiki/F-statistics). 


In anvi'o, the fixation index is calculated in accordance with [Schloissnig et al. (2013)](https://doi.org/10.1038/nature11711)'s methodology to accommodate variant positions with multiple competing alleles.


## Using anvi-gen-fixation-index-matrix

There are two ways to execute this program.  

### Input 1: Variability Profile 

The simplest approach is the one shown [in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst): provide a %(variability-profile)s, as follows: 

{{ codestart }}
anvi-gen-fixation-index-matrix --variability-profile %(variability-profile)s \
                               --output-file my_matrix.txt
{{ codestop }}

This will use the information in your %(variability-profile-txt)s to generate the fixation index for each pairwise sample comparison and store the results in a %(fixation-index-matrix)s named `my_matrix.txt`.  

### Input 2: Anvi'o databases

Instead of providing a %(variability-profile)s, you can provide the inputs to %(anvi-gen-variability-profile)s and allow anvi'o to perform all the work for you. Specifically, this means providing a %(contigs-db)s and %(profile-db)s pair to identify your variability positions and a specific subset to focus on in any of these ways: 

- Provide a list of gene caller IDs (as a parameter with the flag `--gene-caller-ids` or in a file with one ID per line using the flag `--genes-of-interest`)
- Provide a list of splits (in a %(splits-txt)s)
- Provide a %(collection)s and %(bin)s

Additionally, you can add structural annotations by providing a %(structure-db)s (and focus only on genes with structural annotations using the flag `--only-if-structure`) or choose to focus on only a subset of your samples by providing a file of samples of interest.  

When using this approach, you can also set the variability engine to calculate the fixation index for single codon variants (SCVs) (`--engine CDN`) or single amino acid variants (SAAVs) (`--engine AA`). 

You can find more information about these parameters on the page for %(anvi-gen-variability-profile)s. 

### Additional Parameters

While a fixation index is usually between 0 and 1, it is possible for an index to be negative (usually because of outbreeding). By default, anvi'o sets these negative values to 0, but you can choose to preserve the negative values using the flag `--keep-negatives`. 

