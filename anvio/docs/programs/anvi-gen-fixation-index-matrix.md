
This program generates a matrix of the pairwise fixation indices (F<sub>ST</sub>) between your samples.

### What's a fixation index?

As described [in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst), the fixation index is a measure of the distance between two populations, based on their sequence variants (usually SNVs). Specifically, the fixation index is the ratio between the variance in allele frequency between subpopulations and the variance in the total population. 


The fixation index has its own [Wikipedia page](https://en.wikipedia.org/wiki/Fixation_index) and is a special case of [F-statistics](https://en.wikipedia.org/wiki/F-statistics). 


In anvi'o, the fixation index is calculated in accordance with [Schloissnig et al.  (2013)](https://doi.org/10.1038/nature11711)'s work to allow variant positions with multiple competing alleles.


## Anvi-gen-fixation-index 

There are two ways to run this program.  

### Input 1: Variability Profile 

The simplest one is the one shown [in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#measuring-distances-between-metagenomes-with-fst): just provide a %(variability-profile)s, like so: 

{{ codestart }}
anvi-gen-fixation-index-matrix --variability-profile %(variability-profile)s \
                               --output-file my_matrix.txt
{{ codestop }}

This will use the information in your %(variability-profile-txt)s to generate the fixation index for each of the pairwise sample comparisons, and store the results in a %(fixation-index-matrix)s named `my_matrix.txt`.  

### Input 2: Anvi'o databases

Instead of providing a %(variability-profile)s, you can instead provide the inputs to %(anvi-gen-variability-profile)s and let anvi'o do all of the work for you. Specifically, this means providing a %(contigs-db)s and %(profile-db)s pair to find your variability positions and a specific subset to focus on in any of these ways: 

- Provide a list of gene caller IDs (as a parameter with the flag `--gene-caller-ids` or in a file with one ID per line with the flag `--genes-of-interest`)
- Provide a list of splits (in a %(splits-txt)s)
- Provide a %(collection)s and %(bin)s

Additionally, you can add structural annotations by inputting a %(structure-db)s (and focus only on genes with structural annotations with the flag `--only-if-structure`) or choose to focus on only a subset of your samples by providing a file of samples of interest.  

When doing this, you can also set the variability engine to get the fixation index for SCVs (`--engine CDN`) or SAAVs (`--engine AA`). 

You can find more information about these parameters on the page for %(anvi-gen-variability-profile)s. 

### Additional Parameters

While a fixation index is usually between 0 and 1, it is possible for an index to be negative (usually because of out-breeding). By default, anvi'o sets these negative values to 0, but you can choose to keep the negative values with the flag `--keep-negatives` 

