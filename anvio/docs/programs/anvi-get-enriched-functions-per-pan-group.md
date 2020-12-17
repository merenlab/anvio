This program returns a **matrix of functions that are enriched within specific groups in your pangenome**.

{:.notice}
As of `v7`, there is a more general version of this program that can compute enrichment scores for other things. See %(anvi-compute-functional-enrichment)s for details.

You provide a %(pan-db)s and %(genomes-storage-db)s pair, as well as a %(misc-data-layers)s that stores categorical data, and the program will consider each of the categories their own 'pan-group'. It will then find functions that are enriched in that group (i.e., functions that are associated with gene clusters that are characteristic of the genomes in that group). It returns this output as a %(functional-enrichment-txt)s

{:.notice}
Note that your %(genomes-storage-db)s must have at least one functional annotation source for this to work.

This helps you highlight functions or pathways that separate a specific pan-group and determine the functional core of your pangenome. For example, in the *Prochlorococcus* pangenome (the one used in [the pangenomics tutorial, where you can find more info about this program](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)), this program finds that `Exonuclease VII` is enriched in the low-light pan-group. The output file provides various statistics about how confident the program is in making this association.

### How does it work?

What this program does can be broken down into three steps:

1. Determining the pan-groups. Firstly, the program uses a %(misc-data-layers)s (containing categorical, not numerical, data) to split the pangenome into several groups. For example, in the pangenome tutorial, this was the low-light and high-light groups.
2.  Determine the "functional associations" for each of your gene clusters. In short, this is collecting the functional annotations for all of the genes in each cluster and assigning the one that appears most frequently to represent the entire cluster.
3. Looking at the functional associations and their relative levels of abundance across the pan-groups. Specifically, it looks at the level that a particular gene cluster's functional association is unique to a single pan-group and the percent of genomes it appears in in each pan-group. This is what is reported in the output matrix, a %(functional-enrichment-txt)s

If you're still curious, check out [Alon's behind the scenes post](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome), which goes into a lot more detail.

### Okay, cool. Let's talk parameters.

Here is the simplest run of this program:

{{ codestart }}
anvi-get-enriched-functions-per-pan-group -p %(pan-db)s\
                                          -g %(genomes-storage-db)s \
                                          -o %(functional-enrichment-txt)s \
                                          --category-variable CATEGORY
                                          --annotation-source FUNCTION_SOURCE
{{ codestop }}

The parameter `--category-variable` gives the name of the categorical %(misc-data-layers)s that you want to use to define your pan-groups. Note that by default any genomes not in a category will be ignored; you can instead include these in the analysis by using the flag `--include-ungrouped`

You must also provide the `--annotation-source` parameter to indicate which source of functional annotations to use. Use the parameter `--list-annotation-sources` to list the available annotation sources in your %(pan-db)s.

You can choose to not group together gene clusters with the same function by adding the parameter `--include-gc-identity-as-function` and setting the annotation source to `IDENTITY`.

#### Additional Output

You can also output a functional occurrence table, which describes the number of times each of your functional associations occurs in each genome you're looking at. You can interact more with this data by using %(anvi-matrix-to-newick)s.

You can find more information about this [here](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions).
