This program computes functional enrichment within a pangenome and returns a %(functional-enrichment-txt)s file.

For its sister programs, see %(anvi-compute-metabolic-enrichment)s and %(anvi-compute-functional-enrichment-across-genomes)s.

{:.notice}
Please also see %(anvi-display-functions)s which can both calculate functional enrichment, AND give you an interactive interface to display the distribution of functions.

## Enriched functions in a pangenome

For this to run, you must provide a %(pan-db)s and %(genomes-storage-db)s pair, as well as a %(misc-data-layers)s that associates genomes in your pan database with categorical data. The program will then find functions that are enriched in each group (i.e., functions that are associated with gene clusters that are characteristic of the genomes in that group). 

{:.notice}
Note that your %(genomes-storage-db)s must have at least one function annotation source for this to work.

This analysis will help you identify functions that are associated with a specific group of genoes in a pangenome and determine the functional core of your pangenome. For example, in the *Prochlorococcus* pangenome (the one used in [the pangenomics tutorial, where you can find more info about this program](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)), this program finds that `Exonuclease VII` is enriched in the `low-light` genomes and not in `high-light` genomes. The output file provides various statistics about how confident the program is in making this association.

### How does it work?

What this program does can be broken down into three steps:

1. **Determine groups of genomes**. The program uses a %(misc-data-layers)s variable (containing categorical, not numerical, data) to split genomes in a pangenome into two or more groups. For example, in the pangenome tutorial, the categorical variable name was `light` that partitioned genomes into `low-light` and `high-light `groups.

2.  **Determine the "functional associations" of gene clusters**. In short, this is collecting the functional annotations for all of the genes in each cluster and assigning the one that appears most frequently to represent the entire cluster.

3. **Quantify the distribution of functions in each group of genomes**. For this, the program determines to what extent a particular function is enriched in specific groups of genomes and reports it as a %(functional-enrichment-txt)s file.

{:.notice}
Check out [Alon's behind the scenes post](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome), which goes into a lot more detail.

### Basic usage

Here is the simplest way to run of this program:

{{ codestart }}
anvi-compute-functional-enrichment -p %(pan-db)s\
                                   -g %(genomes-storage-db)s \
                                   -o %(functional-enrichment-txt)s \
                                   --category-variable CATEGORY \
                                   --annotation-source FUNCTION_SOURCE
{{ codestop }}

The %(pan-db)s must contain at least one categorical data layer in %(misc-data-layers)s, and you must choose one of these categories to define your pan-groups with the `--category-variable` parameter. You can see available variables with %(anvi-show-misc-data)s program with the parameters `-t layers --debug`.

Note that by default any genomes not in a category will be ignored; you can instead include these in the analysis by using the flag `--include-ungrouped`.

The %(genomes-storage-db)s must have at least one functional annotation source, and you must choose one of these sources with the `--annotation-source`. If you do not know which functional annotation sources are available in your %(genomes-storage-db)s, you can use the `--list-annotation-sources` parameter to find out.

### Additional options

By default, gene clusters with the same functional annotation will be merged. But if you provide the `--include-gc-identity-as-function` parameter and set the annotation source to be 'IDENTITY', anvi'o will treat gene cluster names as functions and enable you to investigate enrichment of each gene cluster independently. This is how you do it:

{{ codestart }}
anvi-compute-functional-enrichment -p %(pan-db)s\
                                   -g %(genomes-storage-db)s \
                                   -o %(functional-enrichment-txt)s \
                                   --category-variable CATEGORY \
                                   --annotation-source IDENTITY \
                                   --include-gc-identity-as-function
{{ codestop }}

To output a functional occurrence table, which describes the number of times each of your functional associations occurs in each genome you're looking at, use the `--functional-occurrence-table-output` parameter, like so:

{{ codestart }}
anvi-compute-functional-enrichment -p %(pan-db)s\
                                   -g %(genomes-storage-db)s \
                                   -o %(functional-enrichment-txt)s \
                                   --category-variable CATEGORY \
                                   --annotation-source FUNCTION_SOURCE \
                                   --functional-occurrence-table-output FUNC_OCCURRENCE.TXT
{{ codestop }}