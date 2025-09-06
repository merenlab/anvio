This program computes functional enrichment within a pangenome and generates a %(functional-enrichment-txt)s file.

{:.warning}
For its sister programs, see %(anvi-compute-metabolic-enrichment)s and %(anvi-compute-functional-enrichment-across-genomes)s.

{:.notice}
Please also see %(anvi-display-functions)s which can both calculate functional enrichment, AND provide an interactive interface to display the distribution of functions.

## Enriched functions in a pangenome

To execute this program, you must provide a %(pan-db)s and %(genomes-storage-db)s pair, along with %(misc-data-layers)s that associates genomes in your pangenome database with categorical data. The program will identify functions that are enriched in each group (i.e., functions associated with gene clusters that are characteristic of the genomes in that group). 

{:.notice}
Note that your %(genomes-storage-db)s must contain at least one functional annotation source for this analysis to work.

This analysis helps identify functions associated with specific groups of genomes in a pangenome and determines the functional core of your pangenome. For example, in the *Prochlorococcus* pangenome (analyzed in [the pangenomics tutorial, where you can find more information about this program](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)), this program identifies that `Exonuclease VII` is enriched in the `low-light` genomes and not in `high-light` genomes. The output file provides various statistics about the confidence of this functional association.

### How does it work?

The analysis performed by this program can be broken down into three steps:

1. **Determine groups of genomes**. The program uses a %(misc-data-layers)s variable (containing categorical, not numerical, data) to partition genomes in a pangenome into two or more groups. For example, in the pangenome tutorial, the categorical variable named `light` partitioned genomes into `low-light` and `high-light` groups.

2.  **Determine the "functional associations" of gene clusters**. This step involves collecting the functional annotations for all genes within each cluster and assigning the most frequently occurring annotation to represent the entire cluster.

3. **Quantify the distribution of functions in each group of genomes**. The program determines the extent to which particular functions are enriched in specific groups of genomes and reports this information as a %(functional-enrichment-txt)s file. This analysis is performed by executing the script `anvi-script-enrichment-stats`. 

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and was first described in [this paper](https://doi.org/10.1186/s13059-020-02195-w).

{:.notice}
Check out [Alon's behind the scenes post](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome), which provides a more detailed explanation of this process.

### Basic usage

The simplest way to execute this program is as follows:

{{ codestart }}
anvi-compute-functional-enrichment-in-pan -p %(pan-db)s\
                                          -g %(genomes-storage-db)s \
                                          -o %(functional-enrichment-txt)s \
                                          --category-variable CATEGORY \
                                          --annotation-source FUNCTION_SOURCE
{{ codestop }}

The %(pan-db)s must contain at least one categorical data layer in %(misc-data-layers)s, and you must select one of these categories to define your pangenome groups using the `--category-variable` parameter. You can view available variables using the %(anvi-show-misc-data)s program with the parameters `-t layers --debug`.

The %(genomes-storage-db)s must contain at least one functional annotation source, and you must specify one of these sources using the `--annotation-source` parameter. If you are unsure which functional annotation sources are available in your %(genomes-storage-db)s, you can use the `--list-annotation-sources` parameter to identify them.

### Additional options

By default, gene clusters with identical functional annotations will be merged. However, if you provide the `--include-gc-identity-as-function` parameter and set the annotation source to 'IDENTITY', anvi'o will treat gene cluster names as functions and enable you to investigate the enrichment of each gene cluster independently. This is accomplished as follows:

{{ codestart }}
anvi-compute-functional-enrichment-in-pan -p %(pan-db)s\
                                          -g %(genomes-storage-db)s \
                                          -o %(functional-enrichment-txt)s \
                                          --category-variable CATEGORY \
                                          --annotation-source IDENTITY \
                                          --include-gc-identity-as-function
{{ codestop }}

To generate a functional occurrence table that describes the number of times each functional association occurs in each genome under analysis, use the `--functional-occurrence-table-output` parameter as shown below:

{{ codestart }}
anvi-compute-functional-enrichment-in-pan -p %(pan-db)s\
                                          -g %(genomes-storage-db)s \
                                          -o %(functional-enrichment-txt)s \
                                          --category-variable CATEGORY \
                                          --annotation-source FUNCTION_SOURCE \
                                          --functional-occurrence-table-output FUNC_OCCURRENCE.TXT
{{ codestop }}
