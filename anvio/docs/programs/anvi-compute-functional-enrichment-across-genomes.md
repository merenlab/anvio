This program computes functional enrichment across groups of genomes and returns a %(functional-enrichment-txt)s file.

{:.warning}
For its sister programs, see %(anvi-compute-functional-enrichment-in-pan)s and %(anvi-compute-metabolic-enrichment)s.

{:.notice}
Please also see %(anvi-display-functions)s which can both calculate functional enrichment, AND give you an interactive interface to display the distribution of functions.

## Functional enrichment

You can use this program by combining genomes described through %(external-genomes)s, %(internal-genomes)s, and/or stored in a %(genomes-storage-db)s. In addition to sources for your genomes, you will need to provide a %(groups-txt)s file to declare which genome belongs to which group for enrichment analysis to consider.

### How does it work?

1. **Aggregate functions from all sources**. Gene calls in each genome are tallied according to their functional annotations from the given annotation source.

2. **Quantify the distribution of functions in each group of genomes**. This information is then used by `anvi-script-enrichment-stats` to fit a GLM to determine (1) the level that a particular functional annotation is unique to a single group and (2) the percent of genomes it appears in in each group. This produces a %(functional-enrichment-txt)s file.

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and described first in [this paper](https://doi.org/10.1186/s13059-020-02195-w).


### Basic usage

You can use it with a single source of genomes:

{{ codestart }}
anvi-compute-functional-enrichment-across-genomes -i %(internal-genomes)s \
                                                  -o %(functional-enrichment-txt)s \
                                                  -G %(groups-txt)s \
                                                  --annotation-source FUNCTION_SOURCE
{{ codestop }}

or many:

{{ codestart }}
anvi-compute-functional-enrichment-across-genomes -i %(internal-genomes)s\
                                                  -e %(external-genomes)s \
                                                  -G %(groups-txt)s \
                                                  -g %(genomes-storage-db)s \
                                                  -o %(functional-enrichment-txt)s \
                                                  --annotation-source FUNCTION_SOURCE
{{ codestop }}

### Additional Parameters

You can get a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter:

{{ codestart }}
anvi-compute-functional-enrichment-across-genomes -i %(internal-genomes)s \
                                                  -G %(groups-txt)s \
                                                  -o %(functional-enrichment-txt)s \
                                                  --annotation-source FUNCTION_SOURCE
                                                  --functional-occurrence-table-output FUNC_OCCURRENCE.TXT
{{ codestop }}
