This program computes functional enrichment across groups of genomes and generates a %(functional-enrichment-txt)s file.

{:.warning}
For its sister programs, see %(anvi-compute-functional-enrichment-in-pan)s and %(anvi-compute-metabolic-enrichment)s.

{:.notice}
Please also see %(anvi-display-functions)s which can both calculate functional enrichment, AND provide an interactive interface to display the distribution of functions.

## Functional enrichment

This program can be executed using genomes described through %(external-genomes)s, %(internal-genomes)s, and/or stored in a %(genomes-storage-db)s. In addition to specifying genome sources, you must provide a %(groups-txt)s file that declares which genome belongs to which group for the enrichment analysis.

### How does it work?

1. **Aggregate functions from all sources**. Gene calls in each genome are tallied according to their functional annotations from the specified annotation source.

2. **Quantify the distribution of functions in each group of genomes**. This information is then processed by `anvi-script-enrichment-stats` to fit a Generalized Linear Model (GLM) that determines (1) the extent to which a particular functional annotation is unique to a single group and (2) the percentage of genomes in which it appears within each group. This analysis produces a %(functional-enrichment-txt)s file.

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and was first described in [this paper](https://doi.org/10.1186/s13059-020-02195-w).


### Basic usage

You can execute this program with a single source of genomes:

{{ codestart }}
anvi-compute-functional-enrichment-across-genomes -i %(internal-genomes)s \
                                                  -o %(functional-enrichment-txt)s \
                                                  -G %(groups-txt)s \
                                                  --annotation-source FUNCTION_SOURCE
{{ codestop }}

or multiple sources:

{{ codestart }}
anvi-compute-functional-enrichment-across-genomes -i %(internal-genomes)s\
                                                  -e %(external-genomes)s \
                                                  -G %(groups-txt)s \
                                                  -g %(genomes-storage-db)s \
                                                  -o %(functional-enrichment-txt)s \
                                                  --annotation-source FUNCTION_SOURCE
{{ codestop }}

### Additional Parameters

You can generate a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter:

{{ codestart }}
anvi-compute-functional-enrichment-across-genomes -i %(internal-genomes)s \
                                                  -G %(groups-txt)s \
                                                  -o %(functional-enrichment-txt)s \
                                                  --annotation-source FUNCTION_SOURCE
                                                  --functional-occurrence-table-output FUNC_OCCURRENCE.TXT
{{ codestop }}
