This program has multiple abilities. It can compute enriched functions across categories in a pangenome, enriched metabolic modules across groups of samples, or enriched functions across groups of genomes. To do this it relies on the script `anvi-script-enrichment-stats` by [Amy Willis](https://github.com/adw96).

Regardless of the situation, it returns a **matrix of things that are enriched within specific groups in your dataset**.

## Enriched functions in a pangenome

This option achieves the same thing as the program %(anvi-get-enriched-functions-per-pan-group)s. You can check that page for additional details and helpful tips, but below you will find some usage information.

### Basic usage

You must provide this program with a %(pan-db)s and its corresponding %(genomes-storage-db)s. The %(pan-db)s must contain at least one categorical data layer in %(misc-data-layers)s, and you must choose one of these categories to group your genomes with the `--category-variable` parameter. The %(genomes-storage-db)s must have at least one functional annotation source, and you must choose one of these sources with the `--annotation-source`. You must also provide an output file name.

{{ codestart }}
anvi-compute-enrichment-scores -p %(pan-db)s\
                               -g %(genomes-storage-db)s \
                               -o %(functional-enrichment-txt)s \
                               --category-variable CATEGORY \
                               --annotation-source FUNCTION_SOURCE
{{ codestop }}

If you do not know which functional annotation sources are available in your %(genomes-storage-db)s, you can use the `--list-annotation-sources` parameter to find out.

### Additional parameters

By default, gene clusters with the same functional annotation will be merged. But if you provide the `--include-gc-identity-as-function` parameter and set the annotation source to be 'IDENTITY', anvi'o will treat gene cluster names as functions and enable you to investigate enrichment of each gene cluster independently. This is how you do it:

{{ codestart }}
anvi-compute-enrichment-scores -p %(pan-db)s\
                               -g %(genomes-storage-db)s \
                               -o %(functional-enrichment-txt)s \
                               --category-variable CATEGORY \
                               --annotation-source IDENTITY \
                               --include-gc-identity-as-function
{{ codestop }}

If you provide the `--exclude-ungrouped` parameter, then genomes without a category in the provided `--category-variable` will be excluded from the analysis. (By default, these genomes go into their own 'ungrouped' category.)

You can get a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter, like so:

{{ codestart }}
anvi-compute-enrichment-scores -p %(pan-db)s\
                               -g %(genomes-storage-db)s \
                               -o %(functional-enrichment-txt)s \
                               --category-variable CATEGORY \
                               --annotation-source FUNCTION_SOURCE \
                               --functional-occurrence-table-output functional_occurrence.txt
{{ codestop }}


## Enriched modules

## Enriched functions in groups of genomes
