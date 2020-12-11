This program has multiple abilities. It can compute enriched functions across categories in a pangenome, enriched metabolic modules across groups of samples, or enriched functions across groups of genomes. To do this it relies on the script `anvi-script-enrichment-stats` by [Amy Willis](https://github.com/adw96).

Regardless of the situation, it returns a **matrix of things that are enriched within specific groups in your dataset**, as a %(functional-enrichment-txt)s file.

## Input option 1: Enriched functions in a pangenome

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

You can get a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter, like so:

{{ codestart }}
anvi-compute-enrichment-scores -p %(pan-db)s\
                               -g %(genomes-storage-db)s \
                               -o %(functional-enrichment-txt)s \
                               --category-variable CATEGORY \
                               --annotation-source FUNCTION_SOURCE \
                               --functional-occurrence-table-output FUNC_OCCURRENCE.TXT
{{ codestop }}

## Input option 2: Enriched modules

This option computes enrichment scores for metabolic modules in groups of samples. In order to do this, you must already have estimated completeness of metabolic modules in your samples using %(anvi-estimate-metabolism)s and obtained a "modules" mode output file (the default). You must provide that file to this program along with a %(groups-txt)s file indicating which samples belong to which groups.

### How it works

1. Determining presence of modules. Each module in the "modules" mode output has a completeness score associated with it in each sample, and any module with a completeness score over a given threshold (set by `--module-completion-threshold`) will be considered to be present in that sample.
2. Examining the distribution of modules in each group of samples to compute an enrichment score for each module. This is done by fitting a generalized linear model (GLM) with a logit linkage function in `anvi-script-enrichment-stats`, and it produces a %(functional-enrichment-txt)s file.

### Basic usage

See %(kegg-metabolism)s for more information on the "modules" mode output format from %(anvi-estimate-metabolism)s, which you must provide with the `-M` flag. The sample names in this file must match those in the %(groups-txt)s file, provided with `-G`. You must also provide the name of the output file.

{{ codestart }}
anvi-compute-enrichment-scores -M MODULES.TXT \
                               -G %(groups-txt)s \
                               -o %(functional-enrichment-txt)s
{{ codestop }}

### Additional parameters

The default completeness threshold for a module to be considered 'present' in a sample is 0.75 (75 percent). If you wish to change this, you can do so by providing a different threshold - as a number in the range (0, 1] - using the `--module-completion-threshold` parameter. For example:

{{ codestart }}
anvi-compute-enrichment-scores -M MODULES.TXT \
                               -G %(groups-txt)s \
                               -o %(functional-enrichment-txt)s \
                               --module-completion-threshold 0.9
{{ codestop }}

By default, the column containing sample names in your MODULES.TXT file will have the header `db_name`, but there are certain cases in which you might have them in a different column - for example, if you did not run %(anvi-estimate-metabolism)s in multi-mode. In those cases, you can specify that a different column contains the sample names by providing its header with `--sample-header`. For example, if you sample names were in the `metagenome_name` column, you would do the following:

{{ codestart }}
anvi-compute-enrichment-scores -M MODULES.TXT \
                               -G %(groups-txt)s \
                               -o %(functional-enrichment-txt)s \
                               --sample-header metagenome_name
{{ codestop }}

If you ran %(anvi-estimate-metabolism)s on a bunch of extra samples but only want to include a subset of those samples in the %(groups-txt)s, that is fine - by default any samples from the MODULES.TXT file that are missing from the %(groups-txt)s will be ignored. However, there is also an option to include those missing samples in the analysis, as one big group called 'UNGROUPED'. To do this, you can use the --include-samples-missing-from-groups-txt parameter. Just be careful that if you are also using the --include-ungrouped flag (see below), any samples without a specified group in the %(groups-txt)s will also be included in the 'UNGROUPED' group.

{{ codestart }}
anvi-compute-enrichment-scores -M MODULES.TXT \
                               -G %(groups-txt)s \
                               -o %(functional-enrichment-txt)s \
                               --include-samples-missing-from-groups-txt
{{ codestop }}


## Input option 3: Enriched functions in groups of genomes

You are not limited to computing functional enrichment in pangenomes, you can do it for regular genomes, too. This option takes either external or internal genomes (or both) which are organized into groups, and computes enrichment scores and associated groups for annotated functions in those genomes.

### How it works

This is similar to computing functional enrichment in pangenomes (as described in %(anvi-get-enriched-functions-per-pan-group)s), but a bit simpler.

1. Counting functions. Gene calls in each genome are tallied according to their functional annotations from the given annotation source.
2. Looking at the functions and their relative levels of abundance across the groups of genomes. This again uses `anvi-script-enrichment-stats` to fit a GLM to determine A) the level that a particular functional annotation is unique to a single group and B) the percent of genomes it appears in in each group. This produces a %(functional-enrichment-txt)s file.

### Basic usage

You can provide either an %(external-genomes)s file or an %(internal-genomes)s file or both, but no matter what these files must contain a `group` column which indicates the group that each genome belongs to. Similar to option 1, you must also provide an annotation source from which to extract the functional annotations of interest. In the example below, we provide both types of input files.

{{ codestart }}
anvi-compute-enrichment-scores -i %(internal-genomes)s\
                               -e %(external-genomes)s \
                               -o %(functional-enrichment-txt)s \
                               --annotation-source FUNCTION_SOURCE
{{ codestop }}

### Additional Parameters

Also similar to option 1, you can get a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter:

{{ codestart }}
anvi-compute-enrichment-scores -i %(internal-genomes)s\
                               -e %(external-genomes)s \
                               -o %(functional-enrichment-txt)s \
                               --annotation-source FUNCTION_SOURCE
                               --functional-occurrence-table-output FUNC_OCCURRENCE.TXT
{{ codestop }}


## Parameters common to all options

If you provide the `--include-ungrouped` parameter, then genomes (or samples) without a group will be included from the analysis. (By default, these genomes/samples are ignored.) For the pangenome case, these genomes are those without a category in the provided `--category-variable`. For metabolic modules or the genomes in groups case, these samples/genomes are those with an empty value in the 'group' column (of either the %(groups-txt)s or the %(external-genomes)s/%(internal-genomes)s files).


## More information on `anvi-script-enrichment-stats`

This program serves as the interface to `anvi-script-enrichment-stats`, an R script which performs an enrichment test on your input. You will find a brief description of how this script works in Alon's "Behind the Scenes" note in [the pangenomics tutorial](https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome). Better yet, check out the methods section of Alon's paper, found as a pre-print [here](https://www.biorxiv.org/content/10.1101/2020.04.29.069278v2).
