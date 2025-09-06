This program computes metabolic module enrichment across groups of genomes or metagenomes and generates a %(functional-enrichment-txt)s file (throughout this text, we will use the term genome to describe both for simplicity).

{:.warning}
For its sister programs, see %(anvi-compute-functional-enrichment-in-pan)s and %(anvi-compute-functional-enrichment-across-genomes)s.

## Module enrichment

To execute this program, you must have previously estimated the completeness of metabolic modules in your genomes using the program %(anvi-estimate-metabolism)s and obtained a "modules" mode output file (which is the default output mode of that program). In addition, you must provide a %(groups-txt)s file that declares which genome belongs to which group for the enrichment analysis.

### How does it work?

1. **Determine the presence of modules**. Each module in the "modules" mode output has an associated completeness score for each genome, and any module with a completeness score exceeding the specified threshold (set by `--module-completion-threshold`) will be considered *present* in that genome.

2. **Quantify the distribution of modules in each group of genomes**. The distribution of a given module across genomes in each group determines its enrichment score. This is accomplished by fitting a generalized linear model (GLM) with a logit linkage function in `anvi-script-enrichment-stats`, which produces a %(functional-enrichment-txt)s file.

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and was first described in [this paper](https://doi.org/10.1186/s13059-020-02195-w).

### Basic usage

See %(kegg-metabolism)s or %(user-metabolism)s for more information on how to generate a "modules" mode output format from %(anvi-estimate-metabolism)s. Please note that the genome names in the modules file must match those specified in the %(groups-txt)s file.

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s
{{ codestop }}

### Additional parameters

The default completeness threshold for a module to be considered 'present' in a genome is 0.75 (75%%). If you wish to modify this threshold, you can specify a different value between (0, 1] using the `--module-completion-threshold` parameter:

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --module-completion-threshold 0.9
{{ codestop }}

By default, this program uses the [pathwise completeness score](https://anvio.org/help/main/programs/anvi-estimate-metabolism/#two-estimation-strategies---pathwise-and-stepwise) to determine which modules are 'present' in a genome, but you can specify stepwise completeness instead by using the `--use-stepwise-completeness` flag.

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --use-stepwise-completeness
{{ codestop }}

By default, the column containing genome names in your MODULES.TXT file will have the header `db_name`, **but there are certain cases where you might have genome names in a different column** (such as cases where you did not run %(anvi-estimate-metabolism)s in multi-mode). In those cases, you can specify a different column name using the `--sample-header` parameter. For example, if your metagenome names are listed under the `metagenome_name` column, you would execute the following:

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --sample-header metagenome_name
{{ codestop }}

If you ran %(anvi-estimate-metabolism)s on additional genomes but only want to include a subset of them in the %(groups-txt)s, this is acceptable. By default, any samples from the `MODULES.TXT` file that are missing from the %(groups-txt)s will be **ignored**. However, you can also include those missing samples in the analysis as one large group called 'UNGROUPED'. To enable this behavior, use the `--include-samples-missing-from-groups-txt` parameter.

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --include-samples-missing-from-groups-txt
{{ codestop }}
