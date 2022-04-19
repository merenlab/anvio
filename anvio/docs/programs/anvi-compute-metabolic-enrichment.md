This program computes metabolic module enrichment across groups of genomes or metagenomes and returns a %(functional-enrichment-txt)s file (throughout this text, we will use the term genome to describe both for simplicity).

{:.warning}
For its sister programs, see %(anvi-compute-functional-enrichment-in-pan)s and %(anvi-compute-functional-enrichment-across-genomes)s.

## Module enrichment

To run this program, you must already have estimated the completeness of metabolic modules in your genomes using the program %(anvi-estimate-metabolism)s and obtained a "modules" mode output file (which is the default output mode of that program). In addition to that, you will need to provide a %(groups-txt)s file to declare which genome belongs to which group for enrichment analysis to consider.

### How does it work?

1. **Determine the presence of modules**. Each module in the "modules" mode output has a completeness score associated with it in each genome, and any module with a completeness score over a given threshold (set by `--module-completion-threshold`) will be considered to be *present* in that genome.

2. **Quantify the distribution of modules in each group of genomes**. The distribution of a given module across genomes in each group will determine its enrichment. This is done by fitting a generalized linear model (GLM) with a logit linkage function in `anvi-script-enrichment-stats`, and it produces a %(functional-enrichment-txt)s file.

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and described first in [this paper](https://doi.org/10.1186/s13059-020-02195-w).

### Basic usage

See %(kegg-metabolism)s or %(user-metabolism)s for more information on how to generate a "modules" mode output format from %(anvi-estimate-metabolism)s. Please note that the genome names in the modules file must match those that you will mention in the %(groups-txt)s file.

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s
{{ codestop }}

### Additional parameters

The default completeness threshold for a module to be considered 'present' in a genome is 0.75 (=75%%). If you wish to change this, you can do so by providing a different threshold between (0, 1], using the `--module-completion-threshold` parameter:

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --module-completion-threshold 0.9
{{ codestop }}

By default, this program uses the [pathwise completeness score](https://anvio.org/help/main/programs/anvi-estimate-metabolism/#two-estimation-strategies-pathwise-and-stepwise) to determine which modules are 'present' in a genome, but you can ask it to use stepwise completeness instead by using the `--use-stepwise-completeness` flag.

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --use-stepwise-completeness
{{ codestop }}

By default, the column containing genome names in your MODULES.TXT file will have the header `db_name`, **but there are certain cases in which you might have them in a different column name for your genomes or metagenomes** (such as those cases where you did not run %(anvi-estimate-metabolism)s in multi-mode). In those cases, you can tell this program to look for a *different* column name to find your genomes or metagenomes using the `--sample-header`. For example, if your metagenome names are listed under the `metagenome_name` column, you would do the following:

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --sample-header metagenome_name
{{ codestop }}

If you ran %(anvi-estimate-metabolism)s on a bunch of extra genomes but only want to include a subset of them in the %(groups-txt)s, that is fine. By default, any samples from the `MODULES.TXT` file that are missing from the %(groups-txt)s will be **ignored**. However, there is also an option to include those missing samples in the analysis, as one big group called 'UNGROUPED'. To do this, you can use the `--include-samples-missing-from-groups-txt` parameter.

{{ codestart }}
anvi-compute-metabolic-enrichment -M MODULES.TXT \
                                  -G %(groups-txt)s \
                                  -o %(functional-enrichment-txt)s \
                                  --include-samples-missing-from-groups-txt
{{ codestop }}
