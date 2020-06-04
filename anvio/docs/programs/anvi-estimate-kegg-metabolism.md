%(anvi-estimate-kegg-metabolism)s predicts the metabolic capabilities of organisms based on their genetic content. It relies upon %(kegg-functions)s and metabolism information from the KEGG resource, which is stored in %(kegg-db)s.

The metabolic pathways that this program currently considers are those defined by KOs in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). Each KO represents a gene function, and a KEGG module is a set of KOs that collectively carry out the steps in a metabolic pathway. Therefore, for this to work, you need to have annotated your %(contigs-db)s with hits to the KEGG KOfam database by running %(anvi-run-kegg-kofams) prior to using this program.

Given a properly annotated %(contigs-db)s, this program determines which KOs are present and from those determines the completeness of each KEGG module. The results are described in a set of output text files, collectively referred to as %(kegg-metabolism)s.

## Running metabolism estimation on a single contigs database

There are several possible inputs to this program. For single genomes - isolate genomes or MAGs - you can provide a %(contigs-db)s. If your %(contigs-db)s describes a metagenome rather than a single genome, you can provide the flag `--metagenome-mode`. In metagenome mode, KOfam hits in the %(contigs-db)s are analyzed as though they belong to one collective genome, despite the fact that the sequences represent multiple different populations. Alternatively, if you have binned your metagenome sequences into separate populations and would like metabolism estimation to be run separately on each bin, you can provide a %(profile-db)s and a %(collection)s.

### Estimation for a single genome

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db
{{ codestop }}

### Estimation for a metagenome

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --metagenome-mode
{{ codestop }}

### Estimation for bins in a metagenome

You can estimate metabolism for each bin in a %(collection)s:

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME
{{ codestop }}

You can also provide a specific %(bin)s in that %(collection)s to run on:

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME -b BIN_NAME
{{ codestop }}

Or, you can provide a specific list of bins in a text file:

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME -B bin_ids.txt
{{ codestop }}

Each line in the `bin_ids.txt` file should be a bin name from the %(collection)s.

## Running metabolism estimation on multiple contigs databases

If you have a set of contigs databases of the same type (ie, all of them are single genomes, or all are binned metagenomes), you can analyze them all at once. What you need to do is put the relevant information for each %(contigs-db)s into a text file and pass that text file to %(anvi-estimate-kegg-metabolism)s. The program will then run estimation individually on each contigs database in the file. The estimation results for each database will be aggregated and printed to the same output file(s).

### Estimation for multiple single genomes

Multiple single genomes (also known as %(external-genomes)s) can be analyzed with the same command by providing an external genomes file to %(anvi-estimate-kegg-metabolism)s. To see the required format for the external genomes file, see %(external-genomes)s.

{{ codestart }}
anvi-estimate-kegg-metabolism -e external-genomes.txt
{{ codestop }}

### Estimation for multiple metagenomes

Multiple metagenomes can be analyzed with the same command by providing a metagenomes input file. Metagenome mode will be used to analyze each contigs database in the file.

{{ codestart }}
anvi-estimate-kegg-metabolism -M metagenomes.txt
{{ codestop }}

### Estimation for multiple bins in different metagenomes

If you have multiple bins (also known as %(internal-genomes)s across different collections or even different metagenomes, they can be analyzed with the same command by providing an internal genomes file to %(anvi-estimate-kegg-metabolism)s. To see the required format for the external genomes file, see %(internal-genomes)s.

{{ codestart }}
anvi-estimate-kegg-metabolism -i internal-genomes.txt
{{ codestop }}

## Adjusting module completion threshold

KEGG module completeness is computed as the percentage of steps in the metabolic pathway that are 'present' based on the KOs found in the contigs database. If this completeness is larger than a certain percentage, then the entire module is considered to be 'present' in the genome or metagenome. By default, this module completion threshold is 75%; that is, 75% of the KOs in a module must have a KOfam hit in the contigs database in order for the module to be considered 'complete' as a whole. This threshold can be adjusted.

### Changing the module completion threshold

In this example, we change the threshold to 50%.

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --module-completion-threshold 0.5
{{ codestop }}

## Controlling output

%(anvi-estimate-kegg-metabolism)s can produce a variety of output files. All will be prefixed with the same string, which by default is "kegg-metabolism".

### Changing the output file prefix

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db -O my-cool-prefix
{{ codestop }}


This program has two major output options - long format (tab-delimited) output files and matrices.

Long format output has several preset "modes" as well as a "custom" mode in which the user can define the contents of the output file. Multiple modes can be used at once, and each requested "mode" will result in a separate output file. The default output mode is "modules" mode.

### Viewing available output modes

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --list-available-modes
{{ codestop }}

### Using a non-default output mode

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --kegg-output-modes kofam_hits
{{ codestop }}

### Using multiple output modes

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --kegg-output-modes kofam_hits,modules
{{ codestop }}

### Viewing available output headers for 'custom' mode

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --list-available-output-headers
{{ codestop }}

### Using custom output mode

{{ codestart }}
anvi-estimate-kegg-metabolism -c CONTIGS.db --kegg-output-modes custom --custom-output-headers kegg_module,module_name,module_is_complete
{{ codestop }}


Matrix format is only available when working with multiple contigs databases. Several output matrices will be generated, each of which describes one KEGG module statistic such as completion score or presence/absence.  

### Obtaining matrix-formatted output

{{ codestart }}
anvi-estimate-kegg-metabolism -i internal-genomes.txt --matrix-format
{{ codestop }}

You can find more details on the output format by looking at %(kegg-metabolism)s.
