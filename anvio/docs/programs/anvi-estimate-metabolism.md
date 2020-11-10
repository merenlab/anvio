%(anvi-estimate-metabolism)s predicts the metabolic capabilities of organisms based on their genetic content. It relies upon %(kegg-functions)s and metabolism information from the KEGG resource, which is stored in a %(kegg-db)s.

The metabolic pathways that this program currently considers are those defined by KOs in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). Each KO represents a gene function, and a KEGG module is a set of KOs that collectively carry out the steps in a metabolic pathway. Therefore, for this to work, you need to have annotated your %(contigs-db)s with hits to the KEGG KOfam database by running %(anvi-run-kegg-kofams)s prior to using this program.

Given a properly annotated %(contigs-db)s, this program determines which KOs are present and from those determines the completeness of each KEGG module. The results are described in a set of output text files, collectively referred to as %(kegg-metabolism)s.

## Running metabolism estimation on a single contigs database

There are several possible inputs to this program. For single genomes (isolate genomes or MAGs, for example) you can provide a %(contigs-db)s. If your %(contigs-db)s describes a metagenome rather than a single genome, you can provide the flag `--metagenome-mode`. In metagenome mode, estimation is run on each contig individually - that is, only KOfam hits within the same contig are allowed to contribute to the completeness score of a given KEGG module. Alternatively, if you have binned your metagenome sequences into separate populations and would like metabolism estimation to be run separately on each bin, you can provide a %(profile-db)s and a %(collection)s.

This program always takes one or more contigs database(s) as input, but what is in those contigs dbs depends on the context (ie, genome, metagenome, bin). In the case of internal genomes (or bins), is possible to have multiple inputs but only one input contigs db. So for clarity's sake, we sometimes refer to the inputs as 'samples' in the descriptions below. If you are getting confused, just try to remember that a 'sample' can be a genome, a metagenome, or a bin.

### Estimation for a single genome

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db
{{ codestop }}

### Estimation for a metagenome

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --metagenome-mode
{{ codestop }}

{: .notice}
In metagenome mode, this program will estimate metabolism for each contig in the metagenome separately. This will tend to underestimate module completeness because it is likely that some modules will be broken up across multiple contigs belonging to the same population.

If you prefer to instead treat all KOs in the metagenome as belonging to one collective genome, you can do so by simply leaving out the `--metagenome-mode` flag (to effectively pretend that you are doing estimation for a single genome, although in your heart you will know that your contigs database really contains a metagenome). Please note that this will result in the opposite tendency to overestimate module completeness (as the KOs will in reality be coming from multiple different populations), and there will be a lot of redundancy.

We are working on improving our estimation algorithm for metagenome mode. In the meantime, if you are worried about the misleading results from either of these situations, we suggest binning your metagenomes first and running estimation for the bins as described below.

### Estimation for bins in a metagenome

You can estimate metabolism for each bin in a %(collection)s:

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME
{{ codestop }}

You can also provide a specific %(bin)s in that %(collection)s to run on:

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME -b BIN_NAME
{{ codestop }}

Or, you can provide a specific list of bins in a text file:

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME -B bin_ids.txt
{{ codestop }}

Each line in the `bin_ids.txt` file should be a bin name from the %(collection)s.

## MULTI-MODE: Running metabolism estimation on multiple contigs databases

If you have a set of contigs databases of the same type (i.e., all of them are single genomes or all are binned metagenomes), you can analyze them all at once. What you need to do is put the relevant information for each %(contigs-db)s into a text file and pass that text file to %(anvi-estimate-metabolism)s. The program will then run estimation individually on each contigs database in the file. The estimation results for each database will be aggregated and printed to the same output file(s).

### Estimation for multiple single genomes

Multiple single genomes (also known as %(external-genomes)s) can be analyzed with the same command by providing an external genomes file to %(anvi-estimate-metabolism)s. To see the required format for the external genomes file, see %(external-genomes)s.

{{ codestart }}
anvi-estimate-metabolism -e external-genomes.txt
{{ codestop }}

### Estimation for multiple metagenomes

Multiple metagenomes can be analyzed with the same command by providing a metagenomes input file. Metagenome mode will be used to analyze each contigs database in the file. To see the required format for the metagenomes file, see %(metagenomes)s.

{{ codestart }}
anvi-estimate-metabolism -M metagenomes.txt
{{ codestop }}

### Estimation for multiple bins in different metagenomes

If you have multiple bins (also known as %(internal-genomes)s) across different collections or even different metagenomes, they can be analyzed with the same command by providing an internal genomes file to %(anvi-estimate-metabolism)s. To see the required format for the internal genomes file, see %(internal-genomes)s.

{{ codestart }}
anvi-estimate-metabolism -i internal-genomes.txt
{{ codestop }}

## Adjusting module completion threshold

KEGG module completeness is computed as the percentage of steps in the metabolic pathway that are 'present' based on the KOs found in the contigs database. If this completeness is larger than a certain percentage, then the entire module is considered to be 'present' in the genome or metagenome. By default, this module completion threshold is 0.75; that is, 75 percent of the KOs in a module must have a KOfam hit in the contigs database in order for the module to be considered 'complete' as a whole. This threshold can be adjusted.

### Changing the module completion threshold

In this example, we change the threshold to 50 percent.

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --module-completion-threshold 0.5
{{ codestop }}

## Working with a non-default KEGG data directory
If you have previously annotated your contigs databases using a non-default KEGG data directory with `--kegg-data-dir` (see %(anvi-run-kegg-kofams)s), or have moved the KEGG data directory that you wish to use to a non-default location, then you will need to specify where to find the KEGG data so that this program can use the right one. In that case, this is how you do it:

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --kegg-data-dir /path/to/directory/KEGG
{{ codestop }}

## Controlling output

%(anvi-estimate-metabolism)s can produce a variety of output files. All will be prefixed with the same string, which by default is "kegg-metabolism".

### Changing the output file prefix

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db -O my-cool-prefix
{{ codestop }}

### Including only complete modules in the output
Remember that module completion threshold? Well, you can use that to control which modules make it into your output files. If you provide the `--only-complete` flag, then any module-related output files will only include modules that have a completeness score at or above the module completion threshold. (This doesn't affect KO-related outputs, for obvious reasons.)

Here is an example of using this flag with long format output (which is the default, as described below, but we are asking for it explicitly here just to be clear):
{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --kegg-output-modes modules --only-complete
{{ codestop }}

And here is an example of using this flag with matrix output. In this case, we are working with multiple input samples, and the behavior of this flag is slightly different: a module will be included in the matrix if it is at or above the module completion threshold in **at least one sample**. If there are any samples in which that module's completeness is below the threshold, its completeness in that sample will be **represented by a 0.0** in the matrix, regardless of its actual completeness score.
{{ codestart }}
anvi-estimate-metabolism -i internal-genomes.txt --matrix-format --only-complete
{{ codestop }}

## Output options
This program has two major output options: long format (tab-delimited) output files and matrices.

**Long Format Output**
Long format output has several preset "modes" as well as a "custom" mode in which the user can define the contents of the output file. Multiple modes can be used at once, and each requested "mode" will result in a separate output file. The default output mode is "modules" mode.

You can find more details on the output format by looking at %(kegg-metabolism)s.

### Viewing available output modes

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --list-available-modes
{{ codestop }}

### Using a non-default output mode

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --kegg-output-modes kofam_hits
{{ codestop }}

### Using multiple output modes

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --kegg-output-modes kofam_hits,modules
{{ codestop }}

### Viewing available output headers for 'custom' mode

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --list-available-output-headers
{{ codestop }}

### Using custom output mode

Here is an example of defining the modules output to contain columns with the module number, the module name, and the completeness score.

{{ codestart }}
anvi-estimate-metabolism -c CONTIGS.db --kegg-output-modes custom --custom-output-headers kegg_module,module_name,module_is_complete
{{ codestop }}

**Matrix Output**
Matrix format is only available when working with multiple contigs databases. Several output matrices will be generated, each of which describes one statistic such as module completion score, module presence/absence, or KO hit counts. Rows will describe modules or KOs, columns will describe your input samples (ie genomes, metagenomes, bins), and each cell of the matrix will be the corresponding statistic for a module in a sample. You can see examples of this output format by viewing %(kegg-metabolism)s.

### Obtaining matrix-formatted output

{{ codestart }}
anvi-estimate-metabolism -i internal-genomes.txt --matrix-format
{{ codestop }}

### Including KEGG metadata in the matrix output
By default, the matrix output is a matrix ready for use in other computational applications, like visualizing as a heatmap or performing clustering. But you may want to instead have a matrix that is annotated with more information, like the names and categories of each module or the functional annotations of each KO. To include this additional information in the matrix output (as columns that occur before the sample columns), use the `--include-metadata` flag.

{{ codestart }}
anvi-estimate-metabolism -i internal-genomes.txt --matrix-format --include-metadata
{{ codestop }}

Note that this flag only works for matrix output because, well, the long-format output inherently includes KEGG metadata.

Note also that you can combine this flag with the `--only-complete` flag, like so:
{{ codestart }}
anvi-estimate-metabolism -i internal-genomes.txt --matrix-format --only-complete --include-metadata
{{ codestop }}


## Testing this program
You can see if this program is working by running the following suite of tests, which will check several common use-cases:

{{ codestart }}
anvi-self-test --suite metabolism
{{ codestop }}
