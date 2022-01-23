Output text files produced by %(anvi-estimate-metabolism)s that describe the presence of metabolic pathways in a %(contigs-db)s.

Depending on the output options used when running %(anvi-estimate-metabolism)s, these files will have different formats. This page describes and provides examples of the various output file types.

### How to get to this output
![A beautiful workflow of metabolism reconstruction in anvi'o](../../images/metabolism_reconstruction.png)

## Long-format output modes

The long-format output option produces tab-delimited files. Different output "modes" will result in output files with different information.

### 'Modules' Mode

The 'modules' mode output file will have the suffix `modules.txt`. Each line in the file will represent information about a metabolic module in a given genome, bin, or contig of a metagenome assembly. Here is one example, produced by running metabolism estimation on the _Enterococcus_ external genomes in the [Infant Gut dataset](http://merenlab.org/tutorials/infant-gut/):

module | genome_name | db_name | module_name | module_class | module_category | module_subcategory | module_definition | module_completeness | module_is_complete | proportion_unique_enzymes_present | enzymes_unique_to_module | unique_enzymes_hit_counts | enzyme_hits_in_module | gene_caller_ids_in_module | warnings
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
M00001 | Enterococcus_faecalis_6240 | E_faecalis_6240 | Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "(K00844,K12407,K00845,K25026,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)" | 1.0 | True | NA | No enzymes unique to module | NA | K00134,K00134,K00850,K00873,K00927,K01624,K01689,K01803,K01803,K01810,K01834,K01834,K25026 | 1044,642,225,226,1043,348,1041,1042,1043,600,2342,2646,1608 | K01624 is present in multiple modules: M00001/M00003/M00165/M00167/M00345/M00344/M00611/M00612,K01810 is present in multiple modules: M00001/M00004/M00892/M00909,K00927 is present in multiple modules: M00001/M00002/M00003/M00308/M00552/M00165/M00166/M00611/M00612,K00873 is present in multiple modules: M00001/M00002,K25026 is present in multiple modules: M00001/M00549/M00909,K01689 is present in multiple modules: M00001/M00002/M00003/M00346,K01834 is present in multiple modules: M00001/M00002/M00003,K01803 is present in multiple modules: M00001/M00002/M00003,K00134 is present in multiple modules: M00001/M00002/M00003/M00308/M00552/M00165/M00166/M00611/M00612,K00850 is present in multiple modules: M00001/M00345
M00002 | Enterococcus_faecalis_6240 | E_faecalis_6240 | Glycolysis, core module involving three-carbon compounds | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)" | 1.0 | True | NA | No enzymes unique to module | NA | K00134,K00134,K00873,K00927,K01689,K01803,K01803,K01834,K018341044,642,226,1043,1041,1042,1043,2342,2646 | K00927 is present in multiple modules: M00001/M00002/M00003/M00308/M00552/M00165/M00166/M00611/M00612,K00873 is present in multiple modules: M00001/M00002,K01689 is present in multiple modules: M00001/M00002/M00003/M00346,K01834 is present in multiple modules: M00001/M00002/M00003,K01803 is present in multiple modules: M00001/M00002/M00003,K00134 is present in multiple modules: M00001/M00002/M00003/M00308/M00552/M00165/M00166/M00611/M00612
M00003 | Enterococcus_faecalis_6240 | E_faecalis_6240 | Gluconeogenesis, oxaloacetate => fructose-6P | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "(K01596,K01610) K01689 (K01834,K15633,K15634,K15635) K00927 (K00134,K00150) K01803 ((K01623,K01624,K11645) (K03841,K02446,K11532,K01086,K04041),K01622)" | 0.875 | True | NA | No enzymes unique to module | NA | K00134,K00134,K00927,K01624,K01689,K01803,K01803,K01834,K01834,K04041 | 1044,642,1043,348,1041,1042,1043,2342,2646,617 | K01624 is present in multiple modules: M00001/M00003/M00165/M00167/M00345/M00344/M00611/M00612,K00927 is present in multiple modules: M00001/M00002/M00003/M00308/M00552/M00165/M00166/M00611/M00612,K04041 is present in multiple modules: M00003/M00611/M00612,K01689 is present in multiple modules: M00001/M00002/M00003/M00346,K01834 is present in multiple modules: M00001/M00002/M00003,K01803 is present in multiple modules: M00001/M00002/M00003,K00134 is present in multiple modules: M00001/M00002/M00003/M00308/M00552/M00165/M00166/M00611/M00612
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

What are the data in each of these columns?

- `module`: the module identifier for a metabolic pathway (from the KEGG MODULE database or from user-defined modules)
- `genome_name`/`bin_name`/`contig_name`: the identifier for the current sample, whether that is a genome, bin, or contig from a metagenome assembly
- `db_name`: the name of the contigs database from which this data comes (only appears in output from multi-mode, in which multiple DBs are processed at once)
- `module_name`/`module_class`/`module_category`/`module_subcategory`/`module_definition`: metabolic pathway information, from the KEGG MODULE database or from user-defined metabolic modules
- `module_completeness`: a fraction between 0 and 1 indicating the proportion of steps in the metabolic pathway that have an associated enzyme annotation. To learn how this number is calculated, see [the anvi-estimate-metabolism help page](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/#how-is-the-module-completeness-score-calculated)
- `module_is_complete`: a boolean value indicating whether the `module_completeness` score is above a certain threshold or not (the default threshold is 0.75)
- `proportion_unique_enzymes_present`: some enzymes only belong to one metabolic pathway, which means that their presence is a better indicator for the presence of a metabolic capacity than other enzymes that are shared across multiple pathways. This column contains the fraction of these unique enzymes that are present in the sample. For instance, if the module has only 1 unique enzyme and it is present, you will see a 1 in this column. You can find out the denominator of this fraction (ie, the number of unique enzymes for this module) by either calculating the length of the list in the `enzymes_unique_to_module` column, or by requesting custom modules mode output with the `unique_enzymes_context_string` header
- `enzymes_unique_to_module`: a comma-separated list of the enzymes that only belong to the current module (ie, are not shared across multiple metabolic pathways)
- `unique_enzymes_hit_counts`: a comma-separated list of how many times each unique enzyme appears in the sample, in the same order as the `enzymes_unique_to_module` list
- `enzyme_hits_in_module`: a comma-separated list of the enzyme annotations that were found in the current sample and contribute to this metabolic pathway (these will be enzymes from the metabolic pathway definition in the `module_definition` column)
- `gene_caller_ids_in_module`: a comma-separated list of the genes with enzyme annotations that contribute to this pathway, in the same order as the annotations in the `enzyme_hits_in_module` column
- `warnings`: miscellaneous caveats to consider when interpreting the `module_completeness` score. For example, a warning like "No KOfam profile for K00172" would indicate that we cannot annotate K00172 because we have no HMM profile for that gene family, which means that any metabolic pathway containing this KO can never be fully complete (even if a gene from that family does exist in your sequences). Seeing many warnings like "K01810 is present in multiple modules: M00001/M00004/M00892/M00909" indicates that the current module shares many enzymes with other metabolic pathways, meaning that it may appear to be complete only because its component enzymes are common. Extra caution should be taken when considering the completeness of modules with warnings

**Coverage and detection values in the output**

If you use the flag `--add-coverage` and provide a profile database, additional columns containing coverage and detection data will be added for each sample in the profile database. Here is a mock example of the additional columns you will see (for a generic sample called 'SAMPLE_1'):

| SAMPLE_1_gene_coverages | SAMPLE_1_avg_coverage | SAMPLE_1_gene_detection | SAMPLE_1_avg_detection |
|:--|:--|:--|:--|
| 3.0,5.0,10.0,2.0 | 5.0 | 1.0,1.0,1.0,1.0 | 1.0 |

In this mock example, the module in this row has four gene calls in it. The `SAMPLE_1_gene_coverages` column lists the mean coverage of each of those genes in SAMPLE_1 (in the same order as the gene calls are listed in the `gene_caller_ids_in_module` column), and the `SAMPLE_1_avg_coverage` column holds the average of these values. As you probably expected, the `detection` columns are similarly defined, except that they contain detection values instead of coverage.

### 'Hits in Modules' Mode

The 'hits_in_modules' output file will have the suffix `hits_in_modules.txt`. Each line in the file will represent information about one enzyme annotation (ie, one gene) in a given genome, metagenome, or bin - but _only_ enzymes that are a part of at least one module are included in this output. Hits are organized according to the module that they belong to, and more specifically the path through the module in which the enzyme appears.

What is a path through a module, you ask? Well. There is a lengthier explanation of this [here](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/#what-data-is-used-for-estimation), but we will go through it briefly below.

KEGG modules are metabolic pathways defined by a set of KOs. For example, here is the definition of module [M00001](https://www.genome.jp/kegg-bin/show_module?M00001), better known as "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate":

(K00844,K12407,K00845,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)

Spaces separate steps (reactions) in the metabolic pathway, and commas separate alternative KOs or alternative sub-pathways that can facilitate the same overall reaction. So a definition such as the one above can be "unrolled" into several different linear sequences of KOs, each of which we consider to be a possible "path" through the module. As an example, we can take the first option for every step in the Embden-Meyerhof pathway definition from above:

(**K00844**,K12407,K00845,K00886,K08074,K00918) (**K01810**,K06859,K13810,K15916) (**K00850**,K16370,K21071,K00918) (**K01623**,K01624,K11645,K16305,K16306) **K01803** ((**K00134**,K00150) **K00927**,K11389) (**K01834**,K15633,K15634,K15635) **K01689** (**K00873**,K12406)

to get the following path of KOs (which happens to be the path shown in the output example below):

K00844 K01810 K00850 K01623 K01803 K00134 K00927 K01834 K01689 K00873

For every KO in the path above that has a hit in the %(contigs-db)s, there will be a corresponding line in the 'hits_in_modules' output file. The same will occur for every possible path in every single KEGG module, resulting in a lot of lines and extremely repetitive but nicely parseable information.

The same principle applies to user-defined metabolic modules, except that the enzymes can be from a variety of different annotation sources (not just KOfam).

Without further ado, here is an example of this output mode (also from the Infant Gut dataset):

module | genome_name | db_name | module_is_complete | module_completeness | path_id | path | path_completeness | enzyme_hit | gene_caller_id | contig
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
M00001 | Enterococcus_faecalis_6240 | E_faecalis_6240 | True | 1.0 | 0 | K00844,K01810,K00850,K01623,K01803,K00134,K00927,K01834,K01689,K00873 | 0.8 | K00134 | 1044 | Enterococcus_faecalis_6240_contig_00003_chromosome
M00001 | Enterococcus_faecalis_6240 | E_faecalis_6240 | True | 1.0 | 0 | K00844,K01810,K00850,K01623,K01803,K00134,K00927,K01834,K01689,K00873 | 0.8 | K25026 | 642 | Enterococcus_faecalis_6240_contig_00003_chromosome
M00001 | Enterococcus_faecalis_6240 | E_faecalis_6240 | True | 1.0 | 0 | K00844,K01810,K00850,K01623,K01803,K00134,K00927,K01834,K01689,K00873 | 0.8 | K01689 | 1041 | Enterococcus_faecalis_6240_contig_00003_chromosome
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Many of the columns in this data overlap with the 'modules' mode columns; you can find descriptions of those in the previous section. Below are the descriptions of new columns in this mode:

- `enzyme_hit`: an enzyme annotation from the contigs database that contributes to a module
- `gene_caller_id`: the ID of the gene that is annotated with this enzyme
- `contig`: the contig on which this gene is present
- `path_id`: a unique identifier of the current path through the module
- `path`: the current path of enzymes through the module (described above), which this enzyme annotation contributes to
- `path_completeness`: a fraction between 0 and 1 indicating the proportion of enzymes in the _current path_ that are annotated. To learn how this number is calculated, see [the anvi-estimate-metabolism help page](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/#how-is-the-module-completeness-score-calculated)

**Coverage and detection values in the output**

The flag `--add-coverage` can also be used for this output mode, but in this case, the columns that are added (for each sample in the profile database) are slightly different. Here is a mock example:

| SAMPLE_1_coverage | SAMPLE_1_detection |
|:--|:--|
| 3.0 | 1.0 |

Since each row is a single gene in this output mode, these columns will contain the coverage/detection values for that gene only.

### Enzyme 'Hits' Mode

The 'hits' output file will have the suffix `hits.txt`. Unlike the previous mode, this output will include ALL enzyme hits (from all annotation sources used for metabolism estimation), regardless of whether the enzyme belongs to a metabolic module or not. Since only a subset of these enzymes belong to modules, this output does not include module-related information like paths and module completeness.

Here is an example of this output mode (also from the Infant Gut dataset):

enzyme | genome_name | db_name | gene_caller_id | contig | modules_with_enzyme | enzyme_definition
|:--|:--|:--|:--|:--|:--|:--|
K25026 | Enterococcus_faecalis_6240 | E_faecalis_6240 | 1608 | Enterococcus_faecalis_6240_contig_00003_chromosome | M00001,M00549,M00909 | glucokinase [EC:2.7.1.2]
K01810 | Enterococcus_faecalis_6240 | E_faecalis_6240 | 600 | Enterococcus_faecalis_6240_contig_00003_chromosome | M00001,M00004,M00892,M00909 | glucose-6-phosphate isomerase [EC:5.3.1.9]
K00850 | Enterococcus_faecalis_6240 | E_faecalis_6240 | 225 | Enterococcus_faecalis_6240_contig_00003_chromosome | M00001,M00345 | 6-phosphofructokinase 1 [EC:2.7.1.11]
(...) |(...)|(...)|(...)|(...)|(...)|(...)|

Here are the descriptions of any new columns not yet discussed in the previous sections:

- `enzyme`: an enzyme that was annotated in the contigs database
- `modules_with_enzyme`: the modules (if any) that this enzyme belongs to
- `enzyme_definition`: the function of this enzyme (often includes the enzyme name and EC number)

**Coverage and detection values in the output**

If you use the flag `--add-coverage` and provide a profile database, you will get the same additional columns per row as described above in for `hits_in_modules` mode. That is, you will get one column per sample for coverage (containing the coverage value of the KO hit in the sample) and one column per sample for detection (containing the detection value of the KO hit in the sample).

### Custom Mode (for module data)

The 'modules_custom' output mode will have user-defined content and the suffix `modules_custom.txt` (we currently only support output customization for modules data). See %(anvi-estimate-metabolism)s for an example command to work with this mode. The output file will look similar to the 'modules' mode output, but with a different (sub)set of columns.

## Matrix format output

Matrix format is an output option when %(anvi-estimate-metabolism)s is working with multiple contigs databases at once. The purpose of this output type is to generate matrices of module statistics for easy visualization and clustering. Currently, the matrix-formatted output includes a module completeness matrix, a matrix of binary module presence/absence values, and a matrix of enzyme annotation counts. In these matrices, each row is a module or enzyme, and each column is an input sample.

Here is an example of a module completeness matrix, for bins in a metagenome:

| module | bin_1 | bin_2 | bin_3 | bin_4 | bin_5 | bin_6 |
|:--|:--|:--|:--|:--|:--|:--|
| M00001 | 1.00 | 0.00 | 1.00 | 1.00 | 1.00 | 0.00 |
| M00002 | 1.00 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| M00003 | 0.88 | 0.00 | 1.00 | 0.75 | 1.00 | 0.88 |
| M00004 | 0.88 | 0.00 | 0.88 | 0.88 | 0.88 | 0.00 |
| M00005 | 1.00 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 |
|(...) | (...) | (...) | (...) | (...) | (...) | (...) |

Each cell of the matrix is the completeness score for the corresponding module in the corresponding sample (which is, in this case, a bin).

While the above is the default matrix format, some users may want to include more annotation information in the matrices so that it is easier to know what is going on when looking at the matrix data manually. You can add this metadata to the matrices by using the `--include-metadata` flag when running %(anvi-estimate-metabolism)s, and the output will look something like the following:

| module | module_name | module_class | module_category | module_subcategory | bin_1 | bin_2 | bin_3 | bin_4 | bin_5 | bin_6 |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| M00001 |Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | 1.00 | 0.00 | 1.00 | 1.00 | 1.00 | 0.00 |
| M00002 | Glycolysis, core module involving three-carbon compounds | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | 1.00 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| M00003 | Gluconeogenesis, oxaloacetate => fructose-6P | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | 0.88 | 0.00 | 1.00 | 0.75 | 1.00 | 0.88 |
|(...) | (...) | (...) | (...) | (...) | (...) | (...) |

The module completeness matrix files will have the suffix `completeness-MATRIX.txt`.

Module presence/absence matrix files will have the suffix `presence-MATRIX.txt`. In these files, each cell of the matrix will have either a 1.0 or a 0.0. A 1.0 indicates that the module has a completeness score above the module completeness threshold in that sample, while a 0.0 indicates that the module's completeness score is not above the threshold.

Finally, KO hit matrix files will have the suffix `enzyme_hits-MATRIX.txt`. Each row of the matrix will be an enzyme, and each column will be an input sample. Cells in this matrix will contain an integer value, representing the number of times the enzyme was annotated in that sample. (Note: you will also add metadata to this matrix type when you use the `--include-metadata` flag).
