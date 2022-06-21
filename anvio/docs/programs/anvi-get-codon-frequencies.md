This program **calculates codon or amino acid frequencies from genes or functions**.

A range of options allow calculation of different frequency statistics. This program is "maximalist," in that it has many options that do the equivalent of a couple extra commands in R or pandas -- because we (not you) tend to be lazy and prone to mistakes.

## Basic commands

### Gene frequencies

This command produces a table of codon frequencies from coding sequences in the contigs database, with the first column of the table before frequency data containing gene caller IDs. The decoded amino acid is included in each codon column name.

{{ codestart }}
anvi-get-codon-frequencies -c path/to/contigs.db -o path/to/output.txt --header-amino-acids
{{ codestop }}

### Function frequencies

This command produces a table of function frequencies rather than gene frequencies. By using `--function-sources` without any arguments, the output will include every function from every source used to annotate the genome, e.g., 'KOfam', 'KEGG_BRITE', 'Pfam'. The first four columns of the table before frequency data contain, respectively, gene caller IDs, function sources, accessions, and names.

{{ codestart }}
anvi-get-codon-frequencies -c path/to/contigs.db --function-sources --function-table-output path/to/function_output.txt
{{ codestop }}

### Gene frequencies with function information

In contrast to the previous example, this command produces a table of gene frequencies, but has an entry for every gene/function pair, allowing statistical interrogation of the gene components of functions. The function table output is derived from this table by grouping rows by function source, retaining only one row per gene caller ID, and summing frequencies across rows of the groups.

{{ codestart }}
anvi-get-codon-frequencies -c path/to/contigs.db --function-sources --gene-table-output path/to/gene_output.txt
{{ codestop }}

### Codon frequencies from multiple internal and external genomes

This command produces a table of codon frequencies from coding sequences in multiple genomes. A column is added at the beginning of the table for genome name.

{{ codestart }}
anvi-get-codon-frequencies -i path/to/internal_genomes.txt -e path/to/external_genomes.txt -o path/to/output.txt
{{ codestop }}

## Option examples

Tables show the options to get the requested statistics.

Options can be mixed and matched, with there being relatively few incompatible options.

### _Different_ frequency statistics

| Get | Options |
| --- | ------- |
| Codon frequencies | |
| Codon relative frequencies | `--relative` |
| [Synonymous (per-amino acid) codon relative frequencies](#synonymous-codon-relative-frequencies) | `--synonymous` |
| Amino acid frequencies | `--amino-acid` |
| Amino acid relative frequencies | `--amino-acid --relative` |
| [Summed frequencies across genes](#frequencies-across-genes) | `--sum` |
| [Synonymous relative summed frequencies across genes](#frequencies-across-genes) | `--sum --synonymous` |
| [Summed frequencies across genes annotated by each function source](#frequencies-across-genes) | `--sum --function-sources` |
| [Relative summed frequencies across genes with KOfam annotations](#frequencies-across-genes) | `--sum --relative --function-sources KOfam` |
| [Average frequencies across all genes](#frequencies-across-genes) | `--average` |

### Frequencies from _sets of genes with shared functions_

| Get | Options |
| --- | ------- |
| All function annotation sources | `--function-sources` |
| [All KEGG BRITE categories](#brite-hierarchies) | `--function-sources KEGG_BRITE` |
| All KEGG KOfams and all Pfams | `--function-sources KOfam Pfam` |
| [Certain KEGG BRITE categories](#brite-hierarchies) | `--function-sources KEGG_BRITE --function-names Ribosome Ribosome>>>Ribosomal proteins` |
| Certain KEGG KOfam accessions | `--function-sources KOfam --function-accessions K00001 K00002` |
| [Certain BRITE categories and KOfam accessions](#inputs) | `--functions-txt path/to/functions_table.txt` |

### Frequencies from _selections of genes_

| Get | Options |
| --- | ------- |
| From contigs database | `--contigs-db path/to/contigs.db` |
| From collection of internal genomes | `--contigs-db path/to/contigs.db --profile-db path/to/profile.db --collection-name my_bins` |
| From certain listed bins | `--contigs-db path/to/contigs.db --profile-db path/to/profile.db --bin-names-list my_first_bin,my_second_bin` |
| From certain bins listed in a file | `--contigs-db path/to/contigs.db --profile-db path/to/profile.db --bin-ids-file path/to/bin_list.txt` |
| From external genomes (contigs databases) listed in a file | `--external-genomes path/to/genomes.txt` |
| With certain gene IDs | `--gene-caller-ids 1,2,500` |
| With certain gene IDs or genes annotated with certain KOfams | `--gene-caller-ids 1,2,500 --function-sources KOfam --function-accessions K00001` |

### _Filtering genes and codons_ that are analyzed and reported

| Get | Options |
| --- | ------- |
| [Exclude genes shorter than 300 codons](#gene-length-and-codon-count) | `--gene-min-codons 300` |
| [Exclude functions with <300 codons](#function-codon-count) | `--function-min-codons 300` |
| [Exclude genes shorter than 300 codons from contributing to function codon frequencies](#function-codon-count) | `--gene-min-codons 300 --function-sources` |
| [Exclude stop codons and single-codon amino acids](#codons) | `--exclude-amino-acids STP Met Trp` |
| [Only include certain codons](#codons) | `--include-amino-acids Leu Ile` |
| [Exclude codons for amino acids with <5 codons in >90% of genes](#codons) | `--exclude-amino-acid-count 5 --exclude-amino-acid-fraction 0.9` |
| [Exclude genes with <300 codons after dynamically filtering codons based on all genes of ANY length](#gene-length-and-codon-count) | `--gene-min-codons 300 --exclude-amino-acid-count 5 --exclude-amino-acid-fraction 0.9 --min-codon-filter remaining` |
| [Exclude genes with lengths <300 codons before but not after dynamically filtering codons](#gene-length-and-codon-count) dynamically filtering codons | `--gene-min-codons 300 --exclude-amino-acid-count 5 --exclude-amino-acid-fraction 0.9 --min-codon-filter length` |

## Option details

### Synonymous codon relative frequencies

This flag returns the relative frequency of each codon among the codons encoding the same amino acid, e.g., 0.4 GCC and 0.6 GCT for Ala. By default, stop codons and single-codon amino acids (Met ATG and Trp TGG) in the standard translation table are excluded, equivalent to using `--exclude-amino-acids STP Met Trp` for other frequency statistics.

### Frequencies across genes

`--sum` and `--average` produce a table with a single row of frequencies from across genes. For example, the following command sums the codon frequencies of each decoded amino acid (and STP) across all genes, and then calculates the relative frequencies of the amino acids.

{{ codestart }}
anvi-get-codon-frequencies -c path/to/contigs.db -o path/to/output_table.txt --sum --amino-acid --relative
{{ codestop }}

These options operate on genes; when used with a function option, `--average` subsets the genes annotated by the functions of interest and calculates the average frequency across genes rather than functions (sums of genes with functional annotation). For example, the following command calculates the average synonymous relative frequency across genes annotated by KEGG KOfams.

{{ codestart }}
anvi-get-codon-frequencies -c path/to/contigs.db -o path/to/output_table.txt --average --synonymous --function-sources KOfam
{{ codestop }}

### Functions

Functions and function annotation sources can be provided to subset genes (as seen in the [last section](#frequencies-across-genes) with `--average`) and to calculate statistics for functions in addition to genes (as seen in a [previous example](#function-codon-frequencies).

Using `--output-file` is equivalent to `--gene-table-output` rather than `--function-table-output`, producing rows containing frequencies for annotated genes rather than summed frequencies for functions.

#### Inputs

There are multiple options to define which functions and sources should be used. `--function-sources` without arguments uses all available sources that had been used to annotate genes. `--function-accessions` and `--function-names` select functions from a single provided source. To use different functions from different sources, a tab-delimited file can be provided to `functions-txt`. This headerless file must have three columns, for source, accession, and name of functions, respectively, with an entry in each row for `source`.

#### BRITE hierarchies

Genes are classified in KEGG BRITE functional hierarchies by %(anvi-run-kegg-kofams). For example, a bacterial SSU ribosomal protein is classified in a hierarchy of ribosomal genes, `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`. Codon frequencies can be calculated for genes classified at each level of the hierarchy from the most general, those genes in the `Ribosome`, to the most specific -- in the example, those genes in `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`. Therefore, the following command returns summed codon frequencies for each annotated hierarchy level -- in the example, the output would include four rows for the genes in each level from `Ribosome` to `Small subunit`.

{{ codestart }}
anvi-get-codon-frequencies -c path/to/contigs.db -o path/to/output_table.txt --function-sources KEGG_BRITE
{{ codestop }}

### Filter genes and codons

#### Codons

It may be useful to restrict codons in the analysis to those encoding certain amino acids. Stop codons and the single codons encoding Met and Trp are excluded by default from calculation of synonymous codon relative frequencies (`--synonymous`). Relative frequencies across codons in a gene (`--relative`) are calculated for the selected amino acids, so the following option would return a table of codon frequencies relative to the codons encoding the selected nonpolar amino acids: `--include-amino-acids Gly Ala Val Leu Met Ile`.

Dynamic exclusion of amino acids can be useful in the calculation of synonymous codon frequencies, e.g., 0.5 GCC and 0.5 GCT for Ala may be statistically insignificant for a gene with 1 GCC and 1 GCT; even more meaningless would be 1.0 GCC and 0.0 GCT for a gene with 1 GCC and 0 GCT. Two options used in tandem remove rarer amino acids across the dataset, `--exclude-amino-acid-count` and `--exclude-amino-acid-fraction`. These set a minimum number of codons in a minimum number of genes to retain the amino acid. For example, amino acids with <5 codons in >90% of genes will be excluded from the analysis with the arguments, `--exclude-amino-acid-count 5 --exclude-amino-acid-fraction 0.9`.

#### Gene length and codon count

Removal of genes with few codons can improve the statistical utility of relative frequencies. `--gene-min-codons` sets the minimum number of codons required in a gene, and this filter can be applied before and/or after the removal of rarer codons. Applied before, `--gene-min-codons` filters genes by length; applied after, it filters genes by codons remaining after removing rarer codons. `--min-codon-filter` can take three possible arguments: `length`, `remaining`, or, by default when codons are removed, `both`, which applies the `--gene-min-codons` filter both before and after codon removal.

It may seem redundant for `remaining` and `both` to both be possibilities, but this is due to the possibility of dynamic amino acid exclusion using `--exclude-amino-acid-count` and `--exclude-amino-acid-fraction`. Amino acids are removed based on their frequency in a proportion of genes, so removing shorter genes by length before removing amino acids can affect which amino acids are dynamically excluded.

#### Function codon count

`--function-min-codons` can be used to filter functions with a minimum number of codons. Function codon count filters occur after gene codon count filters: the set of genes contributing to function codon frequency can be restricted by applying `--gene-min-codons`.
