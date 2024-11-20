This program **calculates codon or amino acid frequencies from genes or functions**.

A range of options enables calculation of different frequency statistics.

## Basic commands

This section shows types of codon frequency data that can be generated.

### Gene frequencies

The following command produces a table of codon frequencies from coding sequences in an input contigs database. The first column of the tab-delimited output file contains anvi'o gene caller IDs. Other columns contain codon frequencies. By default, these columns are headed by codons, e.g., "ACG". The `--header-amino-acids` flag includes the decoded amino acid with the codon in the header, e.g, "ThrACG".

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --header-amino-acids
{{ codestop }}

### Function frequencies

The following command produces a table of codon frequencies for functions rather than genes. By using `--function-sources` without any arguments, the output will include every %(functions)s source available in a given %(contigs-db)s, e.g., `KOfam`, `KEGG_BRITE`, `Pfam`. The complete list of function sources in a contigs database can be checked with %(anvi-db-info)s. The first three columns of the output table contain function sources, accessions, and names, e.g., "COG20_FUNCTION", "COG0001", "Glutamate-1-semialdehyde aminotransferase (HemL) (PDB:2CFB)". Other tables in the column contain codon frequencies.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           --function-sources \
                           --function-table-output function_output.txt
{{ codestop }}

More information on the meaning of function codon frequencies can be found in the [functions section](#inputs).

### Gene frequencies with function information

In contrast to the [previous example](#function-frequencies), the following command produces a table of gene codon frequencies, but has an entry for every gene/function pair. This enables analysis of the genes encoding functions, including genes encoding the same function. The function table output of the previous example can be directly derived from this output by summing gene frequencies with the same function source. The `--gene-table-output` option is interchangeable with `-o`/`--output-file`, if you wish to make it more explicit that the command is producing a table of gene rather than function codon frequencies.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           --function-sources \
                           -o gene_output.txt
{{ codestop }}

### Codon frequencies from multiple internal and external genomes

The following command produces a table of gene codon frequencies from coding sequences in multiple genomes. Here, separate files are provided listing both %(internal-genomes)s and %(external-genomes)s. The first column of the output table contains genome names. Other columns are codon frequencies.

{{ codestart }}
anvi-get-codon-frequencies -i %(internal-genomes)s \
                           -e %(external-genomes)s \
                           -o output.txt
{{ codestop }}

## Option examples

The following tables demonstrate how to use options to get various results. These examples are not comprehensive.

### Normalized frequencies

| Get | Options |
| --- | ------- |
| Codon absolute frequencies | |
| Codon relative frequencies | `--relative` |
| [Synonymous (per-amino acid) codon relative frequencies](#synonymous-codon-relative-frequencies) | `--synonymous` |

### Amino acid frequencies

| Amino acid frequencies | `--amino-acid` |
| Amino acid relative frequencies | `--amino-acid --relative` |

### Aggregate frequencies

| Get | Options |
| --- | ------- |
| [Summed frequencies across genes](#aggregate-frequencies-across-genes) | `--sum` |
| [Average frequencies across genes](#aggregate-frequencies-across-genes) | `--average` |
| [Relative summed frequencies across genes](#aggregate-frequencies-across-genes) | `--sum --relative` |
| [Synonymous relative summed frequencies across genes](#aggregate-frequencies-across-genes) | `--sum --synonymous` |

### Report gene function annotations

| Get | Options |
| --- | ------- |
| [Gene frequencies showing all function annotation sources](#gene-frequencies-with-function-information) | `--function-sources` |
| [Gene frequencies showing annotations from select sources, KOs and COG20](#inputs) | `--function-sources KOfam COG20_FUNCTION` |
| [Genes in certain KEGG BRITE categories](#brite-hierarchies) | `--function-sources KEGG_BRITE --function-names Ribosome Ribosome>>>Ribosomal proteins` |
| [Genes with certain KO accessions](#inputs) | `--function-sources KOfam --function-accessions K00001 K00002` |
| [Genes with certain annotations, e.g., in certain BRITE categories and with certain COG20 accessions](#inputs) | `--select-functions-txt select_functions.txt` |

### Function frequencies

| Get | Options |
| --- | ------- |
| [Function frequencies from all annotation sources](#function-frequencies) | `--function-sources --function-table-output output.txt` |
| [Function frequencies from a single annotation source, Pfam](#inputs) | `--function-sources Pfam --function-table-output output.txt` |
| [Frequencies for KEGG BRITE categories of KOs](#brite-hierarchies) | `--function-sources KEGG_BRITE --function-table-output output.txt` |
| [Summed frequencies across genes with each annotation source](#aggregate-frequencies-across-genes-with-function) | `--sum --function-sources` |
| [Relative summed frequencies across genes with KO annotations](#aggregate-frequencies-across-genes-with-function) | `--sum --relative --function-sources KOfam` |

### Select genes

| Get | Options |
| --- | ------- |
| From contigs database | `--contigs-db contigs.db` |
| From collection of internal genomes | `--contigs-db contigs.db --profile-db profile.db --collection-name my_bins` |
| From internal genome | `--contigs-db contigs.db --profile-db profile.db --collection-name my_bins --bin-id my_bin` |
| From internal genomes listed in a file | `--internal-genomes genomes.txt` |
| From external genomes (contigs databases) listed in a file | `--external-genomes genomes.txt` |
| With certain gene caller IDs | `--gene-caller-ids 0 2 500` |
| With certain gene caller IDs or genes annotated by certain KOs | `--gene-caller-ids 0 2 500 --function-sources KOfam --function-accessions K00001 K00002` |

### Filter analyzed codons, genes, and functions

| Get | Options |
| --- | ------- |
| [Exclude codons ending with A](#select-codons-and-amino-acids-in-analysis) | `--exclude-codons ..A` |
| [Only include codons starting with G or T](#select-codons-and-amino-acids-in-analysis) | `--include-codons [GT]..` |
| [Exclude stop codons and single-codon amino acids](#select-codons-and-amino-acids-in-analysis) | `--exclude-amino-acids STP Met Trp` |
| [Only include codons encoding certain amino acids](#select-codons-and-amino-acids-in-analysis) | `--include-amino-acids Leu Ile` |
| [Exclude codons with counts <3 in ≥90%% of genes](#dynamically-exclude-rarer-codons-and-amino-acids-from-analysis) | `--pansequence-min-codons 3 0.9` |
| [Exclude codons for amino acids with <5 codons in ≥90%% of genes](#dynamically-exclude-rarer-codons-and-amino-acids-from-analysis) | `--pansequence-min-amino-acids 5 0.9` |
| [Exclude genes shorter than 250 codons](#filter-genes-by-codon-count) | `--gene-min-codons 250` |
| [Exclude functions with <250 codons](#filter-functions-by-codon-count) | `--function-min-codons 250` |
| [Exclude genes shorter than 250 codons from contributing to function codon frequencies](#filter-functions-by-codon-count) | `--gene-min-codons 250 --function-sources` |

### Filter reported codons and amino acids

| Get | Options |
| --- | ------- |
| [Exclude codons starting with A and ending C or T from output but not analysis](#select-codons-and-amino-acids-in-output) | `--dont-report-codons A.[CT]` |
| [Only include codons with C or G at all positions in output](#select-codons-and-amino-acids-in-output) | `--report-codons [CG][CG][CG]` |
| [Don't report stop codons and single-codon amino acids in output](#select-codons-and-amino-acids-in-output) | `--dont-report-amino-acids STP Met Trp` |
| [Only report codons encoding certain amino acids in output](#select-codons-and-amino-acids-in-output) | `--report-amino-acids Leu Ile` |
| [Replace codons with counts <3 in the gene or function with NaN](#mask-rarer-codons-and-amino-acids-in-output) | `--sequence-min-codons 3` |
| [Replace codons for amino acids with <5 codons in the gene or function with NaN](#mask-rarer-codons-and-amino-acids-in-output) | `--sequence-min-amino-acids 5` |

## Option details

### Normalized frequency

#### Relative frequency

The `--relative` flag normalizes codon frequencies. Used alone, `--relative` divides codon frequencies by the total for the coding sequence row in the output table.

#### Synonymous frequency

The `--synonymous` flag returns the relative frequency of each codon encoding the same amino acid, e.g., 0.4 GAC and 0.6 GAT for the two Asp codons. This is often called relative synonymous codon usage (RSCU). The `--relative` flag is redundant if `--synonymous` is also used. By default, stop codons and single-codon amino acids (Met ATG and Trp TGG) in the standard translation table are excluded from analysis and output of RSCU, equivalent to using `--exclude-amino-acids STP Met Trp` for other frequency statistics.

### Amino acid frequency

The `--return-amino-acids` flag returns codon frequencies for amino acids, summing frequencies of synonymous codons. Each column of the table represents amino acid rather than codon data, and is headed by amino acid abbreviations, plus "STP" for stop codons, rather than codons. For example, Asp frequencies are the summed frequencies of its two codons, GAC and GAT, and Met represents its single codon, ATG.

### Aggregate gene frequency

#### Summed frequency

The `--sum` flag sums codon frequencies across genes. Used alone, it produces a table with a single row of codon frequencies summed across all coding sequences in the input, which, for example, could output codon counts for a whole genome. This flag interacts in different ways with function options.

#### Average frequency

The `--average` flag averages codon frequencies across genes. Used alone, it produces a table with a single row of average codon frequencies from all coding sequences in the input. This flag interacts in different ways with function options.

#### Standard deviation in frequency

The `--std` flag take the standard deviation of frequencies across genes. Used alone, it produces a table with a single row of standard deviations of codon frequencies from all coding sequences in the input.

### Aggregate gene frequencies

The `--sum` and `--average` flags aggregate gene frequencies.

`--sum` and `--average` produce a table with a single row of frequencies aggregated across genes. For example, the following command sums the codon frequencies of each amino acid (and STP) across all genes, and then calculates the relative frequencies of the amino acids.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --sum \
                           --amino-acid \
                           --relative
{{ codestop }}

The first column of the output table has the header, 'gene_caller_ids', and the value, 'all', indicating that the data is aggregated across genes.

Statistics can also be aggregated across functional subsets of genes, as [detailed below](#aggregate-frequencies-across-functional-subsets-of-genes).

### Functions

Analyses can be extended with functional annotation information beyond the level of individual genes. First, without affecting the analysis itself, additional columns of [gene function annotations can be added](#gene-frequencies-with-function-information) to output tables of per-gene codon statistics. Second, per-function codon statistics can be calculated. In the first case, the output table path must be specified with `--output-file`/`-o` or, equivalently, `--gene-table-output`, and in the second, the output path must be specified with `--function-table-output`. Function codon frequencies treat functions as concatenations of genes, which can be useful for understanding the codon composition of 

Frequency statistics can be calculated for functions, treating them as concatenations of genes. This requires `--function-table-output` to specify the output path of a [table of frequencies per function](#function-frequencies) rather than per gene. Alternatively, `--output-file` or equivalently `--gene-table-output` should be used to produce output tables of per-gene frequencies, and with function options, [gene rows can contain columns of their function annotations](#gene-frequencies-with-function-information).

#### Inputs

There are multiple ways to define which functions and annotation sources should be considered. `--function-sources` can be used as a flag to consider all gene annotation sources. It can also be used with argument values to consider specific sources, such as 'COG20_FUNCTION' and 'KOfam'. `--function-accessions` and `--function-names` select functions from a single source value provided to `--function-sources`. The following example uses both accessions and names to select genes annotated with particular COG20 functions.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --function-sources COG20_FUNCTION \
                           --function-accessions COG0004 COG0005 \
                           --function-names "Xaa-Pro aminopeptidase" "Purine nucleoside phosphorylase"
{{ codestop }}

To use different functions from different sources, a tab-delimited file can be provided with `--functions-txt`. This headerless file must have three columns, for source, accession, and name of functions, respectively. A source entry is required in each row, along with a function accession or name, or both.

By default, selected function accessions or names do not need to be present in the input genomes: the program will return data for any selected function accessions or names that annotated genes. This behavior can be changed with the flag, `--expect-functions`, so that the program will throw an error in the absence of gene annotations with the provided accessions and names.

#### BRITE hierarchies

Genes are classified in KEGG BRITE functional hierarchies of KOs by %(anvi-run-kegg-kofams)s. For example, a bacterial SSU ribosomal protein is classified in the ribosome hierarchy, `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`. Codon frequencies can be calculated for genes classified at each level of the hierarchy, from the most general (all possible genes encoding the `Ribosome`) to the most specific (e.g., bacterial SSU ribosomal proteins in `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`). The following command returns summed codon frequencies for each annotated hierarchy level, e.g., a row each for `Ribosome`, `Ribosome>>>Ribosomal proteins`, `Ribosome>>>Ribosomal proteins>>>Bacteria`, and `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           --function-table-output output.txt \
                           --function-sources KEGG_BRITE
{{ codestop }}

#### Aggregate frequencies across genes in functions

`--sum` and `--average` operate on genes. When a function options are used, aggregate statistics are calculated for genes sharing a functional annotation. For example, the following command calculates average synonymous relative frequencies across genes annotated by each KO.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           --function-table-output output.txt \
                           --average \
                           --synonymous \
                           --function-sources KOfam
{{ codestop }}



`--sum` and `--average` operate on genes. When used with a function option, genes are subset by annotations of interest. For example, the following command calculates the average synonymous relative frequency across all genes annotated by KOs.


### Filter codons in analysis

Codons can be selected for inclusion in the analysis. There are also [options](#filter-codons-and-amino-acids-in-output) that select codons for inclusion in the output table without an effect on the codons analyzed. Removal of codons from the analysis causes relative frequency calculations to be normalized to the reduced set of retained codons.

#### Select codons and amino acids in analysis

`--exclude-codons` removes specified codons from the analysis. `--include-codons` retains select codons, to the exclusion of others, in the analysis. These options accept regular expressions. For example, `--exclude-codons ..A` removes codons ending in A, and `--include-codons G.. T..` retains codons starting with G or T.

`--exclude-amino-acids` removes codons for amino acids specified by their three-letter codes from the analysis. `--include-amino-acids` retains codons for select amino acids, to the exclusion of other amino acids, in the analysis.

The following example would return a table of relative frequencies just among the set of codons for amino acids that can be positively charged: `--include-amino-acids Arg Lys His`.

In calculating synonymous codon frequencies with `--synonymous`, and in the absence of `--exclude-amino-acids` and `--include-amino-acids`, the program behaves as if `--exclude-amino-acids` were passed `STP Met Trp` with the [standard genetic code](#genetic-code): stop codons and the single codons encoding Met and Trp are excluded calculation of synonymous codon relative frequencies. If `--exclude-amino-acids` or `--include-amino-acids` are used, then this behavior does not apply. For example, to also exclude Cys from the analysis, then use `--exclude-amino-acids Cys STP Met Trp`. If for some reason you wish to retain all codons in the synonymous codon analysis, including stop, Met, and Trp codons, then use `--exclude-amino-acids` as a flag without passing any values.

#### Dynamically exclude rarer codons and amino acids from analysis

`--pansequence-min-amino-acids` dynamically excludes codons for amino acids from the analysis, which can be useful in the calculation of synonymous codon relative frequencies. Consider the following case. Gene A has two Asn codons, AAT and AAC, resulting in synonymous relative frequencies of 0.5 AAT and 0.5 AAC. Gene B has one Asn codon, AAT, resulting in 1.0 AAT and 0.0 AAC. Gene C has 20 Asn codons, 10 AAT and 10 AAC, resulting in 0.5 AAT and 0.5 AAC. Asn has more statistical power in gene C than genes A or B for the calculation of codon usage bias or other purposes. If there are many genes like A and B and few like C in the dataset, then Asn can distort codon usage bias calculations, giving this rare amino acid the same weight as common ones. `--pansequence-min-amino-acids` addresses this problem by removing rare amino acids across the dataset, setting a minimum number of codons in a minimum number of genes required to retain the amino acid. For example, amino acids with <5 codons in ≥90%% of genes are excluded from the analysis with `--pansequence-min-amino-acids 5 0.9`. Note that if this option is used with `--relative`, relative codon frequencies are normalized to the summed frequency of retained codons, not all codons.

`--pansequence-min-codons` dynamically excludes select codons from the analysis in the same fashion. The use cases are not as clear for this option. Note the interactions that this option has with `--relative` and, particularly, `--synonymous`, since relative frequencies are normalized to the summed frequency of retained codons, not all codons. In synonymous calculations, retained codons in an amino acid set are considered, causing affected synonymous codon relative frequencies to deviate from the expected meaning.

### Filter by codon count

`--gene-min-codons` and `--function-min-codons` remove genes or functions with fewer than the specified number of codons from the analysis, which can be used to improve the statistical utility of relative frequencies.

#### Filter genes by codon count

`--gene-min-codons` sets the minimum number of codons required in a gene. This filter can be applied before and/or after removing codons with [other options](#filter-codons). The order is controlled by `--min-codon-filter`. Applied before, `--gene-min-codons` filters genes by the full set of codons. Applied after, `--gene-min-codons` filters genes by codons remaining after removing codons. `--min-codon-filter` can take three possible arguments: `length`, `remaining`, or, by default when codons are removed, `both`, which applies the `--gene-min-codons` filter both before and after codon removal.

`remaining` and `both` are not equivalent when `--pansequence-min-amino-acids` or `--pansequence-min-codons` is used for dynamic codon exclusion. Amino acids are removed based on their frequency in a proportion of genes, so removing shorter genes by length before removing amino acids can affect which amino acids are dynamically excluded. For example, `--gene-min-codons 250` removes genes with fewer than 250 codons, and `--pansequence-min-amino-acids 5 0.9` removes amino acids with <5 codons in ≥90%% of genes. Say Cys has <5 codons in 92%% of all genes and 88%% of genes with at least 250 codons. `--min-codon-filter remaining` would cause the Cys codons to first be removed from genes, since 92%% ≥ 90%%, and then the remaining genes to be filtered by codon count minus Cys codons. `--min-codon-filter both` (the default) would cause genes first to be filtered by total codon count. Then, Cys codons would be retained in the analysis, since 88%% < 90%%. No more genes would fall under the threshold of the minimum codon filter in the second pass due to the retention of Cys.

#### Filter functions by codon count

`--function-min-codons` removes reported functions lacking the specified minimum number of codons. Function codon count filters occur after gene codon count filters: the set of genes contributing to function codon frequencies can be restricted by also using `--gene-min-codons`.

### Filter codons and amino acids in output

Codons and amino acids (with `--return-amino-acids`) can be selected for inclusion in the output table. There are also [options](#filter-codons-in-analysis) for selecting codons in the analysis rather than just changing how the output table is displayed.

#### Select codons and amino acids in output

`--dont-report-codons` drops columns for the specified codons from the output table. `--report-codons` retains columns for select codons, to the exclusion of others, in the output table. Again, these options accept regular expressions. For example, `--dont-report-codons A.[CT]` drops columns for codons starting with A and ending with C or T, and `--report-codons [CG][CG][CG]` retains columns for codons with C or G at all positions.

`--dont-report-amino-acids` drops codon columns from the output table for the specified amino acids. `--report-amino-acids` retains codon columns in the output table for select amino acids, to the exclusion of others.

The output table has columns for amino acids rather than codons with `--return-amino-acids`. In this case, `--dont-report-amino-acids` and `--report-amino-acids` select amino acid columns rather than codon columns.

#### Mask rarer codons and amino acids in output

`--sequence-min-codons` and `--sequence-min-amino-acids` can be useful for screening statistically significant data.

`--sequence-min-codons` replaces values for rarer codons within each gene or function row with NaN ([or 0](#replace-nan-with-zero-in-output) using `--infinity-to-zero`). For example, codons with <3 codons in a gene or function are discarded in the table with `--sequence-min-codons 3`.

`--sequence-min-amino-acids` replaces values for rarer amino acids within each gene or function row with NaN ([or 0](#replace-nan-with-zero-in-output) with `--infinity-to-zero`). For example, amino acids with <5 codons in a gene or function are discarded in the table with `--sequence-min-amino-acids 5`.

### Replace NaN with 0 in output

`--infinity-to-zero` replaces null values in the output table with 0. Codons for absent amino acids have null synonymous relative frequency values. Be careful using this option, since NaN replaced by 0 in synonymous analysis does not have the intended meaning, introducing an inaccuracy which can propagate in downstream calculations like codon usage bias.

### Genetic code

`--encodings-txt` takes a tab-delimited file of two columns to make changes to the standard genetic code. Each entry in the first column is a codon. Each entry in the second collumn is the three-letter code for the decoded amino acid. For example, to recode the stop codon, TGA, as Trp, a row would contain 'TGA' and 'Trp'. In calculating synonymous codon frequencies with `--synonymous`, and in the absence of `--exclude-amino-acids` and `--include-amino-acids`, the program excludes stop codons and single-codon amino acids from the analysis.
