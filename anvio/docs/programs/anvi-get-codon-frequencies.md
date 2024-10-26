This program **calculates codon or amino acid frequencies from genes or functions**.

A range of options allows calculation of different frequency statistics.

## Basic commands

### Gene frequencies

The following command produces a table of codon frequencies from coding sequences in the input contigs database. The first column of the tab-delimited output file contains anvi'o gene caller IDs. Other columns contain codon frequencies. By default, these columns are headed by codons, e.g., "ACG". The flag, `--header-amino-acids`, includes the decoded amino acid with the codon in the header, e.g, "ThrACG".

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --header-amino-acids
{{ codestop }}

### Function frequencies

The following command produces a table of function frequencies rather than gene frequencies. By using `--function-sources` without any arguments, the output will include every %(functions)s source available in a given %(contigs-db)s, e.g., `KOfam`, `KEGG_BRITE`, `Pfam` (you can always see the complete list of function sources in your contigs database by running the program %(anvi-db-info)s on it). The first three columns of the output table contain function sources, accessions, and names, e.g., "COG20_FUNCTION", "COG0001", "Glutamate-1-semialdehyde aminotransferase (HemL) (PDB:2CFB)". Other tables in the column contain codon frequencies.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           --function-sources \
                           --function-table-output function_output.txt
{{ codestop }}

### Gene frequencies with function information

In contrast to the [prior example](#function-frequencies), the following command produces a table of gene codon frequencies, but has an entry for every gene/function pair. This enables analysis of the genes encoding functions, including genes encoding the same function. (The function table output in the previous example is directly derived from this data by summing gene frequencies with the same function source.)

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           --function-sources \
                           --gene-table-output gene_output.txt
{{ codestop }}

### Codon frequencies from multiple internal and external genomes

The following command produces a table of gene codon frequencies from coding sequences in multiple genomes, here a combination of %(internal-genomes)s and %(external-genomes)s listed in respective files. The first column of the output table contains genome names. Other columns are codon frequencies.

{{ codestart }}
anvi-get-codon-frequencies -i %(internal-genomes)s \
                           -e %(external-genomes)s \
                           -o output.txt
{{ codestop }}

## Option examples

The following tables show how to use options to get certain results.

### _Different_ frequency statistics

| Get | Options |
| --- | ------- |
| Codon absolute frequencies | |
| Codon relative frequencies | `--relative` |
| [Synonymous (per-amino acid) codon relative frequencies](#synonymous-codon-relative-frequencies) | `--synonymous` |
| Amino acid frequencies | `--amino-acid` |
| Amino acid relative frequencies | `--amino-acid --relative` |
| [Summed frequencies across genes](#frequencies-across-genes) | `--sum` |
| [Average frequencies across all genes](#frequencies-across-genes) | `--average` |
| [Synonymous relative summed frequencies across genes](#frequencies-across-genes) | `--sum --synonymous` |
| [Summed frequencies across genes annotated by each function source](#frequencies-across-genes-with-function) | `--sum --function-sources` |
| [Relative summed frequencies across genes with KOfam annotations](#frequencies-across-genes-with-function) | `--sum --relative --function-sources KOfam` |

### Frequencies from _sets of genes with shared functions_

| Get | Options |
| --- | ------- |
| All function annotation sources | `--function-sources` |
| [All KEGG BRITE categories of KEGG Orthologs (KOs)](#brite-hierarchies) | `--function-sources KEGG_BRITE` |
| All KOs and all Pfam entries | `--function-sources KOfam Pfam` |
| [Certain KEGG BRITE categories](#brite-hierarchies) | `--function-sources KEGG_BRITE --function-names Ribosome Ribosome>>>Ribosomal proteins` |
| [Certain KO accessions](#inputs) | `--function-sources KOfam --function-accessions K00001 K00002` |
| [Certain BRITE categories and KO accessions](#inputs) | `--select-functions-txt path/to/select_functions.txt` |

### Frequencies from _selections of genes_

| Get | Options |
| --- | ------- |
| From contigs database | `--contigs-db path/to/contigs.db` |
| From collection of internal genomes | `--contigs-db path/to/contigs.db --profile-db path/to/profile.db --collection-name my_bins` |
| From internal genome | `--contigs-db path/to/contigs.db --profile-db path/to/profile.db --collection-name my_bins --bin-id my_bin` |
| From internal genomes listed in a file | `--internal-genomes path/to/genomes.txt` |
| From external genomes (contigs databases) listed in a file | `--external-genomes path/to/genomes.txt` |
| With certain gene IDs | `--gene-caller-ids 0 2 500` |
| With certain gene IDs or genes annotated with certain KOs | `--gene-caller-ids 0 2 500 --function-sources KOfam --function-accessions K00001` |

### _Filter codons, genes, and functions_ that are analyzed and reported

| Get | Options |
| --- | ------- |
| [Exclude codons ending with A from analysis](#select-codons) | `--exclude-codons ..A` |
| [Include codons starting with G or T in analysis](#select-codons) | `--include-codons G.. T..` |
| [Exclude codons starting with A and ending C or T from output but not analysis](#select-codons) | `--dont-report-codons A.[CT]` |
| [Only include codons with C or G at all positions in output](#select-codons) | `--report-codons [CG][CG][CG]` |
| [Exclude stop codons and single-codon amino acids from analysis](#select-amino-acids) | `--exclude-amino-acids STP Met Trp` |
| [Only include codons encoding certain amino acids in analysis](#select-amino-acids) | `--include-amino-acids Leu Ile` |
| [Exclude stop codons and single-codon amino acids from output but not analysis](#select-amino-acids) | `--dont-report-amino-acids STP Met Trp` |
| [Only include codons encoding certain amino acids in output](#select-amino-acids) | `--report-amino-acids Leu Ile` |
| [Replace codons for amino acids with <5 codons in the gene or function with NaN](#exclude-rarer-amino-acids-from-output) | `--sequence-min-amino-acids 5` |
| [Exclude codons for amino acids with <5 codons in >90%% of genes](#dynamically-exclude-rarer-amino-acids-from-analysis) | `--pansequence-min-amino-acids 5 0.9` |
| [Exclude genes shorter than 250 codons](#filter-genes-by-codon-count) | `--gene-min-codons 250` |
| [Exclude functions with <250 codons](#filter-functions-by-codon-count) | `--function-min-codons 250` |
| [Exclude genes shorter than 250 codons from contributing to function codon frequencies](#filter-functions-by-codon-count) | `--gene-min-codons 250 --function-sources` |

## Option details

### Synonymous codon relative frequencies

The `--synonymous` flag returns the relative frequency of each codon encoding the same amino acid, e.g., 0.4 GAC and 0.6 GAT for Asp. By default, stop codons and single-codon amino acids (Met ATG and Trp TGG) in the standard translation table are excluded from the analysis and output, equivalent to using `--exclude-amino-acids STP Met Trp` for other frequency statistics.

### Frequencies across genes

`--sum` and `--average` produce a table with a single row of frequencies from across genes. For example, the following command sums the codon frequencies of each decoded amino acid (and STP) across all genes, and then calculates the relative frequencies of the amino acids.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --sum \
                           --amino-acid \
                           --relative
{{ codestop }}

The first column of the output table has the header, 'gene_caller_ids', and the value, 'all', indicating that the data is aggregated across genes.

#### Frequencies across genes with function

`--sum` and `--average` operate on genes. When used with a function option, the program subsets the genes annotated by the functions of interest. With `--average`, it calculates the average frequency across genes rather than functions (sums of genes with functional annotation). For example, the following command calculates the average synonymous relative frequency across genes annotated with KOs.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --average \
                           --synonymous \
                           --function-sources KOfam
{{ codestop }}

### Functions

Functions and function annotation sources can be provided to [subset genes](#frequencies-across-genes-with-function) and to [calculate statistics for functions](#function-codon-frequencies) in addition to genes.

Output tables can report codon frequencies of functionally annotated genes by using `--output-file` or equivalently `--gene-table-output` to specify the output table path. Alternatively, output tables can report codon frequencies summed across genes in functions by using `--function-table-output`.

#### Inputs

There are multiple options to define which functions and sources should be considered. `--function-sources` as can be used with arguments, such as 'COG14_FUNCTION' or 'KOfam', or without arguments as a flag to consider all sources that had been used to annotate genes. `--function-accessions` and `--function-names` select functions from a single source provided by `--function-sources`. The following example uses both accessions and names to select COG functions.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --function-sources COG14_FUNCTION \
                           --function-accessions COG0004 COG0005 \
                           --function-names "Ammonia channel protein AmtB" "Purine nucleoside phosphorylase"
{{ codestop }}

To use different functions from different sources, a tab-delimited file can be provided with `--functions-txt`. This headerless file must have three columns, for source, accession, and name of functions, respectively. A source entry is required in each row, along with a function accession or name, or both.

By default, selected function accessions or names do not need to be present in the input genomes; the program will return data for any selected function accessions or names that annotated genes. This behavior can be changed using the flag, `--expect-functions`, so that the program will throw an error when any of the selected accessions or names are absent.

#### BRITE hierarchies

Genes are classified in KEGG BRITE functional hierarchies of KOs by %(anvi-run-kegg-kofams)s. For example, a bacterial SSU ribosomal protein is classified in the ribosome hierarchy, `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`. Codon frequencies can be calculated for genes classified at each level of the hierarchy, from the most general (all possible genes encoding the `Ribosome`) to the most specific (e.g.,  `Ribosome>>>Ribosomal proteins>>>Bacteria>>>Small subunit`). The following command returns summed codon frequencies for each annotated hierarchy level -- given the example ribosomal protein annotated in the contigs database, the output would include four rows for the genes in each level from `Ribosome` to `Small subunit`.

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \
                           -o output.txt \
                           --function-sources KEGG_BRITE
{{ codestop }}

### Filter codons

Codons can be excluded from or included in the analysis or just the output table. Exclusion from the analysis causes relative frequency calculations to be relative to the restricted set of retained codons.

#### Select codons

`--exclude-codons` removes codons from the analysis. `--include-codons` retains select codons, to the exclusion of others, in the analysis. These options accept regular expressions. For example, `--exclude-codons ..A` removes codons ending in A, and `--include-codons G.. T..` retains codons starting with G or T.

`--dont-report-codons` drops columns for the specified codons from the output table. `--report-codons` retains columns for select codons, to the exclusion of others, in the output table. Again, these options accept regular expressions. For example, `--dont-report-codons A.[CT]` drops columns for codons starting with A and ending with C or T, and `--report-codons [CG][CG][CG]` retains columns for codons with C or G at all positions.

#### Select amino acids

`--exclude-amino-acids` removes codons encoding amino acids specified by their three-letter codes from the analysis. `--include-amino-acids` retains codons encoding specified amino acids, to the exclusion of other amino acids, in the analysis.

In calculating synonymous codon frequencies with `--synonymous`, and in the absence of `--exclude-amino-acids` and `--include-amino-acids`, the program behaves as if `--exclude-amino-acids` were passed `STP Met Trp`: stop codons and the single codons encoding Met and Trp are excluded calculation of synonymous codon relative frequencies. If `--exclude-amino-acids` or `--include-amino-acids` are used, then this behavior does not apply. For example, to also exclude Cys from the analysis, then use `--exclude-amino-acids Cys STP Met Trp`.

The following example would return a table of relative frequencies just among the set of codons for amino acids that can be positively charged: `--include-amino-acids Arg Lys His`.

#### Exclude rarer amino acids from output

`--sequence-min-amino-acids` only affects how the output table is displayed, replacing codon values for rarer amino acids within each gene or function row with NaN (or 0 with `--infinity-to-zero`). For example, amino acids with <5 codons in the gene or function are discarded in the table with `--sequence-min-amino-acids 5`. This may be useful for screening statistically significant data.

#### Dynamically exclude rarer amino acids from analysis

`--pansequence-min-amino-acids` dynamically excludes amino acids from the analysis and can be useful in the calculation of synonymous codon relative frequencies. Consider the following case. Gene A has two Asn codons, AAT and AAC, resulting in synonymous relative frequencies of 0.5 AAT and 0.5 AAC. Gene B has one Asn codon, AAT, resulting in 1.0 AAT and 0.0 AAC. Gene C has 20 Asn codons, 10 AAT and 10 AAC, resulting in 0.5 AAT and 0.5 AAC. Asn has more statistical power in gene C than genes A or B for the calculation of codon usage bias or other purposes. If there are many genes like A and B and few like C in the dataset, then Asn can distort codon usage bias calculations. `--pansequence-min-amino-acids` addresses this problem by removing rare amino acids across the dataset, setting a minimum number of codons in a minimum number of genes required to retain the amino acid. For example, amino acids with <5 codons in >90%% of genes are excluded from the analysis with `--pansequence-min-amino-acids 5 0.9`.

### Filter by codon count

`--gene-min-codons` and `--function-min-codons` remove genes or functions with fewer than the specified number of codons from the analysis, which can be used to improve the statistical utility of relative frequencies.

#### Filter genes by codon count

`--gene-min-codons` sets the minimum number of codons required in a gene. This filter can be applied before and/or after removing codons with [other options](#filter-codons). The order is controlled by `--min-codon-filter`. Applied before, `--gene-min-codons` filters genes by the full set of codons. Applied after, `--gene-min-codons` filters genes by codons remaining after removing codons. `--min-codon-filter` can take three possible arguments: `length`, `remaining`, or, by default when codons are removed, `both`, which applies the `--gene-min-codons` filter both before and after codon removal.

`remaining` and `both` are not equivalent when `--pansequence-min-amino-acids` is used for dynamic amino acid exclusion. Amino acids are removed based on their frequency in a proportion of genes, so removing shorter genes by length before removing amino acids can affect which amino acids are dynamically excluded. For example, `--gene-min-codons 250` removes genes with fewer than 250 codons, and `--pansequence-min-amino-acids 5 0.9` removes amino acids with <5 codons in >90%% of genes. Say Cys has <5 codons in 92%% of all genes and 88%% of genes with at least 250 codons. `--min-codon-filter remaining` would cause the Cys codons to first be removed from genes, since 92%% > 90%%, and then the remaining genes to be filtered by codon count minus Cys codons. `--min-codon-filter both` (the default) would cause genes first to be filtered by total codon count. Then, Cys codons would be retained in the analysis, since 88%% ≤ 90%%. No more genes would fall under the threshold of the minimum codon filter in the second pass due to the retention of Cys.

#### Filter functions by codon count

`--function-min-codons` removes reported functions lacking the specified minimum number of codons. Function codon count filters occur after gene codon count filters: the set of genes contributing to function codon frequencies can be restricted by also using `--gene-min-codons`.
