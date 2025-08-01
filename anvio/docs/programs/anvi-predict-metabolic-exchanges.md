%(anvi-predict-metabolic-exchanges)s uses the metabolic capacity encoded in two genomes to predict which metabolites could be exchanged between two organisms, and which metabolites are unique to one of the organisms. The output files produced by the program are described in %{metabolite-exchange-predictions}s.

By leveraging the %(reaction-network)s infrastructure of anvi'o, the program examines the overlap between the metabolic reactions that each of the two organisms can catalyze and identifies which metabolic compounds fall into one of the following categories:

1. can be produced by only one organism but consumed by the other (or both)
2. can be consumed by only one organism but produced by the other (or both)
3. can be produced and/or consumed by only one organism

It reports the metabolites in categories 1 and 2 as 'potentially-exchanged compounds', and the metabolites in category 3 as 'unique' compounds. Here is a table of what it looks for (everything marked with a star is reported by the program, and everything with a gray background is ignored):

![Table of possible combinations for metabolite production and consumption, some of which can indicate a potential exchange](../../images/anvi-predict-metabolic-exchanges-01.png){:.center-img .width-70}

Briefly, the program works by checking for these conditions for every metabolite in the genomes' (merged) reaction networks. It does so in two orthogonal ways:

1. **Walking over KEGG Pathway Maps.** For each compound that is part of a KEGG Pathway Map, we examine the chains of reactions that produce or consume the compound within the Pathway Map for each organism. This not only allows us to identify the aforementioned situations of potentially-exchanged or unique compounds, but also allows us to compute evidence for a given prediction (such as length and overlap of reaction chains) that can help you filter and interpret the output.
2. **Isolating metabolites in the merged reaction network.** For each compound that is part of the genomes' merged reaction network, we examine the subset of the network centered around it and determine whether it fits into one of the situations described above. This strategy currently doesn't offer any supporting evidence for a given prediction (we only look at the potential 'transfer point' of the compound), but it does allow predictions for compounds that are not included in any KEGG Pathway Maps.

Either of these strategies can be skipped, in case you prefer to use only one method. If you want to know more, check the Technical Details section.

## Prerequisites to using this program

This program uses the %(reaction-network)s stored in each %(contigs-db)s to search for potentially-exchanged or unique compounds. Since the %(reaction-network)s and this program both rely on the [ModelSEED](https://modelseed.org/) and [KEGG](https://www.kegg.jp/) databases (especially KEGG Pathway Maps and KOfam annotation profiles for KEGG Orthlogs), you need to have access to that data on your computer. You also need to have KEGG Ortholog (KO) annotations (%(kegg-functions)s) in your input genomes (which translates to having the functional annotation source 'KOfam' in the %(contigs-db)ss).

Here are the steps you need to run before this program:
1. %(anvi-setup-kegg-data)s to get data from KEGG onto your computer. This step only needs to be done once -- if you've already ran this setup program in the past, you normally don't have to do it again. When in doubt, skip this step and let anvi'o tell you if you are missing something.
2. %(anvi-setup-modelseed-database)s to get data from ModelSEED onto your computer. Like the previous step, this only needs to be done once.
3. %(anvi-run-kegg-kofams)s to annotate _each_ of your %(contigs-db)ss with KEGG Ortholog protein families.
4. %(anvi-reaction-network)s to create a reaction network for _each_ of of your %(contigs-db)ss.

If you've done all that, you are good to go.

## How to run this program

### Running on a single pair of genomes

Provide both contigs databases for your genomes and a prefix for the output files:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX
{{ codestop }}

The %{metabolite-exchange-predictions}s output file names will start using the prefix you provided.

## Adjustable Parameters

### Setting some compound IDs as equivalent

The ModelSEED database sometimes has multiple compound ID numbers for what (we think) should be the same metabolite, at least for the purposes of predicting exchanges. A prime example of this is amino acids, which can have compounds where their chirality is specified (as in `L-Lysine (cpd00039)`) and where the chirality is generic (as in `Lysine (cpd19182)`). This is not so much a problem when we are using KEGG Pathway Maps to predict exchanges, but can lead to missing predictions when using the reaction network.

For amino acids, the program can automatically detect which compound IDs should be considered equivalent, and take that into account when looking for potentially-exchanged compounds. If you use the `--use-equivalent-amino-acids` flag, the program will search through the ModelSEED database for any conventional amino acids (plus Selenocysteine and Pyrrolysine) that have both an 'L-' version and a chiral-unspecific version, and set those two compounds equivalent to each other.

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --use-equivalent-amino-acids
{{ codestop }}

If you do this, you will get an output file listing the amino acid compounds that were deemed equivalent, so you can make sure you agree with them.

If you want to specify a custom set of equivalent compound IDs, you can instead provide an %{equivalent-compounds-txt}s file to the `--custom-equivalent-compounds-file` parameter:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --custom-equivalent-compounds-file %{equivalent-compounds-txt}s
{{ codestop }}

### Using only one prediction method

If you want to skip the first prediction step of walking over KEGG Pathway Maps to find potential exchanges, use the `--no-pathway-walk` flag:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --no-pathway-walk
{{ codestop }}

If you want to skip the second prediction step of examining the local reaction network around each compound, use the `--pathway-walk-only` flag:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --pathway-walk-only
{{ codestop }}

It is hopefully understandable that these two flags are incompatible with each other.

### Changing the number of allowed gaps in the Pathway Map walks

The `--maximum-gaps` parameter applies to the first prediction step of walking over KEGG Pathway Maps, and allows a certain number of missing enzyme annotations in the reaction chains. By default, we don't allow any gaps, but if you think missing annotations in either genome might be throwing off your predictions, you can set this parameter to an integer greater than 0:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --maximum-gaps 1
{{ codestop }}

Changing this parameter will mostly affect the evidence attributed to a given prediction (i.e., length and overlap of reaction chains), but there could be some cases where increasing the gap number enables new predictions to be made.

### Using non-default data directories

If you set up your KEGG or ModelSEED data in a custom directory, you can make sure this program knows where to find it by providing the paths:

{{ codestart }}
anvi-predict-metabolic-exchanges -c1 %(contigs-db)s -c2 %(contigs-db)s \
                                 -O ANY_PREFIX \
                                 --kegg-data-dir /path/to/directory/KEGG \
                                 --modelseed-data-dir /path/to/directory/MODELSEED
{{ codestop }}

The data directories are relevant for loading the %(reaction-network)s in the contigs database, so it is best to use the same data directories that were utilized when running %(anvi-reaction-network)s.
