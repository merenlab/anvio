This is an R script written by [Amy Willis](https://github.com/adw96) that performs an enrichment test on things (functions, metabolic modules) in groups of things (genomes in a %(pan-db)s, samples that have output from %(anvi-estimate-metabolism)s, %(external-genomes)s or %(internal-genomes)s).

## Technical details
You will find a brief description of how this script works in Alon's "Behind the Scenes" note in [the pangenomics tutorial](https://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome). Better yet, check out the methods section of Alon's paper, found as a pre-print [here](https://www.biorxiv.org/content/10.1101/2020.04.29.069278v2).

## Using this script through other programs
Most people will not use this script directly, and instead run it by calling the driver programs %(anvi-compute-enrichment-scores)s or %(anvi-get-enriched-functions-per-pan-group)s. Those programs massage the data from various sources to produce an input file suitable to this script before running it, thereby saving you a lot of work.

{:.notice}
Want to compute enrichment scores for something not supported by those driver programs? Please [contact the developers](https://merenlab.org/2019/10/07/getting-help/#getting-help-from-the-community) to offer your suggestions on how we can extend their functionality. Alternatively, you can try using this script directly - see the section below.

## Using this script directly
For anyone interested in computing enrichment scores for something that we don't already support, you can use this script directly if you get your data into the right input format, which is described below. Please remember to contact [Amy Willis](https://github.com/adw96) for information on how to cite this script! And if you want to use what you've learned to extend anvi'o's functionality in this department, even better. Contact the developers in that case :)

To run the script, you must provide a properly-formatted input file and the name of an output file.

{{ codestart }}
anvi-script-enrichment-stats --input=input.txt --output=output.txt
{{ codestop }}

The input file format is expected to be tab-delimited, with the following columns:

| Column | Header | Contents |
|-----|-----|-----|
| 1 | does not matter (examples: COG_FUNCTION, KEGG_MODULE) | A string which is the annotation or name of the entity (function, module, etc) in this row |
| 2 | 'accession' | A string which is the unique key (or ID) of the entity in this row |
| 3 | does not matter (examples: gene_clusters_ids, sample_ids) | This column doesn't matter for the enrichment analysis and could hold anything you want, even NAs. Historically it has been used for: A comma-separated list of strings describing which gene clusters (in a pangenome) the row's function belongs to, or a comma-separated list of strings describing which (metagenome) samples a KEGG module was found in. |
| 4 | 'associated_groups' | a comma-separated list of strings describing in which groups the row entity was found to be present more than expected. See Alon's note [here](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome) for a description on how these associated groups are calculated. Anvi'o utils.py has get_enriched_groups() function to use when calculating these groups programmatically within anvi'o. |
| Subsequent x columns, where x is the number of groups | 'p_[group]', where [group] is replaced with the name/acronym for the group | A number (float or double) between 0 and 1 which is the proportion of genomes/samples in the group that the row entity was found in. |
| Subsequent x columns, where x is the number of groups | 'N_[group]', where [group] is replaced with the name/acronym for the group | A number (integer) which is the total number of genomes/samples in the group |

NOTE: in practice, it does not matter if the 'p_' and 'N_' columns are out of order/intermixed, as long as there is one pair of them per group

Here is an example input based on the pangenomics tutorial (with NCBI COG functions):

| COG_FUNCTION | accession | gene_clusters_ids | associated_groups | p_HL | p_LL | N_HL | N_LL |
| NADPH-dependent glutamate synthase beta chain or related oxidoreductase | COG0493 | GC 00003904 | LL | 0 | 0.1818 | 20 | 11 |
| Urease beta subunit | COG0832 | GC 00000921 | HL | 0.95 | 0.4545 | 20 | 11 |
| Phosphoribosyl-AMP cyclohydrolase | COG0139 | GC 00000450 | 1 | 1 | 20 | 11 |
| Predicted phage phi-C31 gp36 major capsid-like protein | COG4653 | GC 00005619 | HL | 0.05 | 0 | 20 | 11 |
| Cytochrome c biogenesis protein CcdA | COG0785 | GC 00000956 | HL | 1 | 0.2727 | 20 | 11 |
| FoF1-type ATP synthase, epsilon subunit | COG0355 | GC 00000685 | 1 | 1 | 20 | 11 |
| Predicted mannosyl-3-phosphoglycerate phosphatase, HAD superfamily | COG3769 | GC 00001447,GC 00002275,GC 00003041 | HL | 1 | 0.8182 | 20 | 11 |
| Predicted flavoprotein YhiN | COG2081 | GC 00000967,GC 00002155 | HL | 1 | 0.9091 | 20 | 11 |
| Phosphatidylserine/phosphatidylglycerophosphate/cardiolipin synthase or related enzyme | COG1502 | GC 00003553 | LL | 0 | 0.0909 | 20 | 11 |
(...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) |

Running this script on a file like this will add the following columns into it: 'enrichment_score', 'unadjusted_p_value', 'adjusted_q_value'
