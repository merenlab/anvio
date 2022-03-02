This is a TAB-delimited output file that describes enrichment scores and associated groups for functions or metabolic modules in groups of genomes or samples. It is produced by %(anvi-compute-functional-enrichment)s.

## General format

Each row in the matrix describes an entity (a function, functional association of a gene cluster, or metabolic module) that is associated with one or more groups of samples or genomes. These are listed with the highest enrichment scores displayed first.

The following columns of information are listed in the file:

- the name of the enriched entity, which can be a functional association, metabolic module, or function. The header of this column is either your functional annotation source, OR 'KEGG_MODULE' if you are working with metabolic modules
- enrichment_score: a measure of much this particular entity is enriched in the group it is associated with (i.e., measures how unique this entity [see column 1] is to this group(s) [see column 5])
- unadjusted_p_value: the significance value of the hypothesis test for enrichment, unadjusted for multiple hypothesis testing
- adjusted_q_value: the adjusted p-value after taking into account multiple hypothesis testing
- associated groups: the list of groups that this entity is associated with
- accession: a function accession number or KEGG module number (depends on your input option for %(anvi-compute-functional-enrichment)s)
- a list of gene cluster ids, sample names, or genome names that this entity is found in (also depends on your input option for %(anvi-compute-functional-enrichment)s)
- p values for each group: gives the proportion of the group's member genomes or samples in which this entity was found.
- N values for each group: gives the total number of genomes or samples in each group.

## A specific example - enriched functions in pangenomes

When you run %(anvi-compute-functional-enrichment)s (with input option 1) to compute enrichment scores for functions in a pangenome, the resulting matrix describes the gene cluster-level functional associations that are enriched within specific groups of your pangenome. This is described in more detail [in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome).

Here is a more concrete example (the same example as in the [pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)). Note that that tutorial uses `COG_FUNCTION` as the functional annotation source, and has `LL` (low light) and `HL` (high light) as the two pan-groups.

|COG_FUNCTION | enrichment_score | unadjusted_p_value | adjusted_q_value | associated_groups | accession | gene_clusters_ids | p_LL | p_HL | N_LL | N_HL|
|-- | -- | -- | -- | -- | -- | -- | -- | -- | --| --|
|Proteasome lid subunit RPN8/RPN11, contains Jab1/MPN domain metalloenzyme (JAMM) motif | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG1310 | GC_00002219, GC_00003850, GC_00004483 | 1 | 0 | 11 | 20|
|Adenine-specific DNA glycosylase, acts on AG and A-oxoG pairs | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG1194 | GC_00001711 | 1 | 0 | 11 | 20|
|Periplasmic beta-glucosidase and related glycosidases | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG1472 | GC_00002086, GC_00003909 | 1 | 0 | 11 | 20|
|Single-stranded DNA-specific exonuclease, DHH superfamily, may be involved in archaeal DNA replication intiation | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG0608 | GC_00002752, GC_00003786, GC_00004838, GC_00007241 | 1 | 0 | 11 | 20|
|Ser/Thr protein kinase RdoA involved in Cpx stress response, MazF antagonist | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG2334 | GC_00002783, GC_00003936, GC_00004631, GC_00005468 | 1 | 0 | 11 | 20|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|Signal transduction histidine kinase | -7.34E-41 | 1 | 1 | NA | COG5002 | GC_00000773, GC_00004293 | 1 | 1 | 11 | 20|
|tRNA A37 methylthiotransferase MiaB | -7.34E-41 | 1 | 1 | NA | COG0621 | GC_00000180, GC_00000851 | 1 | 1 | 11 | 20|
