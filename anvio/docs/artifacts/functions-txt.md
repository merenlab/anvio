This artifact is a TAB-delimited file that **associates genes and functions**. 

The user can generate this file to import gene functions into a %(contigs-db)s via %(anvi-import-functions)s or can acquire this file by recovering it from a %(contigs-db)s via %(anvi-export-functions)s. It is also the output of %(anvi-search-functions)s which searches for specific terms in your functional annotations.

In general, this is the simplest way to get gene functions into anvi'o, and all downstream analyses, including pangenomics. For other ways to get gene functions into anvi'o you can take a look at [this page](http://merenlab.org/2016/06/18/importing-functions/). 


## Simple matrix file format


The TAB-delimited file for this artifact has five columns:

1. `gene_callers_id`: The gene caller ID recognized by anvi'o (see the note below).
2. `source`: The name of the functional annotation source (i.e., the database that you got this function data from).
3. `accession`: A unique accession id per function, better if a single word.
4. `function`: Full name / description of the function.
5. `e_value`: The significance score of this annotation, where zero is maximum significance. This information may be used by anvi'o in operations that require filtering of functions based on their significance.

Through this file format **you can import functions from any source** into anvi'o, whether those sources are commonly used programs to annotate genes with functions or your ad hoc manual curations for genes of interest. But **please note while there are many ways to have your genes annotated with functions, there is only one way to make sure the gene caller ids anvi'o knows will match perfectly to the gene caller ids in your input file**. The best way to ensure that linkage is to export your gene DNA or amino acid sequences for your an %(contigs-db)s using the anvi'o program `anvi-get-sequences-for-gene-calls`.

## An example matrix

Here is an example file that matches to this format that can be used with %(anvi-import-functions)s to import functions into a %(contigs-db)s:

|gene_callers_id|source|accession|function|e_value|
|:--|:--:|:--:|:--|:--:|
|1|Pfam|PF01132|Elongation factor P (EF-P) OB domain|4e-23|
|1|Pfam|PF08207|Elongation factor P (EF-P) KOW-like domain|3e-25|
|1|TIGRFAM|TIGR00038|efp: translation elongation factor P|1.5e-75|
|2|Pfam|PF01029|NusB family|2.5e-30|
|2|TIGRFAM|TIGR01951|nusB: transcription antitermination factor NusB|1.5e-36|
|3|Pfam|PF00117|Glutamine amidotransferase class-I|2e-36|
|3|Pfam|PF00988|Carbamoyl-phosphate synthase small chain, CPSase domain|1.2e-48|
|3|TIGRFAM|TIGR01368|CPSaseIIsmall: carbamoyl-phosphate synthase, small subunit|1.5e-132|
|4|Pfam|PF02787|Carbamoyl-phosphate synthetase large chain, oligomerisation domain|1.4e-31|
|4|TIGRFAM|TIGR01369|CPSaseII_lrg: carbamoyl-phosphate synthase, large subunit|0|
|5|TIGRFAM|TIGR02127|pyrF_sub2: orotidine 5'-phosphate decarboxylase|1.9e-59|
|6|Pfam|PF00625|Guanylate kinase|5.7e-39|
|6|TIGRFAM|TIGR03263|guanyl_kin: guanylate kinase|3.5e-62|
|8|Pfam|PF01192|RNA polymerase Rpb6|4.9e-13|
|8|TIGRFAM|TIGR00690|rpoZ: DNA-directed RNA polymerase, omega subunit|1.7e-20|
|9|TIGRFAM|TIGR01034|metK: methionine adenosyltransferase|2.5e-169|
|11|Pfam|PF13419|Haloacid dehalogenase-like hydrolase|2.8e-27|
|11|TIGRFAM|TIGR01509|HAD-SF-IA-v3: HAD hydrolase, family IA, variant 3|1.2e-11|
|12|Pfam|PF00551|Formyl transferase|1.4e-34|
|12|TIGRFAM|TIGR00460|fmt: methionyl-tRNA formyltransferase|2.9e-70|
|13|Pfam|PF12710|haloacid dehalogenase-like hydrolase|2.3e-14|
|13|TIGRFAM|TIGR00338|serB: phosphoserine phosphatase SerB|4.9e-76|
|13|TIGRFAM|TIGR01488|HAD-SF-IB: HAD phosphoserine phosphatase-like hydrolase, family IB|6e-29|
|14|Pfam|PF00004|ATPase family associated with various cellular activities (AAA)|7.7e-45|
|14|Pfam|PF16450|Proteasomal ATPase OB/ID domain|1.8e-34|
|14|TIGRFAM|TIGR03689|pup_AAA: proteasome ATPase|1e-206|
|(...)|(...)|(...)|(...)|(...)|


Please note that,

* Not every gene call has to be present in the matrix,

* It is OK if there are multiple annotations from the same source for a given gene call,

* It is OK if a give gene is annotated only by a single source.

* If the **accession** information is not available to you, it is OK to leave it blank (but it will prevent you from being able to use some toys, such as functional enrichment analyses later for pangenomes).

* If you have no e-values associated with your annotations, it is OK to put `0` for every entry (you should make sure you keep this in mind for your downstream analyses that may require filtering of weak hits).

* If there are multiple annotations from a single source for a single gene call, anvi'o uses e-values in this file to use only the most significant one to show in interfaces.