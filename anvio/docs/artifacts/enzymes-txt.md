This artifact is a TAB-delimited file that describes a set of enzymes.

The user can generate this file to define an arbitrary set of enzymes that they want to estimate metabolism on, using the program %(anvi-estimate-metabolism)s.


## Minimal file format

Each row (besides the header) in this file represents one enzyme in the set. At minimum, the file must contain three columns:
- a `gene_id` column containing a _unique_ value to identify a gene for the enzyme. The value can be either a string (like a gene name) or an integer (like a gene callers id), but it has to be unique because sometimes multiple genes can have the same enzyme annotation.
- an `enzyme_accession` column containing the accession of the enzyme, such as a KEGG Ortholog accession for KOfams, a COG accession for NCBI COGs, a Pfam, etc
- a `source` column containing the name of the database that would be used to annotate the enzyme. For example, "KOfam", "COG20_FUNCTION", "Pfam", etc.

{: .notice}
Ideally, all annotation sources in this column would match to those used to define the metabolic pathways you are estimating completeness for (whether those are [KEGG Modules](https://www.genome.jp/kegg/module.html) or user-defined modules as in %(user-modules-data)s, but in practice, we don't currently check for this. If you include some enzymes that are not part of any metabolic modules, they simply will not contribute to the completeness scores of any pathways, and you would therefore only see them in "hits" mode output files.
So the `source` column is (at this time) mostly for you to make sure you know which database these enzymes are coming from and that at least some (hopefully most) will actually be part of the metabolic pathways you are interested, because otherwise the results from %(anvi-estimate-metabolism)s might not make much sense. However, you do have complete freedom to define the 'source' value arbitrarily, if you want. But please keep in mind that this may change in the very near future - one day these `source` values might actually matter for the functioning of %(anvi-estimate-metabolism)s (in which case this documentation will be updated to reflect that). So it is best to get used to setting them properly. :)

## Minimal file example

Here is an example file with the minimum set of columns:

|**gene_id**|**enzyme_accession**|**source**|
|:--|:--|:--|
|aad:TC41_3038|K02886|KOfam|
|aca:ACP_1744|K02626|KOfam|
|aco:Amico_1604|K00606|KOfam|
|ade:Adeh_0623|K02669|KOfam|

## Adding gene coverage and detection values

If you want downstream programs like %(anvi-estimate-metabolism)s to have access to the coverage and detection data for each enzyme (well, technically, its gene), then you can add two additional columns to this file:
- the  `coverage` column should contain the numerical coverage value for the gene encoding the enzyme
- the `detection` column should contain the numerical detection value for the gene encoding the enzyme

If these columns are included, you can use the `--add-coverage` flag with %(anvi-estimate-metabolism)s so that this data is included in the output for each metabolic pathway and/or enzyme. However, you do need to include _both_ of the columns - that program does not currently support adding just coverage or just detection.

## Example with coverage and detection

|**gene_id**|**enzyme_accession**|**source**|**coverage**|**detection**|
|:--|:--|:--|:--|:--|
|aad:TC41_3038|K02886|KOfam|4.44|0.7862318840579711|
|aca:ACP_1744|K02626|KOfam|4.522875816993464|0.7790055248618785|
|aco:Amico_1604|K00606|KOfam|2.63953488372093|0.8063380281690141|
|ade:Adeh_0623|K02669|KOfam|2.011764705882353|0.6639344262295082|
