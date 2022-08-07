This is a 3-column, tab-delimited file that lists the enzymes to be included in a user-defined metabolic module. It is used as input to %(anvi-script-gen-user-module-file)s, which will create a module file using the enzyme definitions within the file.

The first column in this file must be an enzyme accession, the second column must be the annotation source of the enzyme, and the third column specifies the orthology (or functional definition) of the enzyme. Note that if the annotation source is 'KOfam' and the enzyme is a KEGG Ortholog that is present in the KEGG KOfam profiles in %(kegg-data)s, then the orthology field can be blank. In this case, the orthology field will be filled in automatically with the enzyme's known orthology in the KOfam data.

Here is an example file:

|**enzyme**|**source**|**orthology**|
|:--|:--|:--|
|K01657|KOfam||
|K01658|KOfam||
|PF06603.14|METABOLISM_HMM|UpxZ|
|COG1362|COG20_FUNCTION|Aspartyl aminopeptidase|
|TIGR01709.2|TIGRFAM|type II secretion system protein GspL|
