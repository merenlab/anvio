By default, anvi'o predicts protein structures (with MODELLER or ColabFold) when creating a %(structure-db)s. Yet, if the user provides an external structures file, then anvi'o does not predict anything, and instead uses this file to import pre-computed structures into the %(structure-db)s. This is handy when you already have structures, or when you predicted them with another tool (e.g. AlphaFold3).

External structures is a user-provided TAB-delimited file. Its header determines the format, and the format must match the input you are building the %(structure-db)s from (a %(contigs-db)s, or a pangenome). Anvi'o auto-detects which one you provided.

## File format

### From a contigs database

A two-column file mapping each gene caller id to a structure file:

|**gene_callers_id**|**path**|
|:--|:--|
|1|path/to/gene1/structure.pdb|
|2|path/to/gene2/structure.pdb|
|(...)|(...)|

where,

* **gene_callers_id** is a gene call in the %(contigs-db)s,
* **path** points to a %(protein-structure-txt)s for that gene.

### From a pangenome

When you build the %(structure-db)s from a pangenome (a %(pan-db)s plus its %(genomes-storage-db)s), a gene is identified by the genome it comes from *and* its gene caller id, so the file has a genome column:

|**genome_name**|**gene_callers_id**|**path**|
|:--|:--|:--|
|G_01|1|path/to/structure_A.pdb|
|G_01|42|path/to/structure_B.pdb|
|G_02|17|path/to/structure_C.pdb|
|(...)|(...)|(...)|

where,

* **genome_name** is a genome in the %(genomes-storage-db)s,
* **gene_callers_id** is a gene call in that genome,
* **path** points to a %(protein-structure-txt)s for that gene.

You do not provide the gene cluster: anvi'o derives it from the %(pan-db)s for each gene. Each gene must therefore belong to a gene cluster in your pangenome.

{:.notice}
Anvi'o will test the integrity of each file, and ultimately you may be subject to the strict requirements anvi'o sets forth (for example, if a structure has a missing residue, you will hear about it). In particular, the amino acid sequence in each structure file must match the sequence anvi'o has for that gene (from the %(contigs-db)s, or from the %(genomes-storage-db)s for a pangenome).
