By default, anvi'o predicts protein structures using MODELLER when creating a %(structure-db)s. Yet, if the user provides an external structures file, then anvi'o does not perform template-based homology modelling, and instead uses this file to obtain the structure information for the %(structure-db)s.

External structures is a user-provided TAB-delimited file that should follow this format:

|gene_callers_id|path|
|:---:|:---|
|1|path/to/gene1/structure.pdb|
|2|path/to/gene2/structure.pdb|
|3|path/to/gene3/structure.pdb|
|4|path/to/gene4/structure.pdb|
|7|path/to/gene5/structure.pdb|
|8|path/to/gene6/structure.pdb|
|(...)|(...)|

Please note that anvi'o will try its best to test the integrity of each file, and work with any limitations, however ultimately the user may be subject to the strict requirements set forth by anvi'o. For example, if a structure has a missing residue, you will hear about it.

