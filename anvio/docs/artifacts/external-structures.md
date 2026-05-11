By default, anvi'o predicts protein structures using MODELLER when creating a %(structure-db)s. Yet, if the user provides an external structures file, then anvi'o does not perform template-based homology modelling, and instead uses this file to obtain the structure information for the %(structure-db)s.

External structures is a specialization of the generic %(structures-txt)s artifact: a TAB-delimited file with two columns mapping a gene identifier to a path. For the structure-db workflow the identifier is a `gene_callers_id` from the source %(contigs-db)s:

|gene_callers_id|path|
|:---:|:---|
|1|path/to/gene1/structure.pdb|
|2|path/to/gene2/structure.pdb|
|3|path/to/gene3/structure.pdb|
|4|path/to/gene4/structure.pdb|
|7|path/to/gene5/structure.pdb|
|8|path/to/gene6/structure.pdb|
|(...)|(...)|

The newer canonical header `gene_id` is also accepted (it is the column name used by the generic %(structures-txt)s artifact). The two are interchangeable here.

Each path should point to a %(protein-structure-txt)s.

In addition to the generic %(structures-txt)s sanity checks (file existence, supported extension, no duplicate IDs), this workflow enforces a few extra constraints:

- Only `.pdb` files are accepted. The broader extension list (CIF, mmCIF, gzipped variants, FoldComp) advertised by the generic %(structures-txt)s artifact applies to the structure-informed pangenomics path only — the structure-db workflow runs each file through Biopython's PDB parser, which understands PDB and nothing else.
- Every `gene_callers_id` must parse as an integer.
- Every `gene_callers_id` must exist in the source %(contigs-db)s and be coding (have an amino-acid sequence).
- By default, anvi'o opens each structure file and confirms that its sequence matches the contigs-db record for the corresponding gene. (You can skip this with `--lazy`, where supported.)

{:.notice}
Please note that anvi'o will try its best to test the integrity of each file, and work with any limitations, however ultimately the user may be subject to the strict requirements set forth by anvi'o. For example, if a structure has a missing residue, you will hear about it.
