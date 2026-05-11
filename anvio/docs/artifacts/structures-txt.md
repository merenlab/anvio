A user-provided **TAB-delimited** file that maps gene identifiers to predicted protein structure files. It is the generic anvi'o input artifact for "here are some 3D structures, keyed by something".

The file has exactly two columns:

|gene_id|path|
|:--|:--|
|GC_00000001|path/to/structures/GC_00000001.pdb|
|GC_00000002|path/to/structures/GC_00000002.pdb|
|GC_00000003|path/to/structures/GC_00000003.cif.gz|
|(...)|(...)|

The meaning of `gene_id` is set by whichever anvi'o program is consuming the file:

- For %(anvi-pan-genome)s in structure-informed pangenomics, `gene_id` is a **gene cluster ID** (the FASTA defline of the GC representative emitted by the conventional pangenome).
- For %(anvi-gen-structure-database)s and %(anvi-update-structure-database)s, `gene_id` is a **gene-callers-id** in the source %(contigs-db)s. The legacy column header `gene_callers_id` is also accepted for backward compatibility — see %(external-structures)s for the structure-db-specific notes.

### Format and sanity checks

Anvi'o validates a `structures-txt` file with the following rules:

- The file must be tab-delimited and have exactly two columns.
- The first column header must be `gene_id` (the legacy `gene_callers_id` is still accepted for the structure-db workflow).
- The second column header must be `path`.
- Each `gene_id` must appear at most once.
- Each row must have a non-empty `gene_id` and a non-empty `path`.
- Every `path` must point to a file that exists on disk. Paths can be absolute or relative to the directory of the `structures-txt` file itself.
- Every `path` must end with one of these extensions: `.pdb`, `.cif`, `.mmcif`, `.pdb.gz`, `.cif.gz`, `.mmcif.gz`, `.fcz` (FoldComp).

Programs that consume `structures-txt` may add their own additional checks on top of these.

### A note on file naming

You don't need to name your structure files after the `gene_id`, but it makes life easier when assembling the `structures-txt` file. For instance, if your conventional pangenome produced a representative FASTA, predicting one structure per defline and naming each file `<defline>.pdb` makes the `structures-txt` essentially write itself.
