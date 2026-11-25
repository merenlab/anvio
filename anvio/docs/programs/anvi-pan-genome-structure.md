This program merges the gene clusters of a conventional %(pan-db)s into **protein structure-informed gene clusters** (PSGCs), based on the similarity of their predicted structures. Gene clusters computed from amino acid sequences can split proteins that share a fold but diverged too much to be recognized by sequence alone; comparing structures brings those back together.

The output is a new %(pan-db)s (named `PROJECT-STRUCTURE-PAN.db` by default) whose items are PSGCs instead of gene clusters. You can use it with every program that works with a %(pan-db)s, such as %(anvi-display-pan)s, %(anvi-summarize)s, or %(anvi-compute-functional-enrichment-in-pan)s. The input %(pan-db)s is never modified.

### The workflow

A structure-informed pangenome is built on top of a conventional one, in three steps:

1. Compute a conventional pangenome with %(anvi-pan-genome)s. Keep the alignments (i.e., don't use `--skip-alignments`), as they are needed to pick the representative of each gene cluster.
2. Predict one structure per gene cluster with %(anvi-gen-structure-database)s, using the pangenome as input and the flag `--select-representative`, which predicts the structure of the representative of each gene cluster. Structure prediction is the slow part, and it is the step you would run on a GPU when using ColabFold.
3. Run this program with the %(pan-db)s, its %(genomes-storage-db)s, and the %(structure-db)s.

{{ codestart }}
anvi-pan-genome -g %(genomes-storage-db)s \
                -n PROJECT

anvi-gen-structure-database --pan-db PROJECT-PAN.db \
                            -g %(genomes-storage-db)s \
                            --select-representative \
                            --engine colabfold \
                            -o PROJECT-STRUCTURE.db

anvi-pan-genome-structure -p PROJECT-PAN.db \
                          -g %(genomes-storage-db)s \
                          --structure-db PROJECT-STRUCTURE.db
{{ codestop }}

Behind the scenes, anvi'o exports one structure per gene cluster, compares all of them with each other with [foldseek](https://github.com/steineggerlab/foldseek), keeps the hits with a TM-score of at least `--min-tm-score` (the minimum of the query and target TM-scores, 0.5 by default), and uses MCL to group gene clusters with similar structures into PSGCs. The MCL inflation parameter for this step is `--structure-mcl-inflation`, which is independent of the one used to compute the gene clusters of the input pangenome.

### Which structure represents a gene cluster

The representative of a gene cluster is the protein the %(structure-db)s holds for it. Anvi'o does not pick it again. Before anything else, anvi'o makes sure the %(structure-db)s goes with your %(pan-db)s: both must come from the same %(genomes-storage-db)s, and every protein with a structure must still belong to the gene cluster it was predicted for. A %(structure-db)s generated from an earlier version of your pangenome will not pass this check.

If the %(structure-db)s holds more than one protein with a structure for a gene cluster (for instance, because it was generated without `--select-representative`), anvi'o uses the one that is the representative of that gene cluster in your pangenome. If none of them is, anvi'o stops rather than guessing.

### Gene clusters without a structure

Structure prediction does not always succeed, and you may decide to predict structures for some gene clusters only. By default, gene clusters without a structure are not a part of the structure-informed pangenome (they remain in the input %(pan-db)s, of course). With the flag `--keep-gcs-without-structure`, each of them becomes its own PSGC instead, so the structure-informed pangenome covers every gene of the input pangenome. In that case, the item additional data key `has_structure` tells you which PSGCs come from structures.

### Structures from elsewhere

If you predicted structures outside of anvi'o, you can provide them with a %(structures-txt)s instead of a %(structure-db)s. In that context, `gene_id` is a gene cluster ID of your %(pan-db)s, and every file must be a `.pdb`. You can get one representative sequence per gene cluster to predict structures for with %(anvi-get-sequences-for-gene-clusters)s and its flag `--representative-sequences`.

{{ codestart }}
anvi-pan-genome-structure -p PROJECT-PAN.db \
                          -g %(genomes-storage-db)s \
                          --structures-txt %(structures-txt)s
{{ codestop }}

### Tracing PSGCs back to gene clusters

The structure-informed %(pan-db)s keeps a copy of the gene clusters of the input pangenome, and records which gene cluster went into which PSGC. If you would like to see, for each PSGC, how the gene clusters it merged are distributed (core, accessory, singleton), you can compare the two pan-dbs with %(anvi-compare-pan)s, or use the flag `--compare-with-input-pan` to do it right away. The comparison adds these layers to the structure-informed %(pan-db)s, and writes a TAB-delimited report (`PROJECT-STRUCTURE-VS-SEQUENCE.txt`) next to it.
