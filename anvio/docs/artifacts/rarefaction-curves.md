In general, this artifact will have visualization similar to this one for a given %(pan-db)s:

![Contents of the contigs and profile databases](../../images/anvi-compute-rarefaction-curves-output.png)

The same artifact can also be computed for a %(pan-graph-db)s, in which case the curves count synteny-aware gene clusters (SynGCs) rather than gene clusters, and the figure is labeled accordingly.

Alongside the figure, this artifact includes four TAB-delimited text files: per-genome-count averages and every individual subsampling observation, for the whole pangenome and for its core. Their columns are named `num_genomes`, `avg_num_gene_clusters` and `standard_deviation` (averages), and `num_genomes` and `GeneClusters` (iterations), regardless of which kind of pangenome they were computed from.
