An anvi'o genes database is a %(profile-db)s-like database that contains statistics, such their coverage and detection across samples, rather than contigs in a given %(contigs-db)s.

A gene database for a given %(bin)s stored in a %(collection)s will be automatically generated when %(anvi-interactive)s is run in 'gene mode'. For details, see the [relevant section](../../programs/anvi-interactive/#visualizing-genes-instead-of-contigs) in %(anvi-interactive)s

Alternatively, genes databases can be explicitly generated using the program %(anvi-gen-gene-level-stats-databases)s. By default, this program will generate a gene database for each %(bin)s for a given %(collection)s.

In addition to gene-level coverage stats (i.e., mean coverage, detection, and non-outlier mean coverage and standard deviation of each gene in each sample), a genes database by default also stores five normalizations of non-outlier mean coverage values in a table called `gene_level_normalized_coverages`: `log1p`, `cov_per_million` (mean coverage depth per million mapped reads), and per-gene z-scores across samples for raw, log1p, and cov_per_million values (`zscore_raw`, `zscore_log1p`, and `zscore_cov_per_million`). These become additional views in gene mode. See %(anvi-gen-gene-level-stats-databases)s for details, including how to skip computing them, and why that is required for profile databases with a single sample.

Genes databases generated with earlier versions of anvi'o can be upgraded with %(anvi-migrate)s. Their normalized coverages will be computed the next time they are used by %(anvi-interactive)s in gene mode, or by %(anvi-gen-gene-level-stats-databases)s.

Due to the structural similarities between a %(genes-db)s and a %(profile-db)s, many of the anvi'o programs that operate on profile databases will also run on genes databases. These programs include those that import/export states and import/export misc additional data.
