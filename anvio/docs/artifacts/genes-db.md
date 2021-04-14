An anvi'o genes database is a %(profile-db)s-like database that contains statistics, such their coverage and detection across samples, rather than contigs in a given %(contigs-db)s.

A gene database for a given %(bin)s stored in a %(collection)s will be automatically generated when %(anvi-interactive)s is run in 'gene mode'. For details, see the [relevant section](../programs/anvi-interactive/#visualizing-genes-instead-of-contigs) in %(anvi-interactive)s

Alternatively, genes databases can be explicitly generated using the program %(anvi-gen-gene-level-stats-databases)s. By default, this program will generate a gene database for each %(bin)s for a given %(collection)s. 

Due to the strucutral similarities between a %(genes-db)s and a %(profile-db)s, many of the anvi'o programs that operate on profile databases will also run on genes databases. These programs include those that import/export states and import/export misc additional data.
