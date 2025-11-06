A tRNA-seq profile database is a **%(profile-db)s variant containing tRNA seed coverage information from one or more samples**.

This database is created by the program, %(anvi-merge-trnaseq)s, which is part of the [trnaseq-workflow](../../workflows/trnaseq/). This program also creates a %(trnaseq-contigs-db)s.

## Specific and nonspecific coverage

The coverage of tRNA seeds by tRNA-seq reads is determined differently than the coverage of contigs by metagenomic reads. Metagenomic contigs are constructed by an assembly tool and reads are assigned to contigs by a mapping tool. Reads that map to multiple contigs are randomly assigned to one contig. tRNA-seq seeds and coverages are found simultaneously by %(anvi-trnaseq)s and %(anvi-merge-trnaseq)s. Two types of coverage are tracked: **specific** coverage of reads unique to seeds and **nonspecific** coverage of reads in multiple seeds. tRNA-seq reads are often short fragments found in numerous tRNAs; random assignment of these reads would distort tRNA abundances and coverage patterns.

Separate tRNA-seq profile databases are produced for specific and nonspecific coverages. A "combined" database containing both sets of data is produced by default for convenience, allowing specific and nonspecific coverages to be compared side-by-side in the %(anvi-interactive)s interface. A "summed" database of specific + nonspecific coverage can optionally be produced. The `--nonspecific-output` option of %(anvi-merge-trnaseq)s controls the production of nonspecific, combined, and summed databases.

## Modifications versus SNVs

The other significant difference between a tRNA-seq profile database and a normal %(profile-db)s is that variable nucleotides are restricted to tRNA modification positions predicted from mutation signatures. Single nucleotide variants are purposefully excluded, though they can be mistaken for modifications, especially with permissive parameterization of %(anvi-merge-trnaseq)s (see that artifact for more information).

## Uses

Tabulation of tRNA-seq data by %(anvi-tabulate-trnaseq)s takes a specific and optionally nonspecific profile database in addition to a %(trnaseq-contigs-db)s.

Interactive visualization of tRNA-seq data in %(anvi-interactive)s requires a specific, nonspecific, combined, or summed profile database in addition to a %(trnaseq-contigs-db)s.
