A tRNA-seq contigs database is a **%(contigs-db)s variant containing information on tRNA transcripts identified from tRNA-seq experiments**.

This database is created by the program, %(anvi-merge-trnaseq)s, which is part of the [trnaseq-workflow](../../workflows/trnaseq/). This program also creates %(trnaseq-profile-db)ss. %(anvi-run-trna-taxonomy)s populates the tRNA-seq contigs database with taxonomic annotations.

This database functions in a manner equivalent to the normal metagenomic-style contigs database. As normal contigs databases are associated with a normal %(profile-db)s containing coverage-related data, tRNA-seq contigs databases are associated with %(trnaseq-profile-db)ss. The name can be misleading: tRNA-seq contigs databases do not contain information on assembled contigs as such. Rather, the fundamental type of sequence reconstructed from a tRNA-seq experiment is a **tRNA seed**, representing a mature tRNA sequence (minus the 3'-CCA acceptor) found in one or more samples in the experiment. tRNA seeds are not predicted by assembly at all, but by the specialized software of %(anvi-trnaseq)s and %(anvi-merge-trnaseq)s.

A variety of information on tRNA seeds is contained in a tRNA-seq contigs database, including structural profiles, taxonomic annotations, and user-defined bins.

## Uses

Tabulation of tRNA-seq data by %(anvi-tabulate-trnaseq)s requires a tRNA-seq contigs database and %(trnaseq-profile-db)s.

Interactive visualization of tRNA-seq datasets in %(anvi-interactive)s requires this database and a %(trnaseq-profile-db)s.

Visualization of grouped seeds by %(anvi-plot-trnaseq)s requires this database in addition to files produced by %(anvi-tabulate-trnaseq)s.
