This program **generates tabular files of tRNA-seq seed coverage and modification data that are easily manipulable by the user**.

anvi-tabulate-trnaseq is part of the [trnaseq-workflow](../../workflows/trnaseq/), and is run following the finalization of tRNA seeds by %(anvi-merge-trnaseq)s.

This program generates a table, %(seeds-specific-txt)s, containing the specific coverage of each nucleotide position in each seed in every sample. If a nonspecific %(trnaseq-profile-db)s is also provided, this program generates a table of nonspecific coverages, %(seeds-non-specific-txt)s. The distinction between specific and nonspecific coverage is explained in the %(trnaseq-profile-db)s artifact. These coverage tables have one row per seed per sample. They have three header rows for different ways of describing tRNA nucleotide positions: canonical position name (e.g., "discriminator_1"), canonical position (e.g., "73"), and "ordinal" position relative to all the other **possible** positions (e.g., "95").

anvi-tabulate-trnaseq also generates a table, %(modifications-txt)s, containing information on each predicted modification position in each seed, with one row per modification per seed per sample. This table includes four columns of position coverage counts of the four nucleotides.

All tables include taxonomic annotations of the seeds; annotations are added to the %(trnaseq-contigs-db)s by %(anvi-run-trna-taxonomy)s.
