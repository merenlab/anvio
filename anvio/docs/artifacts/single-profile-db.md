An anvi'o database that contains the same information as a merged %(profile-db)s, namely **key information about the mapping of short reads *in a single sample* to your contigs.** 

You can think of this as a extension of a %(contigs-db)s that contains information about how your contigs align with a single one of your individual samples. If you have more than one sample, you'll probably want to use %(anvi-merge)s to merge your databases into a merged %(profile-db)s. The vast majority of programs that use a profile database will also ask for the contigs database associated with it. 

A single profile database contains information about how the short reads in a single BAM-file (see %(bam-file)s) map to the contigs in a %(contigs-db)s. Specificially, a profile database contains 
* the coverage and abundance per nucleotide position for each contig 
* variants of various kinds (single-nucleotide, single-codon, and single-amino acid)
* structural variants (ex insertions and deletions)

Once created, a single profile database is almost interchangable with a %(profile-db)s (even though the names can be a little confusing. Think of a single-profile-db as a type of profile-db, since it has only a few differences). The main differences between the two are as follows: 
* You cannot run %(anvi-cluster-contigs)s or %(anvi-mcg-classifier)s on a single profile db, since these two programs look at the alignment data in many samples. 
* You can run %(anvi-import-taxonomy-for-layers)s on a single profile database but not a merged one. 
* You can only run %(anvi-merge)s on a single profile database.

If you want to look at the contents of a single profile database, you can do so using %(anvi-interactive)s. 
