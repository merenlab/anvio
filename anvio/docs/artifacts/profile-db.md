An anvi'o database that **contains key information about the mapping of short reads to your contigs.**

You can think of this as a extension of a %(contigs-db)s that contains other information. The vast majority of programs that use the profile database will also ask for the contigs database associated with it. 

A profile database contains information about how short reads map to the contigs in a %(contigs-db)s. Specificially, a profile database contains 
* the coverage and abundance per nucleotide position for each contig 
* varience of various kinds (single-nucleotide, single-codon, and single-amino acid)
* structural variance (ex insertions and deletions)

This information is neccessary to run anvi'o programs like %(anvi-cluster-contigs)s, %(anvi-estimate-metabolism)s, and %(anvi-gen-gene-level-stats-databases)s. You can also interact with a profile database using programs like %(anvi-interactive)s and %(anvi-display-structure)s.

### How to make a profile database
1. Prepare your %(contigs-db)s
2. Run $(anvi-profile)s with an appropriate %(bam-file)s. The output of this will give you a %(single-profile-db)s. You may need to do this step multiple times if you have many BAM files of short reads.
3. Run %(anvi-merge)s on your %(contigs-db)s (from step 1) and your %(single-profile-db)ss (from step 2) 
