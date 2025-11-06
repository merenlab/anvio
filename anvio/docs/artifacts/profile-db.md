An anvi'o database that **contains key information about the mapping of short reads *from multiple samples* to your contigs.**

You can think of this as a extension of a %(contigs-db)s that contains information about how your contigs align with each of your samples. The vast majority of programs that use a profile database will also ask for the contigs database associated with it.

A profile database contains information about how short reads map to the contigs in a %(contigs-db)s. Specifically, for each sample, a profile database contains
* the coverage and abundance per nucleotide position for each contig
* variants of various kinds (single-nucleotide, single-codon, and single-amino acid)
* structural variants (ex. insertions and deletions)
These terms are explained on the [anvi'o vocabulary page.](http://merenlab.org/vocabulary/)

![Contents of the contigs and profile databases](../../images/contigs-profile-db.png)

This information is necessary to run anvi'o programs like %(anvi-cluster-contigs)s, %(anvi-estimate-metabolism)s, and %(anvi-gen-gene-level-stats-databases)s. You can also interact with a profile database using programs like %(anvi-interactive)s.

Technically, "profile-db" refers to a profile database that contains the data from several samples -- in other words, the result of running %(anvi-merge)s on several %(single-profile-db)s. However, since a %(single-profile-db)s has a lot of the functionality of a profile-db, it might be easier to think of a profile database as a header referring to both single-profile-dbs and profile-dbs (which can also be called a merged-profile-dbs). For simplicity's sake, since most users are dealing with multiple samples, the name was shortened to just profile-db. The following are a list of differences in functionality between a single profile database and a merged profile database:
* You can run %(anvi-cluster-contigs)s on only a merged profile database (or profile-db), since they look at the allignment data in many samples
* You cannot run %(anvi-merge)s or %(anvi-import-taxonomy-for-layers)s on a merged profile database, only on a %(single-profile-db)s.

## How to make a profile database

### If you have multiple samples
1. Prepare your %(contigs-db)s
2. Run %(anvi-profile)s with an appropriate %(bam-file)s. The output of this will give you a %(single-profile-db)s. You will need to do this for each of your samples, which have been converted into a %(bam-file)s with your short reads.
3. Run %(anvi-merge)s on your %(contigs-db)s (from step 1) and your %(single-profile-db)ss (from step 2). The output of this is a profile-db.

### If you have a single sample
1. Prepare your %(contigs-db)s
2. Run %(anvi-profile)s with an appropriate %(bam-file)s. The output of this will give you a %(single-profile-db)s. You can see that page for more information, but essentially you can use a single-profile-db instead of a profile database to run most anvi'o functions.

## Variants

Profile databases, like %(contigs-db)ss, are allowed to have different variants, though the only currently implemented variant, the %(trnaseq-profile-db)s, is for tRNA transcripts from tRNA-seq experiments. The default variant stored for "standard" profile databases is `unknown`. Variants should indicate that substantially different information is stored in the database. For instance, single codon variability is applicable to protein-coding genes but not tRNA transcripts, so SCV data is not recorded for the `trnaseq` variant. The $(trnaseq-workflow)s generates %(trnaseq-profile-db)ss using a very different approach to %(anvi-profile)s.
