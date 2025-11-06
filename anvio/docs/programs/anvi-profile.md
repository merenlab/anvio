This program **creates a %(single-profile-db)s from a %(bam-file)s and %(contigs-db)s**.

Once you have a %(single-profile-db)s, you can run programs like %(anvi-cluster-contigs)s, %(anvi-estimate-metabolism)s, and %(anvi-gen-gene-level-stats-databases)s, as well as use the interactive interface with %(anvi-interactive)s. If you want to run these same contigs against multiple BAM files (because you have multiple samples), you'll combine your %(single-profile-db)ss into a %(profile-db)s after you've created them all using %(anvi-merge)s. See the pages for %(single-profile-db)s or %(profile-db)s for more you can do with these artifacts.

In short, this program runs various analyses on the contigs in your %(contigs-db)s and how they relate to the sample information stored in the %(bam-file)s you provided. It then stores this information into a %(single-profile-db)s. Specifically, this program calculates
* coverage per nucleotide position (if you're unsure what coverage refers to, check out [this page](http://merenlab.org/vocabulary/#coverage))
* single-nucleotide, single-codon, and single-amino acid variants (You can find all of those terms on the vocab page linked above, as well as a more detailed explaination [here](http://merenlab.org/2015/07/20/analyzing-variability/#an-intro-to-single-nucleotidecodonamino-acid-variation))
* structural variants such as insertions or deletions

## Basic Usage

### Inputs

This program takes in an [indexed](https://merenlab.org/software/anvio/help/programs/anvi-init-bam) %(bam-file)s and a %(contigs-db)s. The BAM file contains the short reads from a single sample that will be used to create the profile database. Thus, here is a standard run with default parameters:

{{ codestart }}
anvi-profile -i %(bam-file)s \
             -c %(contigs-db)s
{{ codestop }}

Alternatively, if you lack mapping data, you can add the flag `--blank-profile` so that you can still get the functionality of a profile database.

{{ codestart }}
anvi-profile -c %(contigs-db)s  \
            --blank-profile
{{ codestop }}

### Checking your BAM file: always a good idea

If you want to first check your BAM file to see what contigs it contains, just use the flag `--list-contigs` to see a comprehensive list.

### Profiling a subset of contigs

*Note: This describes how to profile a named subset of contigs. To profile a subset of contigs based on their characterists (for example, only contigs of a certain length or that have a certain coverage), see the section below on "contig specifications"*

By default, anvi'o will use every contig in your %(contigs-db)s. However, if you wish to focus specifically on a subset of these contigs, just provide a file that contains only the names of the contigs you want to analyze, one per line, using the tag `--contigs-of-interest`.

For example, you could run

{{ codestart }}
anvi-profile -c Ross_sea_contigs.db  \
             --blank-profile \
             --contigs-of-interest contigs_i_like.txt
{{ codestop }}

Where `contigs_i_like.txt` looks like this:

    SF15-RossSeacontig4922
    SF15-RossSeacontig702

## Analysis Parameters

Changing these will affect the way that your sequences are analyzed.

Keep in mind that if you plan to merge your resulting %(single-profile-db)s with others later in the project, you'll want to keep these parameters consistent.

### Contig Specification

To profile only contigs within a specific length, you can use the flags `--min-contig-length` and `-max-contig-length`. By default, the minimum length for analysis is 1000 and there is no maximum length.

But beyond these flags, you can specify which contigs you would like to profile much more explicitly using the flag `--contigs-of-interest`.

For instance, if you wish to work only with contigs that have more than a certain coverage across your samples, you can first run the program %(anvi-profile-blitz)s on all BAM files, then use the resulting output file %(bam-stats-txt)s to identify contigs of interest based on their coverages across samples, then put their names in a text file, and pass this file to %(anvi-profile)s using the flag `--contigs-of-interest` (the anvi'o profile used to have a flag for this, `--min-mean-coverage`, that allowed users to remove contigs based on their coverage in a given sample, but [we recently removed it](https://github.com/merenlab/anvio/issues/2047) to promote explicit specification of contigs.

### Filter reads

You can also ignore reads in your BAM file with a percent identity to the reference less than some threshold using the flag `--min-percent-identity`.  By default, all reads are used.

For example, the following code will only look at contigs longer than 2000 nts and will ignore BAM file reads with less than 95 percent identity to the reference:

{{ codestart }}
anvi-profile -c Ross_sea_contigs.db  \
            -i bam_file.bam \
            --min-contig-length 2000 \
            --min-percent-identity 95
{{ codestop }}

By default, anvi'o fetches all reads from the bam file. With `--fetch-filter` you can determine which reads from a bam file will be used for profiling. The current filters are:

* `double-forwards`: only paired-end reads with both R1 and R2 with a 'forward' orientation,
* `double-reverses`: only paired-end reads with both R1 and R2 with a 'reverse' orientation,
* `inversions`: only paired-end reads with both R1 and R2 either 'forward' or 'reverse' and a maximum insert size of 2000 nts,
* `single-mapped-reads`: only single mapped reads (mate is unmapped),
* `distant-pairs-1K`: only paired-end reads with a minimum 1000 nts insert size.

For example, the following code only considers 'inversions' reads:

{{ codestart }}
anvi-profile -c Ross_sea_contigs.db \
             -i bam_file.bam \
             --fetch-filter inversions
{{ codestop }}

### Ancient DNA friendly SNV profiling

By default, anvi'o will report variable nucleotides and their allele frequencies from any nucleotide position in a given short read found in the BAM file. Although, there are some applications where the observed variation in short reads will depend on the location of the nucleotide positions in the read. For instance, in ancient DNA sequencing, the start and the end of short reads are often suffer from DNA damage, leading to an increased number of single-nucleotide variants that emerge from the edges of short reads when they are aligned to a reference. The program %(anvi-profile)s has an additional parameter to mitigate this: `--skip-edges`, which can be used like this:

{{ codestart }}
anvi-profile -c %(contigs-db)s \
             -i %(bam-file)s \
             --skip-edges 5
{{ codestop }}


This parameter offers a means to ameliorate the inflation of SNVs due to biases associated with short-read edges by enabling the user to ask anvi'o to ignore a few number of nucleotides from the beginning and the end of short reads as they're being profiled. For instance, a parameter of `5` like shown above, will make sure that the mismatches of a given read to the reference sequence at the first and the last 5 nucleotides will not be reported in the variability table. Please note that the coverage data will not be impacted by the use of this parameter -- only what is reported as variants will be impacted, decreasing the impact of noise in specific applications.

See [https://github.com/merenlab/anvio/pull/2081](https://github.com/merenlab/anvio/pull/2081) for more information.

### Hierarchical Clustering

#### To cluster or not to cluster?

By default, anvi'o will not try to cluster your splits (since it takes quite a bit of runtime) unless you are using the tag `--blank-profile`. If you don't want to run this, use the tag `--skip-hierarchical-clustering`.

If you're planning to later merge this sample with others, it is better to perform clustering while running %(anvi-merge)s than at this stage.

However, if you want to bin this single sample or otherwise want clustering to happen, just use the tag `--cluster-contigs`.

If you do plan to cluster, you can set a custom distance metric or a custom linkage method.

### Variability

Anvi-profile will throw away variability data below certain thresholds to reduce noise. After all, if you have a single C read at a position with a 1000X coverage where all other reads are T, this is probably not a variant position that you want to investigate further. By default, it will not analyze positions with coverage less than 10X, and it will further discard variants based on [this criteria](https://merenlab.org/2015/07/20/analyzing-variability/#de-novo-characterization-and-reporting-of-snvs).

However, you can change the coverage threshold using the  `--min-coverage-for-variability` flag. You can also report every variability position using the flag `--report-variability-full`.

For example, if you wanted to view every variant, you would profile with the following:

{{ codestart }}
anvi-profile -c Ross_sea_contigs.db  \
            -i bam_file.bam \
            --min-coverage-for-variability 1 \
            --report-variability-full
{{ codestop }}

## Other Parameters

You should provide the sample name with the flag `-S` and can provide a description of your project using the `--description` tag followed by a text file. These will help anvi'o name output files and will show up in the anvi'o interfaces down the line.

You can characterize the codon frequencies of genes in your sample at the cost of some runtime. Despite time being money, codon frequency analysis can be helpful downstream. Simply add the tag `--profile-SCVs` and watch the magic happen.

{:.notice}
If you have prior experience with `--profile-SCVs` being slow, you will be surprised how fast it is
since v6.2

Alternatively, you can choose not to store insertion and deletion data or single nucleotide variant data.

If you know the limits of your system, you can also multithread this program. See the program help menu for more information.
