This program allows you to find genomic inversions using short-read mapping information.

An inversion is typically carried out by an invertase. This enzyme recognizes a pair of inverted repeat (IR), which are a special case of palindromic sequence where the repeats are facing inward on different DNA strand. The IRs are distant from each other and the invertase will invert the DNA fragment between the IRs.

In brief, anvi'o leverages paired-read orientation to locate regions of interest in a set of contigs. It screens for IRs whithin these regions, and uses short-reads to confirm which IRs corrrespond to real inversions. Eventually, anvi'o can compute the inversion activity: the relative proportion of an inversion's orientation in each sample.

### Anvi'o's philosophy to find inversions

Much like a T-Rex, anvi'o's vision rely on movement and it cannot see an inversion if it does not move, or in this case, invert. So let's start with what you cannot do with this command: you cannot find inversions in a set of contigs alone.

To find an inversion, you need to have short-reads from at least one sample. If there is even a small fraction of a microbial population that have an inverted sequence compare to your contigs of reference, then anvi-report-inversions is for you!

### Prerequistes to run this program

Anvi'o is able to locate inversion using the paired-end read orientation. Regular paired-end reads are facing inward with a FWD/REV orientation, but when an inversion happens, some reads will be mapping in the opposite orientation regarding the reference. As a consequence, some paired-end reads will have the same orientation: FWD/FWD or REV/REV.

To leverage that information, anvi'o can profile bam files for FWD/FWD and REV/REV reads only with %(anvi-profile)s to make special %(single-profile-db)s.

{{ codestart }}
anvi-profile -i %(bam-file)s \
             -c %(contigs-db)s \
             --fetch-filter inversion
{{ codestop }}

### Inputs

The main input for this command is a %(bams-and-profiles-txt)s, which is a TAB-delimited file composed of at least four columns:

* Sample name,
* %(contigs-db)s,
* %(single-profile-db)s generated with the inversion fetch filter,
* %(bam-file)s.

You can also add two column for the R1 and R2 fastq files so that anvi'o can compute the inversion's activity.

Here is a standard run with default parameters:

{{ codestart }}
anvi-report-inversion -P %(bams-and-profiles-txt)s \
                      -o inversions_output
{{ codestop }}

### Identifying regions of interest

While anvi'o could directly search of inverted repeats in all the contigs, it would be a waste of time as many IRs are actually not related to inversions. Instead, anvi'o uses the FWD/FWD and REV/REV reads to identify region of interest and constrain the seach for IRs only in these regions.

For this step, you can set the minimum coverage of FWD/FWD and REV/REV reads to define 'stretches' with `--min-coverage-to-define-stretches`. Lower threshold yield more stretches, but also more noise.

The parameter `--min-stretch-length` defines the minimun length for a stretch to be considered.

FWD/FWD reads are found on the most left side of an inversion, while the REV/REV reads are found on the right side. When an inversion is quite long, a region of low to no coverage can separate the two group of reads. With the flag `--min-distance-between-independent-stretches`, anvi'o merges fragmented stretches if they are closer than this value.

When the coverage is quite low, a stretch can be a little bit too short and miss potential IRs, so you use `--num-nts-to-pad-a-stretch` to extend the stretch by x bp upstream and downstream.

Here are the default values for these flags:

{{ codestart }}
anvi-report-inversions -P %(bams-and-profiles-txt)s \
                       -o inversions_output \
                       --min-coverage-to-define-stretches 10 \
                       --min-stretch-length 50 \
                       --min-distance-between-independent-stretches 2000 \
                       --num-nts-to-pad-a-stretch 100
{{ codestop }}

### Finding palindromes

There are a few paramters to constrain the search for palindromic sequences of the IRs, like a minimum length that can be set with `--min-palindrome-length`, and a maximum number of mismatches with `--max-num-mismatches`.

You can set the minimum distance between two palindromic sequence with `--min-distance`. A distance of 0 would correspond to a in-place palindrome, though they don't relate to genomic inversions.

When searching for palindromes with mismatches, the algorithm will extend the palindrom length as much as possible, often including mismatches which are outside of the true palindrome sequences. The flag `--min-mismatch-distance-to-first-base` allows you to trim the palindrome when one or more mismatches are n nucleotides away from a palindrome's start or stop. The default value is 1, meaning that a palindrome `MMMMMM(X)M`, where M denotes matching nucleotides and X a mismatch, will be trimmed to the first 6 matches `MMMMMM`.

There are currently two algorithms to find palindromes in anvi'o: numba and BLAST. Numba is very fast when looking for palindromes in short sequences, and BLAST is more efficient for longer stretches. Anvi'o dynamically set the algorithm accoding to each stretch length: numba for stretches under 5,000 bp and BLAST for longer stretches. You can use the flag `--palindrome-search-algorithm` to ask anvi'o to use either of these methods explicitly. Note that results between the two methods may differs.

{{ codestart }}
anvi-report-inversions -P %(bams-and-profiles-txt)s \
                       -o inversions_output \
                       --min-palindrome-length 10 \
                       --max-num-mismatches 1 \
                       --min-distance 50 \
                       --min-mismatch-distance-to-first-base 1 \
                       --palindrome-search-algorithm numba
{{ codestop }}

### Confirming inversions

Multiple palindromes are usualy reported for each stretch and to confirm which one actually relates to an inversions, anvi'o searches short-reads in the bam file for unique constructs that can only occur when a genomic region inverted.

By default, anvi'o reports the first confirmed palindrome and move to the next stretch. This process is very efficient as a strech usually have only one inversion. But in rare cases, you can have multiple inversions happening in a single stretch. Then, you can use the flag `--check-all-palindromes` and anvi'o will look for inversion evidences in the short-reads for every palindrome in a stretch.

Anvi'o looks for inversion evidence in the FWD/FWD and REV/REV reads first. If no evidence are found, then it searches the rest of the reads mapping to the region of interest. If you want to only search the FWD/FWD and REV/REV reads you can use the flag `--process-only-inverted-reads`

### Computing inversion activity

If you provide the short-reads R1 and R2 in the %(bams-and-profiles-txt)s, anvi'o can compute the proportion of the inversion's orientation in each sample.

This is a very time consuming step, and if you have multiple sample, you can use the parameter `--num-threads` to set the maximum of threads for multithreading when possible.

To compute the inversion's ratios, anvi'o design in silco primers based on the palidrome sequence and the upstream/downstream genomic context to search short-reads in the raw fastq files. The variable `--oligo-primer-base-length` is used to control how much of the palindrome should be used to design the primers. The longer, the more specific but if it is too long, less reads will match to the primer.

This step is very computationally intense, but you can test it with the parameter `--end-primer-search-after-x-hits`. Once the total number of reads reach this parameter, anvi'o will stop searching further and will continue with the next sample. This flag is only good for testing.

If you want to skip this step, you can use the flag `--skip-compute-inversion-activity`.

An example command with 12 threads:

{{ codestart }}
anvi-report-inversions -P %(bams-and-profiles-txt)s \
                       -o inversions_output \
                       --num-threads 12 \
                       --oligo-primer-base-length 12
{{ codestop }}

### Reporting genomic context around inversions

For every inversion, anvi'o can report the surrounding genes and their function as additional files.

You can use the flag `--num-genes-to-consider-in-context` to choose how many genes to consider upstream/downstream of the inversion. By default, anvi'o report three genes downstream, and three genes upstream.

To select a specific gene caller, you can use `--gene-caller`. The default is prodigal.

If you want to skip this step, you can use the flag `--skip-recovering-genomic-context`.

### Targeted search

If you are interested in a given contig region you can use the following flags to limit the search:

* `--target-contig`: contig of interest,
* `--target-region-start`: the start position of the region of interest,
* `--target-region-end`: the end position of the region of interest.

### Output
%(anvi-report-inversions)s searches for inversions in every single sample at a time and thus genereates a TAB-delimited table for every sample: `INVERSIONS-IN-SAMPLE_01.txt`, `INVERSIONS-IN-SAMPLE_02`, ...

These tables contains the following information:

* entry ID,
* contig name,
* first palindrome sequence,
* aligment midline,
* second palindrome sequence,
* start and stop position of the first and second palindrome sequence,
* number of mismatches,
* number of gaps,
* length of the palindrome sequence,
* distance between the first and second palindrome seqeuences, i.e. the size of the inversion,
* the number of samples in which it was detected and confirmed,
* the in silico primers used to compute the inversion's activity, for the first and second palindrome,
* the oligo corresponding to the reference sequence.

Anvi'o eventually create a consensus table with all the unique inversions found accross all your samples in a file called `INVERSIONS-CONSENSUS.txt`. This table has the same format as the individual sample outputs, with the 'entry ID' replaced by a unique inversion ID.

Another default output table is named `ALL-STRETCHES-CONSIDERED.txt` and it reports every stretch that passed the ['Identifying regions of interest'](#identifying-regions-of-interest) parameters. It reports the maximum coverage of FWD/FWD and REV/REV in that stretch, per sample. It also reports the number of palindromes found and if a true inversion was confirmed.

If the user enable the reporting of the genomic context, two addition TAB-delimited tables are generated: `INVERSIONS-CONSENSUS-SURROUNDING-GENES.txt` and `INVERSIONS-CONSENSUS-SURROUNDING-FUNCTIONS.txt`.
The first table report the gene calls surrounging every inversion when possible (inversions_id, gene_caller_id, start and stop position, orientation, gene_caller and contig).
The second table report the function associated to every gene call reported in the first file (inversions_id, gene_caller_id, source, accession, function).

Finally, if the user provide R1 and R2 fastq files and enable the reporting of inversion's activity, %(anvi-report-inversions)s will generate a long-format file named `INVERSION-ACTIVITY.txt`. This file reports, for every inversion and sample, the relative proportion and read abundance of unique oligos, which either correspond to the reference contig (no inversion), or to an inversion sequence. The inversion's activity is computed and reported for both side of each inversion.





