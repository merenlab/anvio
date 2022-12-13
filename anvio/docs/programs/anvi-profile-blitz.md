This program **produces a %(bam-stats-txt)s from one or more %(bam-file)s given a %(contigs-db)s**. It is designed to serve people who only need to process read recruitment data stored in a %(bam-file)s to recover coverage and detection statistics as well as the number of reads mapped reads (along with other statistics) for their genes and/or contigs. It will report what's going on nicely with memory usage information and estimated time of completion:

[![anvi-profile-blitz](../../images/anvi-profile-blitz.png){:.center-img}](../../images/anvi-profile-blitz.png)

There are other programs in anvi'o software ecosystem that are similar to this one:

* %(anvi-profile)s also takes a %(bam-file)s and profiles it. **They both require a %(contigs-db)s**. But while %(anvi-profile)s produces a %(single-profile-db)s for downstream analyses in anvi'o, %(anvi-profile-blitz)s produces text files for downstream analyses by the user (via R, Python, or other solutions). In contrast to %(anvi-profile)s, %(anvi-profile-blitz)s is orders of magnitude faster with similar memory usage.

* %(anvi-script-get-coverage-from-bam)s also takes a %(bam-file)s and profiles it. **They both produce text output files.** But while %(anvi-script-get-coverage-from-bam)s does not require a %(contigs-db)s, %(anvi-profile-blitz)s requires one to work. They will both run very rapidly, %(anvi-script-get-coverage-from-bam)s will work with much smaller amount of memory.

## Output files

For output file formats, please see %(bam-stats-txt)s.

## Running

You can use this program with one or more BAM files to recover minimal or extended statistics for contigs or genes in a %(contigs-db)s.

{:.warning}
Since the program will not be able to ensure the %(contigs-db)s was generated from the same %(contigs-fasta)s that was used for read recruitment that resulted in %(bam-file)ss for analysis, you can make serious mistakes unless you mix up your workflow and start profiling BAM files that have nothing to do with a %(contigs-db)s. If you make a mistake like that, in the best case scenario you will get an empty output file because the program will skip all contigs with non-matching name. In the worst case scenario you will get a file if some names in %(contigs-db)s incorrectly matches to some names in the %(bam-file)s. While this warning may be confusing, you can avoid all these if you use the SAME FASTA FILE both as reference for read recruitment and as input for %(anvi-gen-contigs-database)s.

### Contigs mode, default output

Profile contigs, produce a default output:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   -o OUTPUT.txt
{{ codestop }}

This example is with a single BAM file, but you can also have multiple BAM files as a parameter by using wildcards,

{{ codestart }}
anvi-profile-blitz *.bam \
                   -c %(contigs-db)s \
                   -o OUTPUT.txt
{{ codestop }}

or by providing multiple paths:

{{ codestart }}
anvi-profile-blitz /path/to/SAMPLE-01.bam \
                   /path/to/SAMPLE-02.bam \
                   /another/path/to/SAMPLE-03.bam
                   -c %(contigs-db)s \
                   -o OUTPUT.txt
{{ codestop }}

### Contigs mode, minimal output

Profile contigs, produce a minimal output. This is the fastest option:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --report-minimal \
                   -o OUTPUT.txt
{{ codestop }}

### Genes mode, default output

Profile genes, produce a default output:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --gene-mode \
                   -o OUTPUT.txt
{{ codestop }}

### Genes mode, minimal output

Profile genes, produce a default output:

{{ codestart }}
anvi-profile-blitz %(bam-file)s \
                   -c %(contigs-db)s \
                   --gene-mode \
                   --report-minimal \
                   -o OUTPUT.txt
{{ codestop }}


## Performance

The memory use will be correlated linaerly with the size of the %(contigs-db)s, but once everything is loaded, the memory usage will not increase substantially over time.

With the flag `--report-minimal`, %(anvi-profile-blitz)s profiled on a laptop computer 100,000 contigs that contained 1 billion nts in 6 minutes and used  ~300 Mb memory. This contigs database had 1.5 million genes, and memory usage increased to 1.7 Gb when %(anvi-profile-blitz)s run in `--gene-mode`. The flag `--gene-mode` does not change time complexity dramatically.

Anvi'o has this program because [Emile Faure](https://twitter.com/faureemile) presented us with a challenge: Emile had a ~140 Gb anvi'o %(contigs-db)s that contained nearly 70 million contig sequences from over 200 single-assembled metagenomes, and wanted to learn the coverages of each gene in the contigs database in 200 metagenomes individually. Yet the combination of %(anvi-profile)s and %(anvi-summarize)s jobs would take **more than 40 days** to complete. Since all Emile needed was to learn the coverages from BAM files, we implemented %(anvi-profile-blitz)s to skip the profiling step. The run took **8 hours to compute and report coverage values for 175 million genes in 70 million contigs**, and the memory use remained below 200 Gb.
