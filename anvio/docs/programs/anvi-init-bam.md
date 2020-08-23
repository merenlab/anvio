This program sorts and indexes your BAM files, essentially converting a %(raw-bam-file)s into a %(bam-file)s, which are ready to be used in anvi'o. 

If you're unsure what a BAM file is, check out the %(bam-file)s page or [this file](https://samtools.github.io/hts-specs/SAMv1.pdf), written by the developers of samtools. For a description of what indexing a BAM file does, check out the page for %(raw-bam-file)s. 

To run this program, just provide a path to the bam files that you want to index. For example, 

{{ codestart }}
anvi-init-bam %(raw-bam-file)s 
{{ codestop }}

You can also multithread this to shorten runtime with the flag `-T` followed by the desired number of threads if your system is capable of this. 

To see it in action (plus a description on how to run it on an entire folder), check out [this page](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-init-bam). 
