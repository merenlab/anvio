This is a **%(bam-file)s (which contains aligned sequence data) that has not yet been indexed and sorted**. 

### What does being "indexed" mean? 

Think of your BAM file as a long, complex book. In order to get the most out of it when trying to perform analysis, it will be super helpful to have a table of contents. Indexing your BAM file basically creates a second file that serves as an external table of contents, so that anvi'o doesn't have to keep looking through the entire BAM file during analysis. 

You can tell whether or not your BAM file is indexed based on the presence of this second file, which will have the same title as your BAM file, but end with the extension `.bai`. For example, if your directory contained these files:

{{ codestart }}
Lake_Michigan_Sample_1.bam
Lake_Michigan_Sample_1.bam.bai
Lake_Michigan_Sample_2.bam 
{{ codestop }}

then you would still need to index `Lake_Michigan_Sample_2.bam`. 

### How do you index a BAM file?

You can either do this directly using samtools, or you can just run the anvi'o program %(anvi-init-bam)s (which uses samtools for you). 
