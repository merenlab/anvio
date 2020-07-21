Similarly to a %(genes-fasta)s, a short reads fasta is what it sounds like: a %(fasta)s file containing short reads. 

Short reads usually refer to the initial pieces of sequencing data that you had before you assembled them into longer contigs. In other words, these are the kinds of reads you could get out of a technique like Sanger sequencing. Knowing how those short reads align to your contigs is vital for analysis (In fact, that's a lot of the functionality of a %(profile-db)s!). 

In anvi'o, you can get short reads out of two sources: either from a %(bam-file)s by running the program %(anvi-get-short-reads-from-bam)s or from a %(contigs-db)s by running the program %(anvi-get-short-reads-mapping-to-a-gene)s.  
