This script get the short reads (in the form of a %(short-reads-fasta)s) out of a %(bam-file)s.  

A basic run of this program is as follows: 

{{ codestart }}
anvi-get-short-reads-from-bam -o path/to/output \ 
                              BAM_FILE_1.bam BAM_FILE_2.bam
{{ codestop }}

This will get all of the short reads out of the provided bam files (`BAM_FILE_1.bam` and `BAM_FILE_2.bam`) and put them into a single file. 

### Narrowing the input 

You can choose to only return the short reads that are contained within a %(collection)s or %(bin)s, as so:

{{ codestart }}
anvi-get-short-reads-from-bam -o path/to/output \ 
                              -c %(contigs-db)s \
                              -p %(profile-db)s \
                              -C %(collection)s \
                              BAM_FILE_1.bam BAM_FILE_2.bam
{{ codestop }}

### Changing the output format

You can split the output based on the directionality of paired-end reads. Adding the tag `--split-R1-and-R2` causes the program to create three separate output files: one for R1 (sequences in the forward direction), one for R2 (sequences in the reverse direction; i.e. reverse complement of R1 sequences), and one for unparied reads. When doing this, you can name these three files with a prefix by using the flag `-O`.  

{{ codestart }}
anvi-get-short-reads-from-bam -o path/to/output \ 
                              --split-R1-and-R2 \ 
                              -O BAM_1_and_BAM_2 \
                              BAM_FILE_1.bam BAM_FILE_2.bam
{{ codestop }}

You can also compress the output by adding the flag `--gzip-output`
