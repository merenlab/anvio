This program takes in a %(short-reads-fasta)s file and tries to recreate what paired reads for the data in that fasta file might look like. 

An arbitrarily chosen half of the reads will be put into the R1 output, while the other half will be reverse complemented and put into the R2 output. 

For example, if you ran 

{{ codestart }}
anvi-script-gen-pseudo-paired-reads-from-fastq -f %(short-reads-fasta)s \
                                               -O MY_READS 
{{ codestop }}

Then you would end up with two files: 

- `MY_READS_1.fastq` which contains half of the reads straight out of your input file
- `MY_READS_2.fastq` which contains the reverse complement of the other half of the reads. 
