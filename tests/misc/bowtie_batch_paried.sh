#!/bin/bash

# $1: R1.fastq
# $2: R2.fastq
# $3: SAMPLE_NAME
# $4: REF

# build ref:
#       bowtie2-build REF.fa REF

set -e

bowtie2 -x $4 -1 $1 -2 $2 -S $3.sam
samtools view -F 4 -bS $3.sam > $3-RAW.bam
samtools sort $3-RAW.bam $3
samtools index $3.bam
rm $3.sam $3-RAW.bam
