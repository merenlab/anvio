#!/bin/bash

# $1: R1.fa
# $2: SAMPLE_NAME
# $3: REF

# build ref:
#       bowtie2-build REF.fa REF

set -e

bowtie2 -x $3 -f $1 -S $2.sam --no-unal
samtools view -F 4 -bS $2.sam > $2-RAW.bam
samtools sort $2-RAW.bam -o $2.bam
samtools index $2.bam
rm $2.sam $2-RAW.bam
