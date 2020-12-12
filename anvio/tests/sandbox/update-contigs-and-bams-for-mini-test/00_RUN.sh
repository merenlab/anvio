#!/bin/bash

set -e

# if you update this variable, make sure there is also a config file
# associated with it in the same directory
samples="SAMPLE-01 SAMPLE-02 SAMPLE-03"

# restart the output directory
rm -rf output
mkdir output

# generate short reads for each sample
for sample in $samples
do
    ~/github/reads-for-assembly/gen-paired-end-reads samples/$sample.ini
done

# put all 'main' contigs into a single file
cat contigs/*orig.fa > output/contigs.fa

# change work directory
cd output

# generate a bowtie2 reference database
bowtie2-build contigs.fa contigs

# map short reads
for sample in $samples
do
    bowtie2 -x contigs \
            -1 $sample-R1.fastq \
            -2 $sample-R2.fastq \
            -S $sample.sam \
            --threads 4

    samtools view -F 4 -bS $sample.sam -o $sample-RAW.bam

    gzip $sample-R1.fastq
    gzip $sample-R2.fastq
done

# go back to where you came from
cd ../
