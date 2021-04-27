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
cat contigs/*orig.fa > output/contigs-orig.fa

# change work directory
cd output

anvi-script-reformat-fasta contigs-orig.fa --simplify-names -o contigs.fa

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

#############################################
# why not go through a test run.
#############################################
anvi-gen-contigs-database -f contigs.fa -o CONTIGS.db -L 1000 --project-name "Mini test"
anvi-run-hmms -c CONTIGS.db --num-threads 4
for sample in $samples
do
    anvi-init-bam $sample-RAW.bam --output-file $sample.bam

    anvi-profile -i $sample.bam \
                 -o $sample \
                 -c CONTIGS.db \
                 --profile-SCVs
done
anvi-merge SAMPLE-*/PROFILE.db -o SAMPLES-MERGED -c CONTIGS.db
mv SAMPLES-MERGED/*db .
anvi-interactive -p PROFILE.db -c CONTIGS.db
#############################################
# Done with the test run
#############################################

# go back to where you came from
cd ../
