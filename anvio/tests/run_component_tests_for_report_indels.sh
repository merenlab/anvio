#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Skip gracefully if long-read mapping tools are missing.
if ! command -v minimap2 >/dev/null 2>&1; then
    INFO "minimap2 not installed -- skipping anvi-report-indels component test"
    exit 0
fi
if ! command -v samtools >/dev/null 2>&1; then
    INFO "samtools not installed -- skipping anvi-report-indels component test"
    exit 0
fi

INFO "Setting up the report-indels mock-MGE analysis directory"
mkdir -p $output_dir
cp -r $files/mock_data_for_indels/* $output_dir/
cd $output_dir/

INFO "Building per-sample variant FASTAs (shared backbone + MGE at distinct sites)"
python3 build_variants.py

INFO "Generating PacBio-HiFi reads for each variant (20 kb mean, ~Q40)"
mkdir -p 01_FASTQ
for sample in S01 S02 S03
do
    anvi-script-gen-reads -f 01_FASTA/${sample}.fa \
                          --preset pacbio-hifi \
                          --read-length 20000 \
                          --error-rate 0.0001 \
                          --coverage 50 \
                          --seed 42 \
                          -o 01_FASTQ/${sample}
done

INFO "Building reference contigs-db from S01 (transposase at its original location)"
mkdir -p 02_CONTIGS
anvi-gen-contigs-database -f 01_FASTA/S01.fa \
                          -o 02_CONTIGS/CONTIGS.db \
                          -n MGE_TEST \
                          --quiet \
                          $thread_controller

INFO "Mapping HiFi reads of every variant to the S01 reference with minimap2"
mkdir -p 03_MAPPING
for sample in S01 S02 S03
do
    minimap2 -ax map-hifi 01_FASTA/S01.fa 01_FASTQ/${sample}.fastq 2>/dev/null \
        | samtools sort -o 03_MAPPING/${sample}.bam -
    samtools index 03_MAPPING/${sample}.bam
done

INFO "Profiling each BAM (indels are profiled by default)"
for sample in S01 S02 S03
do
    anvi-profile -c 02_CONTIGS/CONTIGS.db \
                 -i 03_MAPPING/${sample}.bam \
                 -o 04_PROFILE/${sample} \
                 --sample-name ${sample} \
                 --min-contig-length 1000 \
                 --cluster-contigs \
                 --report-variability-full \
                 --quiet \
                 $thread_controller
done

INFO "Merging the per-sample profile databases"
anvi-merge 04_PROFILE/S01/PROFILE.db \
           04_PROFILE/S02/PROFILE.db \
           04_PROFILE/S03/PROFILE.db \
           -c 02_CONTIGS/CONTIGS.db \
           -o 04_PROFILE/MERGED \
           -W \
           --quiet

INFO "Running anvi-report-indels with --min-length 400 (the mock transposase is 435 bp) and --min-samples 1 (each INS site is unique to one sample in this 3-variant mock)"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_DEFAULT \
                   $thread_controller

SHOW_FILE INDELS_DEFAULT/INDELS-VIEW-DATA.txt
SHOW_FILE INDELS_DEFAULT/INDELS-ITEMS-ADDITIONAL.txt
SHOW_FILE INDELS_DEFAULT/INDELS-LOCI-LONG.txt

INFO "Running anvi-report-indels against multiple single profile-dbs via -p (repeated)"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/S01/PROFILE.db \
                   -p 04_PROFILE/S02/PROFILE.db \
                   -p 04_PROFILE/S03/PROFILE.db \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_MULTI_P \
                   $thread_controller

INFO "Running anvi-report-indels with --event-type DEL only"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --event-type DEL \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_DEL_ONLY \
                   $thread_controller

INFO "Running anvi-report-indels with --event-type INS only"
anvi-report-indels -c 02_CONTIGS/CONTIGS.db \
                   -p 04_PROFILE/MERGED/PROFILE.db \
                   --event-type INS \
                   --min-samples 1 \
                   --min-length 400 \
                   -o INDELS_INS_ONLY \
                   $thread_controller
