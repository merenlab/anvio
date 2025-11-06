#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the inversions analysis directory"
mkdir -p $output_dir
cp -r $files/mock_data_for_inversions/* $output_dir/
cd $output_dir/

INFO "Migrating the contigs database (quietly)"
anvi-migrate 02_CONTIGS/CONTIGS.db \
             --migrate-quickly \
             --quiet

INFO "Generating single profile databases with '--fetch-filter inversions' (quietly)"
for sample in S01 S02 S03
do
    anvi-profile -c 02_CONTIGS/CONTIGS.db \
                 -i 03_MAPPING/${sample}.bam \
                 --fetch-filter inversions \
                 -o 04_PROFILE/$sample \
                 --quiet \
                 $thread_controller
done

INFO "Running the base analysis (without reporting the activity or genomic context of inversions)"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_BASIC \
                       --skip-compute-inversion-activity \
                       --skip-recovering-genomic-context \
                       --skip-search-for-motifs \
                       $thread_controller

SHOW_FILE INVERSION_BASIC/INVERSIONS-CONSENSUS.txt
SHOW_FILE INVERSION_BASIC/PER_SAMPLE/INVERSIONS-IN-S01.txt
SHOW_FILE INVERSION_BASIC/ALL-STRETCHES-CONSIDERED.txt

INFO "Running the analysis with 1 mismatch allowed in palindromes"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_MISMATCH \
                       -m 1 \
                       --skip-compute-inversion-activity \
                       --skip-recovering-genomic-context \
                       --skip-search-for-motifs \
                       $thread_controller

SHOW_FILE INVERSION_MISMATCH/INVERSIONS-CONSENSUS.txt

INFO "Running the analysis with using BLAST"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_BLAST \
                       -m 1 \
                       --palindrome-search-algorithm BLAST \
                       --skip-compute-inversion-activity \
                       --skip-recovering-genomic-context \
                       --skip-search-for-motifs \
                       $thread_controller

INFO "Running the complete analysis (with --verbose)"
anvi-report-inversions -P bams-and-profiles.txt \
                       -m 1 \
                       -o INVERSION_COMPLETE \
                       --verbose \
                       $thread_controller

SHOW_FILE INVERSION_COMPLETE/INVERSION-ACTIVITY.txt
SHOW_FILE INVERSION_COMPLETE/PER_INV/INV_0001/SURROUNDING-FUNCTIONS.txt
SHOW_FILE INVERSION_COMPLETE/PER_INV/INV_0001/SURROUNDING-GENES.txt

INFO "Re-computing inversion activity using previous results (quietly)"
anvi-report-inversions -P bams-and-profiles.txt \
                       --pre-computed-inversions INVERSION_COMPLETE/INVERSIONS-CONSENSUS.txt \
                       --quiet \
                       $thread_controller

INFO "Running the analysis on a target region (quietly)"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_TARGET \
                       -m 1 \
                       --target-contig c_000000000001 \
                       --target-region-start 5000 \
                       --target-region-end 15000 \
                       --quiet \
                       $thread_controller

INFO "Running the analysis with 'check all palindromes' directive (quietly)"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_ALL_PALINDROMES \
                       -m 1 \
                       --check-all-palindromes \
                       --quiet \
                       $thread_controller

INFO "Running the analysis with 'process only inverted reads' directive (quietly)"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_INVERTED_READS \
                       -m 1 \
                       --process-only-inverted-reads \
                       --quiet \
                       $thread_controller

INFO "Running the analysis and searching for only one motif (quietly)"
anvi-report-inversions -P bams-and-profiles.txt \
                       -o INVERSION_ONE_MOTIF \
                       -m 1 \
                       --num-of-motifs 1 \
                       --quiet \
                       $thread_controller
