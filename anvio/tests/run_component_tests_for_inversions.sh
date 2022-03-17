#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the inversions analysis directory"
mkdir $output_dir/inversions_test
cp -r $files/mock_data_for_inversions/* $output_dir/inversions_test/
cd $output_dir/inversions_test

INFO "Migrating contigs database"
anvi-migrate 02_CONTIGS/CONTIGS.db --migrate-dbs-quickly

INFO "Generating single profile databases"
for sample in S01 S02 S03
do
  anvi-profile -c 02_CONTIGS/CONTIGS.db -i 03_MAPPING/${sample}.bam --fetch-filter inversions -o 04_PROFILE/$sample
done

INFO "Running the base analysis (without inversion's activity and gene context report)"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_BASIC \
                       --skip-compute-inversion-activity \
                       --skip-recovering-genomic-context

INFO "Running the analysis with 1 mismatch allowed"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_MISMATCH \
                       -m 1 \
                       --skip-compute-inversion-activity \
                       --skip-recovering-genomic-context

INFO "Running the analysis with using blast"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_BLAST \
                       -m 1 \
                       --palindrome-search-algorithm BLAST \
                       --skip-compute-inversion-activity \
                       --skip-recovering-genomic-context


INFO "Running the analysis with inversion's activity"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_ACTIVITY \
                       -m 1 \
                       --skip-recovering-genomic-context


INFO "Running the complete analyis"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -m 1 \
                       -o INVERSION_COMPLETE 


INFO "Running the analysis on a target region"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_TARGET \
                       -m 1 \
                       --target-contig c_000000000001 \
                       --target-region-start 5000 \
                       --target-region-end 15000

INFO "Running the analysis and check for all palindromes"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_ALL_PALINDROMES \
                       -m 1 \
                       --check-all-palindromes \
                       --verbose

INFO "Running the analysis and only use inverted-reads"
anvi-report-inversions --my-name-is florian \
                       -P bams-and-profiles.txt \
                       -o INVERSION_INVERTED_READS \
                       -m 1 \
                       --process-only-inverted-reads \
                       --verbose
