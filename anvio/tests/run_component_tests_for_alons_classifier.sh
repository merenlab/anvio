#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Setting the folder where the files are
files=$files/mock_files_for_alons_classifier

INFO "Initializing raw BAM files"
# init raw bam files.
for f in 41 62 74 75 79 94
do
    anvi-init-bam $files/HMP00$f-RAW.bam --output-file $output_dir/HMP00$f.bam
    echo
done

INFO "Generating the TEST contigs database"
anvi-gen-contigs-database -f $files/TEST.fa -o $output_dir/CONTIGS.db $thread_controller

# profiling generates individual directiorues uner $output_dir directory for each sample.
for f in 41 62 74 75 79 94
do
    INFO "Profiling sample SAMPLE-$f"
    anvi-profile -i $output_dir/HMP00$f.bam -o $output_dir/HMP00$f -c $output_dir/CONTIGS.db --skip-SNV-profiling $thread_controller
    echo
done


INFO "Merging profiles"
anvi-merge $output_dir/HMP00*/PROFILE.db -o $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db --skip-concoct-binning

INFO "Importing collection"
anvi-import-collection -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db $files/TEST-COLLECTION.txt -C TEST

# INFO "Run anvi-mcg-classifier on PROFILE database"
# anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-ALL

#INFO "Running anvi-mcg-classifier on TAB-delimited files (no PROFILE database)"
#anvi-mcg-classifier -d $output_dir/TEST-ALL-gene-coverages.txt -D $output_dir/TEST-ALL-gene-detections.txt -O $output_dir/TEST-ALL-TAB-delim

INFO "Running anvi-mcg-classifier on a collection"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-ALL-BINS -C TEST

INFO "Running anvi-mcg-classifier on a bin"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-Bin_1 -C TEST -b Bin_1

INFO "Running anvi-mcg-classifier on a bin with samples to exclude"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-Bin_exclude -C TEST -b Bin_1 --exclude-samples $files/samples_to_exclude.txt

INFO "Running anvi-mcg-classifier on a bin with samples to include"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-Bin_include -C TEST -b Bin_1 --include-samples $files/samples_to_include.txt

INFO "Summarizing collection TEST with --init-gene-coverages"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C TEST -o $output_dir/SAMPLES-MERGED-SUMMARY --init-gene-coverages

INFO "A round of dry run to get the profile db created for Bin_1"
anvi-interactive -p $output_dir/mock_profile_Bin_1.db \
                 -d $output_dir/SAMPLES-MERGED-SUMMARY/bin_by_bin/Bin_1/Bin_1-gene_coverages.txt \
                 --manual \
                 --dry-run

INFO "Importing a default state into newly generated profile database"
anvi-import-state -p $output_dir/mock_profile_Bin_1.db \
                  --state $files/default.json \
                  --name default

INFO "Importing MCG data for Bin 1 into the ad hoc profile database"
anvi-import-misc-data $output_dir/TEST-Bin_1-additional-layers.txt \
                      -p $output_dir/mock_profile_Bin_1.db \
                      -t items
anvi-import-misc-data $output_dir/TEST-Bin_1-samples-information.txt \
                      -p $output_dir/mock_profile_Bin_1.db \
                      -t layers

INFO "Firing up the interactive interface for genes in Bin 1"
anvi-interactive -p $output_dir/mock_profile_Bin_1.db \
                 -d $output_dir/SAMPLES-MERGED-SUMMARY/bin_by_bin/Bin_1/Bin_1-gene_coverages.txt \
                 --title "Alon's gene classifier" \
                 --manual
