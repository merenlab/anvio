#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
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
anvi-gen-contigs-database -f $files/TEST.fa -o $output_dir/CONTIGS.db

# profiling generates individual directiorues uner $output_dir directory for each sample.
for f in 41 62 74 75 79 94
do
    INFO "Profiling sample SAMPLE-$f"
    anvi-profile -i $output_dir/HMP00$f.bam -o $output_dir/HMP00$f -c $output_dir/CONTIGS.db --skip-SNV-profiling
    echo
done


INFO "Merging profiles"
# merge samples
anvi-merge $output_dir/HMP00*/PROFILE.db -o $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db --skip-concoct-binning

INFO "Importing collection"
anvi-import-collection -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db $files/TEST-COLLECTION.txt -C TEST

# INFO "Run anvi-mcg-classifier on PROFILE database"
# anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-ALL

##
#INFO "Running anvi-mcg-classifier on TAB-delimited files (no PROFILE database)"
#anvi-mcg-classifier -d $output_dir/TEST-ALL-gene-coverages.txt -D $output_dir/TEST-ALL-gene-detections.txt -O $output_dir/TEST-ALL-TAB-delim
#
# INFO "Generating a samples information database with samples information"
# anvi-gen-samples-info-database -D $output_dir/TEST-ALL-samples-information.txt -o $output_dir/TEST-ALL-SAMPLES.db
# #
INFO "Running anvi-mcg-classifier on a collection"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-ALL-BINS -C TEST

# #
INFO "Running anvi-mcg-classifier on a bin"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-Bin_1 -C TEST -b Bin_1
# #
INFO "Running anvi-mcg-classifier on a bin with samples to exclude"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-Bin_exclude -C TEST -b Bin_1 --exclude-samples $files/samples_to_exclude.txt

INFO "Running anvi-mcg-classifier on a bin with samples to include"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/TEST-Bin_include -C TEST -b Bin_1 --include-samples $files/samples_to_include.txt
# INFO "A round of dry run to get the profile db created"
# ## a dry-run of the interactive so it creates a profile database
# anvi-interactive -d $output_dir/TEST-ALL-gene-coverages.txt \
#                  -A $output_dir/TEST-ALL-additional-layers.txt \
#                  -p $output_dir/TEST-ALL-manual-profile.db \
#                  -s $output_dir/TEST-ALL-SAMPLES.db \
#                  --manual \
#                  --dry-run
# #
# INFO "Importing a default state into newly generated profile database"
# anvi-import-state -p $output_dir/TEST-ALL-manual-profile.db --state $files/default.json --name default
# #
# INFO "Firing up the interactive interface"
# ## fire up the browser to show how does the merged samples look like.
# anvi-interactive -d $output_dir/TEST-ALL-gene-coverages.txt \
#                  -A $output_dir/TEST-ALL-additional-layers.txt \
#                  -p $output_dir/TEST-ALL-manual-profile.db \
#                  -s $output_dir/TEST-ALL-SAMPLES.db \
#                  --title "Alon's gene classifier" \
#                  --manual
