#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

# Setting the folder where the files are
files=$files/example_files_for_alons_classifier/
INFO "Run anvi-alons-classifier on PROFILE database"
anvi-alons-classifier -p $files/PROFILE.db -c $files/CONTIGS.db -O $output_dir/TEST-ALL --store-gene-detections-and-gene-coverages-tables

INFO "Running anvi-alons-classifier on TAB-delimited files (no PROFILE database)"
anvi-alons-classifier -d $output_dir/TEST-ALL-gene-coverages.txt -D $output_dir/TEST-ALL-gene-detections.txt -O $output_dir/TEST-ALL-TAB-delim

INFO "Generating a samples information database with samples information"
anvi-gen-samples-info-database -D $output_dir/TEST-ALL-samples-information.txt -o $output_dir/TEST-ALL-SAMPLES.db

#INFO "Importing a state file into the merged profile"
#touch $output_dir/TEST-ALL-manual-profile.db
#anvi-import-state -p $output_dir/TEST-ALL-manual-profile.db --state $files/TEST-state.json --name default

INFO "Running anvi-alons-classifier on a collection"
anvi-alons-classifier -p $files/PROFILE.db -c $files/CONTIGS.db -O $output_dir/TEST-ALL-BINS -C TEST

INFO "Running anvi-alons-classifier on a bin"
anvi-alons-classifier -p $files/PROFILE.db -c $files/CONTIGS.db -O $output_dir/TEST-Bin_1 -C TEST -b Bin_1

INFO "Firing up the interactive interface"
# fire up the browser to show how does the merged samples look like.
anvi-interactive -d $output_dir/TEST-ALL-gene-coverages.txt \
                 -A $output_dir/TEST-ALL-additional-layers.txt \
                 --manual \
                 -p $output_dir/TEST-ALL-manual-profile.db \
                 -s $output_dir/TEST-ALL-SAMPLES.db
