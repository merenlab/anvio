#!/bin/bash
source 00.sh
set -e

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the pan analysis directory"
cp $files/mock_data_for_pangenomics/*.db $output_dir

INFO "Generating an external genomes file for DBs"
anvi-script-gen-genomes-file --input-dir $output_dir \
                             -o $output_dir/external-genomes.txt

INFO "Running genome view"
anvi-display-genomes -e $output_dir/external-genomes.txt \
                     -V $output_dir/genome-view.db
