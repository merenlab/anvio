#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the anvi-export-locus test directory"
cp $files/data/genomes/bacteria/*.db $output_dir

INFO "Running anvi-run-cazymes help menu"
anvi-run-cazymes -h

INFO "Running anvi-run-cazymes"
anvi-run-cazymes -c $output_dir/B_thetaiotamicron_VPI-5482.db

INFO "Running anvi-run-cazymes with --noise-cutoff-terms"
anvi-run-cazymes -c $output_dir/B_thetaiotamicron_VPI-5482.db --noise-cutoff-terms "-E 1e-12"
