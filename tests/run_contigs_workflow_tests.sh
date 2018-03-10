#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/workflow_test
cp $files/mock_data_for_pangenomics/*.fa $output_dir/workflow_test/
cp $files/workflows/contigs/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "Creating a default config for contigs workflow"
anvi-run-snakemake-workflow -w contigs --get-default-config default-config.json

INFO "Listing dependencies for contigs workflow"
anvi-run-snakemake-workflow -w contigs -c default-config.json --list-dependencies

INFO "Running contigs workflow"
anvi-run-snakemake-workflow -w contigs -c default-config.json
