#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the metagenomics analysis directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/metagenomics/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "unzipping fasta files"
gzip -d three_samples_example/*.fa.gz

INFO "Creating a default config for metagenomics workflow"
anvi-run-snakemake-workflow -w metagenomics --get-default-config default-config.json

INFO "List dependencies for metagenomics workflow"
anvi-run-snakemake-workflow -w metagenomics --config config.json --list-dependencies

INFO "Runnind a dry run"
anvi-run-snakemake-workflow -w metagenomics --config config.json --dry-run
