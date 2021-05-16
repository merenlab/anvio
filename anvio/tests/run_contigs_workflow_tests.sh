#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/workflow_test
cp $files/mock_data_for_pangenomics/*.fa $output_dir/workflow_test/
cp $files/workflows/contigs/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "compressing sample 1"
gzip 01.fa

INFO "Creating a default config for contigs workflow"
anvi-run-workflow -w contigs --get-default-config default-config.json

INFO "Listing dependencies for contigs workflow"
anvi-run-workflow -w contigs -c default-config.json --list-dependencies

INFO "Saving a workflow graph"
anvi-run-workflow -w contigs -c default-config.json --save-workflow-graph

INFO "Running contigs workflow with a dry-run"
anvi-run-workflow -w contigs -c default-config.json --dry-run

INFO "Running contigs workflow"
anvi-run-workflow -w contigs -c default-config.json

INFO "Examine contigs databases with anvi-display-contigs-stats"
anvi-display-contigs-stats 02_CONTIGS/*db
