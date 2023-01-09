#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################
INFO "Setting up the ecophylo workflow test directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/sra_download/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "Creating a default config for sra_download workflow"
anvi-run-workflow -w sra_download --get-default-config sra_download_config.json

INFO "Listing dependencies for sra_download workflow"
anvi-run-workflow -w sra_download -c sra_download_config.json --list-dependencies

INFO "Saving a workflow graph sra_download workflow"
anvi-run-workflow -w sra_download -c sra_download_config.json --save-workflow-graph

INFO "Running workflow graph sra_download workflow"
anvi-run-workflow -w sra_download -c sra_download_config.json