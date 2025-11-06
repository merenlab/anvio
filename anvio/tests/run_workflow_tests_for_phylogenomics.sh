#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/workflow_test
cp $files/workflows/phylogenomics/* $output_dir/workflow_test/
cp -r $files/workflows/pangenomics/five* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "Creating a default config for contigs workflow"
anvi-run-workflow -w phylogenomics --get-default-config default-config.json

INFO "Listing dependencies for phylogenomics workflow"
anvi-run-workflow -w phylogenomics -c config.json --list-dependencies

INFO "Saving a workflow graph"
anvi-run-workflow -w phylogenomics -c config.json --save-workflow-graph

INFO "Running phylogenomics workflow with a dry-run"
anvi-run-workflow -w phylogenomics -c config.json --dry-run

INFO "Running phylogenomics workflow"
anvi-run-workflow -w phylogenomics -c config.json

INFO "Display phylogenomic tree"
anvi-interactive -t 01_PHYLOGENOMICS/TEST-proteins_GAPS_REMOVED.fa.contree -p manual.db --manual
