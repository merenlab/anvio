#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/workflow_test
cp $files/mock_data_for_pangenomics/*.db $output_dir/workflow_test/
cp $files/workflows/contigs/* $output_dir/workflow_test/
cd $output_dir/workflow_test

# all we need a FASTA file here to get things to work, so we
# will start by migrating the anvi'o contigs-db files
anvi-migrate *.db --migrate-quickly --quiet

# then generate some files as expected by the test:
anvi-export-contigs -c E_faecalis_6240.db -o E_faecalis_6240.fa
anvi-export-contigs -c E_faecalis_6255.db -o E_faecalis_6255.fa
anvi-export-contigs -c E_faecalis_6512.db -o E_faecalis_6512.fa

# compressing one of the FASTA files to match the `fasta.txt` coming
# from `$files/workflows/contigs/fasta.txt`
gzip E_faecalis_6240.fa

# and get rid of the contigs-db files, and resume normal operations
rm -rf E_fa*.db

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
