#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/workflow_test
cp $files/mock_data_for_pangenomics/*.fa $output_dir/workflow_test/
cp ../anvio/workflows/contigs/Snakefile  $output_dir/workflow_test/
cd $output_dir/workflow_test

