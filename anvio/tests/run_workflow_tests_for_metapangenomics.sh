#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/workflow_test
cp $files/mock_data_for_pangenomics/*.fa $output_dir/workflow_test/
cp $files/mock_data_for_pangenomics/default-state.json        $output_dir/workflow_test/
cp $files/workflows/contigs/fasta.txt $output_dir/workflow_test/
cp $files/workflows/pangenomics/config.json $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "Creating a default config for metapangenomics workflow"
anvi-run-workflow -w metapangenomics --get-default-config default-config.json
exit
INFO "Listing dependencies for pangenomics workflow"
anvi-run-workflow -w pangenomics -c config.json --list-dependencies

INFO "Saving a workflow graph"
anvi-run-workflow -w pangenomics -c config.json --save-workflow-graph

INFO "Running pangenomics workflow"
anvi-run-workflow -w pangenomics -c config.json

INFO "Importing default state for pangenome"
anvi-import-state -p 03_PAN/TEST-PAN.db -s default-state.json -n default

INFO "Vizualize pangenomic results"
anvi-display-pan -g 03_PAN/TEST-GENOMES.db -p 03_PAN/TEST-PAN.db
