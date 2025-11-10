#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the metagenomics analysis directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/metagenomics/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "Creating a default config for metagenomics workflow"
anvi-run-workflow -w metagenomics --get-default-config default-config.json

INFO "Generate a graph from the megahit workflow"
anvi-run-workflow -w metagenomics --config config-megahit.json --save-workflow-graph

INFO "Run with megahit - groups - QC - import Kraken"
anvi-run-workflow -w metagenomics --config config-megahit.json

INFO "Run with metaspades - no groups - no QC - no gzip - use scaffolds"
anvi-run-workflow -w metagenomics --config config-metaspades.json

INFO "Run with idba_ud - all-vs-all - remove SR matching ref - QC"
anvi-run-workflow -w metagenomics --config config-idba_ud.json

INFO "Run with metaflye - LR only"
anvi-run-workflow -w metagenomics --config config-metaflye.json

INFO "Run with metaflye and megahit - all-vs-all - SR QC"
anvi-run-workflow -w metagenomics --config config-mixed-assembly.json

INFO "Run reference-mode with SR and LR - no QC - add collection"
anvi-run-workflow -w metagenomics --config config-references-mode.json

INFO "Run reference-mode with one sample - no QC"
anvi-run-workflow -w metagenomics --config config-references-mode-one-sample.json
