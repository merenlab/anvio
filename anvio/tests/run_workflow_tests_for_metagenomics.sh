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

INFO "List dependencies for metagenomics workflow with megahit"
anvi-run-workflow -w metagenomics --config config-megahit.json --list-dependencies

INFO "Running a dry run in references mode with no qc and no gzip"
anvi-run-workflow -w metagenomics --config config-references-mode-no-qc-no-gzip-no-groups.json
anvi-run-workflow -w metagenomics --config config-references-mode-no-qc-no-gzip-no-groups.json --save-workflow-graph

INFO "Running a dry run for metagenomics workflow with megahit"
anvi-run-workflow -w metagenomics --config config-megahit.json
anvi-run-workflow -w metagenomics --config config-megahit.json --save-workflow-graph

INFO "Running a dry run for metagenomics workflow with megahit with no qc"
anvi-run-workflow -w metagenomics --config config-megahit-no-qc-all-against-all.json
anvi-run-workflow -w metagenomics --config config-megahit-no-qc-all-against-all.json --save-workflow-graph

INFO "List dependencies for metagenomics workflow in references mode"
anvi-run-workflow -w metagenomics --config config-references-mode.json --list-dependencies

INFO "Running a dry run in references mode"
anvi-run-workflow -w metagenomics --config config-references-mode.json
anvi-run-workflow -w metagenomics --config config-references-mode.json --save-workflow-graph

INFO "List dependencies for metagenomics workflow with idba_ud"
anvi-run-workflow -w metagenomics --config config-idba_ud.json --list-dependencies

INFO "Running a dry run with idba_ud"
anvi-run-workflow -w metagenomics --config config-idba_ud.json
anvi-run-workflow -w metagenomics --config config-idba_ud.json --save-workflow-graph

INFO "Running a dry run with idba_ud with no qc"
anvi-run-workflow -w metagenomics --config config-idba_ud-no-qc.json --list-dependencies

INFO "Running a dry run with metaspades using scaffolds"
anvi-run-workflow -w metagenomics --config config-metaspades-no-qc-use-scaffolds.json
anvi-run-workflow -w metagenomics --config config-metaspades-no-qc-use-scaffolds.json --save-workflow-graph
anvi-run-workflow -w metagenomics --config config-metaspades-no-qc-use-scaffolds.json --list-dependencies

INFO "Running a dry run with metaspades"
anvi-run-workflow -w metagenomics --config config-metaspades.json
anvi-run-workflow -w metagenomics --config config-metaspades.json --save-workflow-graph
anvi-run-workflow -w metagenomics --config config-metaspades.json --list-dependencies
