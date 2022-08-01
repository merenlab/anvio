#!/bin/bash

set -e

./run_contigs_workflow_tests.sh
./run_phylogenomics_workflow_tests.sh
./run_pangenomics_workflow_tests.sh
./run_metagenomics_workflow_tests.sh
