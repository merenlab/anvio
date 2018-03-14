#!/bin/bash
source 00.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd sandbox/files_for_manual_interactive
rm -rf test-output
mkdir test-output

INFO "Anvi'o version ..."
anvi-profile --version

INFO "Generating a newick file from the data ..."
anvi-matrix-to-newick view_data.txt -o test-output/tree.txt

INFO "Running a dry run to generate a profile database (to populate it with additional misc data)"
anvi-interactive -d view_data.txt --manual-mode -p test-output/test.db --dry-run

INFO "Importing items additional data"
anvi-import-misc-data -p test-output/test.db additional_view_data.txt -t items

INFO "Importing layers additional data"
anvi-import-misc-data -p test-output/test.db samples-information.txt -t layers

INFO "Importing layer order data"
anvi-import-misc-data -p test-output/test.db samples-order.txt -t layer_orders

INFO "Running the interactive interface on files"
anvi-interactive -f fasta.fa \
                 -d view_data.txt \
                 -p test-output/test.db \
                 -t test-output/tree.txt \
                 --title 'Interactive Tree For User Provided Files' \
                 --manual-mode
