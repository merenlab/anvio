#!/bin/bash
source 00.sh
set -e

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

for f in fasta.fa view_data.txt additional_view_data.txt samples-information.txt samples-order.txt
do
    cp $files/files_for_manual_interactive/$f $output_dir/
done


INFO "Generating a newick file from the data ..."
anvi-matrix-to-newick $output_dir/view_data.txt \
                      -o $output_dir/tree.txt

INFO "Running a dry run to generate a profile database (to populate it with additional misc data)"
anvi-interactive -d $output_dir/view_data.txt \
                 --manual-mode \
                 -p $output_dir/test.db \
                 --dry-run

INFO "Importing items additional data"
anvi-import-misc-data -p $output_dir/test.db \
                      $output_dir/additional_view_data.txt \
                      -t items

INFO "Importing layers additional data"
anvi-import-misc-data -p $output_dir/test.db \
                      $output_dir/samples-information.txt \
                      -t layers

INFO "Importing layer order data"
anvi-import-misc-data -p $output_dir/test.db \
                      $output_dir/samples-order.txt \
                      -t layer_orders

INFO "Running the interactive interface on files"
anvi-interactive -f $output_dir/fasta.fa \
                 -d $output_dir/view_data.txt \
                 -p $output_dir/test.db \
                 -t $output_dir/tree.txt \
                 --title "Interactive Interface Without Anvi'o Files"\
                 --manual-mode \
                 $dry_run_controller
