#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the metabolism test directory"
mkdir $output_dir/functions_display
cp $files/data/genomes/bacteria/*.db                    $output_dir/functions_display
cp $files/data/genomes/archaea/*.db                     $output_dir/functions_display
cp $files/data/input_files/external-genomes.txt         $output_dir/functions_display
cp $files/data/input_files/groups.txt                   $output_dir/functions_display
cd $output_dir/functions_display

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Running anvi-display-functions without groups"
anvi-display-functions -e external-genomes.txt \
                       --annotation-source COG20_FUNCTION \
                       -p COG20_FUNCTION.db \
                       --dry

INFO "Running anvi-display-functions with groups"
anvi-display-functions -e external-genomes.txt \
                       --groups groups.txt \
                       --annotation-source COG20_FUNCTION \
                       -p COG20_FUNCTION_W_GROUPS.db
