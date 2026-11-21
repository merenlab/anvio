#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Set the path for example files
example_files=$files/example_files_for_gene_level_stats

cp $example_files/CONTIGS.db $output_dir/CONTIGS.db
cp $example_files/PROFILE.db $output_dir/PROFILE.db
cp $example_files/AUXILIARY-DATA.db $output_dir/AUXILIARY-DATA.db

INFO "Migrate all dbs"
anvi-migrate $output_dir/*db --migrate-quickly

# Generate gene-level-stats database
INFO "Computing gene level stats database"
anvi-gen-gene-level-stats-databases -p $output_dir/PROFILE.db \
                                    -c $output_dir/CONTIGS.db \
                                    -C DEFAULT -b EVERYTHING

# Visualize the bin in gene-mode:
anvi-interactive -p $output_dir/PROFILE.db \
                 -c $output_dir/CONTIGS.db  \
                 -C DEFAULT \
                 -b EVERYTHING \
                 --gene-mode
