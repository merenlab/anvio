#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# Set INSEQ files path
inseq_files=$files/example_files_for_inseq_tnseq

cp $inseq_files/CONTIGS.db $output_dir/CONTIGS.db
cp $inseq_files/PROFILE.db $output_dir/PROFILE.db
cp $inseq_files/AUXILIARY-DATA.db $output_dir/AUXILIARY-DATA.db

INFO "Migrate all dbs"
anvi-migrate $output_dir/*db --migrate-quickly

# Generate gene-level-stats database
INFO "Computing gene level stats database"
anvi-gen-gene-level-stats-databases -p $output_dir/PROFILE.db \
                                    -c $output_dir/CONTIGS.db \
                                    -C DEFAULT -b EVERYTHING

# Delete the GENE database
rm -r $output_dir/GENES/

# Compute INSeq stats database
INFO "Computing INSeq stats database"
anvi-gen-gene-level-stats-databases -p $output_dir/PROFILE.db \
                                    -c $output_dir/CONTIGS.db  \
                                    -C DEFAULT \
                                    -b EVERYTHING \
                                    --inseq-stats \
                                    --just-do-it

# Visualize the bin in gene-mode:
anvi-interactive -p $output_dir/PROFILE.db \
                 -c $output_dir/CONTIGS.db  \
                 -C DEFAULT \
                 -b EVERYTHING \
                 --gene-mode
