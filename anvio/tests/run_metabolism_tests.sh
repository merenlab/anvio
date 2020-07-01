#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the metabolism test directory"
mkdir $output_dir/metabolism_test
cp $files/real_data_for_testing/*.fa                $output_dir/metabolism_test
cp $files/real_data_for_testing/IGD_SUBSET/*.db     $output_dir/metabolism_test
cp $files/real_data_for_testing/KOfams/*KOfams.txt  $output_dir/metabolism_test
cp $files/real_data_for_testing/*-genomes.txt       $output_dir/metabolism_test
cd $output_dir/metabolism_test

INFO "Generating contigs databases for external genomes"
anvi-script-FASTA-to-contigs-db B_thetaiotamicron_VPI-5482.fa
anvi-script-FASTA-to-contigs-db P_marinus_CCMP1375.fa
anvi-script-FASTA-to-contigs-db S_islandicus_LS215.fa

INFO "Importing KOfams into the contigs databases"
anvi-import-functions -c B_thetaiotamicron_VPI-5482.db -i B_thetaiotamicron_VPI-5482_KOfams.txt
anvi-import-functions -c P_marinus_CCMP1375.db -i P_marinus_CCMP1375_KOfams.txt
anvi-import-functions -c S_islandicus_LS215.db -i S_islandicus_LS215_KOfams.txt

# IGD_SUBSET database is already annotated, but just in case that breaks, here is how to fix it:
# anvi-import-functions -c CONTIGS.db -i IGD_SUBSET_KOfams.txt

INFO "Estimating metabolism on a single contigs database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db -O single_contigs_db
