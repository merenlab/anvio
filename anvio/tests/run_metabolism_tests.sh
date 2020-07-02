#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

# Change this whenever re-annotation of the test databases is done
CURRENT_MODULES_DB_HASH="ba0d0dda2796"

INFO "Setting up the metabolism test directory"
mkdir $output_dir/metabolism_test
cp $files/real_data_for_testing/*.fa                $output_dir/metabolism_test
cp $files/real_data_for_testing/IGD_SUBSET/*.db     $output_dir/metabolism_test
cp $files/real_data_for_testing/KOfams/*KOfams.txt  $output_dir/metabolism_test
cp $files/real_data_for_testing/*genomes.txt       $output_dir/metabolism_test
cp $files/real_data_for_testing/bin-ids.txt         $output_dir/metabolism_test
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

INFO "Setting modules db hash in databases to $CURRENT_MODULES_DB_HASH"
anvi-db-info --self-key modules_db_hash --self-value $CURRENT_MODULES_DB_HASH --just-do-it B_thetaiotamicron_VPI-5482.db
anvi-db-info --self-key modules_db_hash --self-value $CURRENT_MODULES_DB_HASH --just-do-it P_marinus_CCMP1375.db
anvi-db-info --self-key modules_db_hash --self-value $CURRENT_MODULES_DB_HASH --just-do-it S_islandicus_LS215.db

INFO "Estimating metabolism on a single contigs database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db -O single_contigs_db

INFO "Estimating metabolism using metagenome mode"
anvi-estimate-metabolism -c CONTIGS.db --metagenome-mode -O metagenome_mode

INFO "Estimating metabolism on a collection"
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C bins_for_testing -O collection

INFO "Estimating metabolism on a single bin in a collection"
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C bins_for_testing -b E_facealis -O single_bin

INFO "Estimating metabolism using a bin IDs file"
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C bins_for_testing -B bin-ids.txt -O bin_ids_file

INFO "Estimating metabolism on external genomes"
anvi-estimate-metabolism -e external-genomes.txt -O external

INFO "Estimating metabolism on internal genomes"
anvi-estimate-metabolism -i internal-genomes.txt -O internal

INFO "Estimating metabolism on metagenomes file"
anvi-estimate-metabolism -M metagenomes.txt -O metagenomes
