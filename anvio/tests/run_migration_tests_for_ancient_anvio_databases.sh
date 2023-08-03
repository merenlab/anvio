#!/bin/bash
source 00.sh
set -e

# this is necessary for any system that will run these tests
# due to some historical crap:

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

cd $output_dir

#########################################################################################

TEST="POUCHITIS_METAGENOMES_FROM_2015"
mkdir $TEST && cd $TEST
INFO "[$TEST] Downloading the data pack"
curl -L https://ndownloader.figshare.com/files/6035580 -o P214-MERGED.tar.gz
INFO "[$TEST] Unpacking"
tar -zxvf P214-MERGED.tar.gz && cd P214-MERGED/
INFO "[$TEST] Migrating databases"
anvi-migrate *db --migrate-safely
INFO "[$TEST] Downloading state"
curl -L http://merenlab.org/files/P-214-state.json -o P-214-state.json
INFO "[$TEST] Importing state"
anvi-import-state -p PROFILE.db -n default -s P-214-state.json
INFO "[$TEST] Running the interactive interface (--dry)"
anvi-interactive -p PROFILE.db -c CONTIGS.db --split-hmm-layers --dry
INFO "[$TEST] Setting up SCG taxonomy databases"
anvi-setup-scg-taxonomy --redo-databases --num-threads 4
INFO "[$TEST] Running HMMs for Bacterial SCGs"
anvi-run-hmms -c CONTIGS.db -I Bacteria_71 --num-threads 4
INFO "[$TEST] Running SCG taxonomy"
anvi-run-scg-taxonomy -c CONTIGS.db --num-threads 4
INFO "[$TEST] Estimating SCG taxonomy"
anvi-estimate-scg-taxonomy -p PROFILE.db -c CONTIGS.db -C FINAL
INFO "[$TEST] Done!"
cd ..

#########################################################################################

TEST="A_TARA_MAG_FROM_2017"
mkdir $TEST && cd $TEST
INFO "[$TEST] Downloading the data pack"
curl -L https://ndownloader.figshare.com/files/8248433 -o TARA_ANW_MAG_00006.tar.gz
INFO "[$TEST] Unpacking"
tar -zxvf TARA_ANW_MAG_00006.tar.gz && cd TARA_ANW_MAG_00006
gzip -d AUXILIARY-DATA.h5.gz
INFO "[$TEST] Migrating databases"
anvi-migrate PROFILE.db CONTIGS.db --migrate-quickly
INFO "[$TEST] Running the interactive interface (--dry)"
anvi-interactive -p PROFILE.db -c CONTIGS.db --dry
INFO "[$TEST] Setting up SCG taxonomy databases"
anvi-setup-scg-taxonomy --redo-databases --num-threads 4
INFO "[$TEST] Running HMMs"
anvi-run-hmms -c CONTIGS.db -T 4
INFO "[$TEST] Estimating genome completeness"
anvi-estimate-genome-completeness -c CONTIGS.db
INFO "[$TEST] Running SCG taxonomy"
anvi-run-scg-taxonomy -c CONTIGS.db --num-threads 4
INFO "[$TEST] Estimating SCG taxonomy"
anvi-estimate-scg-taxonomy -c CONTIGS.db
INFO "[$TEST] Done!"
cd ..

#########################################################################################
