#!/bin/bash
source 00.sh
set -e

# this is necessary for any system that will run these tests
# due to some historical crap:
pip install h5py==3.9.0

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

#########################################################################################
# POUCHITIS DATA FROM 2015
#########################################################################################

cd $output_dir
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
anvi-setup-scg-taxonomy $thread_controller
INFO "[$TEST] Running HMMs for Bacterial SCGs"
anvi-run-hmms -c CONTIGS.db -I Bacteria_71 $thread_controller
INFO "[$TEST] Running SCG taxonomy"
anvi-run-scg-taxonomy -c CONTIGS.db $thread_controller
INFO "[$TEST] Estimating SCG taxonomy"
anvi-estimate-scg-taxonomy -p PROFILE.db -c CONTIGS.db -C FINAL
INFO "[$TEST] Done!"

#########################################################################################
# TARA GENOMES POUCHITIS DATA FROM 2027
#########################################################################################

cd $output_dir
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
anvi-setup-scg-taxonomy $thread_controller
INFO "[$TEST] Running HMMs"
anvi-run-hmms -c CONTIGS.db $thread_controller
INFO "[$TEST] Estimating genome completeness"
anvi-estimate-genome-completeness -c CONTIGS.db
INFO "[$TEST] Running SCG taxonomy"
anvi-run-scg-taxonomy -c CONTIGS.db $thread_controller
INFO "[$TEST] Estimating SCG taxonomy"
anvi-estimate-scg-taxonomy -c CONTIGS.db
INFO "[$TEST] Done!"

#########################################################################################
# PHROCLOROCOCCUS METAPANGENOME FROM 2018
#########################################################################################

cd $output_dir
TEST="PHROCLOROCOCCUS_METAPANGENOME_FROM_2018"
mkdir $TEST && cd $TEST
INFO "[$TEST] Downloading the data pack"
curl -L https://ndownloader.figshare.com/files/9416623 -o ANVIO-METAPANGENOME-FOR-PROCHLOROCOCCUS-ISOLATES.tar.gz
INFO "[$TEST] Unpacking"
tar -zxvf ANVIO-METAPANGENOME-FOR-PROCHLOROCOCCUS-ISOLATES.tar.gz && cd ANVIO-METAPANGENOME-FOR-PROCHLOROCOCCUS-ISOLATES
INFO "[$TEST] Migrating the pan-db"
ANVIO_SAMPLES_DB=Prochlorococcus-METAPAN-SAMPLES.db anvi-migrate Prochlorococcus-PAN-PAN.db --migrate-quickly
INFO "[$TEST] Migrating the genomes-storage-db"
anvi-migrate Prochlorococcus-GENOMES.h5 --migrate-quickly
INFO "[$TEST] Importing additional data layers into the pan-db"
# Converting additional data layers to the new format -- this will
# only be necessary for this particular project
echo -e "pc_name\tDetection#EDG\tDetection#ECG\tDetection#NA" | sed 's/#/!/g' > ENVIRONMENTAL-CORE-FIXED.txt
grep -v 'pc_name' ENVIRONMENTAL-CORE.txt | awk 'BEGIN{FS=";"}{print $1 "\t" $2 "\t" $3}' >> ENVIRONMENTAL-CORE-FIXED.txt
anvi-import-misc-data -t items ENVIRONMENTAL-CORE-FIXED.txt -p Prochlorococcus-PAN-PAN.db
# let's remove unnecessary files while we're at it
rm -rf *SAMPLES* ENVIRONMENTAL-CORE.txt 00_run.sh ENVIRONMENTAL-CORE-GENES.txt
INFO "[$TEST] Running anvi-db-info on the pan-db"
anvi-db-info -p Prochlorococcus-PAN-PAN.db
INFO "[$TEST] Done!"
