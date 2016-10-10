#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/pan_test
cp $files/mock_data_for_pangenomics/*.fa                 $output_dir/pan_test/
cp $files/mock_data_for_pangenomics/external-genomes.txt $output_dir/pan_test/
cp $files/mock_data_for_pangenomics/default-state.json   $output_dir/pan_test/
cd $output_dir/pan_test

INFO "Generating contigs databases for external genomes"
anvi-script-FASTA-to-contigs-db 01.fa
anvi-script-FASTA-to-contigs-db 02.fa
anvi-script-FASTA-to-contigs-db 03.fa

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt -o TEST-GENOMES.h5

INFO "Running the pangenome anaysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.h5 -o TEST/ -J TEST

INFO "Importing the default state for pretty outputs"
anvi-import-state -p TEST/TEST-PAN.db -s default-state.json -n default

INFO "Displaying the pangenome analysis results"
anvi-display-pan -p TEST/TEST-PAN.db -s TEST/TEST-SAMPLES.db -g TEST-GENOMES.h5 --title "A mock pangenome analysis"
