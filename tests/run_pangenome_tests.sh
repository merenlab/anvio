#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the pan analysis directory"
mkdir $output_dir/pan_test
cp $files/mock_data_for_pangenomics/*.fa                      $output_dir/pan_test/
cp $files/mock_data_for_pangenomics/emapper/*.annotations     $output_dir/pan_test/
cp $files/mock_data_for_pangenomics/external-genomes.txt      $output_dir/pan_test/
cp $files/mock_data_for_pangenomics/example-PC-collection.txt $output_dir/pan_test/
cp $files/mock_data_for_pangenomics/default-state.json        $output_dir/pan_test/
cd $output_dir/pan_test

INFO "Generating contigs databases for external genomes"
anvi-script-FASTA-to-contigs-db 01.fa
anvi-script-FASTA-to-contigs-db 02.fa
anvi-script-FASTA-to-contigs-db 03.fa

# INFO "Importing functions into the contigs database"
# anvi-script-run-eggnog-mapper -c 01.db --annotation aa_sequences_01.emapper.annotations --use-version 0.12.6
# anvi-script-run-eggnog-mapper -c 02.db --annotation aa_sequences_02.emapper.annotations --use-version 0.12.6
# anvi-script-run-eggnog-mapper -c 03.db --annotation aa_sequences_03.emapper.annotations --use-version 0.12.6

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt -o TEST-GENOMES.h5

INFO "Running the pangenome anaysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.h5 -o TEST/ -J TEST --use-ncbi-blast

INFO "Running the pangenome analysis again utilizing previous search results"
anvi-pan-genome -g TEST-GENOMES.h5 -o TEST/ -J ANOTHER_TEST --use-ncbi-blast --min-occurrence 2

INFO "Importing the default state for pretty outputs"
anvi-import-state -p TEST/TEST-PAN.db -s default-state.json -n default
anvi-import-state -p TEST/ANOTHER_TEST-PAN.db -s default-state.json -n default

INFO "Importing an example collection of protein clusters"
anvi-import-collection -p TEST/TEST-PAN.db -C test_collection example-PC-collection.txt

INFO "Exporting the collection 'test_collection'"
anvi-export-collection -p TEST/TEST-PAN.db -C test_collection -O exported_collection --include-unbinned

INFO "Summarizing the pan, using the test collection (in quick mode)"
anvi-summarize -p TEST/TEST-PAN.db -g TEST-GENOMES.h5 -C test_collection -o TEST_SUMMARY_QUICK --quick

INFO "Summarizing the pan, using the test collection"
anvi-summarize -p TEST/TEST-PAN.db -g TEST-GENOMES.h5 -C test_collection -o TEST_SUMMARY

INFO "Displaying the initial pangenome analysis results"
anvi-display-pan -p TEST/TEST-PAN.db -s TEST/TEST-SAMPLES.db -g TEST-GENOMES.h5 --title "A mock pangenome analysis"

INFO "Displaying the second pangenome analysis results"
anvi-display-pan -p TEST/ANOTHER_TEST-PAN.db -s TEST/ANOTHER_TEST-SAMPLES.db -g TEST-GENOMES.h5 --title "A mock pangenome analysis (with --min-occurrence 2)"
