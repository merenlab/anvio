#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the structure-informed pangenome analysis directory"
mkdir -p $output_dir/
cp $files/mock_data_for_pangenomics/*.db                      $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt      $output_dir/
cp $files/mock_data_for_pangenomics/example-gene-clusters-collection.txt $output_dir/
cp $files/mock_data_for_pangenomics/scg-gene-clusters-for-phylogenomics.txt $output_dir/
cp $files/mock_data_for_pangenomics/default-state.json        $output_dir/
cp $files/example_description.md                              $output_dir/
cp $files/mock_data_for_pangenomics/group-information.txt     $output_dir/
cd $output_dir/

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o TEST-GENOMES.db \
                         --no-progress

INFO "Running the structure-informed pangenome analysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.db \
                -o TEST/ \
                -n TEST \
                --pan-mode structure-informed \
                --description example_description.md \
                --no-progress \
                $thread_controller

INFO "Importing group information as misc data for layers"
anvi-import-misc-data -p TEST/TEST-STRUCTURE-PAN.db \
                      -t layers \
                      group-information.txt \
                      --no-progress

INFO "Estimating enriched functions per pan group"
anvi-compute-functional-enrichment-in-pan -p TEST/TEST-STRUCTURE-PAN.db \
                                          -g TEST-GENOMES.db \
                                          --category group \
                                          --annotation-source COG20_FUNCTION \
                                          -o functions-enrichment.txt \
                                          -F functional-occurence.txt \
                                          --no-progress
SHOW_FILE functions-enrichment.txt
SHOW_FILE functional-occurence.txt

INFO "Exporting concatenated amino acid sequences for some SCG gene clusters for phylogenomics"
anvi-get-sequences-for-gene-clusters -p TEST/TEST-STRUCTURE-PAN.db \
                                     -g TEST-GENOMES.db \
                                     -C collection_for_phylogenomics \
                                     -b SCGs \
                                     -o SOME_VARIABLE_SCGs.fa \
                                     --concatenate-gene-clusters \
                                     --no-progress

INFO "Summarizing the pan, using the test collection (in quick mode)"
anvi-summarize -p TEST/TEST-STRUCTURE-PAN.db \
               -g TEST-GENOMES.db \
               -C test_collection \
               -o TEST_SUMMARY_QUICK \
               --quick \
               --no-progress

INFO "Summarizing the pan, using the test collection"
anvi-summarize -p TEST/TEST-STRUCTURE-PAN.db \
               -g TEST-GENOMES.db \
               -C test_collection \
               -o TEST_SUMMARY \
               --no-progress

INFO "Splitting bins in the pan genome into smaller, self-contained pan databases"
anvi-split -p TEST/TEST-STRUCTURE-PAN.db \
           -g TEST-GENOMES.db \
           -C test_collection \
           -o TEST_SPLIT_PAN

INFO "Resulting split pans"
ls -l TEST_SPLIT_PAN/*/*db

INFO "Taking a look at the make up one of the split pans"
anvi-db-info TEST_SPLIT_PAN/GENE_CLUSTER_BIN_1_CORE/PAN.db

INFO "Listing collections available"
anvi-show-collections-and-bins -p TEST/TEST-STRUCTURE-PAN.db \
                               --no-progress

INFO "Importing the default state for pretty outputs"
anvi-import-state -p TEST/TEST-STRUCTURE-PAN.db -s default-state.json -n default
anvi-import-state -p TEST/ANOTHER_TEST-PAN.db -s default-state.json -n default

INFO "Displaying the initial structure informed pangenome analysis results"
anvi-display-pan -p TEST/TEST-STRUCTURE-PAN.db \
                 -g TEST-GENOMES.db \
                 --title "A mock structure informed pangenome analysis" \
                 --no-progress \
                 $dry_run_controller
