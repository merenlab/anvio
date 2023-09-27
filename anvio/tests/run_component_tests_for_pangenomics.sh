#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the pan analysis directory"
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

INFO "Running the pangenome anaysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.db \
                -o TEST/ \
                -n TEST \
                --use-ncbi-blast \
                --description example_description.md \
                --no-progress \
                $thread_controller

INFO "Running ANI on genomes and storing results in the PAN database"
anvi-compute-genome-similarity -e external-genomes.txt \
                               --program pyANI \
                               -o ANI_TEST \
                               --log-file ANI_LOG.txt \
                               -p TEST/TEST-PAN.db \
                               --no-progress \
                               $thread_controller

INFO "Running the pangenome analysis again utilizing previous search results"
anvi-pan-genome -g TEST-GENOMES.db \
                -o TEST/ \
                -n ANOTHER_TEST \
                --use-ncbi-blast \
                --min-occurrence 2 \
                --description example_description.md \
                --no-progress \
                $thread_controller

INFO "Importing collections of gene clusters"
anvi-import-collection -p TEST/TEST-PAN.db \
                       -C test_collection example-gene-clusters-collection.txt \
                       --no-progress

anvi-import-collection -p TEST/TEST-PAN.db \
                       -C collection_for_phylogenomics scg-gene-clusters-for-phylogenomics.txt \
                       --no-progress

INFO "Exporting the collection 'test_collection'"
anvi-export-collection -p TEST/TEST-PAN.db \
                       -C test_collection \
                       -O exported_collection \
                       --include-unbinned \
                       --no-progress

INFO "List available aligners for aligning sequences in gene clusters"
anvi-get-sequences-for-gene-clusters --list-aligners

INFO "Exporting aligned amino acid sequences for some gene clusters"
anvi-get-sequences-for-gene-clusters -p TEST/TEST-PAN.db \
                                     -g TEST-GENOMES.db \
                                     -C test_collection \
                                     -b GENE_CLUSTER_BIN_1_CORE \
                                     -o aligned_gene_sequences_in_GENE_CLUSTER_BIN_1_CORE_AA.fa \
                                     --no-progress

mkdir SPLIT_GENE_CLUSTER_FILES
INFO "Exporting aligned amino acid sequences for some gene clusters as individual files"
anvi-get-sequences-for-gene-clusters -p TEST/TEST-PAN.db \
                                     -g TEST-GENOMES.db \
                                     -C test_collection \
                                     -b GENE_CLUSTER_BIN_1_CORE \
                                     --split-output-per-gene-cluster \
                                     -O SPLIT_GENE_CLUSTER_FILES/MY_PROJECT \
                                     --no-progress

INFO "Exporting aligned DNA sequences for some gene clusters"
anvi-get-sequences-for-gene-clusters -p TEST/TEST-PAN.db \
                                     -g TEST-GENOMES.db \
                                     -C test_collection \
                                     -b GENE_CLUSTER_BIN_1_CORE \
                                     -o aligned_gene_sequences_in_GENE_CLUSTER_BIN_1_CORE_DNA.fa \
                                     --report-DNA-sequences \
                                     --no-progress

INFO "First five line from the AA output"
head -n 5 aligned_gene_sequences_in_GENE_CLUSTER_BIN_1_CORE_AA.fa

INFO "First five line from the DNA output"
head -n 5 aligned_gene_sequences_in_GENE_CLUSTER_BIN_1_CORE_DNA.fa

INFO "Importing group information as misc data for layers"
anvi-import-misc-data -p TEST/TEST-PAN.db \
                      -t layers \
                      group-information.txt \
                      --no-progress

INFO "Estimating enriched functions per pan group"
anvi-compute-functional-enrichment-in-pan -p TEST/TEST-PAN.db \
                                          -g TEST-GENOMES.db \
                                          --category group \
                                          --annotation-source COG20_FUNCTION \
                                          -o functions-enrichment.txt \
                                          -F functional-occurence.txt \
                                          --no-progress
SHOW_FILE functions-enrichment.txt
SHOW_FILE functional-occurence.txt

INFO "Exporting concatenated amino acid sequences for some SCG gene clusters for phylogenomics"
anvi-get-sequences-for-gene-clusters -p TEST/TEST-PAN.db \
                                     -g TEST-GENOMES.db \
                                     -C collection_for_phylogenomics \
                                     -b SCGs \
                                     -o SOME_VARIABLE_SCGs.fa \
                                     --concatenate-gene-clusters \
                                     --no-progress

INFO "Summarizing the pan, using the test collection (in quick mode)"
anvi-summarize -p TEST/TEST-PAN.db \
               -g TEST-GENOMES.db \
               -C test_collection \
               -o TEST_SUMMARY_QUICK \
               --quick \
               --no-progress

INFO "Summarizing the pan, using the test collection"
anvi-summarize -p TEST/TEST-PAN.db \
               -g TEST-GENOMES.db \
               -C test_collection \
               -o TEST_SUMMARY \
               --no-progress

INFO "Splitting bins in the pan genome into smaller, self-contained pan databases"
anvi-split -p TEST/TEST-PAN.db \
           -g TEST-GENOMES.db \
           -C test_collection \
           -o TEST_SPLIT_PAN

INFO "Resulting split pans"
ls -l TEST_SPLIT_PAN/*/*db

INFO "Taking a look at the make up one of the split pans"
anvi-db-info TEST_SPLIT_PAN/GENE_CLUSTER_BIN_1_CORE/PAN.db

INFO "Listing collections available"
anvi-show-collections-and-bins -p TEST/TEST-PAN.db \
                               --no-progress

INFO "Computing homogeneity for a single gene cluster"
anvi-compute-gene-cluster-homogeneity -p TEST/TEST-PAN.db \
                                      -g TEST-GENOMES.db \
                                      --gene-cluster-id GC_00000001 \
                                      -o gene_cluster_homogeneity_results.txt \
                                      --no-progress \
                                      $thread_controller
SHOW_FILE gene_cluster_homogeneity_results.txt

INFO "Computing homogeneity for a list of gene clusters"
echo -e "GC_00000001\nGC_00000003" > gene_clusters_for_homogeneity.txt
anvi-compute-gene-cluster-homogeneity -p TEST/TEST-PAN.db \
                                      -g TEST-GENOMES.db \
                                      --gene-cluster-ids gene_clusters_for_homogeneity.txt \
                                      -o gene_cluster_homogeneity_results.txt \
                                      --no-progress \
                                      $thread_controller
SHOW_FILE gene_cluster_homogeneity_results.txt

INFO "Computing homogeneity for gene clusters in a bin"
anvi-compute-gene-cluster-homogeneity -p TEST/TEST-PAN.db \
                                      -g TEST-GENOMES.db \
                                      -C test_collection \
                                      -b GENE_CLUSTER_BIN_2 \
                                      -o gene_cluster_homogeneity_results.txt \
                                      --no-progress \
                                      $thread_controller
SHOW_FILE gene_cluster_homogeneity_results.txt

INFO "Importing the default state for pretty outputs"
anvi-import-state -p TEST/TEST-PAN.db -s default-state.json -n default
anvi-import-state -p TEST/ANOTHER_TEST-PAN.db -s default-state.json -n default

INFO "Displaying the initial pangenome analysis results"
anvi-display-pan -p TEST/TEST-PAN.db \
                 -g TEST-GENOMES.db \
                 --title "A mock pangenome analysis" \
                 --no-progress \
                 $dry_run_controller

INFO "Displaying the second pangenome analysis results"
anvi-display-pan -p TEST/ANOTHER_TEST-PAN.db \
                 -g TEST-GENOMES.db \
                 --title "A mock pangenome analysis (with --min-occurrence 2)" \
                 --no-progress \
                 $dry_run_controller
