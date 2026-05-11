#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# This test requires a `structures-txt` fixture pointing at predicted protein
# structure files for the gene-cluster representatives produced from the mock
# E. faecalis genomes. To regenerate the fixture:
#
#   1. Run anvi-pan-genome in default (sequence) mode against TEST-GENOMES.db.
#   2. Predict structures (AlphaFold/ColabFold/ESMFold) for every defline of
#      gene-cluster-representatives-aa.fa, naming each file <defline>.pdb.
#   3. Drop them into mock_data_for_pangenomics/structures/ and write
#      mock_data_for_pangenomics/structures-txt.txt with `gene_id\tpath`
#      pointing at each PDB.

INFO "Setting up the structure-informed pangenome analysis directory"
mkdir -p $output_dir/
cp $files/mock_data_for_pangenomics/*.db                      $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt      $output_dir/
cp $files/mock_data_for_pangenomics/example-gene-clusters-collection.txt $output_dir/
cp $files/mock_data_for_pangenomics/scg-gene-clusters-for-phylogenomics.txt $output_dir/
cp $files/mock_data_for_pangenomics/default-state.json        $output_dir/
cp $files/example_description.md                              $output_dir/
cp $files/mock_data_for_pangenomics/group-information.txt     $output_dir/
if [ -f $files/mock_data_for_pangenomics/structures-txt.txt ]; then
    cp $files/mock_data_for_pangenomics/structures-txt.txt    $output_dir/
    cp -r $files/mock_data_for_pangenomics/structures         $output_dir/
fi
cd $output_dir/

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o TEST-GENOMES.db \
                         --no-progress

if [ ! -f structures-txt.txt ]; then
    echo
    echo "ERROR: the structure-informed pangenomics component test cannot run because the"
    echo "structures-txt fixture is missing. Anvi'o is failing loudly on purpose: a silent"
    echo "skip would leave this test green in CI even though no new code is being exercised."
    echo
    echo "To regenerate the fixture, follow the recipe at the top of this script:"
    echo "  1. Run anvi-pan-genome (default mode) against TEST-GENOMES.db."
    echo "  2. Predict structures (AlphaFold/ColabFold/ESMFold) for every defline of"
    echo "     gene-cluster-representatives-aa.fa, naming each file <defline>.pdb."
    echo "  3. Drop them into mock_data_for_pangenomics/structures/ and write"
    echo "     mock_data_for_pangenomics/structures-txt.txt with 'gene_id<TAB>path' rows."
    echo
    echo "Until that lands, this code path is covered by unit tests in"
    echo "  anvio/tests/unit/test_structures_txt.py"
    echo "which exercise StructuresTxt and the structure-mode branch of gen_mcl_input"
    echo "without needing a heavy structure-prediction fixture."
    echo
    exit 1
fi

INFO "Running the structure-informed pangenome analysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.db \
                -o TEST/ \
                -n TEST \
                --pan-mode structure-informed \
                --structures-txt structures-txt.txt \
                --min-tm-score 0.5 \
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

INFO "Importing collections of gene clusters"
anvi-import-collection -p TEST/TEST-STRUCTURE-PAN.db \
                       -C test_collection example-protein-structure-informed-gene-clusters-collection.txt \
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
