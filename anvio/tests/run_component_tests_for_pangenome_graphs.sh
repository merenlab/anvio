#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the pan analysis directory"
mkdir -p $output_dir/
cp $files/mock_data_for_pangenomics/*.db                      $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt      $output_dir/
cp $files/example_description.md                              $output_dir/
cp $files/example-pangenome-graph.yaml                        $output_dir/
cp $files/example-pangenome-graph.json.gz                     $output_dir/
cd $output_dir/

# uncompress the JSON file
gzip -d example-pangenome-graph.json.gz

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o TEST-GENOMES.db \
                         --no-progress

INFO "Running the pangenome analysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.db \
                -n TEST \
                --use-ncbi-blast \
                --description example_description.md \
                --no-progress \
                $thread_controller

INFO "Running ANI on genomes and storing results in the PAN database"
anvi-compute-genome-similarity -e external-genomes.txt \
                               --program fastANI \
                               --fragment-length 250 \
                               --min-num-fragments 1 \
                               -o ANI_TEST \
                               --log-file ANI_LOG.txt \
                               -p TEST-PAN.db \
                               --no-progress \
                               $thread_controller

INFO "Generating a pangenome graph from a YAML file"
anvi-pan-genome-graph --pan-graph-yaml example-pangenome-graph.yaml \
                      --project-name FROM-YAML

INFO "Generating a pangenome graph from a pan-db"
anvi-pan-genome-graph -p TEST-PAN.db \
                      -g TEST-GENOMES.db \
                      --project-name TEST \
                      -e external-genomes.txt \
                      $thread_controller

INFO "Generating a summary output for the pangenome graph"
anvi-summarize -p TEST-PAN-GRAPH.db \
               -g TEST-GENOMES.db \
               -o TEST-PAN-GRAPH-SUMMARY
SHOW_FILE TEST-PAN-GRAPH-SUMMARY/REGIONS.txt
SHOW_FILE TEST-PAN-GRAPH-SUMMARY/SYNGCs.txt

INFO "Displaying bona fide pangenome graph from a PAN database"
anvi-display-pan-graph -p TEST-PAN-GRAPH.db \
                       -g TEST-GENOMES.db \
                       $dry_run_controller

INFO "Displaying pangenome graph from a YAML file"
anvi-display-pan-graph -p FROM-YAML-PAN-GRAPH.db
