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

INFO "Generating a pangenome graph from a YAML file"
anvi-pan-genome-graph --pan-graph-yaml example-pangenome-graph.yaml \
                      -o FROM-YAML

INFO "Generating a pangenome graph from a JSON file"
anvi-pan-genome-graph -i example-pangenome-graph.json \
                      -o FROM-JSON

INFO "Generating a pangenome graph from a pan-db"
anvi-pan-genome-graph -p TEST/TEST-PAN.db \
                      -g TEST-GENOMES.db \
                      -o FROM-PAN \
                      --project-name TEST \
                      -e external-genomes.txt \
                      $thread_controller

INFO "Displaying pangenome graph"
anvi-display-pan-graph -i FROM-PAN/TEST-JSON.json \
                       -p TEST/TEST-PAN.db \
                       -g TEST-GENOMES.db \
                       $dry_run_controller
