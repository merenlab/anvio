#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Setting up the metabolism test directory"
mkdir $output_dir/metabolism_test
cp $files/data/genomes/bacteria/*.db                    $output_dir/metabolism_test
cp $files/data/genomes/archaea/*.db                     $output_dir/metabolism_test
cp $files/data/metagenomes/human_gut/IGD_SUBSET/*.db    $output_dir/metabolism_test
cp $files/data/input_files/*.txt                        $output_dir/metabolism_test
cd $output_dir/metabolism_test

INFO "Migrating all databases"
anvi-migrate *db --migrate-dbs-quickly


INFO "Estimating metabolism on a single contigs database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         -O single_contigs_db
SHOW_FILE single_contigs_db_modules.txt

INFO "Estimating metabolism using metagenome mode"
anvi-estimate-metabolism -c CONTIGS.db \
                         --metagenome-mode \
                         -O metagenome_mode
SHOW_FILE metagenome_mode_modules.txt

INFO "Estimating metabolism on a collection"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -O collection
SHOW_FILE collection_modules.txt

INFO "Estimating metabolism on a single bin in a collection"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -b E_facealis \
                         -O single_bin
SHOW_FILE single_bin_modules.txt

INFO "Estimating metabolism using a bin IDs file"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -B bin-ids.txt \
                         -O bin_ids_file
SHOW_FILE bin_ids_file_modules.txt

INFO "Estimating metabolism on external genomes"
anvi-estimate-metabolism -e external-genomes.txt \
                         -O external
SHOW_FILE external_modules.txt

INFO "Estimating metabolism on internal genomes"
anvi-estimate-metabolism -i internal-genomes.txt \
                         -O internal
SHOW_FILE internal_modules.txt

INFO "Estimating metabolism on metagenomes file"
anvi-estimate-metabolism -M metagenomes.txt \
                         -O metagenomes
SHOW_FILE metagenomes_modules.txt

INFO "Trying a different module completeness threshold"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --module-completion-threshold 0 \
                         -O nondefault_threshold
SHOW_FILE nondefault_threshold_modules.txt

INFO "Listing output modes for single database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         --list-available-modes


INFO "Listing output modes for multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --list-available-modes

INFO "Listing custom output headers for single database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         --list-available-output-headers

INFO "Listing custom output headers for multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --list-available-output-headers

INFO "Generating long format output files on single database"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --kegg-output-modes kofam_hits,modules,kofam_hits_in_modules,modules_custom \
                         --custom-output-headers kegg_module,module_is_complete,module_name \
                         -O long_format_single
SHOW_FILE long_format_single_kofam_hits.txt
SHOW_FILE long_format_single_modules.txt
SHOW_FILE long_format_single_kofam_hits_in_modules.txt
SHOW_FILE long_format_single_modules_custom.txt

INFO "Generating long format output files in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --kegg-output-modes kofam_hits,modules,kofam_hits_in_modules \
                         -O long_format_multi
SHOW_FILE long_format_multi_kofam_hits.txt
SHOW_FILE long_format_multi_modules.txt
SHOW_FILE long_format_multi_kofam_hits_in_modules.txt

# below we try every single header for modules_custom mode
INFO "Generating modules custom output file on single database"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --kegg-output-modes modules_custom \
                         -O modules_custom_single \
                         --custom-output-headers kegg_module,module_is_complete,module_completeness,module_name,module_class,module_category,module_subcategory,module_definition,gene_caller_ids_in_module,gene_caller_id,kofam_hits_in_module,kofam_hit,contig,path_id,path,path_completeness,genome_name
SHOW_FILE modules_custom_single_modules_custom.txt

INFO "Generating modules custom output file in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --kegg-output-modes modules_custom \
                         -O modules_custom_multi \
                         --custom-output-headers kegg_module,module_is_complete,module_completeness,module_name,module_class,module_category,module_subcategory,module_definition,gene_caller_ids_in_module,gene_caller_id,kofam_hits_in_module,kofam_hit,contig,path_id,path,path_completeness,genome_name,db_name
SHOW_FILE modules_custom_multi_modules_custom.txt

INFO "Generating matrix output files in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --matrix-format \
                         -O matrix_format_multi
SHOW_FILE matrix_format_multi-completeness-MATRIX.txt
SHOW_FILE matrix_format_multi-presence-MATRIX.txt
SHOW_FILE matrix_format_multi-ko_hits-MATRIX.txt

INFO "Generating JSON output (debug option)"
anvi-estimate-metabolism -c S_islandicus_LS215.db \
                         --get-raw-data-as-json estimation_data \
                         --store-json-without-estimation

INFO "Estimating from JSON output (debug option)"
anvi-estimate-metabolism --estimate-from-json estimation_data.json \
                         -O from_json

INFO "Testing --add-coverage for genome mode"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -O genome_coverage \
                         --add-coverage \
                         --kegg-output-modes modules,kofam_hits_in_modules
SHOW_FILE genome_coverage_modules.txt
SHOW_FILE genome_coverage_kofam_hits_in_modules.txt

INFO "Testing --add-coverage flag for a collection"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -O collection_coverage \
                         --add-coverage \
                         --kegg-output-modes modules,kofam_hits_in_modules
SHOW_FILE collection_coverage_modules.txt

INFO "Testing --add-coverage for metagenome mode"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         --metagenome-mode \
                         -O metagenome_coverage \
                         --add-coverage \
                         --kegg-output-modes modules,kofam_hits_in_modules

INFO "Listing custom output headers with --add-coverage enabled"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         --add-coverage \
                         --list-available-output-headers

INFO "Generating custom output with --add-coverage enabled"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         --add-coverage \
                         --kegg-output-modes modules_custom \
                         -O modules_custom_coverage \
                         --custom-output-headers kegg_module,module_is_complete,DAY_15A_gene_coverages,DAY_15A_avg_coverage,DAY_15A_gene_detection,DAY_15A_avg_detection
