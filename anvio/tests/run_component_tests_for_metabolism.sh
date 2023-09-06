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
anvi-migrate *db --migrate-quickly

# generate a temporary directory to store anvi-setup-kegg-data output,
# and remove it immediately to make sure it doesn't exist:
kegg_data_dir=`mktemp -d`
rm -rf $kegg_data_dir

INFO "Setting up KEGG data"
anvi-setup-kegg-data   --mode modules \
                       --kegg-data-dir $kegg_data_dir \
                       --no-progress

## BASIC TESTS
INFO "Estimating metabolism on a single contigs database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         -O single_contigs_db \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE single_contigs_db_modules.txt

INFO "Estimating metabolism using metagenome mode"
anvi-estimate-metabolism -c CONTIGS.db \
                         --metagenome-mode \
                         -O metagenome_mode \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE metagenome_mode_modules.txt

INFO "Estimating metabolism on a collection"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -O collection \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE collection_modules.txt

INFO "Estimating metabolism on a single bin in a collection"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -b E_facealis \
                         -O single_bin \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE single_bin_modules.txt

INFO "Estimating metabolism using a bin IDs file"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -B bin-ids.txt \
                         -O bin_ids_file \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE bin_ids_file_modules.txt

INFO "Estimating metabolism on external genomes"
anvi-estimate-metabolism -e external-genomes.txt \
                         -O external \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE external_modules.txt

INFO "Estimating metabolism on internal genomes"
anvi-estimate-metabolism -i internal-genomes.txt \
                         -O internal \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE internal_modules.txt

INFO "Estimating metabolism on metagenomes file"
anvi-estimate-metabolism -M metagenomes.txt \
                         -O metagenomes \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE metagenomes_modules.txt

INFO "Trying a different module completeness threshold"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --module-completion-threshold 0 \
                         -O nondefault_threshold \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE nondefault_threshold_modules.txt

INFO "Listing output modes for single database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         --list-available-modes \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir


INFO "Listing output modes for multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --list-available-modes \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir

INFO "Listing custom output headers for single database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         --list-available-output-headers \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir

INFO "Listing custom output headers for multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --list-available-output-headers \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir

INFO "Generating long format output files on single database"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --output-modes hits,modules,module_paths,module_steps,modules_custom \
                         --custom-output-headers module,pathwise_module_is_complete,stepwise_module_is_complete,module_name \
                         -O long_format_single \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE long_format_single_hits.txt
SHOW_FILE long_format_single_modules.txt
SHOW_FILE long_format_single_module_paths.txt
SHOW_FILE long_format_single_module_steps.txt
SHOW_FILE long_format_single_modules_custom.txt

INFO "Generating long format output files in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --output-modes hits,modules,module_paths,module_steps,modules_custom \
                         --custom-output-headers module,pathwise_module_is_complete,stepwise_module_is_complete,module_name \
                         -O long_format_multi \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE long_format_multi_hits.txt
SHOW_FILE long_format_multi_modules.txt
SHOW_FILE long_format_multi_module_paths.txt
SHOW_FILE long_format_multi_module_steps.txt
SHOW_FILE long_format_multi_modules_custom.txt

# below we try every single header for modules_custom mode
INFO "Generating modules custom output file with ALL headers (including path-level) on single database"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --output-modes modules_custom \
                         -O modules_custom_single_path \
                         --custom-output-headers module,stepwise_module_is_complete,stepwise_module_completeness,pathwise_module_is_complete,pathwise_module_completeness,enzymes_unique_to_module,unique_enzymes_hit_counts,proportion_unique_enzymes_present,unique_enzymes_context_string,module_name,module_class,module_category,module_subcategory,module_definition,module_substrates,module_products,module_intermediates,gene_caller_ids_in_module,enzyme_hits_in_module,path_id,path,path_completeness,warnings,genome_name \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_single_path_modules_custom.txt

INFO "Generating modules custom output file with ALL headers (including step-level) on single database"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db \
                         --output-modes modules_custom \
                         -O modules_custom_single_step \
                         --custom-output-headers module,stepwise_module_is_complete,stepwise_module_completeness,pathwise_module_is_complete,pathwise_module_completeness,enzymes_unique_to_module,unique_enzymes_hit_counts,proportion_unique_enzymes_present,unique_enzymes_context_string,module_name,module_class,module_category,module_subcategory,module_definition,module_substrates,module_products,module_intermediates,gene_caller_ids_in_module,enzyme_hits_in_module,step_id,step,step_completeness,warnings,genome_name \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_single_step_modules_custom.txt

INFO "Generating modules custom output file with ALL headers (including path-level) in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --output-modes modules_custom \
                         -O modules_custom_multi_path \
                         --custom-output-headers module,stepwise_module_is_complete,stepwise_module_completeness,pathwise_module_is_complete,pathwise_module_completeness,enzymes_unique_to_module,unique_enzymes_hit_counts,proportion_unique_enzymes_present,unique_enzymes_context_string,module_name,module_class,module_category,module_subcategory,module_definition,module_substrates,module_products,module_intermediates,gene_caller_ids_in_module,enzyme_hits_in_module,path_id,path,path_completeness,warnings,genome_name,db_name \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_multi_path_modules_custom.txt

INFO "Generating modules custom output file with ALL headers (including step-level) in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --output-modes modules_custom \
                         -O modules_custom_multi_step \
                         --custom-output-headers module,stepwise_module_is_complete,stepwise_module_completeness,pathwise_module_is_complete,pathwise_module_completeness,enzymes_unique_to_module,unique_enzymes_hit_counts,proportion_unique_enzymes_present,unique_enzymes_context_string,module_name,module_class,module_category,module_subcategory,module_definition,module_substrates,module_products,module_intermediates,gene_caller_ids_in_module,enzyme_hits_in_module,step_id,step,step_completeness,warnings,genome_name,db_name \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_multi_step_modules_custom.txt

INFO "Generating matrix output files in multi mode"
anvi-estimate-metabolism -e external-genomes.txt \
                         --matrix-format \
                         -O matrix_format_multi \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE matrix_format_multi-module_pathwise_completeness-MATRIX.txt
SHOW_FILE matrix_format_multi-module_pathwise_presence-MATRIX.txt
SHOW_FILE matrix_format_multi-module_stepwise_completeness-MATRIX.txt
SHOW_FILE matrix_format_multi-module_stepwise_presence-MATRIX.txt
SHOW_FILE matrix_format_multi-step_completeness-MATRIX.txt
SHOW_FILE matrix_format_multi-enzyme_hits-MATRIX.txt

INFO "Generating JSON output (debug option)"
anvi-estimate-metabolism -c S_islandicus_LS215.db \
                         --get-raw-data-as-json estimation_data \
                         --store-json-without-estimation \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir

INFO "Estimating from JSON output (debug option)"
anvi-estimate-metabolism --estimate-from-json estimation_data.json \
                         -O from_json \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir


## COVERAGE TESTS
INFO "Testing --add-coverage for genome mode"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -O genome_coverage \
                         --add-coverage \
                         --output-modes modules,hits \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE genome_coverage_modules.txt
SHOW_FILE genome_coverage_hits.txt

INFO "Testing --add-coverage flag for a collection"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         -C bins_for_testing \
                         -O collection_coverage \
                         --add-coverage \
                         --output-modes modules,hits \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE collection_coverage_modules.txt
SHOW_FILE collection_coverage_hits.txt

INFO "Testing --add-coverage for metagenome mode"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         --metagenome-mode \
                         -O metagenome_coverage \
                         --add-coverage \
                         --output-modes modules,hits \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE metagenome_coverage_modules.txt
SHOW_FILE metagenome_coverage_hits.txt

INFO "Listing custom output headers with --add-coverage enabled"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         --add-coverage \
                         --list-available-output-headers \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir

INFO "Generating custom output with --add-coverage enabled"
anvi-estimate-metabolism -c CONTIGS.db \
                         -p PROFILE.db \
                         --add-coverage \
                         --output-modes modules_custom \
                         -O modules_custom_coverage \
                         --custom-output-headers module,pathwise_module_is_complete,DAY_15A_gene_coverages,DAY_15A_avg_coverage,DAY_15A_gene_detection,DAY_15A_avg_detection \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_coverage_modules_custom.txt


## COPY NUMBER TESTS
INFO "Testing --add-copy-number in long-format output"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         -O copy_num \
                         --add-copy-number \
                         --output-modes modules,module_paths,module_steps \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE copy_num_modules.txt
SHOW_FILE copy_num_module_paths.txt
SHOW_FILE copy_num_module_steps.txt

INFO "Testing --add-copy-number in matrix output"
anvi-estimate-metabolism -e external-genomes.txt \
                         -O copy_num \
                         --add-copy-number \
                         --matrix-format \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE copy_num-module_pathwise_copy_number-MATRIX.txt
SHOW_FILE copy_num-module_stepwise_copy_number-MATRIX.txt
SHOW_FILE copy_num-step_copy_number-MATRIX.txt

INFO "Listing custom output headers with --add-copy-number enabled"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                         --add-copy-number \
                         --list-available-output-headers \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir

INFO "Generating custom output with --add-copy-number enabled (including path-level headers)"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                          --add-copy-number \
                          --output-modes modules_custom \
                          -O modules_custom_copynum_path \
                          --custom-output-headers module,path,pathwise_copy_number,num_complete_copies_of_path,stepwise_copy_number,per_step_copy_numbers \
                          --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_copynum_path_modules_custom.txt

INFO "Generating custom output with --add-copy-number enabled (including step-level headers)"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db \
                          --add-copy-number \
                          --output-modes modules_custom \
                          -O modules_custom_copynum_step \
                          --custom-output-headers module,step,pathwise_copy_number,stepwise_copy_number,per_step_copy_numbers,step_copy_number \
                          --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE modules_custom_copynum_step_modules_custom.txt


## ENZYMES TXT TESTS
INFO "Testing enzymes txt input (no coverage/detection)"
anvi-estimate-metabolism --enzymes-txt minimal_enzymes_input.txt \
                         -O enzymes_txt \
                         --output-modes hits,module_paths,module_steps,modules \
                         --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE enzymes_txt_modules.txt
SHOW_FILE enzymes_txt_module_paths.txt
SHOW_FILE enzymes_txt_module_steps.txt
SHOW_FILE enzymes_txt_hits.txt

INFO "Testing enzymes txt input (--add-coverage enabled)"
anvi-estimate-metabolism --enzymes-txt cov_det_enzymes_input.txt \
                          -O enzymes_txt_cov \
                          --add-coverage \
                          --output-modes hits,module_paths,module_steps,modules \
                          --no-progress \
                         --kegg-data-dir $kegg_data_dir
SHOW_FILE enzymes_txt_cov_modules.txt
SHOW_FILE enzymes_txt_cov_module_paths.txt
SHOW_FILE enzymes_txt_cov_module_steps.txt
SHOW_FILE enzymes_txt_cov_hits.txt


## ENRICHMENT TESTS
INFO "Testing metabolic enrichment script"
anvi-compute-metabolic-enrichment -M long_format_multi_modules.txt \
                                  -G groups.txt \
                                  -o enrichment.txt \
                                  --no-progress
SHOW_FILE enrichment.txt

INFO "Testing metabolic enrichment with --include-samples-missing-from-groups-txt"
anvi-compute-metabolic-enrichment -M long_format_multi_modules.txt \
                                  -G groups_with_missing_sample.txt \
                                  --include-samples-missing-from-groups-txt \
                                  -o enrichment_ungrouped.txt \
                                  --no-progress
SHOW_FILE enrichment_ungrouped.txt

# clean up
rm -rf $kegg_data_dir
