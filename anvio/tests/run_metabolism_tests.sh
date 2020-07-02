#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the metabolism test directory"
mkdir $output_dir/metabolism_test
cp $files/data/genomes/bacteria/*.db                    $output_dir/metabolism_test
cp $files/data/genomes/archaea/*.db                     $output_dir/metabolism_test
cp $files/data/metagenomes/human_gut/IGD_SUBSET/*.db    $output_dir/metabolism_test
cp $files/data/input_files/*.txt                        $output_dir/metabolism_test
cd $output_dir/metabolism_test

INFO "Estimating metabolism on a single contigs database"
anvi-estimate-metabolism -c B_thetaiotamicron_VPI-5482.db -O single_contigs_db

INFO "Estimating metabolism using metagenome mode"
anvi-estimate-metabolism -c CONTIGS.db --metagenome-mode -O metagenome_mode

INFO "Estimating metabolism on a collection"
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C bins_for_testing -O collection

INFO "Estimating metabolism on a single bin in a collection"
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C bins_for_testing -b E_facealis -O single_bin

INFO "Estimating metabolism using a bin IDs file"
anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C bins_for_testing -B bin-ids.txt -O bin_ids_file

INFO "Estimating metabolism on external genomes"
anvi-estimate-metabolism -e external-genomes.txt -O external

INFO "Estimating metabolism on internal genomes"
anvi-estimate-metabolism -i internal-genomes.txt -O internal

INFO "Estimating metabolism on metagenomes file"
anvi-estimate-metabolism -M metagenomes.txt -O metagenomes

INFO "Trying a different module completeness threshold"
anvi-estimate-metabolism -c P_marinus_CCMP1375.db --module-completion-threshold 0 -O nondefault_threshold

INFO "Generating long format output files"
anvi-estimate-metabolism -e external-genomes.txt --kegg-output-modes kofam_hits,modules,custom --custom-output-headers kegg_module,module_is_complete,module_name -O long_format

INFO "Generating matrix output files"
anvi-estimate-metabolism -i internal-genomes.txt --matrix-format -O matrix_format

INFO "Generating JSON output (debug option)"
anvi-estimate-metabolism -c S_islandicus_LS215.db --get-raw-data-as-json estimation_data --store-json-without-estimation

INFO "Estimating from JSON output (debug option)"
anvi-estimate-metabolism -c S_islandicus_LS215.db --estimate-from-json estimation_data.json -O from_json
