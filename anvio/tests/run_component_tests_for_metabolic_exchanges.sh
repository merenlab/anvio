#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

# If you don't want to use the default KEGG data directory for testing, you should
# run `export kegg_data_dir=/path/to/data/dir/you/want` in your terminal before starting the self-test
# same thing goes for the default MODELSEED data directory (run `export modelseed_dir=/path/to/dir/you/want` before starting)

if [ x"${kegg_data_dir}" == "x" ]; then
 	INFO "Checking for KEGG database in default location"
    # Here we use Sam's clever function to check for default KEGG data dir
    rn_python_script=`readlink -f run_component_tests_for_reaction_network`
    ${rn_python_script} --check-default-kegg-database
    source_dir=$(dirname -- "$( readlink -f -- "$0"; )";)
    kegg_data_dir=${source_dir%/tests}/data/misc/KEGG
    INFO "Using default KEGG data directory: $kegg_data_dir"
else
 	INFO "Using manually-provided KEGG data directory: $kegg_data_dir"
fi

if [ x"${modelseed_dir}" == "x" ]; then
 	INFO "Checking for ModelSEED database in default location"
    rn_python_script=`readlink -f run_component_tests_for_reaction_network`
    ${rn_python_script} --check-default-modelseed-database
    source_dir=$(dirname -- "$( readlink -f -- "$0"; )";)
    modelseed_dir=${source_dir%/tests}/data/misc/MODELSEED
    INFO "Using default MODELSEED data directory: $modelseed_dir"
else
 	INFO "Using manually-provided MODELSEED data directory: $modelseed_dir"
fi

INFO "Setting up the metabolic-exchanges test directory"
mkdir -p $output_dir/metabolic-exchanges_test
cp $files/data/genomes/bacteria/*.db                    $output_dir/metabolic-exchanges_test
cp $files/data/genomes/archaea/*.db                     $output_dir/metabolic-exchanges_test
cp $files/data/input_files/external-genomes.txt         $output_dir/metabolic-exchanges_test
cp $files/data/input_files/genome-pairs.txt             $output_dir/metabolic-exchanges_test

cd $output_dir/metabolic-exchanges_test

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Annotating all databases with KOfams"
for db in *.db
do
    annotation_hash=$(sqlite3 $db "select value from self where key='modules_db_hash'")
    if [ -z "$annotation_hash" ]; then 
        INFO "Annotating $db"
        anvi-run-kegg-kofams    -c $db \
                                --kegg-data-dir $kegg_data_dir \
                                $thread_controller \
                                --just-do-it
    else
        INFO "KOfam annotations already found in $db"
    fi
done

INFO "Creating reaction networks for all databases"
for db in *.db
do
    anvi-reaction-network -c $db \
                          --kegg-dir $kegg_data_dir \
                          --modelseed-dir $modelseed_dir
done

INFO "Single genome pair comparison with default parameters"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt
SHOW_FILE B_vs_P-evidence.txt
SHOW_FILE B_vs_P-unique-compounds.txt

INFO "Single genome pair comparison with --no-pathway-walk"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P \
                                 --force-overwrite \
                                 --no-pathway-walk
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt

INFO "Single genome pair comparison with --pathway-walk-only"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P \
                                 --force-overwrite \
                                 --pathway-walk-only
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt

INFO "Single genome pair comparison with --exclude-pathway-maps"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P \
                                 --force-overwrite \
                                 --exclude-pathway-maps 00470,00195,00542,00190,00541,00543,00511
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt

INFO "Single genome pair comparison with --add-reactions-to-output"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P \
                                 --force-overwrite \
                                 --add-reactions-to-output
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt
SHOW_FILE B_vs_P-unique-compounds.txt

INFO "Single genome pair comparison with --maximum-gaps 2 (might take a while to finish pathway map walks)"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P \
                                 --force-overwrite \
                                 --maximum-gaps 2
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt

INFO "Single genome pair comparison with --use-equivalent-amino-acids"
anvi-predict-metabolic-exchanges -c1 B_thetaiotamicron_VPI-5482.db \
                                 -c2 P_marinus_CCMP1375.db \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O B_vs_P \
                                 --force-overwrite \
                                 --use-equivalent-amino-acids
SHOW_FILE B_vs_P-potentially-exchanged-compounds.txt

INFO "MULTI-mode prediction: genome pairs file"
anvi-predict-metabolic-exchanges -e external-genomes.txt \
                                 --genome-pairs genome-pairs.txt \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O genome_pairs
SHOW_FILE genome_pairs-potentially-exchanged-compounds.txt
SHOW_FILE genome_pairs-evidence.txt
SHOW_FILE genome_pairs-unique-compounds.txt

INFO "MULTI-mode prediction: all-vs-all"
anvi-predict-metabolic-exchanges -e external-genomes.txt \
                                 --kegg-data-dir $kegg_data_dir \
                                 --modelseed-data-dir $modelseed_dir \
                                 $thread_controller \
                                 -O all_vs_all
SHOW_FILE all_vs_all-potentially-exchanged-compounds.txt
SHOW_FILE all_vs_all-evidence.txt
SHOW_FILE all_vs_all-unique-compounds.txt