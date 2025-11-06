#!/bin/bash
source 00.sh

SETUP_WITH_OUTPUT_DIR $1 $2 $3

rn_python_script=`readlink -f run_component_tests_for_reaction_network`

INFO "Checking for the required KEGG database set up by anvi'o in a default location"
${rn_python_script} --check-default-kegg-database

INFO "Setting up the KEGG mapping analysis directory"
mkdir -p ${output_dir}/
# These databases should already contain KO annotations.
cp ${files}/mock_data_for_pangenomics/*.db ${output_dir}/
cp ${files}/mock_data_for_pangenomics/external-genomes.txt ${output_dir}/
cp ${files}/example_description.md ${output_dir}/
cp ${files}/mock_data_for_pangenomics/group-information.txt ${output_dir}/pan-group-information.txt
awk -F'\t' 'NR==1 {print $0} NR>1 {$1=$1".db"; print $0}' OFS='\t' \
${output_dir}/pan-group-information.txt > ${output_dir}/contigs-db-group-information.txt
cd ${output_dir}/
mkdir contigs_db_kos
mkdir contigs_dbs_kos_count
mkdir pan_db_kos_genome_count_emphasize_shared
mkdir pan_db_kos_genome_count_emphasize_unshared
mkdir pan_db_kos_presence_absence

INFO "Migrating all databases"
anvi-migrate *db --migrate-quickly

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e external-genomes.txt -o TEST-GENOMES.db --no-progress

INFO "Running the pangenome analysis with default parameters"
anvi-pan-genome -g TEST-GENOMES.db \
                -o TEST/ \
                -n TEST \
                --use-ncbi-blast \
                --description example_description.md \
                --no-progress \
                ${thread_controller}

use_default_modelseed_db=`${rn_python_script} --check-default-modelseed-database`
if [ "${use_default_modelseed_db}" == "True" ]
then
    INFO "Using the ModelSEED Biochemistry database already set up by anvi'o in a default location"
else
    INFO "Setting up the ModelSEED Biochemistry database in a temporary directory (a permanent \
ModelSEED database can be installed in the default location with \
'anvi-setup-modelseed-database')"
    data_dir=`mktemp -d`
    anvi-setup-modelseed-database --dir ${data_dir}
    modelseed_data_dir=${data_dir}/MODELSEED
fi

INFO "Generating a pangenomic reaction network"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
if [ ${use_default_modelseed_db} == "False" ]
then
    args+=( "--modelseed-dir" ${modelseed_data_dir} )
fi
args+=( "--no-progress" )
anvi-reaction-network "${args[@]}"

pathway_numbers=( "00010" "01100" "01200" )

INFO "Testing mapping KOs from a genomic contigs database"
args=()
args+=( "--contigs-db" "E_faecalis_6240.db" )
args+=( "--output-dir" ${output_dir}/contigs_db_kos )
args+=( "--name-files" )
args+=( "--categorize-files" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" "05310" )
args+=( "--draw-bare-maps" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying database counts emphasizing \
shared reactions, drawing map grids"
args=()
args+=( "--external-genomes" "external-genomes.txt" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_count )
args+=( "--categorize-files" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--draw-grid" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying database membership"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db"
)
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_membership )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying group membership where KOs \
are in any databases in the group, emphasizing shared reactions, drawing map grids and map files \
for each group"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db" \
"E_faecalis_6557.db" \
"E_faecalis_6563.db" \
)
args+=( "--groups-txt" ${output_dir}/contigs-db-group-information.txt)
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_group_membership_min_threshold )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying group membership where KOs \
are in all databases in the group, emphasizing shared reactions, drawing map grids and map files \
for a subset of the groups"
args=()
args+=( "--contigs-dbs" \
"E_faecalis_6240.db" \
"E_faecalis_6255.db" \
"E_faecalis_6512.db" \
"E_faecalis_6557.db" \
"E_faecalis_6563.db" \
)
args+=( "--groups-txt" ${output_dir}/contigs-db-group-information.txt)
args+=( "--group-threshold" "1" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_group_membership_max_threshold )
args+=( "--draw-individual-files" "G2" )
args+=( "--draw-grid" "G2" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying genome counts emphasizing shared \
reactions"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_genome_count_emphasize_shared )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying genome counts emphasizing \
unshared reactions, drawing map grids and map files for each genome"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_genome_count_emphasize_unshared )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--colormap" "plasma" "0.1" "0.9")
args+=( "--reverse-overlay" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying presence/absence"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_presence_absence )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--set-color" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying group count where KOs are in \
any genomes in the group, emphasizing shared reactions, drawing map grids and map files for a \
subset of groups"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--groups-txt" "pan-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_group_count_emphasize_shared )
args+=( "--categorize-files" )
args+=( "--draw-individual-files" "G1" )
args+=( "--draw-grid" "G1" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying group membership where KOs are \
in any genomes in the group, emphasizing unshared reactions, drawing map grids for each group"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--groups-txt" "pan-group-information.txt" )
args+=( "--group-threshold" "0" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_group_membership_emphasize_unshared )
args+=( "--draw-grid" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--colormap" "plasma" "0.1" "0.9")
args+=( "--colormap-scheme" "by_count" )
args+=( "--reverse-overlay" )
args+=( "--group-colormap" "plasma" "0.1" "0.9" )
args+=( "--group-reverse-overlay" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"
