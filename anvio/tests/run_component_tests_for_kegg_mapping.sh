#!/bin/bash
source 00.sh

SETUP_WITH_OUTPUT_DIR $1 $2 $3

rn_python_script=`readlink -f run_component_tests_for_reaction_network`

INFO "Setting up the KEGG mapping analysis directory"
mkdir -p ${output_dir}/
# These databases should already contain KO annotations.
cp ${files}/mock_data_for_pangenomics/*.db ${output_dir}/
cp ${files}/mock_data_for_pangenomics/external-genomes.txt ${output_dir}/
cp ${files}/example_description.md ${output_dir}/
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
    INFO "Setting up the ModelSEED Biochemistry database in a temporary directory (a permanent ModelSEED database can be installed in the default location with 'anvi-setup-modelseed-database')"
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
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from multiple contigs databases, displaying database counts \
emphasizing shared reactions, drawing grid maps"
args=()
args+=( "--external-genomes" "external-genomes.txt" )
args+=( "--output-dir" ${output_dir}/contigs_dbs_kos_count )
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

INFO "Testing mapping KOs from a pangenomic database, displaying genome counts \
emphasizing shared reactions"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_genome_count_emphasize_shared )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--no-progress" )
anvi-draw-kegg-pathways "${args[@]}"

INFO "Testing mapping KOs from a pangenomic database, displaying genome counts \
emphasizing unshared reactions, drawing grid maps and maps for each genome"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--output-dir" ${output_dir}/pan_db_kos_genome_count_emphasize_unshared )
args+=( "--draw-individual-files" )
args+=( "--draw-grid" )
args+=( "--ko" )
args+=( "--pathway-numbers" "${pathway_numbers[@]}" )
args+=( "--colormap" "magma" )
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
