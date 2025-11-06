#!/bin/bash
source 00.sh

SETUP_WITH_OUTPUT_DIR $1 $2 $3

python_script=`readlink -f run_component_tests_for_reaction_network`

INFO "Checking for the required KEGG database set up by anvi'o in a default location"
${python_script} --check-default-kegg-database

INFO "Setting up the reaction network analysis directory"
mkdir -p ${output_dir}/
# These databases should already contain KO annotations.
cp ${files}/mock_data_for_pangenomics/*.db ${output_dir}/
cp ${files}/mock_data_for_pangenomics/external-genomes.txt ${output_dir}/
cp ${files}/example_description.md ${output_dir}/
cd ${output_dir}/

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

use_default_modelseed_db=`${python_script} --check-default-modelseed-database`
if [ "${use_default_modelseed_db}" == "True" ]
then
    INFO "Using the ModelSEED Biochemistry database already set up by anvi'o in a default location"
else
    INFO "Setting up the ModelSEED Biochemistry database in a temporary directory (a permanent ModelSEED database can be installed in the default location with 'anvi-setup-modelseed-database')"
    data_dir=`mktemp -d`
    anvi-setup-modelseed-database --dir ${data_dir}
    modelseed_data_dir=${data_dir}/MODELSEED
fi

INFO "Testing a genomic reaction network generated from a contigs database"
args=()
args+=( "--contigs-db" "E_faecalis_6240.db" )
args+=( "--test-dir" ${output_dir} )
if [ "${use_default_modelseed_db}" == "False" ]
then
    args+=( "--modelseed-dir" ${modelseed_data_dir} )
fi
args+=( "--no-progress" )
${python_script} "${args[@]}"

INFO "Exporting the genomic reaction network to a file"
anvi-get-metabolic-model-file --contigs-db E_faecalis_6240.db \
                              --output-file E_faecalis_6240-network.json

INFO "Testing a pangenomic reaction network generated from the pan and genomes storage databases"
args=()
args+=( "--pan-db" "TEST/TEST-PAN.db" )
args+=( "--genomes-storage" "TEST-GENOMES.db" )
args+=( "--test-dir" ${output_dir} )
if [ ${use_default_modelseed_db} == "False" ]
then
    args+=( "--modelseed-dir" ${modelseed_data_dir} )
fi
args+=( "--no-progress" )
${python_script} "${args[@]}"

INFO "Exporting the pangenomic reaction network to a file"
anvi-get-metabolic-model-file --pan-db TEST/TEST-PAN.db \
                              --genomes-storage TEST-GENOMES.db \
                              --record-genomes \
                              --output-file TEST-PAN-network.json

if [ ${use_default_modelseed_db} == "False" ]
then
    INFO "Removing the temporary ModelSEED Biochemistry database directory"
    rm -rf ${modelseed_data_dir}
fi
