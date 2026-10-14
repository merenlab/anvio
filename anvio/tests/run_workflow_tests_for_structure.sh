#!/bin/bash
# Workflow test for structure database generation in the contigs workflow.
#
# The contigs workflow first runs without structure prediction. Then, as a user would, the test picks
# a few genes of interest from the resulting contigs databases, adds them to the fasta.txt under the
# `structure_genes_of_interest` column, turns on `anvi_gen_structure_database`, and runs the workflow
# again, which should only run the structure steps.
#
# The default MODELLER engine needs internet (it downloads its templates from the RCSB PDB) and may find
# no template for a given gene, so the `proteins` table is asserted strictly while whether a structure
# was actually produced is treated softly.
#
# ColabFold is opt-in: set COLABFOLD_CONDA_ENV to the conda environment ColabFold is installed in (a real
# run also needs a GPU). With COLABFOLD_DB also set to a local ColabFold database, the test runs the
# split MSA (CPU) / prediction (GPU) steps; without it, it runs a single step with the public MSA server.
#
#     COLABFOLD_CONDA_ENV=colabfold COLABFOLD_DB=/path/to/db bash run_workflow_tests_for_structure.sh

source 00.sh
set -e

assert_eq() {
    # assert_eq <actual> <expected> <message>
    if [ "$1" != "$2" ]; then
        echo ""
        echo "ASSERTION FAILED: $3"
        echo "  expected: '$2'"
        echo "  actual:   '$1'"
        exit 1
    fi
    echo "  OK: $3 (== '$2')"
}

proteins_count() { sqlite3 "$1" "SELECT COUNT(*) FROM proteins;"; }
structures_count() { sqlite3 "$1" "SELECT COUNT(*) FROM proteins WHERE has_structure = 1;"; }

# update_config <config> <rule> <JSON object of params>
update_config() {
    python -c "
import json, sys
config = json.load(open('$1'))
config['$2'].update(json.loads(sys.argv[1]))
json.dump(config, open('$1', 'w'), indent=4)" "$3"
}

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the structure workflow directory"
mkdir $output_dir/workflow_test
cp $files/mock_data_for_pangenomics/E_faecalis_6240.db $output_dir/workflow_test/
cp $files/mock_data_for_pangenomics/E_faecalis_6255.db $output_dir/workflow_test/
cd $output_dir/workflow_test

anvi-migrate *.db --migrate-quickly --quiet
anvi-export-contigs -c E_faecalis_6240.db -o E_faecalis_6240.fa
anvi-export-contigs -c E_faecalis_6255.db -o E_faecalis_6255.fa
rm -rf E_fa*.db

echo -e "name\tpath"                        >  fasta.txt
echo -e "E_faecalis_6240\tE_faecalis_6240.fa" >> fasta.txt
echo -e "E_faecalis_6255\tE_faecalis_6255.fa" >> fasta.txt

INFO "Creating a default config for contigs workflow (annotation off, to keep things fast)"
anvi-run-workflow -w contigs --get-default-config config.json
for rule in anvi_run_hmms anvi_run_ncbi_cogs anvi_run_scg_taxonomy anvi_run_kegg_kofams; do
    update_config config.json $rule '{"run": false}'
done

INFO "Running contigs workflow without structure prediction"
anvi-run-workflow -w contigs -c config.json

INFO "Picking two complete genes of interest per contigs database"
for genome in E_faecalis_6240 E_faecalis_6255; do
    sqlite3 02_CONTIGS/$genome.db "SELECT gene_callers_id FROM genes_in_contigs WHERE partial = 0 \
                                   ORDER BY gene_callers_id LIMIT 2;" > $genome-genes-of-interest.txt
done

echo -e "name\tpath\tstructure_genes_of_interest"                                       >  fasta.txt
echo -e "E_faecalis_6240\tE_faecalis_6240.fa\tE_faecalis_6240-genes-of-interest.txt" >> fasta.txt
echo -e "E_faecalis_6255\tE_faecalis_6255.fa\tE_faecalis_6255-genes-of-interest.txt" >> fasta.txt

INFO "Turning on structure prediction with MODELLER"
cp config.json config-modeller.json
update_config config-modeller.json anvi_gen_structure_database '{"run": true, "--very-fast": true, "--num-models": 1}'

INFO "Listing dependencies for the structure steps"
anvi-run-workflow -w contigs -c config-modeller.json --list-dependencies

INFO "Running contigs workflow with structure prediction (only the structure steps should run)"
anvi-run-workflow -w contigs -c config-modeller.json

for genome in E_faecalis_6240 E_faecalis_6255; do
    assert_eq "$(proteins_count 03_STRUCTURE/$genome-STRUCTURE.db)" "2" "$genome proteins rows (MODELLER)"
    echo "  INFO: $genome has $(structures_count 03_STRUCTURE/$genome-STRUCTURE.db) structure(s) from MODELLER"
done

INFO "Annotating the contigs databases must not trigger a new structure prediction"
mtime_before=$(stat -c %Y 03_STRUCTURE/E_faecalis_6240-STRUCTURE.db 2>/dev/null || stat -f %m 03_STRUCTURE/E_faecalis_6240-STRUCTURE.db)
update_config config-modeller.json anvi_run_hmms '{"run": true}'
anvi-run-workflow -w contigs -c config-modeller.json
mtime_after=$(stat -c %Y 03_STRUCTURE/E_faecalis_6240-STRUCTURE.db 2>/dev/null || stat -f %m 03_STRUCTURE/E_faecalis_6240-STRUCTURE.db)
assert_eq "$mtime_after" "$mtime_before" "structure database untouched by annotation"

if [ -z "$COLABFOLD_CONDA_ENV" ]; then
    INFO "COLABFOLD_CONDA_ENV is not set, skipping the ColabFold workflow tests"
    exit 0
fi

cp config.json config-colabfold.json
update_config config-colabfold.json output_dirs '{"STRUCTURE_DIR": "03_STRUCTURE_COLABFOLD"}'

if [ -n "$COLABFOLD_DB" ]; then
    INFO "Turning on structure prediction with ColabFold, split into MSA and prediction steps"
    update_config config-colabfold.json anvi_gen_structure_database "{\"run\": true, \"--engine\": \"colabfold\", \
        \"--colabfold-conda-env\": \"$COLABFOLD_CONDA_ENV\", \"--colabfold-db\": \"$COLABFOLD_DB\", \
        \"split_msa_and_predict\": true}"
else
    INFO "Turning on structure prediction with ColabFold, using the public MSA server"
    update_config config-colabfold.json anvi_gen_structure_database "{\"run\": true, \"--engine\": \"colabfold\", \
        \"--colabfold-conda-env\": \"$COLABFOLD_CONDA_ENV\", \"--colabfold-msa-server\": true}"
fi

INFO "Running contigs workflow with ColabFold"
anvi-run-workflow -w contigs -c config-colabfold.json

for genome in E_faecalis_6240 E_faecalis_6255; do
    assert_eq "$(proteins_count 03_STRUCTURE_COLABFOLD/$genome-STRUCTURE.db)" "2" "$genome proteins rows (ColabFold)"
    assert_eq "$(structures_count 03_STRUCTURE_COLABFOLD/$genome-STRUCTURE.db)" "2" "$genome structures (ColabFold)"
done

if [ -n "$COLABFOLD_DB" ]; then
    assert_eq "$(ls -d 03_STRUCTURE_COLABFOLD/*-COLABFOLD-MSA 2>/dev/null | wc -l | tr -d ' ')" "0" "temporary MSA directories removed"
fi
