#!/bin/bash
source 00.sh
set -e

make_structure_db() {
    anvi-gen-structure-database -c test-output/one_contig_five_genes.db \
                                --gene-caller-ids 2,4 \
                                --dump-dir test-output/RAW_MODELLER_OUTPUT \
                                --output-db-path test-output/STRUCTURE.db \
                                --very-fast \
                                --debug \
                                --alignment-fraction-cutoff 0.7 \
                                $thread_controller \
                                --num-models 1
}
# Predict structures with ColabFold instead of MODELLER. This is opt-in: it only runs when the
# COLABFOLD_CONDA_ENV environment variable names a conda environment in which ColabFold is installed
# (anvi'o runs every ColabFold command via `conda run -n $COLABFOLD_CONDA_ENV`). A real run also needs
# internet access (we use --colabfold-msa-server) and is slow without a GPU (~10 min per protein on a
# CPU). When COLABFOLD_CONDA_ENV is unset, the ColabFold parts of the routines are simply skipped.
#
# For example, to exercise the ColabFold engine end-to-end:
#     COLABFOLD_CONDA_ENV=colabfold bash run_component_tests_for_SCVs_SAAVs_structure.sh new
make_structure_db_colabfold() {
    anvi-gen-structure-database -c test-output/one_contig_five_genes.db \
                                --engine colabfold \
                                --colabfold-conda-env "$COLABFOLD_CONDA_ENV" \
                                --colabfold-msa-server \
                                --skip-DSSP \
                                --gene-caller-ids 2 \
                                --num-models 1 \
                                --dump-dir test-output/RAW_COLABFOLD_OUTPUT \
                                --output-db-path test-output/STRUCTURE_COLABFOLD.db \
                                $thread_controller \
                                --debug
}
# Exercise the --only-msa / --only-predict checkpoint that splits a ColabFold run into its (local,
# CPU-heavy) MSA step and its (GPU-heavy) prediction step. This is a further opt-in on top of
# COLABFOLD_CONDA_ENV: it also needs a *local* ColabFold database, because the MSA step cannot be split
# out when using the public MSA server. Point COLABFOLD_DB at the directory you set up with ColabFold's
# `setup_databases.sh`. When COLABFOLD_DB is unset, this routine is skipped.
#
# For example:
#     COLABFOLD_CONDA_ENV=colabfold COLABFOLD_DB=/path/to/colabfold_db \
#         bash run_component_tests_for_SCVs_SAAVs_structure.sh new
make_structure_db_colabfold_checkpoint() {
    # step 1: MSA only. Produces the .a3m files + a checkpoint manifest in --dump-dir, and no db.
    anvi-gen-structure-database -c test-output/one_contig_five_genes.db \
                                --engine colabfold \
                                --colabfold-conda-env "$COLABFOLD_CONDA_ENV" \
                                --colabfold-db "$COLABFOLD_DB" \
                                --skip-DSSP \
                                --gene-caller-ids 2 \
                                --num-models 1 \
                                --dump-dir test-output/COLABFOLD_CHECKPOINT \
                                --only-msa \
                                $thread_controller \
                                --debug

    for expected in test-output/COLABFOLD_CHECKPOINT/genes_of_interest.fa \
                    test-output/COLABFOLD_CHECKPOINT/colabfold_checkpoint.json
    do
        if [ ! -f "$expected" ]; then echo "FAIL: --only-msa did not produce $expected"; exit 1; fi
    done
    if ! ls test-output/COLABFOLD_CHECKPOINT/msas/*.a3m >/dev/null 2>&1
    then echo "FAIL: --only-msa did not produce any MSA (.a3m) files"; exit 1; fi
    if [ -f "test-output/STRUCTURE_COLABFOLD_CHECKPOINT.db" ]
    then echo "FAIL: --only-msa should not have produced a structure database"; exit 1; fi

    # step 2: predict only, resuming from the checkpoint written above. Same contigs-db + genes.
    anvi-gen-structure-database -c test-output/one_contig_five_genes.db \
                                --engine colabfold \
                                --colabfold-conda-env "$COLABFOLD_CONDA_ENV" \
                                --skip-DSSP \
                                --gene-caller-ids 2 \
                                --num-models 1 \
                                --dump-dir test-output/COLABFOLD_CHECKPOINT \
                                --only-predict \
                                --output-db-path test-output/STRUCTURE_COLABFOLD_CHECKPOINT.db \
                                $thread_controller \
                                --debug

    if [ ! -f "test-output/STRUCTURE_COLABFOLD_CHECKPOINT.db" ]
    then echo "FAIL: --only-predict did not produce the structure database"; exit 1; fi
}
# Add a new gene to a ColabFold structure database with anvi-update-structure-database. The scientific
# parameters are read back from the database; only the machine-specific ones are re-supplied. This db
# was built with the public MSA server (see make_structure_db_colabfold), so we re-supply
# --colabfold-msa-server to match. Gene 4 is not yet in the db, so this exercises the 'add' path.
update_structure_db_colabfold() {
    anvi-update-structure-database -c test-output/one_contig_five_genes.db \
                                   -s test-output/STRUCTURE_COLABFOLD.db \
                                   --colabfold-conda-env "$COLABFOLD_CONDA_ENV" \
                                   --colabfold-msa-server \
                                   --gene-caller-ids 4 \
                                   $thread_controller \
                                   --debug
}
# Add a new gene to the local-sourced ColabFold checkpoint database using the same --only-msa /
# --only-predict split as creation. The db was built against a local ColabFold database, so we re-supply
# --colabfold-db for the (local) MSA step. Gene 4 is not yet in the db, so this exercises the 'add' path.
update_structure_db_colabfold_checkpoint() {
    # step 1: MSA only for the new gene. Produces .a3m files + a checkpoint manifest, adds nothing to the db.
    anvi-update-structure-database -c test-output/one_contig_five_genes.db \
                                   -s test-output/STRUCTURE_COLABFOLD_CHECKPOINT.db \
                                   --colabfold-conda-env "$COLABFOLD_CONDA_ENV" \
                                   --colabfold-db "$COLABFOLD_DB" \
                                   --gene-caller-ids 4 \
                                   --dump-dir test-output/COLABFOLD_CHECKPOINT_UPDATE \
                                   --only-msa \
                                   $thread_controller \
                                   --debug

    if ! ls test-output/COLABFOLD_CHECKPOINT_UPDATE/msas/*.a3m >/dev/null 2>&1
    then echo "FAIL: --only-msa (update) did not produce any MSA (.a3m) files"; exit 1; fi

    # step 2: predict only, resuming from the checkpoint above, and add the new gene to the db.
    anvi-update-structure-database -c test-output/one_contig_five_genes.db \
                                   -s test-output/STRUCTURE_COLABFOLD_CHECKPOINT.db \
                                   --colabfold-conda-env "$COLABFOLD_CONDA_ENV" \
                                   --gene-caller-ids 4 \
                                   --dump-dir test-output/COLABFOLD_CHECKPOINT_UPDATE \
                                   --only-predict \
                                   $thread_controller \
                                   --debug
}
gen_var_profile1() {
    anvi-gen-variability-profile -p test-output/SAMPLES-MERGED/PROFILE.db \
                                 -c test-output/one_contig_five_genes.db \
                                 -s test-output/STRUCTURE.db \
                                 -C default \
                                 -b bin1 \
                                 --engine AA \
                                 -o test-output/variability_AA.txt
}
gen_var_profile2() {
    anvi-gen-variability-profile -p test-output/SAMPLES-MERGED/PROFILE.db \
                                 -c test-output/one_contig_five_genes.db \
                                 -s test-output/STRUCTURE.db \
                                 -C default \
                                 -b bin1 \
                                 --engine CDN \
                                 -o test-output/variability_CDN.txt
}
display_structure1() {
    anvi-display-structure -p test-output/SAMPLES-MERGED/PROFILE.db \
              -c test-output/one_contig_five_genes.db \
              -s test-output/STRUCTURE.db \
              --gene-caller-ids 2,4 \
              --debug
}
display_structure2() {
    anvi-display-structure -V test-output/variability_AA.txt \
              -c test-output/one_contig_five_genes.db \
              -s test-output/STRUCTURE.db \
              --debug
}
display_structure3() {
    anvi-display-structure -p test-output/SAMPLES-MERGED/PROFILE.db \
              -c test-output/one_contig_five_genes.db \
              -s test-output/EXTERNAL_STRUCTURE.db \
              --gene-caller-ids 2,4 \
              --debug
}
display_structure_colabfold() {
    anvi-display-structure -p test-output/SAMPLES-MERGED/PROFILE.db \
              -c test-output/one_contig_five_genes.db \
              -s test-output/STRUCTURE_COLABFOLD.db \
              --gene-caller-ids 2 \
              --debug
}

make_routine() {
    INFO "anvi-gen-structure-database with DSSP"
    make_structure_db

    INFO "anvi-update-structure-database"
    anvi-update-structure-database -c test-output/one_contig_five_genes.db -s test-output/STRUCTURE.db --gene-caller-ids 2 --rerun

    INFO "anvi-export-structures"
    anvi-export-structures -o test-output/exported_pdbs -s test-output/STRUCTURE.db

    INFO "make db with external structures"
    cat <<EOF > external_structures
gene_callers_id	path
2	test-output/exported_pdbs/gene_2.pdb
EOF
    anvi-gen-structure-database -o test-output/EXTERNAL_STRUCTURE.db --external-structures external_structures -c test-output/one_contig_five_genes.db --debug

    INFO "update db with external structures"
    cat <<EOF > external_structures
gene_callers_id	path
2	test-output/exported_pdbs/gene_2.pdb
4	test-output/exported_pdbs/gene_4.pdb
EOF
    anvi-update-structure-database -s test-output/EXTERNAL_STRUCTURE.db --external-structures external_structures -c test-output/one_contig_five_genes.db --rerun --debug

    # ColabFold is only exercised when the user points the test at a conda environment in which
    # ColabFold is installed (see make_structure_db_colabfold above). Otherwise it is skipped.
    if [ -n "$COLABFOLD_CONDA_ENV" ]
    then
        INFO "anvi-gen-structure-database with ColabFold (conda env: $COLABFOLD_CONDA_ENV)"
        make_structure_db_colabfold

        INFO "anvi-update-structure-database with ColabFold"
        update_structure_db_colabfold

        # exercise the --only-msa / --only-predict checkpoint only when a local ColabFold database is
        # also available (the MSA step cannot be split out when using the public server)
        if [ -n "$COLABFOLD_DB" ]
        then
            INFO "anvi-gen-structure-database ColabFold --only-msa / --only-predict checkpoint (local db: $COLABFOLD_DB)"
            make_structure_db_colabfold_checkpoint

            INFO "anvi-update-structure-database ColabFold --only-msa / --only-predict checkpoint (local db: $COLABFOLD_DB)"
            update_structure_db_colabfold_checkpoint
        fi
    fi
}

display_routine() {
    INFO "anvi-gen-variability-profile --engine AA"
    gen_var_profile1

    INFO "anvi-gen-variability-profile --engine CDN"
    gen_var_profile2

    INFO "anvi-display-structure with profile database"
    display_structure1

    INFO "anvi-display-structure with variability"
    display_structure2

    INFO "anvi-display-structure with external structure database"
    display_structure3

    # display the ColabFold structures whenever they were generated (i.e. the ColabFold engine was
    # exercised in a `new`/`make` run). Displaying needs neither ColabFold nor a conda env, so we gate
    # on the database existing rather than on COLABFOLD_CONDA_ENV.
    if [ -f "test-output/STRUCTURE_COLABFOLD.db" ]
    then
        INFO "anvi-display-structure with ColabFold structure database"
        display_structure_colabfold
    fi
}


if [ $# -eq 0  ]
then
      echo "\

        No arguments supplied. You must call this script with one argument: 'new', 'make',
        or 'display':

        'new'     : generate short reads, map them to the contig, gen contigs.db, make structure db, profile varability, open interactive.
        'make'    : make structure db, profile varability, open interactive.
        'display' : profile varability, open interactive.

        You can pass the number of threads as an optional second argument (defaults to 1), e.g.:

            bash run_component_tests_for_SCVs_SAAVs_structure.sh new 8

        To additionally exercise the ColabFold engine, set the COLABFOLD_CONDA_ENV environment
        variable to the name of a conda environment in which ColabFold is installed, e.g.:

            COLABFOLD_CONDA_ENV=colabfold bash run_component_tests_for_SCVs_SAAVs_structure.sh new

        This requires internet access and is slow without a GPU (~10 min per protein on a CPU).

        To also exercise the --only-msa / --only-predict checkpoint, additionally set COLABFOLD_DB to
        the local ColabFold database directory you set up with ColabFold's setup_databases.sh, e.g.:

            COLABFOLD_CONDA_ENV=colabfold COLABFOLD_DB=/path/to/colabfold_db \\
                bash run_component_tests_for_SCVs_SAAVs_structure.sh new
        "
        exit -1
fi

if [ $# -gt 2  ]
then
      echo "
        This script expects at most two arguments: the mode ('new', 'make', or 'display') and,
        optionally, the number of threads to use (defaults to 1).
        "
        exit -1
fi

# optional second argument: how many threads the anvi-gen-structure-database steps should use. When it
# is omitted the commands run single-threaded (see 00.sh for the same pattern in other component
# tests). Threading mostly helps the MODELLER runs and ColabFold's (CPU-heavy) MSA step.
if [ -z "$2" ]
then
    thread_controller=""
else
    thread_controller="--num-threads $2"
fi

if [ $1 = "new"  ]
then
    INFO "Creating the output directory"
    # change directory and clean the old mess if it exists
    cd sandbox
    rm -rf test-output
    mkdir test-output
    
    INFO "Anvo'o version"
    anvi-profile --version
    
    INFO "Generating a Bowtie2 ref"
    bowtie2-build mock_data_for_structure/one_contig_five_genes.fa test-output/one_contig_five_genes.build

    INFO "Generating anvi'o contigs database"
    anvi-gen-contigs-database -f mock_data_for_structure/one_contig_five_genes.fa -o test-output/one_contig_five_genes.db -n "5 genes concatenated"

    for sample in 01 02 03 04 05
    do
        INFO "Generating short reads for sample $sample"
        anvi-script-gen-short-reads mock_data_for_structure/$sample.ini --output-file-path test-output/$sample.fa

        INFO "Mapping short reads to the ref"
        ../misc/bowtie_batch_single_fasta.sh test-output/$sample.fa test-output/$sample test-output/one_contig_five_genes.build

        INFO "Profiling the BAM file"
        anvi-profile -i test-output/$sample.bam -c test-output/one_contig_five_genes.db -o test-output/$sample-PROFILE -M 0 --profile-SCVs
    done

    INFO "Merging all"
    anvi-merge test-output/*PROFILE/PROFILE.db -c test-output/one_contig_five_genes.db -o test-output/SAMPLES-MERGED

    INFO "Defining a collection and bin"
    anvi-import-collection mock_data_for_structure/collection.txt -c test-output/one_contig_five_genes.db -p test-output/SAMPLES-MERGED/PROFILE.db -C default

    INFO "Importing additional layers data"
    anvi-import-misc-data mock_data_for_structure/additional_layers_data.txt -p test-output/SAMPLES-MERGED/PROFILE.db --target-data-table layers

    make_routine
    display_routine

####################################################################################

elif [ $1 = "make"  ]
then
    cd sandbox

    if [ ! -f "test-output/one_contig_five_genes.db"  ]
    then
      echo "
        You asked to continue with the previously generated files,
        but the contigs database is not there. Maybe you should start
        from scratch by re-running this script with the parameter 'new'.
        "
        exit -1
    else
        INFO "Attempting to continue with the previously generated files"
    fi

    rm -rf test-output/STRUCTURE.db
    rm -rf test-output/EXTERNAL_STRUCTURE.db
    rm -rf test-output/RAW_MODELLER_OUTPUT
    rm -rf test-output/exported_pdbs
    rm -rf test-output/STRUCTURE_COLABFOLD.db
    rm -rf test-output/RAW_COLABFOLD_OUTPUT
    rm -rf test-output/STRUCTURE_COLABFOLD_CHECKPOINT.db
    rm -rf test-output/COLABFOLD_CHECKPOINT
    rm -rf test-output/COLABFOLD_CHECKPOINT_UPDATE

    make_routine
    display_routine

    echo
    echo

elif [ $1 = "display"  ]
then
    cd sandbox

    if [ ! -f "test-output/STRUCTURE.db"  ]
    then
      echo "
        You asked to continue with the previously generated files,
        but the structure database is not there. Maybe you should start
        from scratch by re-running this script with the parameter 'new'.
        "
        exit -1
    else
        INFO "Attempting to continue with the previously generated files"
    fi

    rm -rf test-output/variability_CDN.txt
    rm -rf test-output/variability_AA.txt

    display_routine

    echo
    echo

else
      echo "
        Unknown parameter $1 :/ Try 'new', or 'make', or 'display'.
        "
        exit -1
fi

