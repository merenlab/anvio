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
                                --num-threads 2 \
                                --num-models 1
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


if [ $# -eq 0  ]
then
      echo "\

        No arguments supplied. You must call this script with one argument: 'new', 'make',
        or 'display':

        'new'     : generate short reads, map them to the contig, gen contigs.db, make structure db, profile varability, open interactive.
        'make'    : make structure db, profile varability, open interactive.
        'display' : profile varability, open interactive.
        "
        exit -1
fi

if [ $# -gt 1  ]
then
      echo "
        This scripts expect only one argument ('new', or 'continue').
        "
        exit -1
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

    INFO "anvi-gen-structure-database with DSSP"
    make_structure_db

    INFO "anvi-gen-variability-profile --engine AA"
    gen_var_profile1

    INFO "anvi-gen-variability-profile --engine CDN"
    gen_var_profile2

    INFO "anvi-display-structure with profile database"
    display_structure1

    INFO "anvi-display-structure with variability"
    display_structure2


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
    rm -rf test-output/RAW_MODELLER_OUTPUT
    rm -rf test-output/exported_pdbs

    INFO "anvi-gen-structure-database with DSSP"
    make_structure_db

    INFO "anvi-update-structure-database"
    anvi-update-structure-database -c test-output/one_contig_five_genes.db -s test-output/STRUCTURE.db --gene-caller-ids 2 --rerun

    INFO "anvi-export-structures"
    anvi-export-structures -o test-output/exported_pdbs -s test-output/STRUCTURE.db

    INFO "anvi-gen-variability-profile --engine AA"
    gen_var_profile1

    INFO "anvi-gen-variability-profile --engine CDN"
    gen_var_profile2

    INFO "anvi-display-structure with profile database"
    display_structure1

    INFO "anvi-display-structure with variability"
    display_structure2

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

    INFO "anvi-gen-variability-profile --engine AA"
    gen_var_profile1

    INFO "anvi-gen-variability-profile --engine CDN"
    gen_var_profile2

    INFO "anvi-display-structure with profile and contigs databases"
    display_structure1

    INFO "anvi-display-structure with variability"
    display_structure2

    echo
    echo

else
      echo "
        Unknown parameter $1 :/ Try 'new', or 'make', or 'display'.
        "
        exit -1
fi

