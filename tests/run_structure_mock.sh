#!/bin/bash
source 00.sh
set -e

if [ $# -eq 0  ]
then
      echo "\

        No arguments supplied. You must call this script with one argument: 'new',
        or 'continue':

        'new'     : generate short reads, map them to the contig, gen contigs.db, run profiling.
        'continue': re-run profiling.
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

    for sample in 01 02 03
    do
        INFO "Generating short reads for sample $sample"
        anvi-script-gen-short-reads mock_data_for_structure/$sample.ini --output-file-path test-output/$sample.fa

        INFO "Mapping short reads to the ref"
        ../misc/bowtie_batch_single_fasta.sh test-output/$sample.fa test-output/$sample test-output/one_contig_five_genes.build

        INFO "Profiling the BAM file"
        anvi-profile -i test-output/$sample.bam -c test-output/one_contig_five_genes.db -o test-output/$sample-PROFILE -M 0 --profile-AA-frequencies
    done

    INFO "Merging all"
    anvi-merge test-output/*PROFILE/PROFILE.db -c test-output/one_contig_five_genes.db -o test-output/SAMPLES-MERGED --skip-concoct-binning

    INFO "Defining a collection and bin"
    anvi-import-collection mock_data_for_structure/collection.txt -c test-output/one_contig_five_genes.db -p test-output/SAMPLES-MERGED/PROFILE.db -C default

####################################################################################

INFO "anvi-model-protein-structures with DSSP"
anvi-model-protein-structures -c test-output/one_contig_five_genes.db \
                              --gene-caller-ids 2,4 \
                              --dump-dir test-output/STRUCTURES_full \
                              --output-db-path test-output/STRUCTURE.db

elif [ $1 = "continue"  ]
then
    cd sandbox
    
    #rm -rf test-output/STRUCTURES
    #rm -rf test-output/STRUCTURES_full

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
else
      echo "
        Unknown parameter $1 :/ Try 'new', or 'continue'.
        "
        exit -1
fi


INFO "anvi-gen-variability-profile"
anvi-gen-variability-profile -p test-output/SAMPLES-MERGED/PROFILE.db \
                             -c test-output/one_contig_five_genes.db \
                             -s test-output/STRUCTURE.db \
                             -C default \
                             -b bin1 \
                             --engine AA \
                             --quince-mode

#INFO "anvi-model-protein-structures"
#anvi-model-protein-structures -c test-output/one_contig_five_genes.db \
#                              --gene-caller-ids 2,3,4 \
#                              --output-db-path test-output/STRUCTURE.db \
#                              --skip-DSSP


INFO "Do you want the interactive interface? Run the following:"

echo "anvi-interactive -p `pwd`/test-output/SAMPLES-MERGED/PROFILE.db -c `pwd`/test-output/one_contig_five_genes.db"
echo
echo
