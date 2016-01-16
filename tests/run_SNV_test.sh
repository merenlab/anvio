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
    INFO "Creating the output directory ..."
    # change directory and clean the old mess if it exists
    cd sandbox
    rm -rf test-output
    mkdir test-output
    
    INFO "Anvo'o version ..."
    anvi-profile --version
    
    INFO "Generating short reads for SNV testing ..."
    anvi-script-gen-short-reads short_reads_for_SNV_testing.ini
    
    INFO "Generating a Bowtie2 ref ..."
    bowtie2-build single_contig.fa test-output/single_contig
    
    INFO "Mapping short reads to the ref ..."
    ../misc/bowtie_batch_single_fasta.sh test-output/short_reads.fa test-output/SNV-testing test-output/single_contig
    
    INFO "Generating anvi'o contigs database ..."
    anvi-gen-contigs-database -f single_contig.fa -o test-output/single_contig.db -L 60
    
    INFO "Profiling the BAM file ..."
    anvi-profile -i test-output/SNV-testing.bam -c test-output/single_contig.db -o test-output/SNV-testing -M 0 --cluster
    
    INFO "Calling the interactive interface ..."
    anvi-interactive -p test-output/SNV-testing/PROFILE.db -c test-output/single_contig.db
elif [ $1 = "continue"  ]
then
    cd sandbox

    if [ ! -f "test-output/single_contig.db"  ]
    then
      echo "
        You asked to continue with the previously generated files,
        but the contigs database is not there. Maybe you should start
        from scratch by re-running this script with the parameter 'new'.
        "
        exit -1
    else
        INFO "Attempting to continue with the previously generated files ..."

        rm -rf test-output/SNV-testing

        INFO "Re-profiling the BAM file ..."
        anvi-profile -i test-output/SNV-testing.bam -c test-output/single_contig.db -o test-output/SNV-testing -M 0 --cluster

        INFO "Calling the interactive interface ..."
        anvi-interactive -p test-output/SNV-testing/PROFILE.db -c test-output/single_contig.db

    fi
else
      echo "
        Unknown parameter $1 :/ Try 'new', or 'continue'.
        "
        exit -1
fi

