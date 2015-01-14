#!/bin/bash

C() {
    echo -e "\033[0;30m\033[46m$1\033[0m"
}

INFO() { 
    echo
    C "###############################################################"
    C "#"
    C "# $1"
    C "#"
    C "###############################################################"
    echo
} 

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd mini_test
rm -rf test-output
mkdir test-output


INFO "Initializing raw BAM files ..."
# init raw bam files.
for f in 6M 7M 9M
do
    papi-init-bam 204_3contigs_"$f".bam -o test-output/204-$f
    echo
done


INFO "Generating the annotation database ..."
# generate an annotation db using files obtained from myrast_gui using contigs.fa (contigs.fa
# is the original file all samples were mapped to) using split size 1000 (the default split
# size is better for most projects, small split size here is for testing purposes) (following
# two lines are generating the ANNOTATION database using gui and cmdline outputs, obviously
# the second one overwrites the result of the first one. they are both here for testing
# purposes, but only the result of the second command is used for later steps)
papi-gen-annotation --contigs contigs.fa -p myrast_gui myrast_gui/* -o test-output/ -L 1000 
papi-gen-annotation --contigs contigs.fa -p myrast_cmdline_dont_use myrast_cmdline/svr_assign_to_dna_using_figfams.txt -o test-output/ -L 1000 
papi-gen-annotation --contigs contigs.fa -p myrast_cmdline myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt -o test-output/ -L 1000 


INFO "Profiling samples ..."
# for each sample, run the profiling using the same split size used for the annotation database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 6M 7M 9M
do
    papi-profile -i test-output/204-$f.bam -o test-output/204-$f -a test-output/ANNOTATION.db -L 1000
    echo
done


INFO "Merging profiles ..."
# merge samples
papi-merge-multiple-runs test-output/204*/RUNINFO.cp -o test-output/204-MERGED


INFO "Generating network descriptions for samples based on ORFs and functions ..."
# generate gene and function networks for the merge
papi-gen-network test-output/204-MERGED/RUNINFO.mcp test-output/ANNOTATION.db


INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
papi-interactive-binning -r test-output/204-MERGED/RUNINFO.mcp -a test-output/ANNOTATION.db
