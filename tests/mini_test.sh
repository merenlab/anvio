#!/bin/bash
source 00.sh
set -e

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


# generate an annotation db using files obtained from myrast_gui using contigs.fa (contigs.fa
# is the original file all samples were mapped to) using split size 1000 (the default split
# size is better for most projects, small split size here is for testing purposes) (following
# two lines are generating the ANNOTATION database using gui and cmdline outputs, obviously
# the second one overwrites the result of the first one. they are both here for testing
# purposes, but only the result of the second command is used for later steps)


INFO "Generating an annotation database using 'myrast_gui' parser ..."
papi-gen-annotation myrast_gui/* -p myrast_gui --contigs contigs.fa -o test-output/ANNOTATION-myrast_gui.db -L 1000 --skip-search-tables 

INFO "Generating an annotation database using 'myrast_cmdline_dont_use' parser ..."
papi-gen-annotation myrast_cmdline/svr_assign_to_dna_using_figfams.txt -p myrast_cmdline_dont_use --contigs contigs.fa -o test-output/ANNOTATION-myrast_cmdline_dont_use.db -L 1000 --skip-search-tables

INFO "Generating an annotation database using 'myrast_cmdline' parser ..."
papi-gen-annotation myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt -p myrast_cmdline --contigs contigs.fa -o test-output/ANNOTATION-myrast_cmdline.db -L 1000 --skip-search-tables

INFO "Recovering a standart matrix file from the annotation database generated using 'myrast_cmdline' parser ..." --skip-search-tables
papi-export-annotation-table test-output/ANNOTATION-myrast_cmdline.db -o test-output/ANNOTATION_recovered.txt

INFO "Re-generating an annotation database using the recovered matrix file with 'default_matrix' parser ..."
papi-gen-annotation test-output/ANNOTATION_recovered.txt -p default_matrix --contigs contigs.fa -o test-output/ANNOTATION.db -L 1000 --skip-search-tables


# in the lines below we used '--skip-search-tables' flag, and PaPi skipped generating HMM scan
# results in those databases. Now we will do that step alone specifically:
INFO "Populating search tables in the latest annotation database using default HMM profiles ..."
papi-populate-search-table contigs.fa test-output/ANNOTATION.db

INFO "Populating collections tables using mock clustering results for CONCOCT ..."
papi-populate-collections contigs.fa test-output/ANNOTATION.db --parser concoct -i concoct.txt


INFO "Annotation DB is ready; here are the tables in it:"
sqlite3 test-output/ANNOTATION.db '.tables'


# for each sample, run the profiling using the same split size used for the annotation database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 6M 7M 9M
do
    INFO "Profiling sample 204-$f ..."
    papi-profile -i test-output/204-$f.bam -o test-output/204-$f -a test-output/ANNOTATION.db -L 1000
    echo
done


INFO "Merging profiles ..."
# merge samples
papi-merge test-output/204*/RUNINFO.cp -o test-output/204-MERGED


INFO "Generating network descriptions for samples based on ORFs and functions ..."
# generate gene and function networks for the merge
papi-gen-network test-output/204-MERGED/RUNINFO.mcp test-output/ANNOTATION.db


INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
papi-interactive -r test-output/204-MERGED/RUNINFO.mcp -a test-output/ANNOTATION.db -A additional_metadata.txt
