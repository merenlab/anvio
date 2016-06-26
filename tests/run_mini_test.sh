#!/bin/bash
source 00.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd sandbox
rm -rf test-output
mkdir test-output

INFO "Anvo'o version ..."
anvi-profile --version

INFO "Initializing raw BAM files ..."
# init raw bam files.
for f in 01 02 03
do
    anvi-init-bam SAMPLE-RAW-$f.bam --output-file test-output/SAMPLE-$f.bam
    echo
done


INFO "Generating an EMPTY contigs database ..."
anvi-gen-contigs-database -f contigs.fa -o test-output/CONTIGS.db -L 1000

INFO "Populating taxonomy for splits table in the database using 'centrifuge' parser ..."
anvi-import-taxonomy -c test-output/CONTIGS.db -p centrifuge -i example_files_for_centrifuge_taxonomy/*

INFO "Populating search tables in the latest contigs database using default HMM profiles ..."
anvi-run-hmms -c test-output/CONTIGS.db --num-threads 2

INFO "Importing gene function calls using 'interproscan' parser ..."
anvi-import-functions -c test-output/CONTIGS.db -i example_interpro_output.tsv -p interproscan

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory ..."
anvi-run-hmms -c test-output/CONTIGS.db -H external_hmm_profile

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 test-output/CONTIGS.db '.tables'

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 01 02 03
do
    INFO "Profiling sample SAMPLE-$f ..."
    anvi-profile -i test-output/SAMPLE-$f.bam -o test-output/SAMPLE-$f -c test-output/CONTIGS.db
    echo
done


INFO "Merging profiles ..."
# merge samples
anvi-merge test-output/SAMPLE-*/RUNINFO.cp -o test-output/SAMPLES-MERGED -c test-output/CONTIGS.db

INFO "Generating a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o test-output/SAMPLES.db

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p test-output/SAMPLES-MERGED/PROFILE.db \
                 -c test-output/CONTIGS.db \
                 -s test-output/SAMPLES.db \
                 --split-hmm-layers

INFO "Summarizing CONCOCT results ..."
anvi-summarize -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -o test-output/SAMPLES-MERGED-SUMMARY -C 'CONCOCT'


