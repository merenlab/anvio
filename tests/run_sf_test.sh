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
for f in SF02 SF03 SF15
do
    anvi-init-bam BAMs_SF/"$f".bam --output-file-prefix test-output/$f
    echo
done


INFO "Generating an EMPTY contigs database ..."
anvi-gen-contigs-database -f BAMs_SF/contigs.fa -o test-output/CONTIGS.db -L 100

INFO "Populating search tables in the latest contigs database using default HMM profiles ..."
anvi-run-hmms -c test-output/CONTIGS.db

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 test-output/CONTIGS.db '.tables'

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in SF02 SF03 SF15
do
    INFO "Profiling sample 204-$f ..."
    anvi-profile -i test-output/$f.bam -o test-output/$f -c test-output/CONTIGS.db -M 0
    echo
done


INFO "Merging profiles ..."
# merge samples
anvi-merge test-output/*/RUNINFO.cp -o test-output/SF-MERGED -c test-output/CONTIGS.db --skip-concoct

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p test-output/SF-MERGED/PROFILE.db \
                 -c test-output/CONTIGS.db
