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
for f in 6M 7M 9M
do
    anvi-init-bam 204_3contigs_"$f".bam -o test-output/204-$f
    echo
done


INFO "Generating an EMPTY annotation database ..."
anvi-gen-annotation-database -f contigs.fa -o test-output/ANNOTATION.db -L 1000

INFO "Populating the genes tables in the database using 'myrast_cmdline' parser ..."
anvi-populate-genes-table test-output/ANNOTATION.db -p myrast_cmdline -i myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt

INFO "Populating search tables in the latest annotation database using default HMM profiles ..."
anvi-populate-search-table test-output/ANNOTATION.db

INFO "Populating search tables in the latest annotation database using a mock HMM collection from an external directory ..."
anvi-populate-search-table test-output/ANNOTATION.db -H external_hmm_profile

INFO "Annotation DB is ready; here are the tables in it:"
sqlite3 test-output/ANNOTATION.db '.tables'

# for each sample, run the profiling using the same split size used for the annotation database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 6M 7M 9M
do
    INFO "Profiling sample 204-$f ..."
    anvi-profile -i test-output/204-$f.bam -o test-output/204-$f -a test-output/ANNOTATION.db
    echo
done


INFO "Merging profiles ..."
# merge samples
anvi-merge test-output/204*/RUNINFO.cp -o test-output/204-MERGED -a test-output/ANNOTATION.db

INFO "Summarizing CONCOCT results ..."
anvi-summarize -p test-output/204-MERGED/PROFILE.db -a test-output/ANNOTATION.db -o test-output/204-MERGED-SUMMARY -c 'CONCOCT'

INFO "Generating a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o test-output/SAMPLES.db

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p test-output/204-MERGED/PROFILE.db \
                 -a test-output/ANNOTATION.db \
                 -s test-output/SAMPLES.db \
                 --split-hmm-layers
