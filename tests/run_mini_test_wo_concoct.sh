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
    anvi-init-bam 204_3contigs_"$f".bam --output-file-prefix test-output/204-$f
    echo
done


INFO "Generating an EMPTY contigs database ..."
anvi-gen-contigs-database -f contigs.fa -o test-output/CONTIGS.db -L 1000

INFO "Populating taxonomy for splits table in the database using 'myrast_cmdline' parser ..."
anvi-import-taxonomy-from-gene-annotations -c test-output/CONTIGS.db -p myrast_cmdline -i myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt

INFO "Populating search tables in the latest contigs database using default HMM profiles ..."
anvi-populate-search-table -c test-output/CONTIGS.db

INFO "Importing gene function calls using 'interproscan' parser ..."
anvi-import-functions -c test-output/CONTIGS.db -i example_interpro_output.tsv -p interproscan

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory ..."
anvi-populate-search-table -c test-output/CONTIGS.db -H external_hmm_profile

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 test-output/CONTIGS.db '.tables'

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 6M 7M 9M
do
    INFO "Profiling sample 204-$f ..."
    anvi-profile -i test-output/204-$f.bam -o test-output/204-$f -c test-output/CONTIGS.db
    echo
done


INFO "Merging profiles ..."
# merge samples
anvi-merge test-output/204*/RUNINFO.cp -o test-output/204-MERGED -c test-output/CONTIGS.db --skip-concoct-binning

INFO "Importing external binning results for splits into the profile database as 'SPLITS_IMPORTED'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -p test-output/204-MERGED/PROFILE.db \
                       -c test-output/CONTIGS.db \
                       --source-identifier 'SPLITS_IMPORTED' \
                       --colors example_files_for_external_binning_results/example_colors_file.txt

INFO "Summarizing 'SPLITS_IMPORTED' results ..."
anvi-summarize -p test-output/204-MERGED/PROFILE.db -c test-output/CONTIGS.db -o test-output/204-MERGED-SUMMARY -C 'SPLITS_IMPORTED'

INFO "Generating a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o test-output/SAMPLES.db

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p test-output/204-MERGED/PROFILE.db \
                 -c test-output/CONTIGS.db \
                 -s test-output/SAMPLES.db \
                 --split-hmm-layers
