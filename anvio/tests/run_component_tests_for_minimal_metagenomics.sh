#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Initializing raw BAM files"
# init raw bam files.
for f in 01 02 03
do
    anvi-init-bam $files/SAMPLE-$f-RAW.bam --output-file $output_dir/SAMPLE-$f.bam
    echo
done

INFO "Generating an EMPTY contigs database"
anvi-gen-contigs-database -f $files/contigs.fa -o $output_dir/CONTIGS.db -L 1000 --project-name "Contigs DB for anvi'o mini self-test"

INFO "Populating search tables in the latest contigs database using default HMM profiles"
anvi-run-hmms -c $output_dir/CONTIGS.db $thread_controller

INFO "Importing gene function calls using 'interproscan' parser"
anvi-import-functions -c $output_dir/CONTIGS.db -i $files/example_interpro_output.tsv -p interproscan

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory"
anvi-run-hmms -c $output_dir/CONTIGS.db -H $files/external_hmm_profile

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 $output_dir/CONTIGS.db '.tables'

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner $output_dir directory for each sample.
for f in 01 02 03
do
    INFO "Profiling sample SAMPLE-$f"
    anvi-profile -i $output_dir/SAMPLE-$f.bam \
                 -o $output_dir/SAMPLE-$f \
                 -c $output_dir/CONTIGS.db \
                 --cluster \
                 --profile-SCVs \
                 --display-db-calls

    INFO "Importing short-read-level taxonomy for SAMPLE-$f"
    anvi-import-taxonomy-for-layers -p $output_dir/SAMPLE-$f/PROFILE.db \
                                    -i $files/example_files_for_kraken_hll_taxonomy/SAMPLE-$f.mpa \
                                    --parser krakenuniq
    echo
done

INFO "Merging profiles"
# merge samples
anvi-merge $output_dir/SAMPLE-*/PROFILE.db -o $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db --description $files/example_description.md

INFO "Import layer additional data from file"
anvi-import-misc-data $files/samples-information.txt \
                      -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table layers

INFO "Import layer orders from file"
anvi-import-misc-data $files/samples-order.txt \
                      -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table layer_orders

INFO "Importing a state file into the merged profile"
anvi-import-state -p $output_dir/SAMPLES-MERGED/PROFILE.db --state $files/example_state.json --name default

INFO "Importing a collection file into the merged profile"
anvi-import-collection -c $output_dir/CONTIGS.db \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -C CONCOCT \
                       $files/concoct_mini_test.txt

INFO "Listing collections available"
anvi-show-collections-and-bins -p $output_dir/SAMPLES-MERGED/PROFILE.db

INFO "Firing up the interactive interface"
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 $dry_run_controller

INFO "Summarizing CONCOCT results"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -o $output_dir/SAMPLES-MERGED-SUMMARY -C 'CONCOCT' --init-gene-coverages
