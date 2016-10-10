#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Initializing raw BAM files"
# init raw bam files.
for f in 01 02 03
do
    anvi-init-bam $files/SAMPLE-RAW-$f.bam --output-file $output_dir/SAMPLE-$f.bam
    echo
done

INFO "Generating an EMPTY contigs database"
anvi-gen-contigs-database -f $files/contigs.fa -o $output_dir/CONTIGS.db -L 1000

INFO "Populating taxonomy for splits table in the database using 'centrifuge' parser"
anvi-import-taxonomy -c $output_dir/CONTIGS.db -p centrifuge -i $files/example_files_for_centrifuge_taxonomy/*

INFO "Populating search tables in the latest contigs database using default HMM profiles"
anvi-run-hmms -c $output_dir/CONTIGS.db --num-threads 2

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
    anvi-profile -i $output_dir/SAMPLE-$f.bam -o $output_dir/SAMPLE-$f -c $output_dir/CONTIGS.db
    echo
done


INFO "Merging profiles"
# merge samples
anvi-merge $output_dir/SAMPLE-*/RUNINFO.cp -o $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db

INFO "Generating a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D $files/samples-information.txt -R $files/samples-order.txt -o $output_dir/SAMPLES.db

INFO "Importing a state file into the merged profile"
anvi-import-state -p $output_dir/SAMPLES-MERGED/PROFILE.db --state $files/example_state.json --name default

INFO "Firing up the interactive interface"
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -s $output_dir/SAMPLES.db \
                 --split-hmm-layers

INFO "Summarizing CONCOCT results"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -o $output_dir/SAMPLES-MERGED-SUMMARY -C 'CONCOCT'
