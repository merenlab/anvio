#!/bin/bash
source 00.sh
set -e

INFO "Creating the output directory ..."
# change directory and clean the old mess if it exists
cd sandbox
rm -rf test-output
mkdir test-output


#
# RUN EVERYTHING
#


INFO "Anvo'o version ..."
anvi-profile --version

INFO "Initializing raw BAM files ..."
# init raw bam files.
for f in 6M 7M 9M
do
    anvi-init-bam 204_3contigs_"$f".bam -o test-output/204-$f
    echo
done


# we first generate an empty annotation database using contigs.fa (keep in mind that 'contigs.fa'
# is the original file all samples were mapped to). here we use split size of 1000 (the default split
# size is much better for most projects. the small split size used here is simply for testing purposes)
INFO "Generating an EMPTY annotation database ..."
anvi-gen-annotation-database -f contigs.fa -o test-output/ANNOTATION.db -L 1000

INFO "Populating the genes tables in the annotation database using 'myrast_gui' parser ..."
anvi-populate-genes-table test-output/ANNOTATION.db -i myrast_gui/* -p myrast_gui

INFO "Populating the genes tables in the annotation database using 'myrast_cmdline_dont_use' parser ..."
anvi-populate-genes-table test-output/ANNOTATION.db -i myrast_cmdline/svr_assign_to_dna_using_figfams.txt -p myrast_cmdline_dont_use

INFO "Populating the genes tables in the database using 'myrast_cmdline' parser ..."
anvi-populate-genes-table test-output/ANNOTATION.db -p myrast_cmdline -i myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt

INFO "Exporting a standart matrix file from genes tables that were populated by 'myrast_cmdline' parser ..."
anvi-export-genes-table test-output/ANNOTATION.db -o test-output/ANNOTATION_recovered.txt

INFO "Re-populating the genes tables in the annotation database using the recovered matrix file with 'default_matrix' parser ..."
anvi-populate-genes-table test-output/ANNOTATION.db -p default_matrix -i test-output/ANNOTATION_recovered.txt

INFO "Populating search tables in the latest annotation database using default HMM profiles ..."
anvi-populate-search-table test-output/ANNOTATION.db

INFO "Populating search tables in the latest annotation database using a mock HMM collection from an external directory ..."
anvi-populate-search-table test-output/ANNOTATION.db -H external_hmm_profile

INFO "Populating collections tables using mock clustering results for CONCOCT ..."
anvi-populate-collections-table test-output/ANNOTATION.db --parser concoct -i concoct.txt

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

INFO "Generating coverages and sequences files for splits (for external binning) ..."
anvi-export-splits-and-coverages -a test-output/ANNOTATION.db -p test-output/204-MERGED/PROFILE.db

INFO "Cluster contigs in the newly generated coverages file ..."
anvi-matrix-to-newick test-output/204-MERGED/s204_MERGED-COVs.txt

INFO "Generating network descriptions for samples based on ORFs and functions ..."
# generate gene and function networks for the merge
anvi-gen-network test-output/204-MERGED/RUNINFO.mcp test-output/ANNOTATION.db

INFO "Use anvi-experimental-organization to generate another tree"
# this is meaningless here, but it is an example to show how one could generate new trees
anvi-experimental-organization ../../anvio/data/clusterconfigs/merged/tnf-cov -i test-output/204-MERGED -o test-output/204-MERGED/experimental-tree.txt -a test-output/ANNOTATION.db

INFO "Importing collections from external files into the profile database"
anvi-import-collection example_external_collections/adhoc_collections.txt -S 'C_IMPORTED' -p test-output/204-MERGED/PROFILE.db --colors example_external_collections/adhoc_colors.txt

INFO "Use CONCOCT to cluster splits in the merged profile and export as a text file..."
anvi-cluster-with-concoct -p test-output/204-MERGED/PROFILE.db -a test-output/ANNOTATION.db -o test-output/anvio_concoct_clusters.txt --source-identifier 'cmdline_concoct'

INFO "Summarizing CONCOCT results ..."
anvi-summarize -p test-output/204-MERGED/PROFILE.db -a test-output/ANNOTATION.db -o test-output/204-MERGED-SUMMARY -c 'cmdline_concoct'

INFO "Generate a variabilty profile for Bin_1 using a collection id"
anvi-gen-variability-profile -a test-output/ANNOTATION.db -p test-output/204-MERGED/PROFILE.db -c cmdline_concoct -b Bin_1 -o test-output/variability_Bin_1.txt

INFO "Get sequences for HMM hits for a bin in a collection ..."
anvi-get-sequences-for-hmm-hits -p test-output/204-MERGED/PROFILE.db -a test-output/ANNOTATION.db -c CONCOCT -b Bin_1 -o test-output/hmm_hits_sequences_in_Bin_1.txt

INFO "Generate a variabilty profile for Bin_1 using split ids stored in a file (after summary)"
anvi-gen-variability-profile -a test-output/ANNOTATION.db -p test-output/204-MERGED/PROFILE.db -s test-output/204-MERGED-SUMMARY/bin_by_bin/Bin_1/Bin_1-original_split_names.txt -o test-output/variability_Bin_1_ALT.txt

INFO "Generate a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o test-output/SAMPLES.db

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p test-output/204-MERGED/PROFILE.db \
                 -a test-output/ANNOTATION.db \
                 -s test-output/SAMPLES.db \
                 -A additional_view_data.txt \
                 -t test-output/204-MERGED/experimental-tree.txt \
                 -V additional_view.txt \
                 --split-hmm-layers

INFO "Firing up the interactive interface to refine a bin ..."
anvi-refine -p test-output/204-MERGED/PROFILE.db -a test-output/ANNOTATION.db -s test-output/SAMPLES.db -c CONCOCT -b Bin_1
