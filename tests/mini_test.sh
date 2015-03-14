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


# we first generate an empty annotation database using contigs.fa (keep in mind that 'contigs.fa'
# is the original file all samples were mapped to). here we use split size of 1000 (the default split
# size is much better for most projects. the small split size used here is simply for testing purposes)
INFO "Generating an EMPTY annotation database ..."
papi-gen-annotation-database -f contigs.fa -o test-output/ANNOTATION.db -L 1000

INFO "Populating the genes tables in the annotation database using 'myrast_gui' parser ..."
papi-populate-genes-table test-output/ANNOTATION.db -i myrast_gui/* -p myrast_gui

INFO "Populating the genes tables in the annotation database using 'myrast_cmdline_dont_use' parser ..."
papi-populate-genes-table test-output/ANNOTATION.db -i myrast_cmdline/svr_assign_to_dna_using_figfams.txt -p myrast_cmdline_dont_use

INFO "Populating the genes tables in the database using 'myrast_cmdline' parser ..."
papi-populate-genes-table test-output/ANNOTATION.db -p myrast_cmdline -i myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt

INFO "Exporting a standart matrix file from genes tables that were populated by 'myrast_cmdline' parser ..."
papi-export-genes-table test-output/ANNOTATION.db -o test-output/ANNOTATION_recovered.txt

INFO "Re-populating the genes tables in the annotation database using the recovered matrix file with 'default_matrix' parser ..."
papi-populate-genes-table test-output/ANNOTATION.db -p default_matrix -i test-output/ANNOTATION_recovered.txt

INFO "Populating search tables in the latest annotation database using default HMM profiles ..."
papi-populate-search-table test-output/ANNOTATION.db

INFO "Populating collections tables using mock clustering results for CONCOCT ..."
papi-populate-collections-table test-output/ANNOTATION.db --parser concoct -i concoct.txt

INFO "Annotation DB is ready; here are the tables in it:"
sqlite3 test-output/ANNOTATION.db '.tables'


# for each sample, run the profiling using the same split size used for the annotation database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 6M 7M 9M
do
    INFO "Profiling sample 204-$f ..."
    papi-profile -i test-output/204-$f.bam -o test-output/204-$f -a test-output/ANNOTATION.db
    echo
done


INFO "Merging profiles ..."
# merge samples
papi-merge test-output/204*/RUNINFO.cp -o test-output/204-MERGED

INFO "Generating coverages and sequences files for splits (for external binning) ..."
papi-export-splits-and-coverages test-output/ANNOTATION.db test-output/204-MERGED/PROFILE.db

INFO "Cluster contigs in the newly generated coverages file ..."
papi-matrix-to-newick test-output/204-MERGED/s204_MERGED-COVs.txt

INFO "Generating network descriptions for samples based on ORFs and functions ..."
# generate gene and function networks for the merge
papi-gen-network test-output/204-MERGED/RUNINFO.mcp test-output/ANNOTATION.db

INFO "Use papi-experimental-organization to generate another tree"
# this is meaningless here, but it is an example to show how one could generate new trees
papi-experimental-organization ../../PaPi/data/clusterconfigs/merged/tnf-cov -i test-output/204-MERGED -o test-output/204-MERGED/experimental-tree.txt

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
papi-interactive -r test-output/204-MERGED/RUNINFO.mcp -a test-output/ANNOTATION.db -A additional_metadata.txt -t test-output/204-MERGED/experimental-tree.txt
