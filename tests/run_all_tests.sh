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
    anvi-init-bam 204_3contigs_"$f".bam -O test-output/204-$f
    echo
done


# we first generate an empty contigs database using contigs.fa (keep in mind that 'contigs.fa'
# is the original file all samples were mapped to). here we use split size of 1000 (the default split
# size is much better for most projects. the small split size used here is simply for testing purposes)
INFO "Generating an EMPTY contigs database ..."
anvi-gen-contigs-database -f contigs.fa -o test-output/CONTIGS.db -L 1000

INFO "Populating taxonomy for splits table in the database using 'myrast_gui' parser ..."
anvi-import-taxonomy-from-gene-annotations -c test-output/CONTIGS.db -i myrast_gui/* -p myrast_gui

INFO "Re-populating taxonomy for splits table in the database using 'myrast_cmdline_dont_use' parser ..."
anvi-import-taxonomy-from-gene-annotations -c test-output/CONTIGS.db -i myrast_cmdline/svr_assign_to_dna_using_figfams.txt -p myrast_cmdline_dont_use

INFO "Re-populating taxonomy for splits table in the database using 'myrast_cmdline' parser ..."
anvi-import-taxonomy-from-gene-annotations -c test-output/CONTIGS.db -p myrast_cmdline -i myrast_cmdline/svr_call_pegs.txt myrast_cmdline/svr_assign_using_figfams.txt

INFO "Re-populating taxonomy for splits table in the database using the recovered matrix file with 'default_matrix' parser ..."
anvi-import-taxonomy-from-gene-annotations -c test-output/CONTIGS.db -i gene_calls_sample_matrix.txt

INFO "Populating HMM hits tables in the latest contigs database using default HMM profiles ..."
anvi-populate-search-table -c test-output/CONTIGS.db --num-threads 2

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory ..."
anvi-populate-search-table -c test-output/CONTIGS.db -H external_hmm_profile

INFO "Importing gene function calls using 'interproscan' parser ..."
anvi-import-functional-annotation-of-genes -c test-output/CONTIGS.db -i example_interpro_output.tsv -p interproscan

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
anvi-merge test-output/204*/RUNINFO.cp -o test-output/204-MERGED -c test-output/CONTIGS.db

INFO "Generating coverages and sequences files for splits (for external binning) ..."
anvi-export-splits-and-coverages -c test-output/CONTIGS.db -p test-output/204-MERGED/PROFILE.db

INFO "Cluster contigs in the newly generated coverages file ..."
anvi-matrix-to-newick test-output/204-MERGED/s204_MERGED-COVs.txt

INFO "Generating network descriptions for samples based on ORFs and functions ..."
# generate gene and function networks for the merge
anvi-gen-network -r test-output/204-MERGED/RUNINFO.mcp -c test-output/CONTIGS.db

INFO "Use anvi-experimental-organization to generate a tree from a new configuration to store it in a file (not in the database)"
anvi-experimental-organization example_clustering_configuration.ini -i test-output/204-MERGED -c test-output/CONTIGS.db -o test-output/204-MERGED/EXP-ORG-FILE.txt --skip-store-in-db

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, and add the resulting tree into the database as 'EXP-ORG-DB'"
anvi-experimental-organization example_clustering_configuration.ini -i test-output/204-MERGED -c test-output/CONTIGS.db -p test-output/204-MERGED/PROFILE.db --name EXP-ORG-DB

INFO "Importing external binning results for splits into the contigs database as 'SPLITS_IMPORTED_INTO_CONTIGS_DB'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -c test-output/CONTIGS.db \
                       --collection-name 'SPLITS_IMPORTED_INTO_CONTIGS_DB' \
                       --bins-info example_files_for_external_binning_results/example_bins_info_file.txt

INFO "Importing external binning results for splits into the profile database as 'SPLITS_IMPORTED'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -p test-output/204-MERGED/PROFILE.db \
                       -c test-output/CONTIGS.db \
                       --collection-name 'SPLITS_IMPORTED' \
                       --bins-info example_files_for_external_binning_results/example_bins_info_file.txt

INFO "Importing external binning results for splits into the profile database as 'CONTIGS_IMPORTED'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_contigs.txt \
                       -c test-output/CONTIGS.db \
                       -p test-output/204-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_IMPORTED' \
                       --bins-info example_files_for_external_binning_results/example_bins_info_file.txt \
                       --contigs-mode

INFO "Exporting the 'CONTIGS_IMPORTED' collection that was just imported ..."
anvi-export-collection -p test-output/204-MERGED/PROFILE.db -C CONTIGS_IMPORTED --output-file-prefix test-output/exported-collection

INFO "Re-importing a collection from files just exported for CONTIGS_IMPORTED collection (just to confuse you, and to see if we can import stuff we export) ..."
anvi-import-collection test-output/exported-collection.txt \
                       -c test-output/CONTIGS.db \
                       -p test-output/204-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_RE_IMPORTED' \
                       --bins-info test-output/exported-collection-info.txt

INFO "Use CONCOCT to cluster splits in the merged profile and export as a text file..."
anvi-cluster-with-concoct -p test-output/204-MERGED/PROFILE.db -c test-output/CONTIGS.db -o test-output/anvio_concoct_clusters.txt --collection-name 'cmdline_concoct'

INFO "Recover short reads for Bin_2 in CONCOCT collection and store them in a FASTA file ..."
anvi-get-short-reads-from-bam -p test-output/204-MERGED/PROFILE.db -c test-output/CONTIGS.db -C CONCOCT -b Bin_2 -o test-output/short_reads_for_Bin_2.fasta test-output/*bam

INFO "Summarizing CONCOCT results ..."
anvi-summarize -p test-output/204-MERGED/PROFILE.db -c test-output/CONTIGS.db -o test-output/204-MERGED-SUMMARY -C 'cmdline_concoct'

INFO "Generate a variabilty profile for Bin_1 using a collection id"
anvi-gen-variability-profile -c test-output/CONTIGS.db -p test-output/204-MERGED/PROFILE.db -C cmdline_concoct -b Bin_1 -o test-output/variability_Bin_1.txt --quince-mode

INFO "Generate a variabilty profile for Bin_1 using split ids stored in a file (after summary)"
anvi-gen-variability-profile -c test-output/CONTIGS.db \
                             -p test-output/204-MERGED/PROFILE.db \
                             --splits-of-interest test-output/204-MERGED-SUMMARY/bin_by_bin/Bin_1/Bin_1-original_split_names.txt \
                             -o test-output/variability_Bin_1_ALT.txt

INFO "Generating amino acid frequencies for gene caller id 3 in 204-6M.bam ..."
anvi-get-aa-frequencies -i test-output/204-6M.bam -c test-output/CONTIGS.db --gene-caller-id 3 -o test-output/AA_frequencies_for_gene_caller_id_3.txt

INFO "Getting back the sequence for gene call 3 ..."
anvi-get-sequences-for-gene-calls -c test-output/CONTIGS.db --gene-caller-ids 3 -o test-output/Sequence_for_gene_caller_id_3.fa

INFO "Get sequences for HMM hits for a bin in a collection ..."
anvi-get-sequences-for-hmm-hits -p test-output/204-MERGED/PROFILE.db -c test-output/CONTIGS.db -C CONCOCT -b Bin_1 -o test-output/hmm_hits_sequences_in_Bin_1.txt

INFO "Generate a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D samples-information.txt -R samples-order.txt -o test-output/SAMPLES.db

INFO "Get linkmers from all BAM files for some distant positions"
anvi-report-linkmers --contigs-and-positions distant_positions_for_linkmers.txt -i test-output/*.bam -o test-output/distant_linkmers.txt

INFO "Get linkmers from all BAM files for some adjacent positions connected with --only-complete-links"
anvi-report-linkmers --contigs-and-positions adjacent_positions_for_linkmers.txt -i test-output/*.bam -o test-output/adjacent_linkmers.txt --only-complete-links

INFO "Oligotype linkmers report generated for adjacent nucleotide positions ..."
anvi-oligotype-linkmers -i test-output/adjacent_linkmers.txt -o test-output/

INFO "Search for functions to get split names with matching genes ..."
anvi-search-functions-in-splits -c test-output/CONTIGS.db --search transporter,kinase -o test-output/transporter-hits.txt --verbose

INFO "Get all short reads that map to the gene ID 38 (which is a Zinc transpoprter)"
anvi-get-short-reads-mapping-to-a-gene -c test-output/CONTIGS.db --gene-caller-id 38 --leeway 100 -i test-output/*bam -o test-output/reads-mapping-to-gene-id-38.fa

INFO "Get AA counts for the entire contigs database ..."
anvi-get-aa-counts -c test-output/CONTIGS.db -o test-output/aa_counts_for_contigs_db.txt
column -t test-output/aa_counts_for_contigs_db.txt

INFO "Get AA counts for bins in collection CONCOCT stored in the merged profile ..."
anvi-get-aa-counts -c test-output/CONTIGS.db -p test-output/204-MERGED/PROFILE.db -C CONCOCT -o test-output/aa_counts_for_bins_in_collection_CONCOCT.txt
column -t test-output/aa_counts_for_bins_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile ..."
anvi-get-aa-counts -c test-output/CONTIGS.db -p test-output/204-MERGED/PROFILE.db -C CONCOCT -o test-output/aa_counts_for_bin_3_in_collection_CONCOCT.txt -B sample_CONCOCT_bin_id.txt
column -t test-output/aa_counts_for_bin_3_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile ..."
anvi-get-aa-counts -c test-output/CONTIGS.db --contigs-of-interest sample_contig_ids.txt -o test-output/aa_counts_for_two_contigs.txt
column -t test-output/aa_counts_for_two_contigs.txt

INFO "Get AA counts for five genes ..."
anvi-get-aa-counts -c test-output/CONTIGS.db --gene-caller-ids sample_gene_call_ids.txt -o test-output/aa_counts_for_five_genes.txt
column -t test-output/aa_counts_for_five_genes.txt

INFO "Firing up the interactive interface ..."
# fire up the browser to show how does the merged samples look like.
anvi-interactive -p test-output/204-MERGED/PROFILE.db \
                 -c test-output/CONTIGS.db \
                 -s test-output/SAMPLES.db \
                 -A additional_view_data.txt \
                 -t test-output/204-MERGED/EXP-ORG-FILE.txt \
                 -V additional_view.txt \
                 --split-hmm-layers

INFO "Firing up the interactive interface to refine a bin ..."
anvi-refine -p test-output/204-MERGED/PROFILE.db -c test-output/CONTIGS.db -s test-output/SAMPLES.db -C CONCOCT -b Bin_1
