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
for f in 01 02 03
do
    anvi-init-bam SAMPLE-RAW-$f.bam --output-file test-output/SAMPLE-$f.bam
    echo
done

INFO "Reformat the contigs FASTA"
anvi-script-reformat-fasta contigs.fa -o test-output/contigs.fa -l 0 --simplify-names --prefix test_prefix --report test-output/contigs-reformat-report.txt
echo
column -t test-output/contigs-reformat-report.txt 

# we first generate an empty contigs database using contigs.fa (keep in mind that 'contigs.fa'
# is the original file all samples were mapped to). here we use split size of 1000 (the default split
# size is much better for most projects. the small split size used here is simply for testing purposes)
INFO "Generating an EMPTY contigs database ..."
anvi-gen-contigs-database -f contigs.fa -o test-output/CONTIGS.db -L 1000

INFO "Populating taxonomy for splits table in the database using 'centrifuge' parser ..."
anvi-import-taxonomy -c test-output/CONTIGS.db -p centrifuge -i example_files_for_centrifuge_taxonomy/*

INFO "Populating HMM hits tables in the latest contigs database using default HMM profiles ..."
anvi-run-hmms -c test-output/CONTIGS.db --num-threads 2

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory ..."
anvi-run-hmms -c test-output/CONTIGS.db -H external_hmm_profile

INFO "Importing gene function calls using 'interproscan' parser ..."
anvi-import-functions -c test-output/CONTIGS.db -i example_interpro_output.tsv -p interproscan

INFO "Importing gene function calls INCREMENTALLY using a TAB-delimited default input matrix ..."
anvi-import-functions -c test-output/CONTIGS.db -i example_gene_functions_input_matrix.txt

INFO "REPLACING gene function calls using a TAB-delimited default input matrix ..."
anvi-import-functions -c test-output/CONTIGS.db -i example_gene_functions_input_matrix.txt --drop-previous-annotations

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 test-output/CONTIGS.db '.tables'

INFO "Generating a 'blank profile' with the newly generated contigs database ..."
anvi-profile -c test-output/CONTIGS.db -o test-output/BLANK-PROFILE -S BLANK --blank-profile

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner test-output directory for each sample.
for f in 01 02 03
do
    INFO "Profiling sample SAMPLE-$f ..."
    anvi-profile -i test-output/SAMPLE-$f.bam -o test-output/SAMPLE-$f -c test-output/CONTIGS.db --profile-AA-frequencies
    echo
done


INFO "Merging profiles ..."
# merge samples
anvi-merge test-output/SAMPLE-*/RUNINFO.cp -o test-output/SAMPLES-MERGED -c test-output/CONTIGS.db

INFO "Add a new variable into the RUNINFO file of the merged profile ..."
anvi-script-update-runinfo-variable test-output/SAMPLES-MERGED/RUNINFO.mcp --variable TEST-VARIABLE --set-bool true

INFO "Generating coverages and sequences files for splits (for external binning) ..."
anvi-export-splits-and-coverages -c test-output/CONTIGS.db -p test-output/SAMPLES-MERGED/PROFILE.db

INFO "Cluster contigs in the newly generated coverages file ..."
anvi-matrix-to-newick test-output/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt

INFO "Generating network descriptions for samples based on ORFs and functions ..."
# generate gene and function networks for the merge
anvi-gen-network -r test-output/SAMPLES-MERGED/RUNINFO.mcp -c test-output/CONTIGS.db

INFO "Use anvi-experimental-organization to generate a tree from a new configuration to store it in a file (not in the database)"
anvi-experimental-organization example_clustering_configuration.ini -i test-output/SAMPLES-MERGED -c test-output/CONTIGS.db -o test-output/SAMPLES-MERGED/EXP-ORG-FILE.txt --skip-store-in-db

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, and add the resulting tree into the database as 'EXP-ORG-DB'"
anvi-experimental-organization example_clustering_configuration.ini -i test-output/SAMPLES-MERGED -c test-output/CONTIGS.db -p test-output/SAMPLES-MERGED/PROFILE.db --name EXP-ORG-DB

INFO "Importing external binning results for splits into the contigs database as 'SPLITS_IMPORTED_INTO_CONTIGS_DB'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -c test-output/CONTIGS.db \
                       --collection-name 'SPLITS_IMPORTED_INTO_CONTIGS_DB' \
                       --bins-info example_files_for_external_binning_results/example_bins_info_file.txt

INFO "Importing external binning results for splits into the profile database as 'SPLITS_IMPORTED'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -p test-output/SAMPLES-MERGED/PROFILE.db \
                       -c test-output/CONTIGS.db \
                       --collection-name 'SPLITS_IMPORTED' \
                       --bins-info example_files_for_external_binning_results/example_bins_info_file.txt

INFO "Importing external binning results for splits into the profile database as 'CONTIGS_IMPORTED'"
anvi-import-collection example_files_for_external_binning_results/external_binning_of_contigs.txt \
                       -c test-output/CONTIGS.db \
                       -p test-output/SAMPLES-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_IMPORTED' \
                       --bins-info example_files_for_external_binning_results/example_bins_info_file.txt \
                       --contigs-mode

INFO "Exporting the 'CONTIGS_IMPORTED' collection that was just imported ..."
anvi-export-collection -p test-output/SAMPLES-MERGED/PROFILE.db -C CONTIGS_IMPORTED --output-file-prefix test-output/exported-collection

INFO "Re-importing a collection from files just exported for CONTIGS_IMPORTED collection (just to confuse you, and to see if we can import stuff we export) ..."
anvi-import-collection test-output/exported-collection.txt \
                       -c test-output/CONTIGS.db \
                       -p test-output/SAMPLES-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_RE_IMPORTED' \
                       --bins-info test-output/exported-collection-info.txt

INFO "Use CONCOCT to cluster splits in the merged profile and export as a text file..."
anvi-cluster-with-concoct -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -o test-output/anvio_concoct_clusters.txt --collection-name 'cmdline_concoct'

INFO "Recover short reads for Bin_2 in CONCOCT collection and store them in a FASTA file ..."
anvi-get-short-reads-from-bam -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -C CONCOCT -b Bin_2 -o test-output/short_reads_for_Bin_2.fasta test-output/*bam

INFO "Rename bins in collection 'cmdline_concoct'"
anvi-rename-bins -c test-output/CONTIGS.db -p test-output/SAMPLES-MERGED/PROFILE.db --prefix 'PSAMPLES' -C cmdline_concoct --report-file test-output/renaming-report.txt

echo
column -t test-output/renaming-report.txt
echo

INFO "Summarizing CONCOCT results ..."
anvi-summarize -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -o test-output/SAMPLES-MERGED-SUMMARY -C 'cmdline_concoct'

INFO "Generate a variabilty profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c test-output/CONTIGS.db -p test-output/SAMPLES-MERGED/PROFILE.db -C cmdline_concoct -b PSAMPLES_Bin_00001 -o test-output/variability_PSAMPLES_Bin_00001.txt --quince-mode

INFO "Generate a variabilty profile for PSAMPLES_Bin_00001 using split ids and gene ids of interest (after summary)"
anvi-gen-variability-profile -c test-output/CONTIGS.db \
                             -p test-output/SAMPLES-MERGED/PROFILE.db \
                             --splits-of-interest test-output/SAMPLES-MERGED-SUMMARY/bin_by_bin/PSAMPLES_Bin_00001/PSAMPLES_Bin_00001-original_split_names.txt \
                             --genes-of-interest example_genes_of_interest.txt \
                             -o test-output/variability_PSAMPLES_Bin_00001_ALT.txt

INFO "Generating amino acid frequencies for gene caller id 3 in SAMPLE-01.bam ..."
anvi-get-aa-frequencies -i test-output/SAMPLE-01.bam -c test-output/CONTIGS.db --gene-caller-id 3 -o test-output/AA_frequencies_for_gene_caller_id_3.txt

INFO "Getting back the sequence for gene call 3 ..."
anvi-get-dna-sequences-for-gene-calls -c test-output/CONTIGS.db --gene-caller-ids 3 -o test-output/Sequence_for_gene_caller_id_3.fa

INFO "Show all available HMM sources ..."
anvi-get-dna-sequences-for-hmm-hits -c test-output/CONTIGS.db -o /dev/null -l

INFO "Show available gene names in HMM sources ..."
anvi-get-dna-sequences-for-hmm-hits -c test-output/CONTIGS.db -o /dev/null -L

INFO "Get sequences for HMM hits for a bin in a collection ..."
anvi-get-dna-sequences-for-hmm-hits -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -C CONCOCT -b Bin_1 -o test-output/hmm_hits_sequences_in_Bin_1.txt

INFO "Get all ABC transporter hits defined in an HMM source ..."
anvi-get-dna-sequences-for-hmm-hits -c test-output/CONTIGS.db -o test-output/ABC_transporter_hits_in_external_hmm_profile.txt --gene-names ABC_tran --hmm-source external_hmm_profile

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
anvi-get-aa-counts -c test-output/CONTIGS.db -p test-output/SAMPLES-MERGED/PROFILE.db -C CONCOCT -o test-output/aa_counts_for_bins_in_collection_CONCOCT.txt
column -t test-output/aa_counts_for_bins_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile ..."
anvi-get-aa-counts -c test-output/CONTIGS.db -p test-output/SAMPLES-MERGED/PROFILE.db -C CONCOCT -o test-output/aa_counts_for_bin_3_in_collection_CONCOCT.txt -B sample_CONCOCT_bin_id.txt
column -t test-output/aa_counts_for_bin_3_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile ..."
anvi-get-aa-counts -c test-output/CONTIGS.db --contigs-of-interest sample_contig_ids.txt -o test-output/aa_counts_for_two_contigs.txt
column -t test-output/aa_counts_for_two_contigs.txt

INFO "Get AA counts for five genes ..."
anvi-get-aa-counts -c test-output/CONTIGS.db --gene-caller-ids sample_gene_call_ids.txt -o test-output/aa_counts_for_five_genes.txt
column -t test-output/aa_counts_for_five_genes.txt

INFO "Firing up the interactive interface for merged samples ..."
anvi-interactive -p test-output/SAMPLES-MERGED/PROFILE.db \
                 -c test-output/CONTIGS.db \
                 -s test-output/SAMPLES.db \
                 -A additional_view_data.txt \
                 -t test-output/SAMPLES-MERGED/EXP-ORG-FILE.txt \
                 -V additional_view.txt \
                 --split-hmm-layers

INFO "Firing up the interactive interface with the blank profile ..."
anvi-interactive -c test-output/CONTIGS.db -p test-output/BLANK-PROFILE/PROFILE.db

INFO "Firing up the interactive interface in 'COLLECTION' mode ..."
anvi-interactive -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -C CONCOCT

INFO "Firing up the interactive interface to refine a bin ..."
anvi-refine -p test-output/SAMPLES-MERGED/PROFILE.db -c test-output/CONTIGS.db -s test-output/SAMPLES.db -C CONCOCT -b Bin_1
