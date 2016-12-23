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

INFO "Reformat the contigs FASTA"
anvi-script-reformat-fasta $files/contigs.fa -o $output_dir/contigs.fa -l 0 --simplify-names --prefix test_prefix --report $output_dir/contigs-reformat-report.txt
echo
column -t $output_dir/contigs-reformat-report.txt 

# we first generate an empty contigs database using contigs.fa (keep in mind that 'contigs.fa'
# is the original file all samples were mapped to). here we use split size of 1000 (the default split
# size is much better for most projects. the small split size used here is simply for testing purposes)
INFO "Generating an EMPTY contigs database"
anvi-gen-contigs-database -f $files/contigs.fa -o $output_dir/CONTIGS.db -L 1000

INFO "Exporting gene calls from the contigs database"
anvi-export-gene-calls -c $output_dir/CONTIGS.db -o $output_dir/exported_gene_calls.txt

INFO "Populating taxonomy for splits table in the database using 'centrifuge' parser"
anvi-import-taxonomy -c $output_dir/CONTIGS.db -p centrifuge -i $files/example_files_for_centrifuge_taxonomy/*

INFO "Populating HMM hits tables in the latest contigs database using default HMM profiles"
anvi-run-hmms -c $output_dir/CONTIGS.db --num-threads 2

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory"
anvi-run-hmms -c $output_dir/CONTIGS.db -H $files/external_hmm_profile

INFO "Importing gene function calls using 'interproscan' parser"
anvi-import-functions -c $output_dir/CONTIGS.db -i $files/example_interpro_output.tsv -p interproscan

INFO "Importing gene function calls INCREMENTALLY using a TAB-delimited default input matrix"
anvi-import-functions -c $output_dir/CONTIGS.db -i $files/example_gene_functions_input_matrix.txt

INFO "REPLACING gene function calls using a TAB-delimited default input matrix"
anvi-import-functions -c $output_dir/CONTIGS.db -i $files/example_gene_functions_input_matrix.txt --drop-previous-annotations

INFO "Listing all available function call sources in the contigs database"
anvi-export-functions -c $output_dir/CONTIGS.db --list-annotation-sources

INFO "Export all functional annotations"
anvi-export-functions -c $output_dir/CONTIGS.db -o $output_dir/exported_functions_from_all_sources.txt
echo
head $output_dir/exported_functions_from_all_sources.txt | tr ' ' @@ | column -t | tr @@ ' '

INFO "Export only Pfam annotations"
anvi-export-functions -c $output_dir/CONTIGS.db -o $output_dir/exported_functions_from_source_Pfam.txt --annotation-sources Pfam
echo
head $output_dir/exported_functions_from_source_Pfam.txt | tr ' ' @@ | column -t | tr @@ ' '

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 $output_dir/CONTIGS.db '.tables'

INFO "Generating a 'blank profile' with the newly generated contigs database"
anvi-profile -c $output_dir/CONTIGS.db -o $output_dir/BLANK-PROFILE -S BLANK --blank-profile

INFO "Importing a collection into the blank profile"
anvi-import-collection -c $output_dir/CONTIGS.db -p $output_dir/BLANK-PROFILE/PROFILE.db -C collection_blank $files/collection_for_blank_profile.txt

INFO "Summarizing the collection stored in the blank profile"
anvi-summarize -c $output_dir/CONTIGS.db -p $output_dir/BLANK-PROFILE/PROFILE.db -C collection_blank -o $output_dir/BLANK-SUMMARY

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner $output_dir directory for each sample.
for f in 01 02 03
do
    INFO "Profiling sample SAMPLE-$f"
    anvi-profile -i $output_dir/SAMPLE-$f.bam -o $output_dir/SAMPLE-$f -c $output_dir/CONTIGS.db --profile-AA-frequencies
    echo
done


INFO "Merging profiles"
anvi-merge $output_dir/SAMPLE-*/RUNINFO.cp -o $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db

INFO "Add a new variable into the RUNINFO file of the merged profile"
anvi-script-update-runinfo-variable $output_dir/SAMPLES-MERGED/RUNINFO.mcp --variable TEST-VARIABLE --set-bool true

INFO "Generating coverages and sequences files for splits (for external binning)"
anvi-export-splits-and-coverages -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db

INFO "Cluster contigs in the newly generated coverages file"
anvi-matrix-to-newick $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt

INFO "Cluster contigs in the newly generated coverages file using 'canberra' distance, and 'complete' linkage"
anvi-matrix-to-newick $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt --distance canberra --linkage complete -o $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs_CANB_COMP.newick

INFO "Generating network descriptions for samples based on ORFs and functions"
anvi-gen-network -r $output_dir/SAMPLES-MERGED/RUNINFO.mcp -c $output_dir/CONTIGS.db

INFO "Use anvi-experimental-organization to generate a tree from a new configuration to store it in a file (not in the database)"
anvi-experimental-organization $files/example_clustering_configuration.ini -i $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db -o $output_dir/SAMPLES-MERGED/EXP-ORG-FILE.txt --skip-store-in-db

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, and add the resulting tree into the database as 'experimental'"
anvi-experimental-organization $files/example_clustering_configuration.ini -i $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db --name experimental

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, but overwrite linkage method and distance metric"
anvi-experimental-organization $files/example_clustering_configuration.ini -i $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db --name experimental --distance canberra --linkage complete

INFO "Importing external binning results for splits into the contigs database as 'SPLITS_IMPORTED_INTO_CONTIGS_DB'"
anvi-import-collection $files/example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -c $output_dir/CONTIGS.db \
                       --collection-name 'SPLITS_IMPORTED_INTO_CONTIGS_DB' \
                       --bins-info $files/example_files_for_external_binning_results/example_bins_info_file.txt

INFO "Importing external binning results for splits into the profile database as 'SPLITS_IMPORTED'"
anvi-import-collection $files/example_files_for_external_binning_results/external_binning_of_splits.txt \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -c $output_dir/CONTIGS.db \
                       --collection-name 'SPLITS_IMPORTED' \
                       --bins-info $files/example_files_for_external_binning_results/example_bins_info_file.txt

INFO "Importing external binning results for splits into the profile database as 'CONTIGS_IMPORTED'"
anvi-import-collection $files/example_files_for_external_binning_results/external_binning_of_contigs.txt \
                       -c $output_dir/CONTIGS.db \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_IMPORTED' \
                       --bins-info $files/example_files_for_external_binning_results/example_bins_info_file.txt \
                       --contigs-mode

INFO "Exporting the 'CONTIGS_IMPORTED' collection that was just imported"
anvi-export-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db -C CONTIGS_IMPORTED --output-file-prefix $output_dir/exported-collection

INFO "Exporting the 'CONTIGS_IMPORTED' collection that was just imported *with* unbinned items"
anvi-export-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db -C CONTIGS_IMPORTED --output-file-prefix $output_dir/exported-collection --include-unbinned

INFO "Re-importing a collection from files just exported for CONTIGS_IMPORTED collection (just to confuse you, and to see if we can import stuff we export)"
anvi-import-collection $output_dir/exported-collection.txt \
                       -c $output_dir/CONTIGS.db \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_RE_IMPORTED' \
                       --bins-info $output_dir/exported-collection-info.txt

INFO "Deleting the collection 'CONTIGS_RE_IMPORTED'"
anvi-delete-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db -C CONTIGS_RE_IMPORTED

INFO "Use CONCOCT to cluster splits in the merged profile and export as a text file..."
anvi-cluster-with-concoct -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -o $output_dir/anvio_concoct_clusters.txt --collection-name 'cmdline_concoct'

INFO "Recover short reads for Bin_2 in CONCOCT collection and store them in a FASTA file"
anvi-get-short-reads-from-bam -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -b Bin_2 -o $output_dir/short_reads_for_Bin_2.fasta $output_dir/*bam

INFO "Rename bins in collection 'cmdline_concoct' using SCG averages"
anvi-rename-bins -c $output_dir/CONTIGS.db \
                 -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 --prefix 'PSAMPLES' \
                 --collection-to-read "cmdline_concoct" \
                 --collection-to-write "cmdline_concoct_RENAMED" \
                 --use-SCG-averages \
                 --report-file $output_dir/renaming-report.txt

echo
column -t $output_dir/renaming-report.txt
echo

INFO "Requesting collection info"
anvi-script-get-collection-info -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C "cmdline_concoct_RENAMED"

INFO "Summarizing CONCOCT results"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -o $output_dir/SAMPLES-MERGED-SUMMARY -C 'cmdline_concoct_RENAMED'

INFO "Generate a variabilty profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db -C cmdline_concoct_RENAMED -b PSAMPLES_Bin_00001 -o $output_dir/variability_PSAMPLES_Bin_00001.txt --quince-mode

INFO "Generate a variabilty profile for PSAMPLES_Bin_00001 using split ids and gene ids of interest (after summary)"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             --splits-of-interest $output_dir/SAMPLES-MERGED-SUMMARY/bin_by_bin/PSAMPLES_Bin_00001/PSAMPLES_Bin_00001-original_split_names.txt \
                             --genes-of-interest $files/example_genes_of_interest.txt \
                             -o $output_dir/variability_PSAMPLES_Bin_00001_ALT.txt

INFO "Generating amino acid frequencies for gene caller id 3 in SAMPLE-01.bam"
anvi-get-aa-frequencies -i $output_dir/SAMPLE-01.bam -c $output_dir/CONTIGS.db --gene-caller-id 3 -o $output_dir/AA_frequencies_for_gene_caller_id_3.txt

INFO "Getting back the sequence for gene call 3"
anvi-get-dna-sequences-for-gene-calls -c $output_dir/CONTIGS.db --gene-caller-ids 3 -o $output_dir/Sequence_for_gene_caller_id_3.fa

INFO "Show all available HMM sources"
anvi-get-sequences-for-hmm-hits -c $output_dir/CONTIGS.db -o /dev/null -l

INFO "Show available gene names in HMM sources"
anvi-get-sequences-for-hmm-hits -c $output_dir/CONTIGS.db -o /dev/null -L

INFO "Get DNA sequences for HMM hits for a bin in a collection"
anvi-get-sequences-for-hmm-hits -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -b Bin_1 -o $output_dir/hmm_hits_sequences_in_Bin_1.txt

INFO "Get DNA sequences for all ABC transporter hits defined in an HMM source"
anvi-get-sequences-for-hmm-hits -c $output_dir/CONTIGS.db -o $output_dir/ABC_transporter_hits_in_external_hmm_profile.txt --gene-names ABC_tran --hmm-source external_hmm_profile

INFO "Get AA sequences for HMM hits for a bin in a collection"
anvi-get-sequences-for-hmm-hits -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -b Bin_1 -o $output_dir/hmm_hits_sequences_in_Bin_1.txt --get-aa-sequences

INFO "Generate a samples information database with samples information and samples order"
anvi-gen-samples-info-database -D $files/samples-information.txt -R $files/samples-order.txt -o $output_dir/SAMPLES.db

INFO "Get linkmers from all BAM files for some distant positions"
anvi-report-linkmers --contigs-and-positions $files/distant_positions_for_linkmers.txt -i $output_dir/*.bam -o $output_dir/distant_linkmers.txt

INFO "Get linkmers from all BAM files for some adjacent positions connected with --only-complete-links"
anvi-report-linkmers --contigs-and-positions $files/adjacent_positions_for_linkmers.txt -i $output_dir/*.bam -o $output_dir/adjacent_linkmers.txt --only-complete-links

INFO "Oligotype linkmers report generated for adjacent nucleotide positions"
anvi-oligotype-linkmers -i $output_dir/adjacent_linkmers.txt -o $output_dir/

INFO "Search for functions to get split names with matching genes"
anvi-search-functions-in-splits -c $output_dir/CONTIGS.db --search transporter,kinase -o $output_dir/transporter-hits.txt --verbose

INFO "Get all short reads that map to the gene ID 38 (which is a Zinc transpoprter)"
anvi-get-short-reads-mapping-to-a-gene -c $output_dir/CONTIGS.db --gene-caller-id 38 --leeway 100 -i $output_dir/*bam -o $output_dir/reads-mapping-to-gene-id-38.fa

INFO "Get AA counts for the entire contigs database"
anvi-get-aa-counts -c $output_dir/CONTIGS.db -o $output_dir/aa_counts_for_contigs_db.txt
column -t $output_dir/aa_counts_for_contigs_db.txt

INFO "Get AA counts for bins in collection CONCOCT stored in the merged profile"
anvi-get-aa-counts -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db -C CONCOCT -o $output_dir/aa_counts_for_bins_in_collection_CONCOCT.txt
column -t $output_dir/aa_counts_for_bins_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile"
anvi-get-aa-counts -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db -C CONCOCT -o $output_dir/aa_counts_for_bin_3_in_collection_CONCOCT.txt -B $files/sample_CONCOCT_bin_id.txt
column -t $output_dir/aa_counts_for_bin_3_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile"
anvi-get-aa-counts -c $output_dir/CONTIGS.db --contigs-of-interest $files/sample_contig_ids.txt -o $output_dir/aa_counts_for_two_contigs.txt
column -t $output_dir/aa_counts_for_two_contigs.txt

INFO "Get AA counts for five genes"
anvi-get-aa-counts -c $output_dir/CONTIGS.db --gene-caller-ids $files/sample_gene_call_ids.txt -o $output_dir/aa_counts_for_five_genes.txt
column -t $output_dir/aa_counts_for_five_genes.txt

INFO "Importing a state file into the merged profile"
anvi-import-state -p $output_dir/SAMPLES-MERGED/PROFILE.db --state $files/example_state.json --name default

INFO "Exporting the state named 'default' from the merged profile"
anvi-export-state -p $output_dir/SAMPLES-MERGED/PROFILE.db --state default -o $output_dir/SAMPLES-MERGED/default_state.json

INFO "Firing up the interactive interface for merged samples"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -s $output_dir/SAMPLES.db \
                 -A $files/additional_view_data.txt \
                 -t $output_dir/SAMPLES-MERGED/EXP-ORG-FILE.txt \
                 -V $files/additional_view.txt \
                 --split-hmm-layers

INFO "Firing up the interactive interface with the blank profile"
anvi-interactive -c $output_dir/CONTIGS.db -p $output_dir/BLANK-PROFILE/PROFILE.db

INFO "Firing up the interactive interface in 'COLLECTION' mode"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT

INFO "Firing up the interactive interface to refine a bin"
anvi-refine -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -s $output_dir/SAMPLES.db -C CONCOCT -b Bin_1
