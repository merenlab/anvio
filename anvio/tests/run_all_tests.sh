#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Initializing raw BAM files"
# init raw bam files.
for f in 01 02 03
do
    anvi-init-bam $files/SAMPLE-$f-RAW.bam --output-file $output_dir/SAMPLE-$f.bam
    echo
done

INFO "Reformat the contigs FASTA"
anvi-script-reformat-fasta $files/contigs.fa -o $output_dir/contigs.fa -l 0 --simplify-names --prefix test_prefix --report $output_dir/contigs-reformat-report.txt
echo
column -t $output_dir/contigs-reformat-report.txt

# we first generate an empty contigs database using contigs.fa (keep in mind that 'contigs.fa'
# is the original file all samples were mapped to). here we use split size of 1000 (the default split
# size is much better for most projects. the small split size used here is simply for testing purposes)
INFO "Generating a new contigs database with external gene calls"
anvi-gen-contigs-database -f $files/contigs.fa \
                          -o $output_dir/CONTIGS.db \
                          -L 1000 \
                          --external-gene-calls $files/example_external_gene_calls.txt \
                          --project-name 'Contigs DB with external gene calls'
rm -rf $output_dir/CONTIGS.db

INFO "Generating a new contigs database with the default gene caller"
anvi-gen-contigs-database -f $files/contigs.fa \
                          -o $output_dir/CONTIGS.db \
                          -L 1000 \
                          --project-name "Contigs DB for anvi'o self-tests"

INFO "Displaying the info for the contigs databse"
anvi-db-info $output_dir/CONTIGS.db

INFO "Setting a new self value in the self table of the contigs databse"
anvi-db-info $output_dir/CONTIGS.db \
             --self-key 'a_new_test_key' \
             --self-value "a new test value" \
             --just-do-it

INFO "Exporting gene calls from the contigs database"
anvi-export-gene-calls -c $output_dir/CONTIGS.db -o $output_dir/exported_gene_calls.txt

INFO "Exporting contig sequences from the contigs database"
anvi-export-contigs -c $output_dir/CONTIGS.db -o $output_dir/exported_contig_seqeunces.fa

INFO "Exporting contig sequences from the contigs database in 'splits mode'"
anvi-export-contigs -c $output_dir/CONTIGS.db -o $output_dir/exported_split_seqeunces.fa --splits-mode

INFO "Populating taxonomy for splits table in the database using 'centrifuge' parser"
anvi-import-taxonomy-for-genes -c $output_dir/CONTIGS.db -p centrifuge -i $files/example_files_for_centrifuge_taxonomy/centrifuge_report.tsv $files/example_files_for_centrifuge_taxonomy/centrifuge_hits.tsv

INFO "Trying to remove HMM sources from the contigs database (when there are none in it)"
anvi-delete-hmms -c $output_dir/CONTIGS.db --just-do-it

INFO "Populating HMM hits tables in the latest contigs database using default HMM profiles"
anvi-run-hmms -c $output_dir/CONTIGS.db --num-threads 2

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory"
anvi-run-hmms -c $output_dir/CONTIGS.db -H $files/external_hmm_profile

INFO "Rerunning HMMs for a specific installed profile"
anvi-run-hmms -c $output_dir/CONTIGS.db -I Ribosomal_RNAs

INFO "Listing all available HMM sources in the contigs database"
anvi-delete-hmms -c $output_dir/CONTIGS.db --list

INFO "Removing HMM hits for Rinke_et_al from the contigs database"
anvi-delete-hmms -c $output_dir/CONTIGS.db --hmm-source Rinke_et_al

INFO "Export genomic locus using HMM"
anvi-export-locus -c $output_dir/CONTIGS.db -O $output_dir/exported_locus_from_hmm -n 22,22 -s S-AdoMet_synt_C --use-hmm --hmm-sources Campbell_et_al

INFO "Recovering completeness esimates for the contigs db"
anvi-compute-completeness -c $output_dir/CONTIGS.db

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

INFO "Export genomic locus using functional annotation search"
anvi-export-locus -c $output_dir/CONTIGS.db -O $output_dir/exported_locus_from_functions -n 22,22 -s NusB

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
    anvi-profile -i $output_dir/SAMPLE-$f.bam -o $output_dir/SAMPLE-$f -c $output_dir/CONTIGS.db --profile-SCVs

    INFO "Importing short-read-level taxonomy for SAMPLE-$f"
    anvi-import-taxonomy-for-layers -p $output_dir/SAMPLE-$f/PROFILE.db \
                                    -i $files/example_files_for_kraken_hll_taxonomy/SAMPLE-$f.mpa \
                                    --parser kraken_hll
done


INFO "Merging profiles"
anvi-merge $output_dir/*/PROFILE.db -o $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db --description $files/example_description.md

INFO "Merging profiles without any clustering"
anvi-merge $output_dir/*/PROFILE.db -o $output_dir/SAMPLES-MERGED-WO-CLUSTERING \
                              -c $output_dir/CONTIGS.db \
                              --skip-concoct-binning \
                              --skip-hierarchical-clustering

INFO "Update the description in the merged profile"
anvi-update-db-description $output_dir/SAMPLES-MERGED/PROFILE.db --description $files/example_description.md

INFO "Generating coverages and sequences files for splits (for external binning)"
anvi-export-splits-and-coverages -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db

INFO "Generating per-nt position coverage values for a single split across samples"
anvi-get-split-coverages -p $output_dir/SAMPLES-MERGED/PROFILE.db -o $output_dir/contig_1720_split_00001_coverages.txt --split-name 204_10M_MERGED.PERFECT.gz.keep_contig_1720_split_00001

INFO "Generating per-nt position coverage values for splits in a bin across samples"
anvi-get-split-coverages -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -o $output_dir/split_coverages_in_Bin_1.txt -C CONCOCT -b Bin_1

INFO "Cluster contigs in the newly generated coverages file"
anvi-matrix-to-newick $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt

INFO "Cluster contigs in the newly generated coverages file using 'canberra' distance, and 'complete' linkage"
anvi-matrix-to-newick $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt --distance canberra --linkage complete -o $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs_CANB_COMP.newick

INFO "Generating network descriptions for samples based on gene functions"
anvi-gen-network -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db --annotation-source Pfam

INFO "Use anvi-experimental-organization to generate a tree from a new configuration to store it in a file (not in the database)"
anvi-experimental-organization $files/example_clustering_configuration.ini -i $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db -o $output_dir/SAMPLES-MERGED/EXP-ORG-FILE.txt --skip-store-in-db

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, and add the resulting tree into the database as 'experimental'"
anvi-experimental-organization $files/example_clustering_configuration.ini -i $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db --name experimental

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, but overwrite linkage method and distance metric"
anvi-experimental-organization $files/example_clustering_configuration.ini -i $output_dir/SAMPLES-MERGED -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db --name experimental --distance canberra --linkage complete

INFO "Adding a 'DEFAULT' collection that describes all splits in an 'EVERYTHING' bin to the merged profile"
anvi-script-add-default-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db

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

INFO "Merging Bin_2 and Bin_3 into a new bin in the collection 'CONTIGS_IMPORTED'"
anvi-merge-bins -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                --collection-name 'CONTIGS_IMPORTED' \
                --bin-names-list Bin_2,Bin_3 \
                --new-bin-name merged_bins

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
                 --report-file $output_dir/renaming-report.txt

echo
column -t $output_dir/renaming-report.txt
echo

INFO "Requesting collection info"
anvi-script-get-collection-info -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C "cmdline_concoct_RENAMED"

INFO "Summarizing CONCOCT results"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -o $output_dir/SAMPLES-MERGED-SUMMARY -C 'cmdline_concoct_RENAMED' --init-gene-coverages

INFO "Generate a SNV variabilty profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db -p $output_dir/SAMPLES-MERGED/PROFILE.db -C cmdline_concoct_RENAMED -b PSAMPLES_Bin_00001 -o $output_dir/variability_PSAMPLES_Bin_00001.txt --quince-mode

INFO "Generate a SNV profile for PSAMPLES_Bin_00001 using genes of interest (after summary)"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             --genes-of-interest $files/example_genes_of_interest.txt \
                             -o $output_dir/variability_PSAMPLES_Bin_00001_ALT_GENES_of_INTEREST.txt \
                             --engine NT

INFO "Generate a SNV profile for PSAMPLES_Bin_00001 using split names of interest (after summary)"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             --splits-of-interest $output_dir/SAMPLES-MERGED-SUMMARY/bin_by_bin/PSAMPLES_Bin_00001/PSAMPLES_Bin_00001-original_split_names.txt \
                             -o $output_dir/variability_PSAMPLES_Bin_00001_ALT_SPLITS_of_INTEREST.txt \
                             --engine NT

INFO "Generate a SCV profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             -C cmdline_concoct_RENAMED \
                             -b PSAMPLES_Bin_00001 \
                             -o $output_dir/variability_CDN_PSAMPLES_Bin_00001.txt \
                             --quince-mode \
                             --engine CDN

INFO "Generate a SAAV profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             -C cmdline_concoct_RENAMED \
                             -b PSAMPLES_Bin_00001 \
                             -o $output_dir/variability_AA_PSAMPLES_Bin_00001.txt \
                             --quince-mode \
                             --engine AA

INFO "Generating amino acid frequencies for gene caller id 3 in SAMPLE-01.bam"
anvi-get-codon-frequencies -i $output_dir/SAMPLE-01.bam \
                           -c $output_dir/CONTIGS.db \
                           --gene-caller-id 3 \
                           -o $output_dir/CODON_frequencies_for_gene_caller_id_3.txt

INFO "Generating amino codon frequencies for gene caller id 3 in SAMPLE-01.bam"
anvi-get-codon-frequencies -i $output_dir/SAMPLE-01.bam \
                           -c $output_dir/CONTIGS.db \
                           --gene-caller-id 3 \
                           -o $output_dir/AA_frequencies_for_gene_caller_id_3.txt \
                           --return-AA-frequencies-instead

INFO "Getting back the sequence for gene call 3"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db --gene-caller-ids 3 -o $output_dir/Sequence_for_gene_caller_id_3.fa

INFO "Getting back the AA sequence for gene call 3"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db --gene-caller-ids 3 --get-aa-sequences -o $output_dir/AA_sequence_for_gene_caller_id_3.fa

INFO "Getting back the sequence for gene call 3 (export as GFF3)"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db --gene-caller-ids 3 --export-gff3 -o $output_dir/Sequence_for_gene_caller_id_3.gff

INFO "Export gene coverage and detection data"
anvi-export-gene-coverage-and-detection -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/MERGED

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

INFO "Import layer additional data from file"
anvi-import-misc-data $files/samples-information.txt \
                      -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table layers

INFO "Import layer orders from file"
anvi-import-misc-data $files/samples-order.txt \
                      -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table layer_orders

INFO "Import items additional data from file"
anvi-import-misc-data $files/items_additional_data.txt \
                      -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table items

INFO "Show available misc data keys"
anvi-show-misc-data -p $output_dir/SAMPLES-MERGED/PROFILE.db

INFO "Export items additional data"
anvi-export-misc-data -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table items \
                      -o $output_dir/exported_items_additional_data.txt

INFO "Remove a single data key from items additional data"
anvi-delete-misc-data -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table items \
                      --keys "Özcan's_column"

INFO "Remove all data in items additional data"
anvi-delete-misc-data -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table items \
                      --just-do-it

INFO "Get linkmers from all BAM files for some distant positions"
anvi-report-linkmers --contigs-and-positions $files/distant_positions_for_linkmers.txt -i $output_dir/*.bam -o $output_dir/distant_linkmers.txt

INFO "Get linkmers from all BAM files for some adjacent positions connected with --only-complete-links"
anvi-report-linkmers --contigs-and-positions $files/adjacent_positions_for_linkmers.txt -i $output_dir/*.bam -o $output_dir/adjacent_linkmers.txt --only-complete-links

INFO "Oligotype linkmers report generated for adjacent nucleotide positions"
anvi-oligotype-linkmers -i $output_dir/adjacent_linkmers.txt -o $output_dir/

INFO "Search for functions to get split names with matching genes"
anvi-search-functions -c $output_dir/CONTIGS.db --search transporter,kinase -o $output_dir/transporter-hits.txt --verbose

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
anvi-import-state -p $output_dir/SAMPLES-MERGED-WO-CLUSTERING/PROFILE.db --state $files/example_state.json --name default

INFO "Exporting the state named 'default' from the merged profile"
anvi-export-state -p $output_dir/SAMPLES-MERGED/PROFILE.db --state default -o $output_dir/SAMPLES-MERGED/default_state.json

INFO "Splitting all bins in the CONCOCT collection stored in the merged profile"
anvi-split -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -o $output_dir/CONCOCT_BINS_SPLIT

INFO "Splitting only Bin_1 from the merged profile"
anvi-split -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -o $output_dir/CONCOCT_BINS_SPLIT_ONLY_BIN_1 --bin-id Bin_1

INFO "Listing all collections and bins available in the merged profile"
anvi-show-collections-and-bins -p $output_dir/SAMPLES-MERGED/PROFILE.db

mkdir $output_dir/MCG_CLASSIFIER_OUTPUTS
INFO "Running anvi-mcg-classifier on a single bin in a collection"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_Bin_1 -C CONCOCT -b Bin_1

INFO "Running anvi-mcg-classifier including only SAMPLE-01 and SAMPLE-02"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_INCLUDE -C CONCOCT -b Bin_1 --include-samples $files/samples_to_include_for_mcg.txt

INFO "Running anvi-mcg-classifier excluding SAMPLE-01"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_EXCLUDE -C CONCOCT -b Bin_1 --exclude-samples $files/samples_to_exclude_for_mcg.txt
# 
# INFO "Running anvi-mcg-classifier on a collection"
# anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_CONCOCT -C CONCOCT
# 
INFO 'A dry run with an items order file for the merged profile without any clustering'
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 --items-order $files/example_items_order_file.txt \
                 --dry-run

INFO "Firing up the interactive interface to display the contigs db stats"
anvi-display-contigs-stats $output_dir/CONTIGS.db

INFO "Firing up the interactive interface for merged samples"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -A $files/additional_view_data.txt \
                 -t $output_dir/SAMPLES-MERGED/EXP-ORG-FILE.txt \
                 -V $files/additional_view.txt \
                 --split-hmm-layers

INFO "Firing up the interfae to display the split bin, Bin_1"
anvi-interactive -c $output_dir/CONCOCT_BINS_SPLIT/Bin_1/CONTIGS.db -p $output_dir/CONCOCT_BINS_SPLIT/Bin_1/PROFILE.db --title "Split bin, Bin_1"

INFO "Firing up the interactive interface with the blank profile"
anvi-interactive -c $output_dir/CONTIGS.db -p $output_dir/BLANK-PROFILE/PROFILE.db

INFO "Firing up the interactive interface in 'COLLECTION' mode"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT

INFO "Firing up the interactive interface to refine a bin"
anvi-refine -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -b Bin_1

INFO "Firing up the interactive interface in 'gene' mode"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db -c $output_dir/CONTIGS.db -C CONCOCT -b Bin_1 --gene-mode
