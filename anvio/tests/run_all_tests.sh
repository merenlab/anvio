#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2
#####################################

INFO "Initializing raw BAM files"
# init raw bam files.
for f in 01 02 03
do
    anvi-init-bam $files/SAMPLE-$f-RAW.bam \
                  --output-file $output_dir/SAMPLE-$f.bam
    echo
done

INFO "Reformat the contigs FASTA"
anvi-script-reformat-fasta $files/contigs.fa -o $output_dir/contigs.fa \
                                             -l 0 \
                                             --simplify-names \
                                             --prefix test_prefix \
                                             --report $output_dir/contigs-reformat-report.txt
SHOW_FILE $output_dir/contigs-reformat-report.txt

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
anvi-export-gene-calls -c $output_dir/CONTIGS.db \
                       --gene-caller prodigal \
                       -o $output_dir/exported_gene_calls.txt

INFO "Exporting contig sequences from the contigs database"
anvi-export-contigs -c $output_dir/CONTIGS.db \
                    -o $output_dir/exported_contig_sequences.fa

INFO "Exporting contig sequences from the contigs database in 'splits mode'"
anvi-export-contigs -c $output_dir/CONTIGS.db \
                    -o $output_dir/exported_split_sequences.fa \
                    --splits-mode

INFO "Searching for sequence motifs in the contigs database"
anvi-search-sequence-motifs -c $output_dir/CONTIGS.db \
                            --motifs ATCG,TAAAT \
                            --output-file $output_dir/sequence-motifs-in-contigs.txt
SHOW_FILE $output_dir/sequence-motifs-in-contigs.txt

INFO "Populating taxonomy for splits table in the database using 'centrifuge' parser"
anvi-import-taxonomy-for-genes -c $output_dir/CONTIGS.db \
                               -p centrifuge \
                               -i $files/example_files_for_centrifuge_taxonomy/centrifuge_report.tsv $files/example_files_for_centrifuge_taxonomy/centrifuge_hits.tsv

INFO "Trying to remove HMM sources from the contigs database (when there are none in it)"
anvi-delete-hmms -c $output_dir/CONTIGS.db \
                 --just-do-it

INFO "Populating HMM hits tables in the latest contigs database using default HMM profiles"
anvi-run-hmms -c $output_dir/CONTIGS.db \
              --num-threads 2

INFO "Populating HMM hits tables in the latest contigs database using a mock HMM collection from an external directory"
anvi-run-hmms -c $output_dir/CONTIGS.db \
              -H $files/external_hmm_profile

INFO "Generating an ad hoc HMM source from two PFAM accessions (A STEP THAT REQUIRES INTERNET CONNECTION AND RESPONSE FROM XFAM.ORG)"
anvi-script-pfam-accessions-to-hmms-directory --pfam-accessions-list PF00705 PF00706 \
                                              -O $output_dir/ADHOC_HMMs

INFO "Running the HMMs in the ad hoc user directory"
anvi-run-hmms -c $output_dir/CONTIGS.db \
              -H $output_dir/ADHOC_HMMs

INFO "Rerunning HMMs for a specific installed profile"
anvi-run-hmms -c $output_dir/CONTIGS.db \
              -I Ribosomal_RNA_16S \
              --just-do-it

INFO "Rerunning HMMs with hmmsearch"
anvi-run-hmms -c $output_dir/CONTIGS.db \
              -I Bacteria_71 \
              --hmmer-program hmmsearch \
              --hmmer-output-dir $output_dir \
              --get-domtable-output $output_dir \
              --just-do-it

INFO "Filtering hmm_hits using target coverage"
anvi-script-filter-hmm-hits-table -c $output_dir/CONTIGS.db \
                                  --domtblout $output_dir/hmm.domtable \
                                  --hmm-source Bacteria_71 \
                                  --target-coverage 0.9

INFO "Listing all available HMM sources in the contigs database"
anvi-delete-hmms -c $output_dir/CONTIGS.db \
                 --list

INFO "Removing HMM hits for Rinke_et_al from the contigs database"
anvi-delete-hmms -c $output_dir/CONTIGS.db \
                 --hmm-source Rinke_et_al

INFO "Export genomic locus using HMM"
anvi-export-locus -c $output_dir/CONTIGS.db \
                  -O $output_dir/exported_locus_from_hmm \
                  -n 22,22 \
                  -s RNA_pol_Rpb6 \
                  --use-hmm \
                  --hmm-sources Bacteria_71

INFO "Export genomic locus using HMM (multiple, only 1 expected hit)"
anvi-export-locus -c $output_dir/CONTIGS.db \
                  -O $output_dir/exported_locus_from_hmm_multi \
                  -n 22,22 \
                  -s RNA_pol_Rpb6,fake_gene \
                  --use-hmm \
                  --hmm-sources Bacteria_71

INFO "Recovering completeness esimates for the contigs db"
anvi-compute-completeness -c $output_dir/CONTIGS.db

INFO "Importing gene function calls using 'interproscan' parser"
anvi-import-functions -c $output_dir/CONTIGS.db \
                      -i $files/example_interpro_output.tsv \
                      -p interproscan

INFO "Listing all available function call sources in the contigs database"
anvi-export-functions -c $output_dir/CONTIGS.db \
                      --list-annotation-sources

INFO "Deleting gene functions that do not exist along with those that do"
anvi-delete-functions -c $output_dir/CONTIGS.db \
                      --annotation-sources PRINTS,XXX,PIRSF,YYY

INFO "Importing gene function calls INCREMENTALLY using a TAB-delimited default input matrix"
anvi-import-functions -c $output_dir/CONTIGS.db \
                      -i $files/example_gene_functions_input_matrix.txt

INFO "REPLACING gene function calls using a TAB-delimited default input matrix"
anvi-import-functions -c $output_dir/CONTIGS.db \
                      -i $files/example_gene_functions_input_matrix.txt \
                      --drop-previous-annotations

INFO "Export all functional annotations"
anvi-export-functions -c $output_dir/CONTIGS.db \
                      -o $output_dir/exported_functions_from_all_sources.txt
SHOW_FILE $output_dir/exported_functions_from_all_sources.txt

INFO "Export genomic locus using functional annotation search"
anvi-export-locus -c $output_dir/CONTIGS.db \
                  -O $output_dir/exported_locus_from_functions \
                  -n 22,22 \
                  -s NusB \
                  --overwrite-output-destinations

INFO "Export genomic locus using functional annotation search (multiple; 2 expected)"
anvi-export-locus -c $output_dir/CONTIGS.db \
                  -O $output_dir/exported_locus_from_functions \
                  -n 22,22 \
                  -s NusB,'Glutamine amidotransferase class-I' \
                  --overwrite-output-destinations

INFO "Export genomic locus using functional annotation search in flank-mode"
anvi-export-locus -c $output_dir/CONTIGS.db \
                  -O $output_dir/exported_locus_from_functions \
                  --flank-mode  \
                  -s NusB,rpoz \
                  --overwrite-output-destinations

INFO "Export only Pfam annotations"
anvi-export-functions -c $output_dir/CONTIGS.db \
                      -o $output_dir/exported_functions_from_source_Pfam.txt \
                      --annotation-sources Pfam

INFO "Contigs DB is ready; here are the tables in it:"
sqlite3 $output_dir/CONTIGS.db '.tables'

INFO "Generating a 'blank profile' with the newly generated contigs database"
anvi-profile -c $output_dir/CONTIGS.db \
             -o $output_dir/BLANK-PROFILE \
             -S BLANK \
             --blank-profile

INFO "Adding a default collection to the blank profile"
anvi-script-add-default-collection -c $output_dir/CONTIGS.db \
                                   -p $output_dir/BLANK-PROFILE/PROFILE.db \

INFO "Importing a collection into the blank profile"
anvi-import-collection -c $output_dir/CONTIGS.db \
                       -p $output_dir/BLANK-PROFILE/PROFILE.db \
                       -C collection_blank $files/collection_for_blank_profile.txt

INFO "Summarizing the collection stored in the blank profile"
anvi-summarize -c $output_dir/CONTIGS.db \
               -p $output_dir/BLANK-PROFILE/PROFILE.db \
               -C collection_blank \
               -o $output_dir/BLANK-SUMMARY

# for each sample, run the profiling using the same split size used for the contigs database.
# profiling generates individual directiorues uner $output_dir directory for each sample.
for f in 01 02 03
do
    INFO "Profiling sample SAMPLE-$f"
    anvi-profile -i $output_dir/SAMPLE-$f.bam \
                 -o $output_dir/SAMPLE-$f \
                 -c $output_dir/CONTIGS.db \
                 --profile-SCVs

    INFO "Importing short-read-level taxonomy for SAMPLE-$f"
    anvi-import-taxonomy-for-layers -p $output_dir/SAMPLE-$f/PROFILE.db \
                                    -i $files/example_files_for_kraken_hll_taxonomy/SAMPLE-$f.mpa \
                                    --parser krakenuniq
done

# Run anvi-profile on one of the samples using the multi-process routine, just to make sure it does
# not crash. FIXME Ideally, this step would compare the identicalness of
# MULTI-THREAD-SAMPLE-01/PROFILE.db and SAMPLE-01/PROFILE.db 
INFO "Profiling sample SAMPLE-01 with --force-multi"
anvi-profile -i $output_dir/SAMPLE-01.bam \
             -o $output_dir/MULTI-THREAD-SAMPLE-01 \
             -c $output_dir/CONTIGS.db \
             --cluster \
             --profile-SCVs \
             --force-multi
rm -rf $output_dir/MULTI-THREAD-SAMPLE-01

INFO "Merging profiles"
anvi-merge $output_dir/*/PROFILE.db -o $output_dir/SAMPLES-MERGED \
                                    -c $output_dir/CONTIGS.db \
                                    --description $files/example_description.md

INFO "Merging profiles without any clustering"
anvi-merge $output_dir/*/PROFILE.db -o $output_dir/SAMPLES-MERGED-WO-CLUSTERING \
                              -c $output_dir/CONTIGS.db \
                              --skip-hierarchical-clustering

INFO "Importing a collection file into the merged profile"
anvi-import-collection -c $output_dir/CONTIGS.db \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -C CONCOCT \
                       $files/concoct_mini_test.txt

INFO "Update the description in the merged profile"
anvi-update-db-description $output_dir/SAMPLES-MERGED/PROFILE.db \
                           --description $files/example_description.md

INFO "Generating coverages and sequences files for splits (for external binning)"
anvi-export-splits-and-coverages -c $output_dir/CONTIGS.db \
                                 -p $output_dir/SAMPLES-MERGED/PROFILE.db

INFO "Generating per-nt position coverage values for a single split across samples"
anvi-get-split-coverages -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                         -o $output_dir/contig_1720_split_00001_coverages.txt \
                         --split-name 204_10M_MERGED.PERFECT.gz.keep_contig_1720_split_00001

INFO "Generating per-nt position coverage values for splits in a bin across samples"
anvi-get-split-coverages -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                         -c $output_dir/CONTIGS.db \
                         -o $output_dir/split_coverages_in_Bin_1.txt \
                         -C CONCOCT \
                         -b Bin_1
SHOW_FILE $output_dir/split_coverages_in_Bin_1.txt

INFO "Generating per-nt position coverage values for a single gene with its 20nt flanks across samples"
anvi-get-split-coverages -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                         -c $output_dir/CONTIGS.db \
                         -o $output_dir/gene_caller_id_5_coverages.txt \
                         --gene-caller-id 5 \
                         --flank-length 20
SHOW_FILE $output_dir/gene_caller_id_5_coverages.txt

INFO "Cluster contigs in the newly generated coverages file"
anvi-matrix-to-newick $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt

INFO "Cluster contigs in the newly generated coverages file using 'canberra' distance, and 'complete' linkage, and report and additional items order file"
anvi-matrix-to-newick $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs.txt \
                      --distance canberra \
                      --linkage complete \
                      -o $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs_CANB_COMP.newick \
                      --items-order-file $output_dir/SAMPLES-MERGED/SAMPLES_MERGED-COVs_CANB_COMP-ITEMS-ORDER.txt

INFO "Generating network descriptions for samples based on gene functions"
anvi-gen-network -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 --annotation-source Pfam

INFO "Use anvi-experimental-organization to generate a tree from a new configuration to store it in a file (not in the database)"
anvi-experimental-organization $files/example_clustering_configuration.ini \
                               -i $output_dir/SAMPLES-MERGED \
                               -c $output_dir/CONTIGS.db \
                               -o $output_dir/SAMPLES-MERGED/EXP-ORG-FILE.txt \
                               --skip-store-in-db

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, and add the resulting tree into the database as 'experimental'"
anvi-experimental-organization $files/example_clustering_configuration.ini \
                               -i $output_dir/SAMPLES-MERGED \
                               -c $output_dir/CONTIGS.db \
                               -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                               --name experimental

INFO "Use anvi-experimental-organization to generate a tree from a non-default configuration, but overwrite linkage method and distance metric"
anvi-experimental-organization $files/example_clustering_configuration.ini \
                               -i $output_dir/SAMPLES-MERGED \
                               -c $output_dir/CONTIGS.db \
                               -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                               --name experimental \
                               --distance canberra \
                               --linkage complete

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
anvi-export-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -C CONTIGS_IMPORTED \
                       --output-file-prefix $output_dir/exported-collection

INFO "Exporting the 'CONTIGS_IMPORTED' collection that was just imported *with* unbinned items"
anvi-export-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -C CONTIGS_IMPORTED \
                       --output-file-prefix $output_dir/exported-collection \
                       --include-unbinned

INFO "Re-importing a collection from files just exported for CONTIGS_IMPORTED collection (just to confuse you, and to see if we can import stuff we export)"
anvi-import-collection $output_dir/exported-collection.txt \
                       -c $output_dir/CONTIGS.db \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       --collection-name 'CONTIGS_RE_IMPORTED' \
                       --bins-info $output_dir/exported-collection-info.txt

INFO "Deleting the collection 'CONTIGS_RE_IMPORTED'"
anvi-delete-collection -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -C CONTIGS_RE_IMPORTED


INFO "Importing 'cmdline_concoct' collection file into the merged profile"
anvi-import-collection -c $output_dir/CONTIGS.db \
                       -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                       -C cmdline_concoct \
                       $files/concoct_mini_test.txt

INFO "Rename bins in collection 'cmdline_concoct' using SCG averages"
anvi-rename-bins -c $output_dir/CONTIGS.db \
                 -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 --prefix 'PSAMPLES' \
                 --collection-to-read "cmdline_concoct" \
                 --collection-to-write "cmdline_concoct_RENAMED" \
                 --report-file $output_dir/renaming-report.txt
SHOW_FILE $output_dir/renaming-report.txt

INFO "Requesting collection info"
anvi-estimate-genome-completeness -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                  -c $output_dir/CONTIGS.db \
                                  -C "cmdline_concoct_RENAMED"

INFO "Summarizing CONCOCT results"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db \
               -c $output_dir/CONTIGS.db \
               -o $output_dir/SAMPLES-MERGED-SUMMARY \
               -C 'cmdline_concoct_RENAMED' \
               --init-gene-coverages \
               --reformat-contig-names

INFO "Run quick summary of CONCOCT results"
anvi-summarize -p $output_dir/SAMPLES-MERGED/PROFILE.db \
               -c $output_dir/CONTIGS.db \
               -o $output_dir/SAMPLES-MERGED-SUMMARY-QUICK \
               -C 'cmdline_concoct_RENAMED' \
               --reformat-contig-names \
               --quick-summary

INFO "Generate a SNV variabilty profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             -C cmdline_concoct_RENAMED \
                             -b PSAMPLES_Bin_00001 \
                             -o $output_dir/variability_PSAMPLES_Bin_00001.txt \
                             --quince-mode
SHOW_FILE $output_dir/variability_PSAMPLES_Bin_00001.txt

INFO "Generate a SNV profile for PSAMPLES_Bin_00001 using genes of interest (after summary)"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             --genes-of-interest $files/example_genes_of_interest.txt \
                             -o $output_dir/variability_PSAMPLES_Bin_00001_ALT_GENES_of_INTEREST.txt \
                             --engine NT
SHOW_FILE $output_dir/variability_PSAMPLES_Bin_00001_ALT_GENES_of_INTEREST.txt

INFO "Generate a SNV profile for PSAMPLES_Bin_00001 using split names of interest (after summary)"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             --splits-of-interest $output_dir/SAMPLES-MERGED-SUMMARY/bin_by_bin/PSAMPLES_Bin_00001/PSAMPLES_Bin_00001-original_split_names.txt \
                             -o $output_dir/variability_PSAMPLES_Bin_00001_ALT_SPLITS_of_INTEREST.txt \
                             --engine NT
SHOW_FILE $output_dir/variability_PSAMPLES_Bin_00001_ALT_SPLITS_of_INTEREST.txt

INFO "Generate a SCV profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             -C cmdline_concoct_RENAMED \
                             -b PSAMPLES_Bin_00001 \
                             -o $output_dir/variability_CDN_PSAMPLES_Bin_00001.txt \
                             --quince-mode \
                             --engine CDN
SHOW_FILE $output_dir/variability_CDN_PSAMPLES_Bin_00001.txt

INFO "Generate a SAAV profile for PSAMPLES_Bin_00001 using a collection id"
anvi-gen-variability-profile -c $output_dir/CONTIGS.db \
                             -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                             -C cmdline_concoct_RENAMED \
                             -b PSAMPLES_Bin_00001 \
                             -o $output_dir/variability_AA_PSAMPLES_Bin_00001.txt \
                             --quince-mode \
                             --engine AA

INFO "Computing INSeq stats database"
anvi-gen-gene-level-stats-databases -c $output_dir/CONTIGS.db \
                                    -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                    -C DEFAULT \
                                    -b EVERYTHING \
                                    --inseq-stats

INFO "Searching for sequence motifs in the profile database"
anvi-search-sequence-motifs -c $output_dir/CONTIGS.db \
                            -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                            --motifs ATCG,TAAAT \
                            --output-file $output_dir/sequence-motifs-in-profile.txt
SHOW_FILE $output_dir/sequence-motifs-in-profile.txt

INFO "Generating normalized codon frequencies for all genes in the contigs database"
anvi-get-codon-frequencies -c $output_dir/CONTIGS.db \
                           -o $output_dir/CODON_frequencies_for_the_contigs_db.txt \
                           --merens-codon-normalization
SHOW_FILE $output_dir/CODON_frequencies_for_the_contigs_db.txt

INFO "Getting back the sequence for gene call 3"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db \
                                  --gene-caller-ids 3 \
                                  -o $output_dir/Sequence_for_gene_caller_id_3.fa

INFO "Getting back the AA sequence for gene call 3"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db \
                                  --gene-caller-ids 3 \
                                  --get-aa-sequences \
                                  -o $output_dir/AA_sequence_for_gene_caller_id_3.fa

INFO "Getting back the sequence for gene call 3 (export as GFF3)"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db \
                                  --gene-caller-ids 3 \
                                  --export-gff3 \
                                  -o $output_dir/Sequence_for_gene_caller_id_3.gff

INFO "Getting back the sequence for gene call 3 with flank-length of 50 nucleotides"
anvi-get-sequences-for-gene-calls -c $output_dir/CONTIGS.db \
                                  --gene-caller-ids 3 \
                                  --flank-length 50 \
                                  -o $output_dir/Sequence_for_gene_caller_id_3_50nt.fa

INFO "Export gene coverage and detection data"
anvi-export-gene-coverage-and-detection -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                        -c $output_dir/CONTIGS.db \
                                        -O $output_dir/MERGED

INFO "Export gene coverage and detection data for a single gene"
anvi-export-gene-coverage-and-detection -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                        -c $output_dir/CONTIGS.db \
                                        --gene-caller-id 1 \
                                        -O $output_dir/MERGED-GENE-01
SHOW_FILE $output_dir/MERGED-GENE-01-GENE-COVERAGES.txt
SHOW_FILE $output_dir/MERGED-GENE-01-GENE-DETECTION.txt

INFO "Show all available HMM sources"
anvi-get-sequences-for-hmm-hits -c $output_dir/CONTIGS.db \
                                -o /dev/null \
                                -l

INFO "Show available gene names in HMM sources"
anvi-get-sequences-for-hmm-hits -c $output_dir/CONTIGS.db \
                                -o /dev/null \
                                -L

INFO "Get DNA sequences for HMM hits for a bin in a collection"
anvi-get-sequences-for-hmm-hits -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                -c $output_dir/CONTIGS.db \
                                -C CONCOCT \
                                -b Bin_1 \
                                -o $output_dir/hmm_hits_sequences_in_Bin_1.txt

INFO "Get DNA sequences for all ABC transporter hits defined in an HMM source"
anvi-get-sequences-for-hmm-hits -c $output_dir/CONTIGS.db \
                                -o $output_dir/ABC_transporter_hits_in_external_hmm_profile.txt \
                                --gene-names ABC_tran \
                                --hmm-source external_hmm_profile

INFO "Get AA sequences for HMM hits for a bin in a collection"
anvi-get-sequences-for-hmm-hits -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                                -c $output_dir/CONTIGS.db \
                                -C CONCOCT \
                                -b Bin_1 \
                                -o $output_dir/hmm_hits_sequences_in_Bin_1.txt \
                                --get-aa-sequences

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
                      --keys "test_column"

INFO "Remove all data in items additional data"
anvi-delete-misc-data -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                      --target-data-table items \
                      --just-do-it

INFO "Get linkmers from all BAM files for some distant positions"
anvi-report-linkmers --contigs-and-positions $files/distant_positions_for_linkmers.txt \
                     -i $output_dir/*.bam \
                     -o $output_dir/distant_linkmers.txt

INFO "Get linkmers from all BAM files for some adjacent positions connected with --only-complete-links"
anvi-report-linkmers --contigs-and-positions $files/adjacent_positions_for_linkmers.txt \
                     -i $output_dir/*.bam \
                     -o $output_dir/adjacent_linkmers.txt \
                     --only-complete-links

INFO "Oligotype linkmers report generated for adjacent nucleotide positions"
anvi-oligotype-linkmers -i $output_dir/adjacent_linkmers.txt \
                        -o $output_dir/

INFO "Running anvi'o script to correct homopolymer INDELs in --test-run mode"
anvi-script-fix-homopolymer-indels --test-run

INFO "Running anvi'o script to correct homopolymer INDELs with real files"
anvi-script-fix-homopolymer-indels --input $files/single_contig_with_INDEL_errors.fa \
                                   --reference $files/single_contig.fa \
                                   --output $output_dir/single_contig_INDEL_errors_FIXED.fa \
                                   --verbose

INFO "Search for functions to get split names with matching genes"
anvi-search-functions -c $output_dir/CONTIGS.db \
                      --search transporter,kinase \
                      -o $output_dir/transporter-hits.txt \
                      --verbose

INFO "Get all short reads that map to the gene ID 38 (which is a Zinc transpoprter)"
anvi-get-short-reads-mapping-to-a-gene -c $output_dir/CONTIGS.db \
                                       --gene-caller-id 38 \
                                       --leeway 100 \
                                       -i $output_dir/*bam \
                                       -O $output_dir/reads-mapping-to

INFO "Get AA counts for the entire contigs database"
anvi-get-aa-counts -c $output_dir/CONTIGS.db \
                   -o $output_dir/aa_counts_for_contigs_db.txt
SHOW_FILE $output_dir/aa_counts_for_contigs_db.txt

INFO "Get AA counts for bins in collection CONCOCT stored in the merged profile"
anvi-get-aa-counts -c $output_dir/CONTIGS.db \
                   -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                   -C CONCOCT \
                   -o $output_dir/aa_counts_for_bins_in_collection_CONCOCT.txt
SHOW_FILE $output_dir/aa_counts_for_bins_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile"
anvi-get-aa-counts -c $output_dir/CONTIGS.db \
                   -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                   -C CONCOCT \
                   -o $output_dir/aa_counts_for_bin_3_in_collection_CONCOCT.txt \
                   -B $files/sample_CONCOCT_bin_id.txt
SHOW_FILE $output_dir/aa_counts_for_bin_3_in_collection_CONCOCT.txt

INFO "Get AA counts for bin 'bin_3' in collection CONCOCT stored in the merged profile"
anvi-get-aa-counts -c $output_dir/CONTIGS.db \
                   --contigs-of-interest $files/sample_contig_ids.txt \
                   -o $output_dir/aa_counts_for_two_contigs.txt
SHOW_FILE $output_dir/aa_counts_for_two_contigs.txt

INFO "Get AA counts for five genes"
anvi-get-aa-counts -c $output_dir/CONTIGS.db \
                   --gene-caller-ids $files/sample_gene_call_ids.txt \
                   -o $output_dir/aa_counts_for_five_genes.txt
SHOW_FILE $output_dir/aa_counts_for_five_genes.txt

INFO "Importing a state file into the merged profile"
anvi-import-state -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                  --state $files/example_state.json \
                  --name default
anvi-import-state -p $output_dir/SAMPLES-MERGED-WO-CLUSTERING/PROFILE.db \
                  --state $files/example_state.json \
                  --name default

INFO "Exporting the state named 'default' from the merged profile"
anvi-export-state -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                  --state default \
                  -o $output_dir/SAMPLES-MERGED/default_state.json

INFO "Splitting all bins in the CONCOCT collection stored in the merged profile"
anvi-split -p $output_dir/SAMPLES-MERGED/PROFILE.db \
           -c $output_dir/CONTIGS.db \
           -C CONCOCT \
           -o $output_dir/CONCOCT_BINS_SPLIT

INFO "Splitting only Bin_1 from the merged profile"
anvi-split -p $output_dir/SAMPLES-MERGED/PROFILE.db \
           -c $output_dir/CONTIGS.db \
           -C CONCOCT \
           -o $output_dir/CONCOCT_BINS_SPLIT_ONLY_BIN_1 \
           --bin-id Bin_1

INFO "Listing all collections and bins available in the merged profile"
anvi-show-collections-and-bins -p $output_dir/SAMPLES-MERGED/PROFILE.db

mkdir $output_dir/MCG_CLASSIFIER_OUTPUTS
INFO "Running anvi-mcg-classifier on a single bin in a collection"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                    -c $output_dir/CONTIGS.db \
                    -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_Bin_1 \
                    -C CONCOCT \
                    -b Bin_1

INFO "Running anvi-mcg-classifier including only SAMPLE-01 and SAMPLE-02"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                    -c $output_dir/CONTIGS.db \
                    -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_INCLUDE \
                    -C CONCOCT \
                    -b Bin_1 \
                    --include-samples $files/samples_to_include_for_mcg.txt

INFO "Running anvi-mcg-classifier excluding SAMPLE-01"
anvi-mcg-classifier -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                    -c $output_dir/CONTIGS.db \
                    -O $output_dir/MCG_CLASSIFIER_OUTPUTS/MCG_EXCLUDE \
                    -C CONCOCT \
                    -b Bin_1 \
                    --exclude-samples $files/samples_to_exclude_for_mcg.txt

INFO "Generating mock external genome data"
cp $files/mock_data_for_pangenomics/{01,02,03}.fa $output_dir/
cp $files/mock_data_for_pangenomics/external-genomes.txt $output_dir/
cp $files/mock_data_for_pangenomics/functions/*-functions.txt $output_dir/
cp $files/example_description.md $output_dir/

for g in 01 02 03
do
    echo -n "$g .. "
    anvi-gen-contigs-database -f $output_dir/$g.fa \
                              -o $output_dir/$g.db \
                              --project-name $g >/dev/null 2>&1

    anvi-import-functions -c $output_dir/$g.db \
                          -i $output_dir/$g-functions.txt >/dev/null 2>&1
done; echo

INFO "Dereplicating genomes using pyANI"
anvi-dereplicate-genomes -o $output_dir/DEREPLICATION_FROM_SCRATCH \
                         -e $output_dir/external-genomes.txt \
                         --similarity 0.99 \
                         --program pyANI
SHOW_FILE $output_dir/DEREPLICATION_FROM_SCRATCH/CLUSTER_REPORT.txt

INFO "Computing genome similarity"
anvi-compute-genome-similarity -e $output_dir/external-genomes.txt \
                               -o $output_dir/GENOME_SIMILARITY_OUTPUT \
                               --fragment-length 250 \
                               --min-num-fragments 1 \
                               --program pyANI
SHOW_FILE $output_dir/GENOME_SIMILARITY_OUTPUT/ANIb_percentage_identity.txt

INFO "Dereplicating genomes using an existing genome similarity analysis directory"
anvi-dereplicate-genomes --ani-dir $output_dir/GENOME_SIMILARITY_OUTPUT \
                         -o $output_dir/DEREPLICATION_FROM_PREVIOUS_RESULTS \
                         --similarity 0.99 \
                         --program pyANI
SHOW_FILE $output_dir/DEREPLICATION_FROM_PREVIOUS_RESULTS/CLUSTER_REPORT.txt

INFO "Generating an anvi'o genomes storage"
anvi-gen-genomes-storage -e $output_dir/external-genomes.txt -o $output_dir/TEST-GENOMES.db

INFO "Running the pangenome analysis with default parameters"
anvi-pan-genome -g $output_dir/TEST-GENOMES.db \
                -o $output_dir/TEST/ \
                -n TEST \
                --use-ncbi-blast \
                --description $output_dir/example_description.md

INFO "Testing anvi-analyze-synteny with default parameters using a pangenome for annotations"
anvi-analyze-synteny -g $output_dir/TEST-GENOMES.db \
                     -p $output_dir/TEST/TEST-PAN.db \
                     --ngram-window-range 2:3 \
                     -o $output_dir/synteny_output_no_unknowns.tsv

INFO "Testing anvi-analyze-synteny now including unannotated genes"
anvi-analyze-synteny -g $output_dir/TEST-GENOMES.db \
                     -p $output_dir/TEST/TEST-PAN.db \
                     --ngram-window-range 2:3 \
                     -o $output_dir/synteny_output_with_unknowns.tsv \
                     --analyze-unknown-functions

INFO "Testing anvi-analyze-synteny now including unannotated genes"
anvi-analyze-synteny -g $output_dir/TEST-GENOMES.db \
                     --annotation-source COG_FUNCTION \
                     --ngram-window-range 2:3 \
                     -o $output_dir/synteny_output_with_COGs.tsv \
                     --analyze-unknown-functions

INFO 'A dry run with an items order file for the merged profile without any clustering'
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 --items-order $files/example_items_order_file.txt \
                 --dry-run

INFO "A dry run in 'gene-mode' to store gene-level coverage stats in a new genes database"
rm -rf $output_dir/SAMPLES-MERGED/GENES/*
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -C CONCOCT \
                 -b Bin_1 \
                 --gene-mode \
                 --dry-run

INFO "A dry run to fill in anvi'o dbs"
curdir=`pwd`
cd $output_dir
anvi-interactive --dry-run
anvi-display-pan --dry-run
cd $curdir

INFO "Firing up the interactive interface to display the contigs db stats"
anvi-display-contigs-stats $output_dir/CONTIGS.db \
                           $dry_run_controller

INFO "Firing up the interactive interface for merged samples"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -A $files/additional_view_data.txt \
                 -t $output_dir/SAMPLES-MERGED/EXP-ORG-FILE.txt \
                 -V $files/additional_view.txt \
                 --split-hmm-layers \
                 $dry_run_controller


INFO "Firing up the interface to display the split bin, Bin_1"
anvi-interactive -c $output_dir/CONCOCT_BINS_SPLIT/Bin_1/CONTIGS.db \
                 -p $output_dir/CONCOCT_BINS_SPLIT/Bin_1/PROFILE.db \
                 --title "Split bin, Bin_1" \
                 $dry_run_controller


INFO "Firing up the interactive interface with the blank profile"
anvi-interactive -c $output_dir/CONTIGS.db \
                 -p $output_dir/BLANK-PROFILE/PROFILE.db \
                 $dry_run_controller


INFO "Firing up the interactive interface in 'COLLECTION' mode"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -C CONCOCT \
                 $dry_run_controller

INFO "Firing up the interactive interface to refine a bin"
anvi-refine -p $output_dir/SAMPLES-MERGED/PROFILE.db \
            -c $output_dir/CONTIGS.db \
            -C CONCOCT \
            -b Bin_1 \
            $dry_run_controller


INFO "Importing items and layers additional data into the genes database for CONCOCT::Bin_1"
anvi-import-misc-data -p $output_dir/SAMPLES-MERGED/GENES/CONCOCT-Bin_1.db \
                     $files/items_addtl_data_gene_mode.txt \
                     -t items
anvi-import-misc-data -p $output_dir/SAMPLES-MERGED/GENES/CONCOCT-Bin_1.db \
                     $files/layers_addtl_data_gene_mode.txt \
                     -t layers

INFO "Firing up the interactive interface in 'gene' mode"
anvi-interactive -p $output_dir/SAMPLES-MERGED/PROFILE.db \
                 -c $output_dir/CONTIGS.db \
                 -C CONCOCT \
                 -b Bin_1 \
                 --gene-mode \
                 $dry_run_controller
