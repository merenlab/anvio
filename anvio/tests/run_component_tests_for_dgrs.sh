#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the dgrs analysis directory"
cp -r $files/mock_data_for_dgrs/* $output_dir/
cd $output_dir/

INFO "Migrating the contigs database (quietly)"
anvi-migrate 02_CONTIGS/CONTIGS.db \
            --migrate-quickly \
            --quiet

INFO "Migrating the profile database (quietly)"
anvi-migrate 05_MERGED/PROFILE.db \
            --migrate-quickly \
            --quiet

INFO "Running the base analysis (without reporting the activity or genomic context of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running debug the base analysis (without reporting the activity or genomic context of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                --debug \
                $thread_controller

INFO "Running debug amd verbose the base analysis (without reporting the activity or genomic context of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                --debug \
                --verbose \
                $thread_controller

INFO "Running the base analysis (without reporting the activity of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                $thread_controller

INFO "Running the base analysis only looking at adenine bases (without reporting the activity of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC \
                --only-a-bases
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                $thread_controller

INFO "Running the base analysis with discovery-mode meaning that SNVs are included when they are in any codon base position not just the 1st and 2nd (without reporting the activity or genomic context of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                --discovery-mode \
                -o DGRS_BASIC_DISCOVERY \
                --parameter-output \
                --skip-compute-DGR-variability-profiling\
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running the analysis with the collections mode where each VR is in a different bin, meaning each DGR has 1 VR (without reporting the activity of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                --collections-mode \
                -collection-name DGR \
                -o DGRS_COLLECTION_ONE_VR \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                $thread_controller

INFO "Running the analysis with the collections mode where the Trichodesmium DGR has 1 VRs in bin Tricho_VR_1 (without reporting the activity of dgrs)"
anvi-report-dgrs -c 02_CONTIGS/CONTIGS.db \
                -p 05_MERGED/PROFILE.db \
                -I Reverse_Transcriptase \
                --collections-mode \
                -collection-name DGR_tricho_with_2_VRs \
                -o DGRS_COLLECTION_ALL_VRs \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                $thread_controller

INFO ""
