#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

INFO "Setting up the dgrs analysis directory"
cp $files/mock_data_for_dgrs/reference.fa $output_dir/
cp $files/mock_data_for_dgrs/mutations.tsv $output_dir/
cd $output_dir/

####################################################################################################
#
#   STEP 1: GENERATE SYNTHETIC READS WITH MULTI-ALLELE SNVs AT DGR MISMATCH POSITIONS
#
#   The reference.fa contains 6 contigs from 3 real DGR systems:
#
#     DGR_001 (activity-only, no RT gene):
#       - SAMEA2619974_000000008595_DGR_no_RT: VR + TR on same contig
#
#     DGR_002 (B. fragilis, with RT gene):
#       - B_fragilis_ARW016_000000000001_VR_001_rev_comp: VR
#       - B_fragilis_ARW016_000000000001_TR_RT: TR + RT gene
#
#     DGR_003 (T. erythraeum, with RT gene, 3 VRs for 1 TR):
#       - T_erythraeum_IMS101_000000000001_RT_TR: TR + RT gene + VR_001 (same contig)
#       - T_erythraeum_IMS101_000000000001_VR_002: VR_002
#       - T_erythraeum_IMS101_000000000001_VR_001: VR_003 (NO SNVs -> homology-only)
#
#   mutations.tsv contains multi-allele SNVs (3-4 bases) at VR mismatch positions where
#   the TR has adenine (A-base mutagenesis). VR_003 has no SNVs to test homology-only detection.
#
####################################################################################################

INFO "Generating synthetic paired-end reads for sample 1 (seed=1)"
anvi-script-gen-reads -f reference.fa \
                      -o sample_01 \
                      --preset illumina-paired \
                      --coverage 50 \
                      --mutations-file mutations.tsv \
                      --seed 1

INFO "Generating synthetic paired-end reads for sample 2 (seed=2)"
anvi-script-gen-reads -f reference.fa \
                      -o sample_02 \
                      --preset illumina-paired \
                      --coverage 50 \
                      --mutations-file mutations.tsv \
                      --seed 2

INFO "Generating synthetic paired-end reads for sample 3 (seed=3)"
anvi-script-gen-reads -f reference.fa \
                      -o sample_03 \
                      --preset illumina-paired \
                      --coverage 50 \
                      --mutations-file mutations.tsv \
                      --seed 3

####################################################################################################
#
#   STEP 2: BUILD CONTIGS DATABASE AND RUN HMMs
#
####################################################################################################

INFO "Building the contigs database"
anvi-gen-contigs-database -f reference.fa \
                          -o CONTIGS.db \
                          -n "DGR test contigs"

INFO "Running HMMs (for Reverse Transcriptase detection)"
anvi-run-hmms -c CONTIGS.db \
              -I Reverse_Transcriptase \
              --quiet

####################################################################################################
#
#   STEP 3: MAP READS, INIT BAM, PROFILE, AND MERGE
#
####################################################################################################

INFO "Building bowtie2 index"
bowtie2-build --quiet reference.fa bt2_index

for sample in sample_01 sample_02 sample_03; do
    INFO "Mapping ${sample} with bowtie2"
    bowtie2 -x bt2_index \
            -1 ${sample}-R1.fastq \
            -2 ${sample}-R2.fastq \
            -S ${sample}.sam \
            --no-unal \
            --quiet

    INFO "Initializing BAM for ${sample}"
    anvi-init-bam ${sample}.sam \
                  -o ${sample}.bam \
                  --quiet

    INFO "Profiling ${sample}"
    anvi-profile -i ${sample}.bam \
                 -c CONTIGS.db \
                 -o ${sample}_profile \
                 --quiet \
                 $thread_controller
done

INFO "Merging profiles"
anvi-merge sample_01_profile/PROFILE.db \
           sample_02_profile/PROFILE.db \
           sample_03_profile/PROFILE.db \
           -c CONTIGS.db \
           -o MERGED_PROFILE \
           --quiet

# Create samples.txt for the variability profiling step
printf "sample\tr1\tr2\n" > samples.txt
printf "sample_01\tsample_01-R1.fastq\tsample_01-R2.fastq\n" >> samples.txt
printf "sample_02\tsample_02-R1.fastq\tsample_02-R2.fastq\n" >> samples.txt
printf "sample_03\tsample_03-R1.fastq\tsample_03-R2.fastq\n" >> samples.txt

####################################################################################################
#
#   STEP 4: RUN anvi-report-dgrs WITH VARIOUS FLAGS
#
####################################################################################################

INFO "Running the base analysis (without reporting the activity or genomic context of dgrs)"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running the base analysis with --debug"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC_DEBUG \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                --debug \
                $thread_controller

INFO "Running the base analysis with --debug and --verbose"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC_DEBUG_VERBOSE \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                --debug \
                --verbose \
                $thread_controller

INFO "Running the base analysis with genomic context recovery"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_GENOMIC_CONTEXT \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                $thread_controller

INFO "Running the base analysis allowing any dominant base, not just adenine"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_BASIC_ALLOW_ANY_BASE \
                --allow-any-base \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                $thread_controller

INFO "Running the base analysis with discovery-mode (SNVs from any codon position)"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                --discovery-mode \
                -o DGRS_BASIC_DISCOVERY \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running with metagenome mode (restricts VR/TR pairs to same contig)"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_METAGENOME_MODE \
                --metagenome-mode \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running with detection-mode 'activity' only"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_ACTIVITY_ONLY \
                --detection-mode activity \
                --parameter-output \
                --skip-compute-DGR-variability-profiling \
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running with detection-mode 'homology' only"
anvi-report-dgrs -c CONTIGS.db \
                -I Reverse_Transcriptase \
                -o DGRS_HOMOLOGY_ONLY \
                --detection-mode homology \
                --parameter-output \
                --skip-recovering-genomic-context \
                $thread_controller

INFO "Running the full analysis with activity reporting and variable primers"
anvi-report-dgrs -c CONTIGS.db \
                -p MERGED_PROFILE/PROFILE.db \
                -I Reverse_Transcriptase \
                -o DGRS_ACTIVITY_VARIABLE_PRIMERS \
                --parameter-output \
                --samples-txt samples.txt \
                $thread_controller
