#!/bin/bash
source 00.sh
set -e

INFO "Creating the output directory"
# change directory and clean the old mess if it exists
cd sandbox
rm -rf test-output
mkdir test-output
mkdir test-output/small
mkdir test-output/large

INFO "Anvo'o version"
anvi-profile --version

#########################################################

INFO "Generating a Bowtie2 ref"
bowtie2-build files_for_indel_testing/small/contig.fa test-output/small/contig

INFO "Generating anvi'o contigs database"
anvi-gen-contigs-database -f files_for_indel_testing/small/contig.fa -o test-output/small/contig.db -L -1 -n "Contig" --skip-gene-calling

INFO "Generating short reads for sample $sample"
anvi-script-gen-short-reads files_for_indel_testing/small/sample.ini --output-file-path test-output/small/sample.fa

INFO "Mapping short reads to the ref"
../misc/bowtie_batch_single_fasta.sh test-output/small/sample.fa test-output/small/sample test-output/small/contig

INFO "Profiling the BAM file"
anvi-profile -i test-output/small/sample.bam -c test-output/small/contig.db -o test-output/small/sample-PROFILE -M 0 --cluster-contigs

INFO "Importing the state"
anvi-import-state -s files_for_indel_testing/small/sample.json -p test-output/small/sample-PROFILE/PROFILE.db -n default

INFO "Visualizing the profile"
anvi-interactive -c test-output/small/contig.db -p test-output/small/sample-PROFILE/PROFILE.db

#########################################################

INFO "Generating a Bowtie2 ref"
bowtie2-build files_for_indel_testing/large/contig.fa test-output/large/contig
bowtie2-build files_for_indel_testing/large/contig_T300del.fa test-output/large/contig_T300del

INFO "Generating anvi'o contigs database"
anvi-gen-contigs-database -f files_for_indel_testing/large/contig.fa -o test-output/large/contig.db -L -1 -n "No del" --skip-gene-calling
anvi-gen-contigs-database -f files_for_indel_testing/large/contig_T300del.fa -o test-output/large/contig_T300del.db -L -1 -n "T300del" --skip-gene-calling

INFO "Generating short reads for sample $sample"
anvi-script-gen-short-reads files_for_indel_testing/large/sample.ini --output-file-path test-output/large/sample.fa
anvi-script-gen-short-reads files_for_indel_testing/large/sample_T300del.ini --output-file-path test-output/large/sample_T300del.fa

INFO "Mapping short reads to the ref"
../misc/bowtie_batch_single_fasta.sh test-output/large/sample.fa test-output/large/sample test-output/large/contig
../misc/bowtie_batch_single_fasta.sh test-output/large/sample_T300del.fa test-output/large/sample_T300del test-output/large/contig_T300del

INFO "Mapping short reads onto opposite ref"
../misc/bowtie_batch_single_fasta.sh test-output/large/sample.fa test-output/large/sample_norm_on_del test-output/large/contig_T300del
../misc/bowtie_batch_single_fasta.sh test-output/large/sample_T300del.fa test-output/large/sample_del_on_norm test-output/large/contig

INFO "Profiling the BAM files"
anvi-profile -i test-output/large/sample.bam             -c test-output/large/contig.db         -o test-output/large/sample-PROFILE -M 0 --cluster-contigs
anvi-profile -i test-output/large/sample_del_on_norm.bam -c test-output/large/contig.db         -o test-output/large/sample-PROFILE_del_on_norm -M 0 --cluster-contigs

anvi-profile -i test-output/large/sample_T300del.bam     -c test-output/large/contig_T300del.db -o test-output/large/sample-PROFILE_T300del -M 0 --cluster-contigs
anvi-profile -i test-output/large/sample_norm_on_del.bam -c test-output/large/contig_T300del.db -o test-output/large/sample-PROFILE_norm_on_del -M 0 --cluster-contigs

INFO "Profiling the BAM files"
anvi-merge test-output/large/sample-PROFILE/PROFILE.db test-output/large/sample-PROFILE_del_on_norm/PROFILE.db -c test-output/large/contig.db -o test-output/large/MERGED --skip-concoct-binning
anvi-merge test-output/large/sample-PROFILE_T300del/PROFILE.db test-output/large/sample-PROFILE_norm_on_del/PROFILE.db -c test-output/large/contig_T300del.db -o test-output/large/MERGED_T300del --skip-concoct-binning

INFO "Visualizing the profile"
anvi-interactive -c test-output/large/contig.db -p test-output/large/MERGED/PROFILE.db
anvi-interactive -c test-output/large/contig_T300del.db -p test-output/large/MERGED_T300del/PROFILE.db

