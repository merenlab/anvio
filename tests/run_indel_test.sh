#!/bin/bash
source 00.sh
set -e

INFO "Creating the output directory"
# change directory and clean the old mess if it exists
cd sandbox
rm -rf test-output
mkdir test-output

INFO "Anvo'o version"
anvi-profile --version

INFO "Generating a Bowtie2 ref"
bowtie2-build files_for_indel_testing/contig.fa test-output/contig

INFO "Generating anvi'o contigs database"
anvi-gen-contigs-database -f files_for_indel_testing/contig.fa -o test-output/contig.db -L -1 -n "Contig" --skip-gene-calling

INFO "Generating short reads for sample $sample"
anvi-script-gen-short-reads files_for_indel_testing/sample.ini --output-file-path test-output/sample.fa

INFO "Mapping short reads to the ref"
../misc/bowtie_batch_single_fasta.sh test-output/sample.fa test-output/sample test-output/contig

INFO "Profiling the BAM file"
anvi-profile -i test-output/sample.bam -c test-output/contig.db -o test-output/sample-PROFILE -M 0 --cluster-contigs

INFO "Importing the state"
anvi-import-state -s files_for_indel_testing/sample.json -p test-output/sample-PROFILE/PROFILE.db -n default

INFO "Visualizing the profile"
anvi-interactive -c test-output/contig.db -p test-output/sample-PROFILE/PROFILE.db
