#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

INFO "Setting up the metagenomics analysis directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/metagenomics/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "unzipping fasta files"
gzip -d three_samples_example/*.fa.gz

#making an uncompressed copy of the fastqs for idba_ud
mkdir three_samples_example_fastq
cp three_samples_example/*fastq.gz three_samples_example_fastq/
gzip -d three_samples_example_fastq/*fastq.gz
for f in `ls three_samples_example_fastq/*fastq`; do
    mv $f ${f%.fastq}-for-idba_ud.fastq
done
mv three_samples_example_fastq/* three_samples_example/

INFO "Creating a default config for metagenomics workflow"
anvi-run-snakemake-workflow -w metagenomics --get-default-config default-config.json

INFO "List dependencies for metagenomics workflow with megahit"
anvi-run-snakemake-workflow -w metagenomics --config config-megahit.json --list-dependencies

INFO "Runnind a dry run for metagenomics workflow with megahit"
anvi-run-snakemake-workflow -w metagenomics --config config-megahit.json --dry-run

INFO "Running metagenomics workflow using megahit"
anvi-run-snakemake-workflow -w metagenomics --config config-megahit.json

INFO "Runnind a dry run for metagenomics workflow with megahit with no qc"
anvi-run-snakemake-workflow -w metagenomics --config config-megahit-no-qc-all-against-all.json --dry-run

INFO "List dependencies for metagenomics workflow in references mode"
anvi-run-snakemake-workflow -w metagenomics --config config-references-mode.json --list-dependencies

INFO "Runnind a dry run in references mode"
anvi-run-snakemake-workflow -w metagenomics --config config-references-mode.json --dry-run

INFO "Runnind a dry run in references mode with no qc and no gzip"
anvi-run-snakemake-workflow -w metagenomics --config config-references-mode-no-qc-no-gzip-no-groups.json --dry-run

INFO "List dependencies for metagenomics workflow with idba_ud"
anvi-run-snakemake-workflow -w metagenomics --config config-idba_ud.json --list-dependencies

INFO "Runnind a dry run with idba_ud"
anvi-run-snakemake-workflow -w metagenomics --config config-idba_ud.json --dry-run
