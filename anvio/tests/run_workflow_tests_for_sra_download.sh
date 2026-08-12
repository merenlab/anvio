#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

ASSERT_FILE_EXISTS() {
    if [ ! -f "$1" ]; then
        echo "FAIL: expected the file '$1' to exist, but it does not."
        exit 1
    fi
}

ASSERT_FILE_MISSING() {
    if [ -f "$1" ]; then
        echo "FAIL: expected the file '$1' to be gone, but it is still there."
        exit 1
    fi
}

INFO "Setting up the sra_download workflow test directory"
mkdir $output_dir/workflow_test
cp -r $files/workflows/sra_download/* $output_dir/workflow_test/
cd $output_dir/workflow_test

INFO "Creating a default config for sra_download workflow"
anvi-run-workflow -w sra_download --get-default-config sra_download_config.json

INFO "Listing dependencies for sra_download workflow"
anvi-run-workflow -w sra_download -c sra_download_config.json --list-dependencies

INFO "Saving a workflow graph sra_download workflow"
anvi-run-workflow -w sra_download -c sra_download_config.json --save-workflow-graph

# A decoy FASTQ whose name starts with the accession SRR5965623 but belongs to no accession in
# the list. Compressing FASTQ files by globbing on the accession would sweep this file up along
# with the real ones, so it must still be sitting here, uncompressed, once the workflow is done.
INFO "Planting a decoy FASTQ file to make sure accessions only claim their own reads"
mkdir -p 02_FASTQ
touch 02_FASTQ/SRR59656231_1.fastq

INFO "Running workflow graph sra_download workflow"
anvi-run-workflow -w sra_download -c sra_download_config.json

INFO "Making sure the decoy FASTQ file was left alone"
ASSERT_FILE_EXISTS 02_FASTQ/SRR59656231_1.fastq
ASSERT_FILE_MISSING 02_FASTQ/SRR59656231_1.fastq.gz

INFO "Making sure every accession got the FASTQ files it was supposed to get"
ASSERT_FILE_EXISTS 02_FASTQ/ERR6450080_1.fastq.gz
ASSERT_FILE_EXISTS 02_FASTQ/ERR6450080_2.fastq.gz
ASSERT_FILE_EXISTS 02_FASTQ/ERR6450081_1.fastq.gz
ASSERT_FILE_EXISTS 02_FASTQ/ERR6450081_2.fastq.gz
ASSERT_FILE_EXISTS 02_FASTQ/SRR5965623.fastq.gz

INFO "Making sure no uncompressed FASTQ file was left behind"
ASSERT_FILE_MISSING 02_FASTQ/ERR6450080_1.fastq
ASSERT_FILE_MISSING 02_FASTQ/ERR6450081_1.fastq
ASSERT_FILE_MISSING 02_FASTQ/SRR5965623.fastq

# `Remove_unzipped_SRA_files` defaults to true, so the prefetched archives must be gone once the
# reads have been extracted from them, while the marker files snakemake tracks in the very same
# directory must survive (otherwise a second run would download everything all over again).
INFO "Making sure the prefetched SRA archives were cleaned up"
if [ -n "$(find 01_NCBI_SRA -name '*.sra' -o -name '*.sralite')" ]; then
    echo "FAIL: prefetched SRA archives are still on disk:"
    find 01_NCBI_SRA -name '*.sra' -o -name '*.sralite'
    exit 1
fi
ASSERT_FILE_EXISTS 01_NCBI_SRA/SRR5965623/SRR5965623.prefetch.done

INFO "Making sure a second run of the workflow does not download anything again"
anvi-run-workflow -w sra_download -c sra_download_config.json
if [ -n "$(find 01_NCBI_SRA -name '*.sra' -o -name '*.sralite')" ]; then
    echo "FAIL: the second run of the workflow downloaded the SRA archives again."
    exit 1
fi

INFO "The samples-txt files the workflow generated"
SHOW_FILE samples_single_reads.txt
SHOW_FILE samples.txt
