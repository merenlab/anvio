#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1 $2 $3
#####################################

files_dir=$(pwd)/$files

INFO "Setting up the process-genbank test directory"
mkdir -p $output_dir/process_genbank
cd $output_dir/process_genbank

# Data
gbk=$files_dir/data/genbank/mock_genome.gbk

INFO "Running anvi-script-process-genbank"
anvi-script-process-genbank -i $gbk -O MOCK

INFO "Checking output files"
for f in MOCK-contigs.fa MOCK-external-gene-calls.txt MOCK-external-functions.txt; do
    if [ ! -f $f ]; then
        echo "ERROR: Output file $f was not created"
        exit 1
    fi
done

INFO "Verifying gene call types and sources"
# MOCK_001 (CDS) -> call_type 1, source NCBI_PGAP
# MOCK_002 (tRNA) -> call_type 2, source Transfer_RNAs
# MOCK_003 (rRNA) -> call_type 2, source Ribosomal_RNAs
# MOCK_004 (pseudogene) -> call_type 2, source NCBI_PGAP

num_type_1=$(grep -c $'\t'1$'\t' MOCK-external-gene-calls.txt)
num_type_2=$(grep -c $'\t'2$'\t' MOCK-external-gene-calls.txt)

if [ $num_type_1 -ne 1 ]; then
    echo "ERROR: Expected 1 CODING gene (type 1), found $num_type_1"
    exit 1
fi

if [ $num_type_2 -ne 3 ]; then
    echo "ERROR: Expected 3 NONCODING genes (tRNA, rRNA, pseudogene), found $num_type_2"
    exit 1
fi

# Check for canonical sources in gene calls
grep -q "Transfer_RNAs" MOCK-external-gene-calls.txt || { echo "ERROR: Gene call source 'Transfer_RNAs' missing"; exit 1; }
grep -q "Ribosomal_RNAs" MOCK-external-gene-calls.txt || { echo "ERROR: Gene call source 'Ribosomal_RNAs' missing"; exit 1; }

INFO "Verifying functional annotations"
# Should have 4 annotations in total
num_functions=$(tail -n +2 MOCK-external-functions.txt | wc -l)
if [ $num_functions -ne 4 ]; then
    echo "ERROR: Expected 4 functional annotations, found $num_functions"
    exit 1
fi

# Check for canonical sources
grep -q "Transfer_RNAs" MOCK-external-functions.txt || { echo "ERROR: Source 'Transfer_RNAs' missing"; exit 1; }
grep -q "Ribosomal_RNAs" MOCK-external-functions.txt || { echo "ERROR: Source 'Ribosomal_RNAs' missing"; exit 1; }

INFO "Creating a contigs database with external gene calls"
anvi-gen-contigs-database -f MOCK-contigs.fa \
                          --external-gene-calls MOCK-external-gene-calls.txt \
                          -o MOCK.db \
                          --project-name "Test for process-genbank" \
                          --no-progress \
                          $thread_controller

INFO "SUCCESS: anvi-script-process-genbank test passed"
