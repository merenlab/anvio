#!/bin/bash
source ../../00.sh

set -e

INFO "Anvo'o version ..."
anvi-profile --version

INFO "Running anvi-pan-genome ..."
anvi-pan-genome -i test-output/contig_dbs.txt -o test-output/pan-output --num-threads 10 --min-percent 80
