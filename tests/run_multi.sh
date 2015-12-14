#!/bin/bash
source 00.sh
set -e
cd sandbox/files_for_manual_interactive
if [ ! -d "userdata" ]; then
  mkdir userdata
fi
anvi-server -f fasta.fa -d view_data.txt -A additional_view_data.txt --manual-mode --title "anvi'o Server" -p test.db -t test-output/tree.txt --port-number 80
