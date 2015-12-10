#!/bin/bash
source 00.sh
set -e
cd sandbox/files_for_manual_interactive
if [ ! -d "userdata" ]; then
  mkdir userdata
fi
anvi-interactive -f fasta.fa -d view_data.txt -A additional_view_data.txt --manual-mode --title 'Manual Interactive' -p test.db -t test-output/tree.txt --port-number 80 --server-only --multiuser
