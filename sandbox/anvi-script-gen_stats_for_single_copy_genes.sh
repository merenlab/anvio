#!/bin/bash

set -e

if [ -z "$1"  ]
then
    echo "Pleaes send a FASTA file as a parameter..";
    exit -1
fi

echo "############################################################"
echo "# Working on $1" 
echo "############################################################"
anvi-gen-contigs-database -f $1 -o $1.db
anvi-run-hmms -c $1.db
anvi-script-gen_stats_for_single_copy_genes.py $1.db
anvi-script-gen_stats_for_single_copy_genes.R $1.db.hits $1.db.genes
echo ""
