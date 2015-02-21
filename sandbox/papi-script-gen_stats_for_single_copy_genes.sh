#!/bin/bash

set -e

echo "############################################################"
echo "# Working on $1" 
echo "############################################################"
papi-gen-annotation-database -f $1 -o $1.db
papi-populate-search-table $1.db
papi-script-gen_stats_for_single_copy_genes.py $1.db
papi-script-gen_stats_for_single_copy_genes.R $1.db.hits $1.db.genes
echo ""
