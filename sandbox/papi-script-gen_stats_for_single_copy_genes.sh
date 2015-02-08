#!/bin/bash

echo "############################################################"
echo "#                  Working on $1                           #" 
echo "############################################################"
papi-populate-search-table $1 $1.db -L 10000
papi-script-gen_stats_for_single_copy_genes.py $1 $1.db
papi-script-gen_stats_for_single_copy_genes.R $1.db.hits $1.db.genes
echo ""
