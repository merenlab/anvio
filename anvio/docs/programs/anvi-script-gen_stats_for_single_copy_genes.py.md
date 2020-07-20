This program **provides information about each of the single copy core genes in your %(contigs-db)s**. 

Simply provide a %(contigs-db)s, and it will create a %(genes-stats)s file containing a variety of information about the single copy core genes in your database. 

The console output will tell you the total number of contigs, splits, and nucleotides, while the output will tell you the source, name, and e value of each single copy core gene. 

You can provide a specific source for your single copy core genes. To see what sources are availible in your %(contigs-db)s, run 

{{ codestart }}
anvi-script-gen_stats_for_single_copy_core_genes.py -c %(contigs-db)s \
                --list-sources
{{ codestop }}
