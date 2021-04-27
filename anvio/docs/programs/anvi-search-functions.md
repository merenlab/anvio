This program **searches for keywords in the function annotations of your database.** 

You can use this program to look for specific functon keywords in a %(contigs-db)s, %(genomes-storage-db)s or %(pan-db)s. For example, say you wanted your %(contigs-db)s to search for genes that encoded some type of kinase. You could call 

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms kinase
{{ codestop }}

By default, the output will be a fairly barren %(functions-txt)s, only telling you which contigs contain genes that matched your search. This will be most helpful as an additional layer in the anvi'o interactive interface, so you can quickly see where the kinase-encoding genes are in the genome. To do this, run %(anvi-interactive)s with the `--aditional-layer` parameter with the %(functions-txt)s. 

However, you can also request a much more comprehensive output that contains much more information, including the matching genes' caller id, functional annotation source and full function name. 

For example, to run the same search as above, but with a more comprehensive output, you could call 

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms kinase \
                      --full-report kinase_information.txt \
                      --include-sequences \
                      --verbose
{{ codestop }}

Following this run, the file `kinase_information.txt` will contain comprehensive information about the matching genes, including their sequences. 

You can also search for multiple terms at the same time, or for terms from only specific annotation sources. For example, if you only wanted Pfam hits with functions related to kinases or phosphatases, you could call 

{{ codestart }}
anvi-search-functions -c %(contigs-db)s \
                      --search-terms kinase,phosphatase \
                      --annotation-sources Pfam \ 
                      --full-report kinase_phosphatase_information.txt
{{ codestop }}
