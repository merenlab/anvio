This program exports the taxonomy hits for the splits contained in a %(contigs-db)s, outputting them in a %(splits-taxonomy-txt)s. 

To do this, anvi'o examines all of the annotated genes within your splits and returns the taxon ID with the most genes associated with it. For example, a split with 3 genes identified as E. coli, 2 genes identified as Staphylococcus aureus, and 1 as Streptococcus pneumoniae would be annotated as E. coli. 

To run this program, just provide a %(contigs-db)s:

{{ codestart }}
anvi-export-splits-taxonomy -c %(contigs-db)s \
                            -o PATH/TO/%(splits-taxonomy-txt)s

{{ codestop }}
