This program estimates the completeness and redundancy of genomes provided to it, based on domain-level single-copy core genes. 

{:.notice}
Wondering what single-copy core genes anvi'o uses? Check out %(hmm-source)s. It uses the tables populated when you ran %(anvi-run-hmms)s on your %(contigs-db)s. 

Genomes provided to this program must be contained in either a %(bin)s (within a %(collection)s) or a %(contigs-db)s (which can be provided alone or as part of an %(external-genomes)s). 

### Running on contigs databases 

For example, calling 

{{ codestart }}
anvi-estimate-genome-completeness -c %(contigs-db)s 
{{ codestop }}

will output to the terminal the completition and redundancy of the single-copy core genes in your %(contigs-db)s, assuming that all of its contigs represent a single genome. To output this information to a file, you can add the flag `-o` and provide an output path. 

To get this information for several contigs databases at once, you can provide them as an %(external-genomes)s, as so:

{{ codestart }}
anvi-estimate-genome-completeness -e %(external-genomes)s \
                                  -o completition.txt
{{ codestop }}

### Running on bins 

To get this data for a series of bins, just provide a %(profile-db)s and %(collection)s. 

{{ codestart }}
anvi-estimate-genome-completeness -c %(contigs-db)s \
                                  -p %(profile-db)s \
                                  -C %(collection)s 
{{ codestop }}

To see what collections are contained in your contigs database, call 

{{ codestart }}
anvi-estimate-genome-completeness -c %(contigs-db)s \
                                  -p %(profile-db)s \
                                  --list-collections
{{ codestop }}

or run %(anvi-show-collections-and-bins)s for a more comprehensive overview. 

If you're looking for a more comprehensive overview of your entire collection and its contents, the completition and redunduncy statistics for your bins are also included when you run %(anvi-summarize)s. 
