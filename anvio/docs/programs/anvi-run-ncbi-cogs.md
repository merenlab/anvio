This program **associates genes in your %(contigs-db)s with functions using NCBI's [Clusters of Orthologus Groups (COGs) database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/).**

Before you run this program, you'll have to set up the COGs database on your computer with the program %(anvi-setup-ncbi-cogs)s.  

As mentioned above, the COGs database is no longer actively added to, so might also want to consider using a separate database. As of yet, anvi'o does not have a program to accesss the eggNOG database (instructions to use this database to get function information are [here](http://merenlab.org/2016/06/18/importing-functions/#eggnog-database--emapper)), but does have the functionality to use the Pfams database (check out %(anvi-run-pfams)s for more information). 

To run, you'll need to provide a %(contigs-db)s. If you stored the %(cogs-data)s that you got from running %(anvi-setup-ncbi-cogs)s in a custom location, you'll need to provide that path as well. The output is a %(functions)s artifact. 

By default, this program uses DIAMOND in the "fast" setting for database searching. To instead run in "sensitive" mode, just call: 

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --cog-data-dir %(cogs-data)s \
            --sensitive
{{ codestop }}

You can also use blastp to search, by running: 

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --search-with blastp
{{ codestop }}

*Note: without the flag `--cog-data-dir`, anvi'o will just search in the default location.*


