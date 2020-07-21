This program **downloads and organizes a local copy of the data from NCBI's [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) for use in function annotation.** This program generates a %(cogs-data)s artifact, which is required to run the program %(anvi-run-ncbi-cogs)s. 

{:notice}
The COGs database is no longer actively added to, so you might also want to consider using a separate database for more comprehensive functional annotation. As of yet, anvi'o does not have a program to accesss the eggNOG database (instructions to use this database to get function information are [here](http://merenlab.org/2016/06/18/importing-functions/#eggnog-database--emapper)), but does have the functionality to use the Pfams database (check out %(anvi-run-pfams)s for more information). 

### Set up COGs data
{{ codestart }}
anvi-setup-ncbi-cogs --just-do-it
{{ codestop }}

If you already have a %(cogs-data)s artifact and are trying to redownload this data, run 

{{ codestart }}
anvi-setup-ncbi-cogs --reset
{{ codestop }}
