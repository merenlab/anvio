This program **downloads and organizes a local copy of the data from NCBI's [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) for use in function annotation.** This program generates a %(cogs-data)s artifact, which is required to run the program %(anvi-run-ncbi-cogs)s. 

### Set up COGs data
{{ codestart }}
anvi-setup-ncbi-cogs --just-do-it
{{ codestop }}

If you already have a %(cogs-data)s artifact and are trying to redownload this data, run 

{{ codestart }}
anvi-setup-ncbi-cogs --reset
{{ codestop }}
