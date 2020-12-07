This program **downloads and organizes a local copy of the data from EBI's [Pfam database](https://pfam.xfam.org/) for use in function annotation.** This program generates a %(pfams-data)s artifact, which is required to run the program %(anvi-run-pfams)s. 

### Set up Pfams data
{{ codestart }}
anvi-setup-pfams 
{{ codestop }}

By default, this data is stored at `anvio/data/misc/Pfam`. To set up this data in a non-default location, run 
{{ codestart }}
anvi-setup-pfams --pfam-data-dir path/to/location
{{ codestop }}

If you already have a %(pfams-data)s artifact and are trying to redownload this data, run 

{{ codestart }}
anvi-setup-pfams --reset
{{ codestop }}
