This program downloads a section of the [Protein Data Bank](https://www.rcsb.org/) which is required for strucutral analyses in anvi'o, such as running %(anvi-gen-strcuture-database)s or %(anvi-3-dev)s. 

{{ codestart }}
anvi-setup-pdb-database --just-do-it
{{ codestop }}

If you already have a %(pdb-db)s artifact and are trying to redownload this data, run 

{{ codestart }}
anvi-setup-pdb-database --reset
{{ codestop }}

Or if you just want to update your database, run 

{{ codestart }}
anvi-setup-pdb-database --update
{{ codestop }}

