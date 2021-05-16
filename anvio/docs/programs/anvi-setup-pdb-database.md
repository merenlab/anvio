
## Basic usage 

This program creates a %(pdb-db)s local database that holds PDB structures from [this sequence database](https://salilab.org/modeller/supplemental.html), which is hosted by the [Sali lab](https://salilab.org/).  Their database comprises all PDB RCSB sequences that have been clustered at 95%% sequence similarity. They seem to update their database every couple of months (thank you guys!).


The purpose of %(anvi-setup-pdb-database)s to have a local copy of reference structures that can be used to, for example, get template structures for homology modelling when %(anvi-gen-structure-database)s is ran.


Running this program is easy:

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

## Notes

The output %(pdb-db)s database is ~20GB and its contents may take several hours to download.

