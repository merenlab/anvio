This program **downloads and organizes a local copy of the dbCAN [CAZyme HMMs database](https://bcb.unl.edu/dbCAN2/download/Databases/) for function annotation.** This program generates a %(cazymes-data)s artifact, which is required to run the program %(anvi-run-cazymes)s. 

### Set up CAZyme data
{{ codestart }}
anvi-setup-cazymes 
{{ codestop }}

By default, this data is stored at `anvio/data/misc/CAZyme/`. To set up this data in a non-default location, run the code code below. This is a good way to download multiple versions of the CAZyme database.
{{ codestart }}
anvi-setup-cazymes --cazyme-data-dir path/to/location
{{ codestop }}

If you already have a %(cazyme-data)s artifact and are trying to re-download this data, run 

{{ codestart }}
anvi-setup-cazymes --reset
{{ codestop }}
