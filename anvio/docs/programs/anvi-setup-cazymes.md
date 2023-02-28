This program **downloads and organizes a local copy of the data from [dbCAN2 CAZyme HMMs](https://bcb.unl.edu/dbCAN2/download/Databases/) for use in function annotation.** This program generates a %(cazyme-data)s artifact, which is required to run the program %(anvi-run-cazymes)s. 

### Set up cazymes data

anvi'o will download the newest version of the database (V11) by default:

{{ codestart }}
anvi-setup-cazymes 
{{ codestop }}

You can use `--cazyme-version`, if you want anvi'o to download a different version of the [dbCAN2 CAZyme HMMs](https://bcb.unl.edu/dbCAN2/download/Databases/) database:

{:.warning}
The following versions have been tested for download: V9, V10, V11

{{ codestart }}
anvi-setup-cazymes --cazyme-version V10
{{ codestop }}

By default, this data is stored at `anvio/data/misc/CAZyme/`. To set up this data in a non-default location, run:

{{ codestart }}
anvi-setup-cazymes --cazyme-data-dir path/to/location
{{ codestop }}

If you already have a %(cazyme-data)s artifact and are trying to re-download this data, run:

{{ codestart }}
anvi-setup-cazymes --reset
{{ codestop }}
