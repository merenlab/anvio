This program **downloads and sets up NCBI's [AMRFinderPlus](https://github.com/ncbi/amr/wiki) database** for use with %(anvi-run-amrfinder)s. It produces an %(amrfinder-data)s artifact.

We have not included AMRFinderPlus as a part of the standard anvi'o installation, so if you want to use it you first need to install it into your anvi'o environment:

{{ codestart }}
conda install -c conda-forge -c bioconda ncbi-amrfinderplus
{{ codestop }}

### Basic usage

anvi'o uses AMRFinderPlus's own program, `amrfinder_update`, to download the newest version of the database from NCBI:

{{ codestart }}
anvi-setup-amrfinder
{{ codestop }}

### Custom data directory

By default, anvi'o stores the database at `anvio/data/misc/AMRFinder/`. If you want to store it in a different location you can use:

{{ codestart }}
anvi-setup-amrfinder --amrfinder-data-dir path/to/directory
{{ codestop }}

In case you already have an AMRFinderPlus database on your computer or your computer cluster, you can direct anvi'o to it using %(anvi-run-amrfinder)s with `--amrfinder-data-dir`.

### Updating the database

There are regular updates by NCBI to the AMRFinderPlus database and if you want to make sure that you are using the latest one you can replace the 'old' database in the default location with the newest version by running:

{{ codestart }}
anvi-setup-amrfinder --reset
{{ codestop }}

Anvi'o wants to avoid deleting valuable data in a non-default directory or to reset a database that may be shared with others, therefore you can only use the `--reset` flag for the default directory location. If you set up the database with `--amrfinder-data-dir`, you need to remove that directory yourself before setting it up again.

{:.notice}
%(anvi-run-amrfinder)s records the database version in every %(contigs-db)s it annotates, so you can always find out which version you used with %(anvi-db-info)s.
