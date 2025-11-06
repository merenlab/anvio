This basically stores **a local copy of the data from the NCBI [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) for function annotation.** 

It is required to run %(anvi-run-ncbi-cogs)s and is set up on your computer by the program %(anvi-setup-ncbi-cogs)s. 

## Database version options

Last we checked, the following versions of NCBI COGs were available:
- `COG20` (the default): the 2020 release of the COGs database, described in [Galperin et al. 2021](https://doi.org/10.1093/nar/gkaa1018)
- `COG14`: the 2014 release of the COGs database, described in [Galperin et al. 2015](https://doi.org/10.1093/nar/gku1223)
- `arCOG14`: the archael COGs database from 2014, described in [Makarova, Wolf, and Koonin 2015](https://doi.org/10.3390/life5010818)

Please see the help page for %(anvi-setup-ncbi-cogs)s if you want to learn how to choose which version to set up on your computer.
