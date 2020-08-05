This program **downloads and sets up the search databases used for the scg-taxonomy workflow** (from [GTDB](https://gtdb.ecogenomic.org/)) so that you can run %(anvi-run-scg-taxonomy)s or %(anvi-estimate-scg-taxonomy)s. This program generates a %(scgs-taxonomy-db)s artifact, which is required to run both of those programs. 

You will only have to run this program once per anvio installation. 

But why is this not done by default? It just makes things easier downstream to build these databases with the DIAMOND installed on your computer to avoid incompatibility issues. Besides, it should take under a minute and is as simple as running

{{ codestart }}
anvi-setup-scg-databases
{{ codestop }}

If you have already already run this program and are trying to redownload this data, run 

{{ codestart }}
anvi-setup-scg-databases --reset
{{ codestop }}

You can also download a specific release of this database by providing its URL with the flag `--scg-taxonomy-remote-database-url`. 
