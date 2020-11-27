This program downloads a local copy of a subset of the databases from [GTDB](https://gtdb.ecogenomic.org/) (stored in a %(trna-taxonomy-db)s), so that tRNA sequences in your dataset can be associated with taxonomy information. It is required to run this program before you can run %(anvi-run-trna-taxonomy)s or %(anvi-estimate-trna-taxonomy)s.

Like other `anvi-setup-` programs, this only needs to be run once per anvi'o version. The default path is `anvio/data/misc/TRNA-TAXONOMY`. You can store the resulting %(trna-taxonomy-db)s in a custom location if desired), but then you'll need to provide the path to it whenever you run %(anvi-run-trna-taxonomy)s. 

To run this program, you can simply run

{{ codestart }}
anvi-setup-trna-taxonomy 
{{ codestop }}

If you are trying to redownload these databases, run: 

{{ codestart }}
anvi-setup-trna-taxonomy --reset
{{ codestop }}

Alternatively, you can use `--redo-databases` if you just want to update the database version without redownloading the data. 
