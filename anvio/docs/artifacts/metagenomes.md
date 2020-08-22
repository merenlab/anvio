A metagenome is any set of sequences that collectively describes multiple different populations (rather than just one genome) and has been converted into a %(contigs-db)s.

The metagenomes file format enables anvi'o to work with one or more metagenomes. A TAB-delimited external genomes file will be composed of at least the following two columns:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

In some cases, (for example when running %(anvi-estimate-scg-taxonomy)s), you may also want to provide the %(profile-db)s that is associated with the %(contigs-db)s. Then the metagenomes file will be composed of three columns:

|name|contigs_db_path|profile_db_path|
|:--|:--|:--|
|Name_01|/path/to/contigs-01.db|/path/to/profile.db|
|Name_02|/path/to/contigs-02.db|/path/to/profile.db|
|Name_03|/path/to/contigs-03.db|/path/to/profile.db|
|(...)|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **%(internal-genomes)s** and **%(external-genomes)s**.
