An internal genome is any %(bin)s described in an anvi'o %(collection)s stored in an anvi'o %(profile-db)s. You can obtain one of these by binning a metagenome assembly (stored in an anvi'o %(contigs-db)s), which you can do either manually in the interactive interface or automatically with a binning software, and saving or importing it into a %(collection)s.

The internal genomes file format enables anvi'o to work with one or more bins from one or more collections that may be defined in different anvi'o %(profile-db)s files. A TAB-delimited internal genomes file will be composed of at least the following five columns:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|Name_01|Bin_id_01|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_02|Bin_id_02|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_03|Bin_id_03|Collection_B|/path/to/another_profile.db|/path/to/another/contigs.db|
|(...)|(...)|(...)|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **%(external-genomes)s** and **%(metagenomes)s**.
