In the anvi'o lingo, an internal genome is any %(bin)s stored in an anvi'o %(collection)s that describes a single genome. You can obtain one of these by binning a metagenome manually in the interactive interface, automatically using a binning software, or by importing a %(collection)s into anvi'o using the program %(anvi-import-collection)s.

The purpose of the external genomes file is to describe one or more internal genomes genomes, so this file can be passed to anvi'o programs that can operate on multiple genomes. The internal genomes file format enables anvi'o programs to work with one or more bins from one or more collections that may be defined in different anvi'o %(profile-db)s files.

The internal-genomes file is a TAB-delimited file with at least the following five columns:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|Name_01|Bin_id_01|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_02|Bin_id_02|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_03|Bin_id_03|Collection_B|/path/to/another_profile.db|/path/to/another/contigs.db|
|(...)|(...)|(...)|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **%(external-genomes)s** and **%(metagenomes)s**.
