In the anvi'o lingo, an external genome is any %(contigs-db)s generated from a FASTA file that describes a single genome for a single microbial population (and not a metagenome).

The purpose of the external genomes file is to describe one or more external genomes, so this file can be passed to anvi'o programs that can operate on multiple genomes.

For a given set of %(contigs-db)s files, you can generate an external-genomes file automatically using the program %(anvi-script-gen-genomes-file)s. Alternatively, you can manually create the file using a text editor, or a program like EXCEL.

The external-genomes file is a TAB-delimited file with at least two columns (you can add more columns to this file, and anvi'o will not mind):

* `name`. The name of the external genome. You can call it anything, but you should keep it to a single word witout any spaces or funny characters.
* `contigs_db_path`. The full path to each %(contigs-db)s file (tip: the command `pwd` will tell you the full path to the directory you are in).

The format of the file should look like this:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **%(internal-genomes)s** and **%(metagenomes)s**.
