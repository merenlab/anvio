An external genome is any genome assembly that was converted into a %(contigs-db)s from its original FASTA file format using the program %(anvi-gen-contigs-database)s. You can obtain one of these in a variety of ways, the most common being 1) downloading a genome from a database such as NCBI and 2) assembling a genome yourself from sequencing reads. The key thing is that the sequences in the %(contigs-db)s represent a _single_ microbial population (or species, if you are not working with microbes) - ie, it is not a metagenome.

The external genomes file format enables anvi'o to work with one or more external genomes. A TAB-delimited external genomes file will be composed of at least the following two columns:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **%(internal-genomes)s** and **%(metagenomes)s**.
