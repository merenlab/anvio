An external genome is anything you have in a FASTA file format (i.e., a genome you have downloaded from NCBI, or obtained through any other way) that was converted into a %(contigs-db)s.

The external genomes file format enables anvi'o to work with one or more external genomes. A TAB-delimited external genomes file will be composed of at least the following two columns:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

## Additional columns

In some cases additional columns may be required to be in this file. Below is a table of the possible columns you may need.

| header | description | required for |
|----|----|----|
| group | which group the genome belongs to (can be empty) | %(anvi-compute-enrichment-scores)s |

Also see **%(internal-genomes)s** and **%(metagenomes)s**.
