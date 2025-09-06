Displays information about an anvi'o database and allows users to modify that information when necessary.

This program is particularly useful for debugging and for quickly verifying database properties - to answer questions such as "Have I run HMMs on this %(contigs-db)s yet?" or "Is this a merged %(profile-db)s?" This program can also be potentially problematic when used to inappropriately modify database information, so if you need to change something, please proceed with caution.

### What information will I see?

All anvi'o databases contain a table of self-describing information known as the "self" table. This table helps anvi'o track critical facts such as the database type, version number, and creation date. It also stores information about how the database was generated, what types of data it contains, which programs have been executed on it, and other relevant metadata. In general, this table exists to ensure that anvi'o can verify you are performing appropriate operations with your data and prevent errors. `anvi-db-info` displays the contents of the self table when you execute this program on an anvi'o database.

The information in the self table varies depending on the type of database you are examining. For example, a %(contigs-db)s self table will indicate the number of contigs (and splits) in the database, whether gene calling has been performed (and with which gene callers), and which functional annotation sources have been used to annotate the genes. A %(profile-db)s self table will list which samples contain mapping information, how many reads were mapped from each sample, and whether SNVs have been profiled. A %(modules-db)s (see also %(kegg-data)s) self table will indicate how many KEGG modules are stored in the database and provide the hash value of the database contents. We could continue with additional examples, but the pattern should be clear.

### View information about a database

This is the primary way most users will interact with this program, and it is straightforward. Simply provide the path to any anvi'o database to this program and review the output displayed on your terminal:

{{ codestart }}
anvi-db-info path-to-DB.db
{{ codestop }}

For a more specific example, if you have a %(contigs-db)s called `CONTIGS.db`, you would examine its self table by executing:
{{ codestart }}
anvi-db-info CONTIGS.db
{{ codestop }}

That's all there is to it! Simple and straightforward.

### Example output

Here is an example of what you might see for a %(contigs-db)s:

```
DB Info (no touch)
===============================================
Database Path ................................: CONTIGS.db
Description ..................................: No description is given
Type .........................................: contigs
Variant ......................................: None
Version ......................................: 20


DB Info (no touch also)
===============================================
contigs_db_hash ..............................: d51abf0a
split_length .................................: 20000
kmer_size ....................................: 4
num_contigs ..................................: 4189
total_length .................................: 35766167
num_splits ...................................: 4784
genes_are_called .............................: 1
splits_consider_gene_calls ...................: 1
creation_date ................................: 1466453807.46107
project_name .................................: Infant Gut Contigs from Sharon et al.
gene_level_taxonomy_source ...................:
scg_taxonomy_was_run .........................: 0
external_gene_calls ..........................: 0
external_gene_amino_acid_seqs ................: 0
skip_predict_frame ...........................: 0
scg_taxonomy_database_version ................: None
trna_taxonomy_was_run ........................: 0
trna_taxonomy_database_version ...............: None
modules_db_hash ..............................: 72700e4db2bc
gene_function_sources ........................: KEGG_Module,COG14_CATEGORY,COG14_FUNCTION,KEGG_Class,KOfam

* Please remember that it is never a good idea to change these values. But in some
cases it may be absolutely necessary to update something here, and a programmer
may ask you to run this program and do it. But even then, you should be
extremely careful.

AVAILABLE GENE CALLERS
===============================================
* 'prodigal' (32,265 gene calls)
* 'Ribosomal_RNAs' (9 gene calls)


AVAILABLE FUNCTIONAL ANNOTATION SOURCES
===============================================
* COG14_CATEGORY (21,121 annotations)
* COG14_FUNCTION (21,121 annotations)
* KEGG_Class (2,760 annotations)
* KEGG_Module (2,760 annotations)
* KOfam (14,391 annotations)


AVAILABLE HMM SOURCES
===============================================
* 'Archaea_76' (type 'singlecopy' with 76 models and 404 hits)
* 'Bacteria_71' (type 'singlecopy' with 71 models and 674 hits)
* 'Protista_83' (type 'singlecopy' with 83 models and 100 hits)
* 'Ribosomal_RNAs' (type 'Ribosomal_RNAs' with 12 models and 9 hits)
```

Most of this output is self-explanatory. However, one aspect that may not be immediately obvious is that in many cases we use `0` to indicate 'False' and `1` to indicate 'True'. For this example, you can see that SCG taxonomy has been run on this database, but tRNA taxonomy has not.

### Modifying database information
We must emphasize - you probably should not do this. Manually changing values in the self table has the potential to cause downstream problems because it allows you to bypass some of anvi'o's internal sanity checks that prevent inappropriate operations. If you modify these values and subsequently encounter unexpected errors, this should not be surprising.

That said, sometimes advanced users need to make modifications, and `anvi-db-info` provides this capability. If a programmer has directed you to update a value in the self table or if you are proceeding independently, this is the process. Let's modify the `project_name` value as an example because it is primarily descriptive and relatively safe:

{{ codestart }}
anvi-db-info --self-key project_name --self-value "test" CONTIGS.db
{{ codestop }}

If you execute this command, you will see a warning indicating the current value of `project_name` and what it will be changed to, but the value will not actually be modified yet. If you are certain you want to proceed, you must then execute:

{{ codestart }}
anvi-db-info --self-key project_name --self-value "test" CONTIGS.db  --just-do-it
{{ codestop }}

Then you may proceed with your analysis.
