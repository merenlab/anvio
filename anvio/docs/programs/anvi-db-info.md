Displays information about an anvi'o database, and allows users to modify that information when absolutely necessary.

This program is particularly useful for debugging, but also handy in a pinch if you want to check some facts about your database - to answer questions like "did I run HMMs on this %(contigs-db)s yet?" or "is this a merged %(profile-db)s?" This program can also be very dangerous when used to inappropriately modify database information, so if you want to change something, please proceed with caution.

### What information will I see?

All anvi'o databases contain a table of self-describing information known as the "self" table. It helps anvi'o keep track of critical facts such as the type of the database, its version number, and the date it was created. It also saves information about how the database was generated, what sorts of data it contains, what programs have been run on it, and so on. In general, this table exists so that anvi'o can make sure you are doing the right things with your data and that nothing will blow up. `anvi-db-info` will show you the contents of the self table when you run this program on an anvi'o database.

The information in the self table will be different depending on the kind of database you are looking at. For example, a %(contigs-db)s self table will indicate the number of contigs (and splits) in the database, whether or not gene calling was done (and with what gene callers), and which functional annotation sources have been used to annotate the genes. A %(profile-db)s self table will list which samples it contains mapping information for, how many reads where mapped from each sample, and whether or not SNVs have been profiled. A %(modules-db)s (see also %(kegg-data)s) self table will tell you how many KEGG modules are saved in the database and what is the hash value of the database contents. We could go on, but you probably get the picture.

### View information about a database

This is the only way that most people will use this program, and it is very simple. Just provide the path to any anvi'o database to this program, and check the output on your terminal screen:

{{ codestart }}
anvi-db-info path-to-DB.db
{{ codestop }}

Let's be even more specific and say you have a %(contigs-db)s called `CONTIGS.db`. To look at its self table, you would run the following:
{{ codestart }}
anvi-db-info CONTIGS.db
{{ codestop }}

That's it! Easy-peasy lemon-squeezy.

### Example output

Here is an example of what you might see for a %(contigs-db)s.

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

Most of this output is self-explanatory. But one thing that may not be quite obvious to some is that in many cases we use `0` to indicate 'False' and `1` to indicate 'True'. So for this example, you will see that SCG taxonomy was run on this database, but tRNA taxonomy was not.

### Modifying database information
We just need to start by saying - you probably shouldn't do this. Manually changing the values in the self table has the potential to break things downstream because it lets you avoid some of anvi'o's internal sanity checks which prevent you from doing things you shouldn't. If you change things and start running into ugly errors, do not be surprised.

That being said, sometimes you just need to live on the edge and do some hacking, and `anvi-db-info` will let you do that. If a programmer sent you here to update a value in the self table or if you are just foraging ahead on your own, this is how you would do it. Let's change the `project_name` value as an example because it is mostly descriptive and seems fairly safe:

{{ codestart }}
anvi-db-info --self-key project_name --self-value "test" CONTIGS.db
{{ codestop }}

If you run this, you will see a warning telling you what the current value of `project_name` is and what it will be changed to, but the value will not actually be changed just yet. If you are sure you want to do this, you then need to run:

{{ codestart }}
anvi-db-info --self-key project_name --self-value "test" CONTIGS.db  --just-do-it
{{ codestop }}

Then go on your merry adventuring way.
