This program creates (or drops) a SQLite **index** on a column of a table inside an anvi'o database, so that `WHERE column = ...` lookups on that table become fast.

Anvi'o databases ship **without** column indexes on purpose. For the vast majority of databases the tables are small enough that an index would only take up disk space without buying any noticeable speed. But for the rare very large database -- most notably a merged %(profile-db)s from a read-recruitment experiment, whose `variable_nucleotides` (SNV) table can hold hundreds of millions of rows -- a single, well-chosen index is the difference between a subsecond per-split lookup and a tens-of-minutes full table scan. `anvi-index-table` is how you opt into one of those indexes when you need it, and how you reverse it when you don't.

### Listing what is worth indexing

Anvi'o keeps a small curated list of `(table, column)` combinations that are known to be worth indexing, per database type. To see it, run:

{{ codestart }}
anvi-index-table --list
{{ codestop }}

If you also pass a database, the list is scoped to that database's type and marks the indexes that already exist:

{{ codestart }}
anvi-index-table %(profile-db)s \
                 --list
{{ codestop }}

### Creating an index

To build an index, give the database, the table, and the column (the motivating example -- fast per-split SNV access in a large merged profile database):

{{ codestart }}
anvi-index-table %(profile-db)s \
                 --table variable_nucleotides \
                 --column split_name
{{ codestop }}

Anvi'o will warn you before it starts, since on a very large table this can take several minutes and grow the database file by a few gigabytes. Building an index that already exists is a no-op, so the command is safe to re-run.

To index more than one column at once (a compound index), pass a comma-separated list to `--column`.

If you try to index a `(table, column)` that is not on the curated list, anvi'o will stop and tell you so. If you are sure you want it anyway, add `--just-do-it`:

{{ codestart }}
anvi-index-table %(contigs-db)s \
                 --table genes_in_contigs \
                 --column source \
                 --just-do-it
{{ codestop }}

### Dropping an index

To remove an index you built earlier, add `--drop-index`:

{{ codestart }}
anvi-index-table %(profile-db)s \
                 --table variable_nucleotides \
                 --column split_name \
                 --drop-index
{{ codestop }}

By default this removes the index but does **not** shrink the database file on disk: SQLite frees the space inside the file and reuses it later. To physically reclaim the space right away, add `--reclaim-space`, which runs a `VACUUM`:

{{ codestart }}
anvi-index-table %(profile-db)s \
                 --table variable_nucleotides \
                 --column split_name \
                 --drop-index \
                 --reclaim-space
{{ codestop }}

Be aware that `VACUUM` rewrites the entire database file, so it needs free disk space roughly equal to the current size of the database and can take a long time on very large databases. Anvi'o will warn you before it starts.
