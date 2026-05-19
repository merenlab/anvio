This program **removes contig classification data from a %(contigs-db)s**.

You must provide `--just-do-it` to confirm the deletion — this is a destructive operation and cannot be undone.

To see what classification sources exist in a database before deleting:

{{ codestart }}
anvi-delete-contig-classification -c %(contigs-db)s \
                                   --list-contig-classification-sources
{{ codestop }}

To delete data from a specific source:

{{ codestart }}
anvi-delete-contig-classification -c %(contigs-db)s \
                                   --source genomad \
                                   --just-do-it
{{ codestop }}

To delete ALL contig classification data from the database:

{{ codestart }}
anvi-delete-contig-classification -c %(contigs-db)s \
                                   --just-do-it
{{ codestop }}
