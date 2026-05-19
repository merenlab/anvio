This program **takes in one or more %(contig-classification-txt)s and stores the contig-level classification data in a %(contigs-db)s**, producing a %(contig-classification)s artifact.

All contig names in the input must exist in the contigs database — any mismatch will raise an error.

Multiple sources can coexist in the same contigs database, keyed by the `source` column. If you provide a source that is already stored, anvi'o will raise an error to protect existing data. To overwrite a source, re-run with the `--just-do-it` flag, which will delete all existing rows for that source before inserting the new data.

{{ codestart }}
anvi-import-contig-classification -c %(contigs-db)s \
                                   -i classification.tsv
{{ codestop }}

To import from multiple tools at once:

{{ codestart }}
anvi-import-contig-classification -c %(contigs-db)s \
                                   -i genomad_out.tsv tiara_out.tsv
{{ codestop }}

To overwrite an existing source:

{{ codestart }}
anvi-import-contig-classification -c %(contigs-db)s \
                                   -i classification.tsv \
                                   --just-do-it
{{ codestop }}
