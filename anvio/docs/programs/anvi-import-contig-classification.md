This program **takes in one or more %(contig-classification-txt)s and stores the contig-level classification data in a %(contigs-db)s**, producing a %(contig-classification)s artifact.

There are many software tools available for classifying contigs according to their predicted domain of origin, and the point of this program is to put those classifications within %(contigs-db)s for use in downstream programs such as %(anvi-split)s. In order for this to work, you'll need to convert the output of whichever tool you have used to the standardized tabular format accepted by this program. Once you have a %(contig-classification-txt)s containing the classification data, you can import that data like this:

{{ codestart }}
anvi-import-contig-classification -c %(contigs-db)s \
                                   -i classification.tsv
{{ codestop }}

Multiple classification sources can coexist in the same contigs database, as described by the `source` column in the input file. If you want to import multiple sources, you can put all their classifications into one %(contig-classification-txt)s.

That said, in case you created a different table for each classification tool that you used, you can import from multiple files at once:

{{ codestart }}
anvi-import-contig-classification -c %(contigs-db)s \
                                   -i genomad_out.tsv tiara_out.tsv
{{ codestop }}

If your input file contains classifications from a source that is already stored in the contigs database, anvi'o will raise an error to protect existing data. To overwrite a source, re-run with the `--just-do-it` flag, which will delete all existing rows for that source before inserting the new data:

{{ codestart }}
anvi-import-contig-classification -c %(contigs-db)s \
                                   -i classification.tsv \
                                   --just-do-it
{{ codestop }}
