The primary purpose of this script is to reduce the amount of labor required to generate %(external-genomes)s or %(internal-genomes)s files anvi'o typically uses to learn about your bins and/or genomes.

## Generating an external genomes file

If you provide an input directory and a name for the output file, then every %(contigs-db)s in that directory will get a line in the resulting %(external-genomes)s file:

```
anvi-script-gen-genomes-file --input-dir path/to/dir \
                             --output-file external_genomes.txt
```

Names for genomes in the resulting external genomes file will be set based on the `project_name` variable, and the `contigs_db_path` column will contain absolute paths.

{:.notice}
You can learn the current `project_name` and/or change it for a given %(contigs-db)s using the program %(anvi-db-info)s. This variable is set by the program %(anvi-gen-contigs-database)s.

You can also instruct `anvi-script-gen-genomes-file` to include all subdirectories under a given directory path:

```
anvi-script-gen-genomes-file --input-dir path/to/dir \
                             --output-file external_genomes.txt \
                             --include-subdirs
```

## Generating an internal genomes file

To get an %(internal-genomes)s file containing all bins from a collection, provide a %(profile-db)s, its corresponding %(contigs-db)s, and the %(collection)s name:

{{ codestart }}
anvi-script-gen-genomes-file -c %(contigs-db)s \
                             -p %(profile-db)s \
                             -C %(collection)s \
                             --output-file internal-genomes.txt
{{ codestop }}

The name of each internal genome will be the same as the bin name, and the path columns will contain absolute paths.
