This script can automatically generate an external or internal genomes file.

## Generating an external genomes file
If you provide an input directory and a name for the output file, then every %(contigs-db)s in that directory will get a line in the resulting %(external-genomes)s file:

```
anvi-script-gen-genomes-file --input-dir path/to/dir -e external_genomes.txt
```

The name of each database will be whatever string is in front of the `*.db` extension, and the `contigs_db_path` column will contain absolute paths.

## Generating an internal genomes file
To get an %(internal-genomes)s file containing all bins from a collection, provide a %(profile-db)s, its corresponding %(contigs-db)s, and the %(collection)s name:

```
anvi-script-gen-genomes-file -i internal_genomes.txt -c CONTIGS.db -p PROFILE.db -C default
```

The name of each internal genome will be the same as the bin name, and the path columns will contain absolute paths.
