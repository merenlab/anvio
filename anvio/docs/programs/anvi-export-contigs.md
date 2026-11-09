This program **exports the contig sequences from a %(contigs-db)s**, outputting them as a %(contigs-fasta)s. It also has the ability to output the sequences of your splits instead.

You can run this program as follows, which will return ALL contigs in a given %(contigs-db)s file:

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s
{{ codestop }}

You can also limit the contigs you may be interested in to a subset by providing the list of contig names you wish to export in a file. For example:

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s \
                    --contigs-of-interest my_favorite_contigs.txt
{{ codestop }}

where `my_favorite_contigs.txt` looks like this:

```
contig_0001
contig_0005
contig_0035
```

Alternatively, you may be interested in contigs that include one or more *genes* you are interested. In that case you can use `--genes-of-interest` with a %(genes-of-interest-txt)s, and anvi'o will export only those contigs that contain at least one of the gene calls you listed:

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s \
                    --genes-of-interest my_favorite_genes.txt
{{ codestop }}

where `my_favorite_genes.txt` looks like this:

```
5
13
206
```

Please note that `--contigs-of-interest` and `--genes-of-interest` are mutually exclusive: you can use one or the other in a given command, but not both at the same time.

### Splits mode

Want to look at your splits instead of your contigs? Just run with the flag `splits-mode` attached.

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s \
                    --splits-mode
{{ codestop }}
