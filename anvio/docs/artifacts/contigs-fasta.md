A %(contigs-fasta)s is a %(fasta)s file that is suitable to be used by %(anvi-gen-contigs-database)s to create a %(contigs-db)s.

The most critical requirement for this file is that **it must have simple deflines**. If your %(fasta)s file doesn't have simple deflines, it is not a proper %(contigs-fasta)s. If you intend to use this file with anvi'o, **you must fix your FASTA file prior to mapping**.

Take a look at your deflines prior to mapping, and remove anything that is not a digit, an ASCII letter, an underscore, or a dash character. Here are some example deflines that are not suitable for a %(fasta)s to be considered a %(contigs-fasta)s

``` bash
>Contig-123 length:4567
>Another defline 42
>gi|478446819|gb|JN117275.2|
```

And here are some OK ones:

``` bash
>Contig-123
>Another_defline_42
>gi_478446819_gb_JN117275_2
```

The program %(anvi-script-reformat-fasta)s can do this automatically for you.