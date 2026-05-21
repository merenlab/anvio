Creates individual, self-contained anvi'o projects for one or more %(bin)ss stored in an anvi'o %(collection)s. This program may be useful if you would like to share a subset of an anvi'o project with the community or a collaborator, or focus on a particular aspect of your data without having to initialize very large files. Altogether, %(anvi-split)s promotoes reproducibility, openness, and collaboration.

The program can generate %(split-bins)s from metagenomes or pangenomes. To split bins, you can provide the program %(anvi-split)s with a %(contigs-db)s and %(profile-db)s pair. To split gene clusters, you can provide it with a %(genomes-storage-db)s and %(pan-db)s pair. In both cases you will also need a %(collection)s. If you don't provide any %(bin)s names, the program will create individual directories for each bin that is found in your collection. You can also limit the output to a single bin. Each of the resulting directories in your output folder will contain a stand-alone anvi'o project that can be shared without sharing any of the larger dataset.

### An example run

Assume you have a %(profile-db)s has a %(collection)s with three bins, which are (very creatively) called `BIN_1`, `BIN_2`, and `BIN_3`.

If you ran the following code:

{{ codestart }}
anvi-split -p %(profile-db)s \
           -c %(contigs-db)s \
           -C %(collection)s \
           -o OUTPUT
{{ codestop }}

Alternatively you can specify a bin name to limit the reported bins:

{{ codestart }}
anvi-split -p %(profile-db)s \
           -c %(contigs-db)s \
           -C %(collection)s \
           --bin-id BIN_1
           -o OUTPUT
{{ codestop }}

Similarly, if you provide a %(genomes-storage-db)s and %(pan-db)s pair, the directories will contain their own smaller %(genomes-storage-db)s and %(pan-db)s pairs.

You can always use the program %(anvi-show-collections-and-bins)s to learn available %(collection)s and %(bin)s names in a given %(profile-db)s or %(pan-db)s.

### Performance

For extremely large datasets, splitting bins may be difficult. For metagenomics projets you can,

* Use the flag `--skip-variability-tables` to NOT report single-nucleotide variants or single-amino acid variants in your split bins (which can reach hundreds of millions of lines of information for large and complex metagenomes), and/or,
* Use the flag `--compress-auxiliary-data` to save space. While this is a great option for data that is meant to be stored long-term and shared with the community, the compressed file would need to be manually decompressed by the end-user prior to using the split bin.
