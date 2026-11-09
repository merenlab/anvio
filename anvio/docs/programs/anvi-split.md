Creates individual, self-contained anvi'o projects for one or more %(bin)ss stored in an anvi'o %(collection)s. This program may be useful if you would like to share a subset of an anvi'o project with the community or a collaborator, or focus on a particular aspect of your data without having to initialize very large files. Altogether, %(anvi-split)s promotoes reproducibility, openness, and collaboration.

The program can generate %(split-bins)s from metagenomes, from pangenomes, or from a %(contigs-db)s on its own (without a %(profile-db)s).

Each of the resulting directories in your output folder will contain a stand-alone anvi'o project that can be used or shared without requiring access to any files of the original (larger) dataset.

## Splitting metagenomes and pangenomes

To split bins from a metagenome, you can provide the program %(anvi-split)s with a %(contigs-db)s and %(profile-db)s pair. To split gene clusters from a pangenome, you can provide it with a %(genomes-storage-db)s and %(pan-db)s pair. In both cases you will also need a %(collection)s. If you don't provide any %(bin)s names, the program will create individual directories for each bin that is found in your collection. You can also limit the output to a single bin.

### An example run

Assume you have a %(profile-db)s has a %(collection)s with three bins, which are (very creatively) called `BIN_1`, `BIN_2`, and `BIN_3`.

If you ran the following code:

{{ codestart }}
anvi-split -p %(profile-db)s \
           -c %(contigs-db)s \
           -C %(collection)s \
           -o OUTPUT
{{ codestop }}

You would get 3 new pairs of %(profile-db)s and %(contigs-db)s files, one for each bin, located in `OUTPUT/BIN_1/`, `OUTPUT/BIN_2/`, and `OUTPUT/BIN_3/`.

Alternatively, you can specify a bin name to limit the reported bins:

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

## Splitting a contigs database without a profile database

%(anvi-split)s can split a %(contigs-db)s on its own, without any %(profile-db)s. Each resulting directory will contain a self-contained %(contigs-db)s for that group of contigs. Two input modes are available. You will need either a %(collection-txt)s file mapping contigs to bins, or per-contig domain-level classification data previously imported with %(anvi-import-contig-classification)s.

### Using an external collection file

You can provide a two-column, TAB-delimited file with no header, where column 1 is the contig name and column 2 is the bin name:

{{ codestart }}
anvi-split -c %(contigs-db)s \
           --collection-txt %(collection-txt)s \
           -o OUTPUT
{{ codestop }}

### Using contig classification data

If your %(contigs-db)s has classification data imported with %(anvi-import-contig-classification)s, you can split it by contig class:

{{ codestart }}
anvi-split -c %(contigs-db)s \
           --split-by-contig-classification \
           -o OUTPUT
{{ codestop }}

Each class (e.g., `virus`, `plasmid`, `non-eukaryotic`) will become a separate output database. You can limit the output to specific classes with `--classes-to-keep`:

{{ codestart }}
anvi-split -c %(contigs-db)s \
           --split-by-contig-classification \
           --classes-to-keep virus,plasmid \
           -o OUTPUT
{{ codestop }}

#### Handling classification conflicts

If your %(contigs-db)s has %(contig-classification)s data from multiple sources, the same contig may be assigned different classes by different sources. %(anvi-split)s will raise an error when it encounters such conflicts. For example, the following classification table has data from two sources, `whokaryote` and `alien`. Both agree on `contig1` through `contig3`, but disagree on `contig4` through `contig6` — `whokaryote` assigns them class `1` (eukaryotic) while `alien` assigns them class `2` (virus):

| contig | class | source | tool_classification | confidence |
|--------|-------|--------|---------------------|------------|
| contig1 | 1 | whokaryote | eukaryote | NA |
| contig2 | 1 | whokaryote | eukaryote | NA |
| contig3 | 1 | whokaryote | eukaryote | NA |
| contig4 | 1 | whokaryote | eukaryote | NA |
| contig5 | 1 | whokaryote | eukaryote | NA |
| contig6 | 1 | whokaryote | eukaryote | NA |
| contig1 | 1 | alien | eukaryote | NA |
| contig2 | 1 | alien | eukaryote | NA |
| contig3 | 1 | alien | eukaryote | NA |
| contig4 | 2 | alien | virus | NA |
| contig5 | 2 | alien | virus | NA |
| contig6 | 2 | alien | virus | NA |

%(anvi-split)s will refuse to proceed until you decide how to handle them. You have three options:

* `--only-use-classification-source SOURCE`: only use classifications from one source, ignoring the other sources entirely.
* `--allow-multiple-classifications`: allow conflicting contigs to appear in all output splits they were assigned to.
* `--mark-conflicting-contigs-as-ambiguous`: redirect conflicting contigs into a separate `ambiguous` split and write a report file documenting their original classifications.
