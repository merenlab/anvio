This program **stores a metabolic %(reaction-network)s in a %(contigs-db)s or %(pan-db)s.**

The network consists of data on biochemical reactions predicted to be encoded by the genome or pangenome.

Information on the predicted reactions and the involved metabolites are stored in tables of the %(contigs-db)s or %(pan-db)s. The program, %(anvi-get-metabolic-model-file)s, can be used to export the %(reaction-network)s from the database to a %(reaction-network-json)s file formatted for input into programs for flux balance analysis.

## Setup

%(anvi-setup-kegg-data)s downloads [binary relations files](https://www.genome.jp/brite/br08906) needed to construct a %(reaction-network)s from [KEGG Orthology (KO)](https://www.genome.jp/kegg/ko.html) sequence annotations. Make sure to run that program with the `--kegg-snapshot` option to use the newest snapshot of %(kegg-data)s, [`v2024-08-30`](https://figshare.com/articles/dataset/KEGG_build_2024-08-30/26880559?file=48903154), which includes binary relations files.

{{ codestart }}
anvi-setup-kegg-data --kegg-snapshot v2024-08-30
{{ codestop }}

%(anvi-setup-modelseed-database)s sets up the [ModelSEED Biochemistry database](https://github.com/ModelSEED/ModelSEEDDatabase), which harmonizes biochemical data from various reference databases, including KEGG. The following command sets up the database in a default anvi'o directory.

{{ codestart }}
anvi-setup-modelseed-database
{{ codestop }}

### Download newest available KEGG files

Alternatively, KEGG data including binary relations files can be set up not from a snapshot but by downloading the newest files available from KEGG using the `-D` flag. In the following command, a higher number of download threads than the default of 1 is provided by `-T`, which significantly speeds up downloading.

{{ codestart }}
anvi-setup-kegg-data -D -T 5
{{ codestop }}

### Install in non-default location

At the moment, KEGG data that includes binary relations files does _not_ include "stray" KOs (see %(anvi-setup-kegg-data)s) due to changes in the available model files. To preserve KEGG data that you already have set up, for this reason or another, the new snapshot or download can be placed in a non-default location using the option, `--kegg-data-dir`.

{{ codestart }}
anvi-setup-kegg-data --kegg-snapshot v2024-08-30 --kegg-data-dir path/to/other/directory
{{ codestop }}

`anvi-reaction-network` requires a `--kegg-dir` argument to seek KEGG data in a non-default location.

Likewise, different versions of the ModelSEED Biochemistry database can be set up in non-default locations and used with the `--modelseed-dir` argument.

## Usage

%(anvi-reaction-network)s takes a either a %(contigs-db)s OR a %(pan-db)s and %(genomes-storage-db)s as required input. Genes stored within the %(contigs-db)s or %(genomes-storage-db)s must have KO protein annotations, which can be assigned by %(anvi-run-kegg-kofams)s.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db
{{ codestop }}

If a %(contigs-db)s already contains a %(reaction-network)s from a previous run of this program, the flag `--overwrite-existing-network` can overwrite the existing network with a new one. For example, if %(anvi-run-kegg-kofams)s is run again on a database using a newer version of KEGG, then %(anvi-reaction-network)s should be rerun to update the %(reaction-network)s derived from the KO annotations.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db \
                      --overwrite-existing-network
{{ codestop }}

A %(reaction-network)s can also be generated from consensus KO annotations of gene clusters. This can be used to understand the conservation or divergence of parts of the metabolic network between organisms in the pangenome.

{{ codestart }}
anvi-reaction-network -p /path/to/pan-db \
                      -g /path/to/genomes-storage-db
{{ codestop }}
