This program **stores a metabolic %(reaction-network)s in a %(contigs-db)s or %(pan-db)s.**

The network consists of data on biochemical reactions predicted to be encoded by the genome or pangenome, referencing the [KEGG Orthology (KO)](https://www.genome.jp/kegg/ko.html) and [ModelSEED Biochemistry](https://github.com/ModelSEED/ModelSEEDDatabase) databases.

Information on the predicted reactions and the involved metabolites are stored in two tables of the %(contigs-db)s or %(pan-db)s. The program, %(anvi-get-metabolic-model-file)s, can be used to export the %(reaction-network)s from the database to a %(reaction-network-json)s file formatted for flux balance analysis.

## Usage

%(anvi-reaction-network)s takes a either a %(contigs-db)s OR a %(pan-db)s and %(genomes-storage-db)s as required input. Genes stored within the %(contigs-db)s or %(genomes-storage-db)s must have KO protein annotations, which can be assigned by %(anvi-run-kegg-kofams)s.

The KO and ModelSEED Biochemistry databases must be set up and available to the program. By default, these are expected to be set up in default anvi'o data directories. %(anvi-setup-kegg-data)s and %(anvi-setup-modelseed-database)s must be run to set up these databases.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db
{{ codestop }}

Custom locations for the reference databases can be provided with the flags, `--ko-dir` and `--modelseed-dir`.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db \
                      --ko-dir /path/to/set-up/ko-dir \
                      --modelseed-dir /path/to/set-up/modelseed-dir
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
