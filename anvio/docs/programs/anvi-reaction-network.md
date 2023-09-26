This program **stores a metabolic %(reaction-network)s in a %(contigs-db)s.**

The network consists of data on biochemical reactions predicted to be encoded by the genome, referencing the [KEGG Orthology (KO)](https://www.genome.jp/kegg/ko.html) and [ModelSEED Biochemistry](https://github.com/ModelSEED/ModelSEEDDatabase) databases.

Information on the predicted reactions and the involved metabolites are stored in two tables of the %(contigs-db)s. The program, %(anvi-get-metabolic-model-file)s, can be used to export the %(reaction-network)s from the database to a JSON-formatted file suitable for inspection and flux balance analysis.

## Usage

%(anvi-reaction-network)s takes a %(contigs-db)s as required input. Genes stored within the database must have KO protein annotations, which can be assigned by %(anvi-run-kegg-kofams)s.

The KO and ModelSEED Biochemistry databases must be set up and available to the program. By default, these are expected to be set up in default anvi'o data directories. %(anvi-setup-kegg-data)s and %(anvi-setup-modelseed-database)s must be run to set up these databases.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db
{{ codestop }}

Custom locations for the reference databases can be provided with the flags, `--ko-dir` and `--modelseed-dir`.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db --ko-dir /path/to/set-up/ko-dir --modelseed-dir /path/to/set-up/modelseed-dir
{{ codestop }}

If a %(contigs-db)s already contains a %(reaction-network)s from a previous run of this program, the flag `--overwrite-existing-network` can overwrite the existing network with a new one. For example, if %(anvi-run-kegg-kofams)s is run again on a database using a newer version of KEGG, then %(anvi-reaction-network)s should be rerun to update the %(reaction-network)s derived from the KO annotations.

{{ codestart }}
anvi-reaction-network -c /path/to/contigs-db --overwrite-existing-network
{{ codestop }}
