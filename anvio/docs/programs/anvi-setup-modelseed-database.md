This program **downloads and sets up the latest version of the ModelSEED Biochemistry database.**

[The ModelSEED Biochemistry database](https://github.com/ModelSEED/ModelSEEDDatabase) consists of two tab-delimited files of reaction and compound data, respectively, and is valuable due to harmonization of IDs and properties from multiple reference databases commonly used in metabolic modeling.

%(anvi-reaction-network)s relies upon ModelSEED Biochemistry in conjunction with the KEGG Orthology database. [KEGG Orthology (KO)](https://www.genome.jp/kegg/ko.html) protein annotations of genes are associated with predicted enzymatic reactions. These KEGG reactions are cross-referenced to the ModelSEED Biochemistry database to retrieve information on properties including reaction stoichiometry and reversibility. %(anvi-reaction-network)s stores reactions and metabolites thereby predicted in the %(contigs-db)s for the genome. The program, %(anvi-setup-kegg-data)s, sets up the requisite KO database.

## Usage

The simplest %(anvi-setup-modelseed-database)s command sets up the database in the default anvi'o ModelSEED data directory.

{{ codestart }}
anvi-setup-modelseed-database
{{ codestop }}

A custom directory can be provided instead. Within the provided directory, a subdirectory named `ModelSEED` is created for storage of the database.

{{ codestart }}
anvi-setup-modelseed-database --dir /path/to/dir
{{ codestop }}

Finally, in conjunction with either of the previous commands, the `--reset` flag can be used to delete any existing target database directory and its contents before setting up the latest version of the ModelSEED Biochemistry database there.

{{ codestart }}
anvi-setup-modelseed-database --reset
{{ codestop }}
