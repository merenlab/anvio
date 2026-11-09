A tab-delimited file describing the classification of contigs in a %(contigs-db)s. This file is the input to %(anvi-import-contig-classification)s.

You can create this file based on the output from whichever software you used to classify your contigs (or based on your own classification strategy). Different tools produce different output formats and often use different names for domain-level classes; hence, for consistency, it is up to you to wrangle those outputs so that they fit into the standardized set of classes that anvi'o knows about (as described in %(contig-classification)s).

The file must contain the following columns:

| Column | Type | Description |
|---|---|---|
| `contig` | text | Contig name (must match a contig in the contigs database) |
| `class` | integer | Standardized class: 0=non-eukaryotic, 1=eukaryotic, 2=virus, 3=plasmid, 4=organelle, 5=unclassified |
| `source` | text | Name of the tool that produced the classification (e.g. `genomad`, `tiara`, `custom`) |
| `tool_classification` | text | Raw classification string from the tool; use semicolons for multi-level entries (e.g. `Viruses;Duplodnaviria;Heunggongvirae`) |
| `confidence` | text | Confidence score from the tool, or `NA` if not available |

A few important notes:
- all contig names must match contigs in the target %(contigs-db)s
- a single file may contain classifications from multiple tools (multiple `source` values)
- only the `class` column is used downstream (i.e., by %(anvi-split)s), while the `tool_classification` column is just for you to be able to remember the original classification assignment coming from whichever tool(s) you used. It's therefore important for you to make sure the `class` column accurately reflects the original classification (_note that anvi'o cannot automatically sanity check this for you_, because the classifications from different software can be vastly different in format)

Here is an example file:

|**`contig`**|**`class`**|**`source`**|**`tool_classification`**|**`confidence`**|
|:--|:--|:--|:--|:--|
|Day17a_QCcontig1|0|genomad|Bacteria;Firmicutes|0.95|
|Day17a_QCcontig33|0|genomad|Bacteria;Firmicutes|0.92|
|Day17a_QCcontig4|0|genomad|Bacteria;Proteobacteria|0.88|
|Day17a_QCcontig74|0|genomad|Bacteria;Bacteroidetes|0.91|
|Day17a_QCcontig6|1|tiara|eukaryote|NA|
|Day17a_QCcontig86|0|tiara|bacteria|NA|
|Day17a_QCcontig90|0|tiara|bacteria|NA|
|Day17a_QCcontig205|2|tiara|virus|NA|

The standardized classes aren't enough for you? Let us know by [opening an issue in the anvi'o Github repository](https://github.com/merenlab/anvio/issues).
