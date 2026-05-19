A tab-delimited file describing the classification of contigs in a %(contigs-db)s. This file is the input to %(anvi-import-contig-classification)s.

The file must contain the following columns:

| Column | Type | Description |
|---|---|---|
| `contig` | text | Contig name (must match a contig in the contigs database) |
| `class` | integer | Standardized class: 0=non-eukaryotic, 1=eukaryotic, 2=virus, 3=plasmid, 4=organelle, 5=unclassified |
| `source` | text | Name of the tool that produced the classification (e.g. `genomad`, `tiara`) |
| `tool_classification` | text | Raw classification string from the tool; use semicolons for multi-level entries (e.g. `Viruses;Duplodnaviria;Heunggongvirae`) |
| `confidence` | text | Confidence score from the tool, or `NA` if not available |

All contig names must match contigs in the target %(contigs-db)s. A single file may contain classifications from multiple tools (multiple `source` values).
