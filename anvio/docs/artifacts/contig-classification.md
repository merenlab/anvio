Contig-level classification data stored in a %(contigs-db)s.

Each contig is assigned a standardized integer class representing its predicted domain of origin, the name of the tool that produced the classification, the raw classification string from that tool, and an optional confidence score. Multiple sources (tools) can coexist in the same database, keyed by the `source` column value. See %(contig-classification-txt)s to learn how to format classification data into an input file for %(anvi-import-contig-classification)s.

The standardized class vocabulary is:

| Class | Meaning |
|---|---|
| 0 | Non-eukaryotic |
| 1 | Eukaryotic |
| 2 | Virus |
| 3 | Plasmid |
| 4 | Organelle |
| 5 | Unclassified |

These standardized classes aren't enough for you? Let us know by [opening an issue in the anvi'o Github repository](https://github.com/merenlab/anvio/issues).
