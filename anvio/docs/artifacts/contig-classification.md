Contig-level classification data stored in a %(contigs-db)s. This artifact is produced by %(anvi-import-contig-classification)s.

Each contig is assigned a standardized integer class, the name of the tool that produced the classification, the raw classification string from that tool, and an optional confidence score. Multiple sources (tools) can coexist in the same database, keyed by the `source` column value.

The standardized class vocabulary is:

| Class | Meaning |
|---|---|
| 0 | Non-eukaryotic |
| 1 | Eukaryotic |
| 2 | Virus |
| 3 | Plasmid |
| 4 | Organelle |
| 5 | Unclassified |
