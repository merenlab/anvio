A tab-delimited text file that describes one or more custom annotation databases to be set up by %(anvi-setup-user-annotation-db)s.

## Format

Each row defines one database. Required columns:

| Column | Description |
|--------|-------------|
| `name` | Unique identifier for the database (used as source name in anvio) |
| `path` | Absolute or relative path to the HMM file (`.hmm`, `.hmm.gz`) or protein FASTA (`.faa`, `.fasta`, `.fa`) |
| `type` | Either `hmm` or `diamond` |

Optional columns:

| Column | Description |
|--------|-------------|
| `companion_fasta` | Path to a protein FASTA containing the sequences used to build an HMM database. When provided, %(anvi-setup-user-annotation-db)s also runs DIAMOND search and produces a second source `<name>_DIAMOND` in the `gene_functions` table. |
| `db` | Scope override for `--cut-tc` per-database TC thresholds |

## Example

```
name	path	type
MyHMMs	/path/to/models.hmm	hmm
MyProteins	/path/to/proteins.faa	diamond
HydDB	/path/to/HydDB.hmm	hmm	/path/to/HydDB.faa
```

## Produced by

User (manually created)

## Consumed by

%(anvi-setup-user-annotation-db)s
