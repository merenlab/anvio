Runs functional annotation against a %(contigs-db)s using custom databases prepared by
%(anvi-setup-user-annotation-db)s and stores all results in the `gene_functions` table.

## Overview

`anvi-run-user-annotation` reads the `manifest.json` from a directory created by
%(anvi-setup-user-annotation-db)s and runs the appropriate search for each registered database.
All results end up in the `gene_functions` table so they are accessible through
`anvi-search-functions`, `anvi-export-functions`, and the interactive interface.

The annotation source name always encodes the search method used:

| Source type | Source name in gene_functions | Example |
|-------------|-------------------------------|---------|
| HMM profile | `{name}_HMM` | `MyHMMs_HMM` |
| Protein FASTA | `{name}_DIAMOND` | `MyProteins_DIAMOND` |

When an HMM database was set up with a companion FASTA, both `{name}_HMM` and `{name}_DIAMOND`
are populated in a single run.

## Re-annotation behaviour

By default the program **skips** any database whose annotation source already exists in the
contigs database and prints a warning. The remaining databases in the same run are still
processed normally.

To drop an existing annotation and re-annotate, pass `--force-overwrite`. Each database is
evaluated independently, so you can safely run a mix of new and already-annotated databases
without the flag — only the ones that need overwriting require it.

## What is stored in gene_functions

### HMM sources

| Column | Content |
|--------|---------|
| `source` | `{db_name}_HMM` |
| `accession` | HMM model accession (ACC field, version-stripped); model name used when ACC is absent |
| `function` | `[HMM] {model_name}` — model NAME field prefixed with method tag |
| `e_value` | Full-sequence E-value from HMMER |

### DIAMOND sources

| Column | Content |
|--------|---------|
| `source` | `{db_name}_DIAMOND` |
| `accession` | Auto-normalized ID (gene name, accession without version suffix, etc.) |
| `function` | `[DMND] {description} [pident: {pct:.1f}%, aln_len: {len} aa, bitscore: {score:.1f}, qcov: {qcov:.1f}%]` |
| `e_value` | BLASTP E-value |

The `description` in the DIAMOND `function` field is taken from the FASTA header line built at
setup time (everything after the first space). If no description was present in the header the
normalized ID is used instead. The `[HMM]` / `[DMND]` method tags make the search method
immediately visible in any tabular export. Query coverage (`qcov`) is always computed and shown
in the brackets when DIAMOND reports query length; `--qcov` is a **filter** threshold that
discards hits below the given percentage but does not affect whether qcov appears in the output.

## HMM cutoff groups

Per-model cutoff annotations (TC, GA, NC) embedded in HMM profiles are used **automatically** —
no flag is required. Each model is assigned to the best available cutoff group (TC > GA > NC >
evalue fallback `-E 1e-15`). When all models share the same type, a single search is run. When
types are mixed, the program splits models into homogeneous groups, runs each with its optimal
flag, then merges all hits under the single `{name}_HMM` source. This is fully transparent.

## HMM options

| Flag | Default | Description |
|------|---------|-------------|
| `--cut-tc` | none | Path to a TSV file for per-model custom trusted cutoff values (see below). |
| `--hmmer-program` | `hmmscan` | HMMER program to use (`hmmscan` or `hmmsearch`). |

### Custom trusted cutoffs (`--cut-tc`)

Override TC values for specific models without touching the rest of the database. The TSV has
two or three columns: `model_name`, `seq_tc`, and optionally `dom_tc` (defaults to `seq_tc`):

```
# optional header / comments ignored
ModelA	25.0	25.0
ModelB	30.0
```

Models listed here have their TC values injected into the extracted HMM file and are searched
with `--cut_tc`. All other models in the same database continue using their embedded cutoffs or
the evalue fallback. Models named in the file but not found in any loaded database are silently
ignored, so one TSV can be reused across runs even when only a subset of databases is annotated.

## Basic usage

Run all databases in the default annotation directory:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db
{{ codestop }}

Run all databases from a specific directory:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs
{{ codestop }}

Run only selected databases (comma-separated):

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --database MyHMMs,MyProteins
{{ codestop }}

Use multiple threads and a custom HMMER program:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --num-threads 8 \
                          --hmmer-program hmmsearch
{{ codestop }}

Drop existing annotations and re-annotate:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --force-overwrite
{{ codestop }}

## DIAMOND options

| Flag | Default | Description |
|------|---------|-------------|
| `--evalue` | `1e-15` (anvio standard) | E-value cutoff for blastp. |
| `--min-pident` | none | Minimum percent identity (0–100). Hits below this threshold are discarded. |
| `--qcov` | none | Minimum query coverage (0–100 %). Hits where the aligned region covers less than this fraction of the query sequence are discarded. |
| `--max-hsps` | DIAMOND default | Maximum number of HSPs per target per query. |
| `--diamond-sensitivity` | DIAMOND default | Sensitivity mode: `fast`, `mid-sensitive`, `sensitive`, `very-sensitive`, `ultra-sensitive`. |

All DIAMOND flags are independent and can be combined freely. For HMM-based databases the
threshold is fixed at setup time from the profile annotations and cannot be changed here.

## Checking what databases are available

{{ codestart }}
anvi-setup-user-annotation-db --list
{{ codestop }}

## Searching for results afterwards

{{ codestart }}
anvi-search-functions --contigs-db CONTIGS.db --search "MyHMMs_HMM"
{{ codestop }}

{{ codestart }}
anvi-export-functions --contigs-db CONTIGS.db --annotation-sources MyHMMs_HMM,MyProteins_DIAMOND
{{ codestop }}

## Notes

* `--annotation-db-dir` defaults to the anvio standard user-annotation data directory. Use
  `--output-dir` / `--annotation-db-dir` consistently between setup and run if you use a custom
  path.
* Moving or renaming files inside the annotation directory after setup will cause failures because
  the manifest stores absolute paths.
* DIAMOND blastp logs are written to `<annotation-db-dir>/logs/` for post-run inspection.
* DIAMOND searches use `--max-target-seqs 1` to keep only the single best hit per gene.
