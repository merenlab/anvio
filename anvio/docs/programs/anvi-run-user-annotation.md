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

### What is stored in gene_functions for HMM sources

| Column | Content |
|--------|---------|
| `source` | `{db_name}_HMM` |
| `accession` | HMM model accession (ACC field), or `"-"` if the model has none |
| `function` | `[HMM] {model_name}` — model NAME field with method tag |
| `e_value` | Full-sequence E-value from HMMER |

### What is stored in gene_functions for DIAMOND sources

| Column | Content |
|--------|---------|
| `source` | `{db_name}_DIAMOND` |
| `accession` | Target sequence ID (best blastp hit) |
| `function` | `[DMND] {target_id} [pident: {pct:.1f}%, aln_len: {len} aa, bitscore: {score:.1f}]` |
| `e_value` | BLASTP E-value |

The `[HMM]` and `[DMND]` method tags in the `function` field make the search method immediately
visible in any tabular export. The enriched DIAMOND `function` string encodes percent identity,
alignment length, and bit score so no search statistics are discarded by the fixed table schema.

## Basic usage

Run all databases in the annotation directory:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs
{{ codestop }}

Run only selected databases (comma-separated):

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs \
                          --database MyHMMs,MyProteins
{{ codestop }}

Use multiple threads and a custom HMMER program:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs \
                          --num-threads 8 \
                          --hmmer-program hmmsearch
{{ codestop }}

Override the DIAMOND e-value cutoff and require at least 50 % identity:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs \
                          --evalue 1e-10 \
                          --min-pident 50
{{ codestop }}

Run DIAMOND in ultra-sensitive mode:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs \
                          --diamond-sensitivity ultra-sensitive
{{ codestop }}

Re-run and overwrite existing annotations from a previous run:

{{ codestart }}
anvi-run-user-annotation --contigs-db CONTIGS.db \
                          --annotation-db-dir /path/to/annotation_dbs \
                          --just-do-it
{{ codestop }}

## DIAMOND options

| Flag | Default | Description |
|------|---------|-------------|
| `--evalue` | `1e-15` | E-value cutoff for blastp. Increase to find more hits, decrease to tighten. |
| `--min-pident` | none | Minimum percent identity (0–100). Hits below this threshold are discarded before being stored. |
| `--diamond-sensitivity` | DIAMOND default | Sensitivity mode: `fast`, `mid-sensitive`, `sensitive`, `very-sensitive`, `ultra-sensitive`. Higher sensitivity finds more distant homologues at the cost of runtime. |

## Checking what databases are available

Use %(anvi-setup-user-annotation-db)s to list or manage the manifest:

{{ codestart }}
anvi-setup-user-annotation-db --output-dir /path/to/annotation_dbs --list
{{ codestop }}

## Searching for results afterwards

Once annotation is complete, use standard anvi'o functions to explore results:

{{ codestart }}
anvi-search-functions --contigs-db CONTIGS.db --search "MyHMMs_HMM"
{{ codestop }}

{{ codestart }}
anvi-export-functions --contigs-db CONTIGS.db --annotation-sources MyHMMs_HMM,MyProteins_DIAMOND
{{ codestop }}

## Notes

* The `--annotation-db-dir` must point to the directory created by
  %(anvi-setup-user-annotation-db)s. Moving files inside that directory after setup will cause
  failures.
* For HMM-based databases the noise cutoff terms are fixed at setup time (TC, GA, NC, or -E 1e-5).
  The `--evalue`, `--min-pident`, and `--diamond-sensitivity` flags apply only to DIAMOND searches.
* DIAMOND searches use `--max-target-seqs 1` to keep only the single best hit per gene.
* DIAMOND blastp logs are written to `<annotation-db-dir>/logs/` and persist after the run so
  you can inspect them if a search fails or produces unexpected results.
* Use `--just-do-it` to overwrite existing annotations from a previous run of the same databases.
