Prepares user-provided HMM profiles or protein FASTA files as custom functional annotation
databases and organises them in a structured output directory that %(anvi-run-user-annotation)s
can consume directly.

## Overview

Functional annotation in anvi'o is normally driven by curated collections such as KEGG KOfams,
Pfams, or COGs. `anvi-setup-user-annotation-db` lets you bring your own database — either a set
of HMM profiles (HMMER3 format) or a collection of protein sequences (FASTA format) — and
prepares everything needed to run annotation without manual intermediate steps.

The program reads a TSV that maps a unique name to a file path, then for each entry:

* **HMM profiles** — parses every model, extracts per-model cutoff annotations (TC, GA, NC),
  creates a standard anvi'o HMM source directory (genes.hmm.gz, genes.txt, kind.txt,
  noise_cutoff_terms.txt, reference.txt, target.txt). Models are grouped by their best available
  cutoff so each group is searched with the optimal flag at run time. Missing `ACC` lines are
  detected at setup time: a single warning lists all affected models and the model name is used as
  the accession so downstream steps are never blocked. Pfam/TIGRFAM version suffixes (e.g.
  `PF00001.23`) are stripped from accessions automatically.

* **Protein FASTA files** — calls `diamond makedb` to build a binary DIAMOND search database.
  FASTA headers are parsed at setup time to build a per-sequence ID-normalization mapping stored
  alongside the `.dmnd` file. Common database formats are auto-detected and normalized without any
  user input (see *ID normalization* below).

All prepared databases and their metadata are recorded in a `manifest.json` file inside the
output directory so that %(anvi-run-user-annotation)s can find them later.

If one database fails during setup the program skips it, continues with the remaining entries,
and reports a summary of succeeded and failed databases at the end.

The program **validates file content** against the file extension. If a file claims to be an HMM
profile by extension but contains FASTA sequences (or vice versa), setup will halt with an
actionable error message rather than silently producing wrong results downstream.

## Input TSV format

The input file is a plain tab-delimited text file with two or three columns:

| Column | Required | Description |
|--------|----------|-------------|
| `name` | yes | A unique identifier for this database. Used as the base of the annotation source name inside the contigs database (see *Source naming* below). |
| `path` | yes | Absolute or relative path to the HMM profile file or protein FASTA file. |
| `companion_fasta` | no | Path to a protein FASTA file containing the sequences that were used to build the HMM profiles. When provided, setup also builds a DIAMOND database from this FASTA so that %(anvi-run-user-annotation)s runs both an HMM search and a DIAMOND search in a single call. Useful for cross-validating HMM hits against their source sequences. Ignored for FASTA-only entries. |

An optional header line and comment lines (starting with `#`) are allowed:

```
# My custom annotation databases
name	path	companion_fasta
MyHMMs	/data/my_models.hmm	/data/my_model_sequences.faa
MyProteins	/data/reference_proteins.faa
```

Accepted extensions:

* HMM profiles: `.hmm`, `.hmm.gz`
* Protein FASTA: `.faa`, `.fasta`, `.fa`, `.fas`, `.fna`

If no recognised extension is found the program sniffs the file content (HMMER3/f header for
HMM; `>` prefix for FASTA). If the extension and content disagree, an error is raised.

## Source naming in the contigs database

%(anvi-run-user-annotation)s stores all results in the `gene_functions` table. To make the
search method immediately visible, the annotation source name carries a suffix:

| File type | Annotation source name |
|-----------|----------------------|
| HMM profile | `{name}_HMM` (e.g. `MyHMMs_HMM`) |
| Protein FASTA | `{name}_DIAMOND` (e.g. `MyProteins_DIAMOND`) |

When a companion FASTA is provided for an HMM entry, both `{name}_HMM` and `{name}_DIAMOND`
appear in the gene_functions table after running %(anvi-run-user-annotation)s.

## Noise cutoff strategy for HMM profiles

For each model in an HMM source the program picks the best available cutoff annotation:

| Per-model annotation | Cutoff used |
|----------------------|-------------|
| `TC` (trusted cutoffs) present | `--cut_tc` |
| `GA` (gathering thresholds) present | `--cut_ga` |
| `NC` (noise cutoffs) present | `--cut_nc` |
| None of the above | `-E 1e-15` (anvio standard evalue) |

If all models share the same cutoff type the whole database is searched in one pass. If models
have **mixed cutoff types**, %(anvi-run-user-annotation)s splits the search into groups and runs
each group with its optimal cutoff flag, then merges all hits under a single source name. This
happens automatically — no manual intervention needed.

If all models lack cutoff annotations a warning is emitted and the evalue fallback is used.

## DIAMOND ID normalization

When setting up a protein FASTA database the program parses every FASTA header and auto-detects
the format to produce a clean, human-readable accession. No user configuration is required.

| Input sseqid format | Normalized accession | Example |
|---------------------|---------------------|---------|
| `sp\|P12345\|RECA_ECOLI` | gene name | `RECA` |
| `tr\|A0A000\|GENE_SPECIES` | gene name | `GENE` |
| `ref\|NP_123456.1\|` | accession without version | `NP_123456` |
| `WP_012345678.1` | accession without version | `WP_012345678` |
| `PF00001.23` | accession without version | `PF00001` |
| `lcl\|gene_name` | inner id | `gene_name` |
| `pdb\|4HHB\|A` | structure + chain | `4HHB_A` |
| plain `recA` | unchanged | `recA` |

The protein description (everything after the first space in the FASTA header) is stored in an
ID-mapping JSON file alongside the `.dmnd` database. At run time this description is used as the
`function` field in the gene_functions table, so results read as human-friendly text rather than
raw database identifiers.

## HMM accession normalization

HMM file extensions are stripped from accessions automatically. Pfam and TIGRFAM version
suffixes are preserved because they carry version information:

* `GT5.hmm` → `GT5`
* `GT5.hmm.gz` → `GT5`
* `PF00001.23` → `PF00001.23` (unchanged)
* `TIGR00001.1` → `TIGR00001.1` (unchanged)

## Basic usage

Set up all databases listed in a TSV (output goes to the default directory):

{{ codestart }}
anvi-setup-user-annotation-db --input-tsv my_databases.tsv
{{ codestop }}

Set up to a custom directory:

{{ codestart }}
anvi-setup-user-annotation-db --input-tsv my_databases.tsv \
                               --output-dir /path/to/annotation_dbs
{{ codestop }}

Reset an existing directory and rebuild from scratch:

{{ codestart }}
anvi-setup-user-annotation-db --input-tsv my_databases.tsv \
                               --output-dir /path/to/annotation_dbs \
                               --reset
{{ codestop }}

Use multiple threads for `diamond makedb`:

{{ codestart }}
anvi-setup-user-annotation-db --input-tsv my_databases.tsv \
                               --num-threads 8
{{ codestop }}

## Adding databases to an existing directory

If the output directory already has a `manifest.json`, running setup again with a new `--input-tsv`
only processes entries whose names are **not yet in the manifest**. Already-registered databases
are skipped with a warning:

```
WARNING
'MyHMMs' is already in the manifest — skipping.
Use `--remove MyHMMs` first if you want to replace it.
```

The summary at the end reports how many databases were newly set up, how many were skipped, and
how many failed, so the outcome is always explicit. Use `--reset` to wipe the directory and
rebuild everything from scratch.

## Managing the manifest

Once an annotation directory exists you can inspect and modify it without re-running setup.
Use `--list` and `--remove` as **standalone** flags — they cannot be combined with `--input-tsv`
or `--reset`.

List all registered databases:

{{ codestart }}
anvi-setup-user-annotation-db --list
{{ codestop }}

Remove a single database and delete its prepared files (other databases are untouched):

{{ codestart }}
anvi-setup-user-annotation-db --remove MyOldHMMs
{{ codestop }}

Both commands use the default output directory unless `--output-dir` points elsewhere.

## Output directory structure

```
annotation_dbs/
├── manifest.json                   # database registry
├── hmm/
│   └── MyHMMs/                    # standard anvi'o HMM source directory
│       ├── genes.hmm.gz           # HMM file (ACC lines injected for models that lacked them)
│       ├── genes.txt              # model name → normalized accession
│       ├── kind.txt
│       ├── noise_cutoff_terms.txt
│       ├── reference.txt
│       └── target.txt
└── diamond/
    ├── MyProteins.dmnd            # binary DIAMOND database
    ├── MyProteins_id_mapping.json # sseqid → {clean_id, description}
    └── MyHMMs.dmnd                # companion DIAMOND database (if companion_fasta provided)
```

## Requirements

* `diamond` must be installed and in `PATH` when setting up FASTA-based databases or companion
  FASTA entries.
* `hmmscan` or `hmmsearch` must be available when later running %(anvi-run-user-annotation)s
  with HMM-based databases.
