Prepares user-provided HMM profiles or protein FASTA files as custom functional annotation
databases and organises them in a structured output directory that %(anvi-run-user-annotation)s
can consume directly.

## Overview

Functional annotation in anvi'o is normally driven by curated collections such as KEGG KOfams,
Pfams, or COGs. `anvi-setup-user-annotation-db` lets you bring your own database — either a set
of HMM profiles (HMMER3 format) or a collection of protein sequences (FASTA format) — and
prepares everything needed to run annotation without manual intermediate steps.

The program reads a TSV that maps a unique name to a file path, then for each entry:

* **HMM profiles** — parses every model, extracts profile-level cutoff annotations (TC, GA, NC),
  creates a standard anvi'o HMM source directory (genes.hmm.gz, genes.txt, kind.txt,
  noise_cutoff_terms.txt, reference.txt, target.txt), and selects the best cutoff strategy
  automatically.

* **Protein FASTA files** — calls `diamond makedb` to build a binary DIAMOND search database.

All prepared databases and their metadata are recorded in a `manifest.json` file inside the
output directory so that %(anvi-run-user-annotation)s can find them later.

The program **validates file content** against the file extension. If a file claims to be an HMM
profile by extension but contains FASTA sequences (or vice versa), setup will halt with an
actionable error message rather than silently producing wrong results downstream.

## Input TSV format

The input file is a plain tab-delimited text file with two or three columns:

| Column | Required | Description |
|--------|----------|-------------|
| `name` | yes | A unique identifier for this database. Used as the base of the annotation source name inside the contigs database (see note on source naming below). |
| `path` | yes | Absolute or relative path to the HMM profile file or protein FASTA file. |
| `companion_fasta` | no | Path to a protein FASTA file that contains the sequences used to build the HMM profiles. When provided, setup also builds a DIAMOND database from this FASTA so that %(anvi-run-user-annotation)s runs both an HMM search and a DIAMOND search in a single call. Useful for cross-validating HMM hits against their source sequences. Ignored for FASTA entries. |

An optional header line (`name<TAB>path`) and comment lines (starting with `#`) are allowed:

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

For each HMM source the program inspects every model and determines the most appropriate
HMMER noise cutoff:

| Condition | Strategy applied |
|-----------|-----------------|
| All models carry `TC` (trusted cutoffs) | `--cut_tc` |
| All models carry `GA` (gathering thresholds) | `--cut_ga` |
| All models carry `NC` (noise cutoffs) | `--cut_nc` |
| Mixed or no annotations | `-E 1e-5` |

If the fallback e-value is applied, the program emits a warning. You can override the choice
by editing `noise_cutoff_terms.txt` inside the HMM source directory before running
%(anvi-run-user-annotation)s.

## Basic usage

Set up all databases in a TSV:

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
                               --output-dir /path/to/annotation_dbs \
                               --num-threads 8
{{ codestop }}

## Managing the manifest

Once an annotation directory exists, you can inspect and modify it without re-running setup.

List all registered databases:

{{ codestart }}
anvi-setup-user-annotation-db --output-dir /path/to/annotation_dbs \
                               --list
{{ codestop }}

Remove a single database and delete its prepared files (other databases are untouched):

{{ codestart }}
anvi-setup-user-annotation-db --output-dir /path/to/annotation_dbs \
                               --remove MyOldHMMs
{{ codestop }}

`--list` and `--remove` do not require `--input-tsv`.

## Output directory structure

```
annotation_dbs/
├── manifest.json            # database registry (read by anvi-run-user-annotation)
├── hmm/
│   └── MyHMMs/              # standard anvi'o HMM source directory
│       ├── genes.hmm.gz
│       ├── genes.txt
│       ├── kind.txt
│       ├── noise_cutoff_terms.txt
│       ├── reference.txt
│       └── target.txt
└── diamond/
    ├── MyProteins.dmnd      # binary DIAMOND database
    └── MyHMMs.dmnd          # companion DIAMOND database (if companion_fasta was provided)
```

## Requirements

* `diamond` must be installed and in `PATH` when setting up FASTA-based databases or companion
  FASTA entries.
* `hmmscan` or `hmmsearch` must be available when later running %(anvi-run-user-annotation)s
  with HMM-based databases.
