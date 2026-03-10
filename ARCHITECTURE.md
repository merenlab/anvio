# Architecture and Codebase Guide for Programmers

The purpose of this document is to serve as a comprehensive reference for anyone — or anything — working with the anvi'o codebase, covering its architecture, major subsystems, coding conventions, and developer workflow.

## What Is Anvi'o

Anvi'o is an interactive analysis and visualization platform for 'omics data (metagenomics, pangenomics, phylogenomics, metabolomics, tRNA-seq, etc.). It is a Python package that processes input files and stores results in SQLite databases, and exposes both a CLI (~180 programs) and a web-based visualization interface for various analyses.

Here is some information where things are currently:

- Version: `9-dev`, codename `eunice`
- Python requirement: **exactly 3.10.x** (enforced at import time in `anvio/__init__.py`)
- License: GPL 3.0
- Help: https://anvio.org
- Repository: https://github.com/merenlab/anvio

---

## Repository Layout

```
anvio/                    # The package
  __init__.py             # Global flags, arg-dict D, helpers; ~4100 lines
  version.py              # All version numbers (DB versions, anvio_version)
  errors.py               # Custom exception hierarchy
  terminal.py             # Run, Progress, Timer output classes
  constants.py            # Global constants (coverage fields, cigar ops, etc.)
  db.py                   # Low-level SQLite3 wrapper
  dbops.py                # High-level DB classes (ContigsDatabase, etc.)
  utils.py                # General-purpose utilities
  filesnpaths.py          # File/path validation helpers
  argparse.py             # Thin wrapper around ArgumentParser
  interactive.py          # Bottle web server backend
  bottleroutes.py         # Web API route handlers
  profiler.py             # BAM → profile database
  bamops.py               # BAM file operations
  variabilityops.py       # SNV / codon variability analysis
  summarizer.py           # Collection summary generation
  reactionnetwork.py      # Metabolic reaction networks (largest file, ~10k lines)
  trnaseq.py              # tRNA-seq specific logic (~8400 lines)
  keggmapping.py          # KEGG pathway mapping
  kgmlnetworkops.py       # KGML XML network parsing
  contigops.py            # Contig/split manipulation
  ccollections.py         # Bin collections
  genomestorage.py        # Genome storage DB for pangenomics
  auxiliarydataops.py     # HDF5/auxiliary coverage data
  homogeneityindex.py     # Gene cluster homogeneity
  fastalib.py             # Fast FASTA I/O (custom, not biopython)
  ttycolors.py            # Terminal color helpers
  dbinfo.py               # DB introspection (DBInfo class)
  cli/                    # One file per anvi-* program (~180 files)
  tables/                 # Database schema definitions
    __init__.py           # All table_name / table_structure / table_types triplets
    genecalls.py          # TablesForGeneCalls class
    miscdata.py           # TableForItemAdditionalData, LayerAdditionalData
    states.py             # TablesForStates
    kmers.py              # KMerTables
    ntpositions.py        # TableForNtPositions
    genelevelcoverages.py # TableForGeneLevelCoverages
    contigsplitinfo.py    # TableForContigsInfo, TableForSplitsInfo
  workflows/              # Snakemake workflow definitions
    __init__.py           # WorkflowSuperClass
    contigs/              # Contigs workflow
    metagenomics/
    pangenomics/
    phylogenomics/
    trnaseq/
  parsers/                # Parsers for external tool outputs
    base.py               # Parser base class + TaxonomyHelper
    hmmer.py, centrifuge.py, kaiju.py, interproscan.py, …
  drivers/                # Wrappers around external executables
  taxonomyops/            # Taxonomy estimation logic
  tests/                  # Integration tests (shell scripts)
    00.sh                 # Main test import
    run_component_tests_for_*.sh
  migrations/             # DB migration scripts (anvi-migrate)
  data/                   # HMM profiles, reference data, static assets
  scripts/                # Non-Python entry points (R scripts, etc.)
pyproject.toml            # Project metadata, entry points, dependencies
```

---

## Major Subsystems

### 1. Databases (SQLite3)

The core data model. Each analysis type has its own versioned SQLite database file. Databases are **never opened directly** — always go through the DB classes.

| DB type | Version | Class | Purpose |
|---|---|---|---|
| `contigs` | 24 | `ContigsDatabase` | Assembled sequences, gene calls, annotations, HMM hits, taxonomy |
| `profile` | 40 | `ProfileDatabase` | Per-sample coverage/variability from BAM files |
| `pan` | 21 | `PanDatabase` | Gene clusters from pangenomics |
| `genes` | 6 | `GenesDatabase` | Gene-level stats split out from profile |
| `structure` | 2 | `StructureDatabase` | Protein 3D structure predictions |
| `trnaseq` | 2 | `TRNASeqDatabase` | tRNA-seq profiling |
| `genomestorage` | 7 | (HDF5-based) | Multi-genome storage for pangenomics |

Version numbers live in `anvio/version.py`. Changing them without a corresponding migration script in `anvio/migrations/` will break existing databases.

### 2. CLI Programs (`anvio/cli/`)

~180 Python modules, one per `anvi-*` command. Each is registered in `pyproject.toml` under `[project.scripts]` mapping `"anvi-foo"` → `"anvio.cli.foo_bar:main"`. Scripts are installed as console entry points.

### 3. Visualization Server

`anvio/interactive.py` + `anvio/bottleroutes.py` implement a Bottle web server. `anvi-interactive` and related display commands start this server and open a browser. The front-end lives in `anvio/data/interactive/`.

### 4. Workflows (`anvio/workflows/`)

Snakemake-based pipelines. Each workflow inherits from `WorkflowSuperClass` in `anvio/workflows/__init__.py` and defines `self.rules`, `self.default_config`, etc. Run via `anvi-run-workflow`.

### 5. Parsers (`anvio/parsers/`)

Import results from external tools (HMMER, Centrifuge, Kaiju, InterProScan, etc.) into anvi'o databases. All inherit from `anvio/parsers/base.py:Parser`.

---

## Key Coding Patterns and Idioms

### The `A` Lambda (Universal arg extraction)

Every class that takes an `args` argparse.Namespace uses this pattern:

```python
A = lambda x: args.__dict__[x] if x in args.__dict__ else None
self.db_path = A('contigs_db')
self.num_threads = A('num_threads') or 1
```

This avoids `AttributeError` when args come from different sources (CLI vs. programmatic use).

### Progress and Run objects

`anvio/terminal.py` provides the two universal output classes. At module level in most files:

```python
run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize
```

**`Run`** prints structured key-value info:

```python
run.info('Database path', '/path/to/db')         # cyan key : yellow value
run.warning('Something odd happened')             # formatted warning block
run.info_single('An indented message', level=1)  # * indented message
run.warning(None, 'Section Header', lc='green')  # prints a header with no message body
```

**`Progress`** prints a live status line to stderr:

```python
progress.new('Processing splits', progress_total_items=1000)
for i, split in enumerate(splits):
    progress.update(f'split {i}')
    progress.increment()
progress.end()
```

Progress automatically goes silent when `--no-progress`, `--quiet`, or output is not a TTY.

### Error Raising

All errors are raised as custom exceptions from `anvio/errors.py`. **Never use bare `raise Exception(...)`**.

```python
from anvio.errors import ConfigError, FilesNPathsError

raise ConfigError(f"The database at '{db_path}' does not exist.")
raise FilesNPathsError(f"File not found: {path}")
```

All exception classes inherit from `AnvioError` and automatically format with colored output. With `--debug`, they also print a full traceback. The most common ones:

- `ConfigError` — bad config, wrong arguments, logical validation failures
- `FilesNPathsError` — missing files, bad paths
- `CommandError` — external command failure
- `TerminalError` — terminal/progress object misuse
- `GenesDBError`, `TRNAIdentifierError`, etc. — domain-specific

### Global Flags from sys.argv

Set at import time in `anvio/__init__.py` by inspecting `sys.argv` directly:

```python
anvio.DEBUG          # --debug
anvio.FORCE          # --force
anvio.QUIET          # --quiet
anvio.NO_PROGRESS    # --no-progress
anvio.AS_MARKDOWN    # --as-markdown
anvio.DISPLAY_DB_CALLS  # --display-db-calls
```

These are read throughout the codebase. Checking `if anvio.DEBUG:` is the standard way to add debug output.

### Lazy Loading via `LazyProperty` descriptor

`dbops.py` defines a `LazyProperty` thread-safe descriptor for expensive data loads:

```python
class ContigsSuperclass:
    @LazyProperty
    def splits_basic_info(self):
        return self.db.get_table_as_dict('splits_basic_info')
    # Loaded on first access, cached thereafter
```

### Module-level singleton `run` and `progress`

Most modules define module-level `run` and `progress` instances used for module-scope logging (not inside classes). Classes receive their own `run`/`progress` as constructor parameters and use those for instance-scope logging.

### `multiprocess` instead of `multiprocessing`

The codebase uses `import multiprocess as multiprocessing` (a fork using `dill` serializer) because Python 3.10 cannot pickle local lambdas with stdlib `multiprocessing`.

### Table Schemas

`anvio/tables/__init__.py` defines every database table as a triplet:

```python
genes_in_contigs_table_name      = 'genes_in_contigs'
genes_in_contigs_table_structure = ['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type', 'source', 'version']
genes_in_contigs_table_types     = ['numeric', 'text', 'numeric', 'numeric', 'text', 'numeric', 'numeric', 'text', 'text']
```

Always import as `import anvio.tables as t` and refer to `t.genes_in_contigs_table_name` etc.

### File-level Metadata

Every CLI module has these module-level variables:

```python
__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__   = "GPL 3.0"
__version__   = anvio.__version__
__authors__   = ['meren', 'ekiefl']          # GitHub handles
__requires__  = ['contigs-db', 'profile-db'] # anvi'o artifact types consumed
__provides__  = ['collection']               # anvi'o artifact types produced
__description__ = "One-sentence description of what this program does"
```

`__requires__` and `__provides__` are used by `anvi-help` and the programs network graph.

---

## Coding Style

Anvi'o does not follow PEP8 strictly. The rules below reflect what is actually in the codebase; new code should match these patterns, not PEP8 defaults.

### Whitespace and Layout

- **Two blank lines** between every top-level definition (classes and module-level functions), and also **between methods inside a class**.
- **Function `def` is always a single line**, no matter how long the signature. Do not break a long parameter list across multiple lines:

```python
# correct — one line even when long
def init_gene_level_coverage_stats_dicts(self, min_cov_for_detection=0, outliers_threshold=1.5, zeros_are_outliers=False, callback=None, callback_interval=100, init_split_coverage_values_per_nt=False, gene_caller_ids_of_interest=set([])):

# wrong — do not wrap the def line
def init_gene_level_coverage_stats_dicts(self,
                                         min_cov_for_detection=0,
                                         ...):
```

- Long **function call** arguments may be wrapped and aligned when it aids readability:

```python
self.genomes_storage = genomestorage.GenomeStorage(self.genomes_storage_path,
                                                   self.p_meta['genomes_storage_hash'],
                                                   genome_names_to_focus=self.p_meta['genome_names'],
                                                   run=self.run,
                                                   progress=self.progress)
```

- The `# pylint: disable=line-too-long` comment appears at the top of many files to silence line-length warnings. Long lines are acceptable.

### Naming Conventions

- `snake_case` for everything: variables, functions, methods, module names.
- Class names use `CamelCase` (`ContigsSuperclass`, `ProfileDatabase`, `DBClassFactory`).
- Constants are `UPPER_SNAKE_CASE` (`MAX_DEPTH`, `CODING`, `QUIET`).
- Single-letter uppercase lambdas are used as local shorthand:

```python
A = lambda x: args.__dict__[x] if x in args.__dict__ else None   # args extraction
D = lambda x: self.__dict__[x] if x in self.__dict__ else None   # self dict check
F = lambda x: '[YES]' if x else '[NO]'                            # bool formatter
G = lambda g: gene_caller_id_conversion_dict[g]                  # ad hoc lookup
```

These are defined at the start of a method, used locally, and never exported.

### Imports

Standard layout within a file, in this order with blank lines between groups:

```python
# 1. stdlib
import os
import sys
import copy

# 2. third-party
import numpy as np
import pandas as pd

# 3. anvio package
import anvio
import anvio.db as db
import anvio.tables as t
import anvio.utils as utils
import anvio.terminal as terminal
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
```

`import anvio.tables as t` is the universal convention — always `t`, never spelled out.

Avoid top-level imports that would create circular dependencies. Import inside a function when necessary, and add a comment explaining why:

```python
def get_metabolism_estimates_for_a_list_of_genes(self, ...):
    # we import here to avoid circular imports since dbops is such a popular hub
    from anvio.metabolism.estimate import KeggMetabolismEstimator
```

### Docstrings

Docstrings use a specific section-header style with `=` underlines. The standard format is:

```python
def get_items_additional_data_for_functions_per_split_summary(self, source, split_names_of_interest, data_dict={}, keys_list=[]):
    """Get items additional data layers to display the frequency of function names
       for each split in a given contigs database so it can be shown as a stacked bar
       chart in the anvi'o interactive interface.

    Parameters
    ==========
    source : str
        A functional annotation source that is in the contigs database.
    split_names_of_interest : list
        Split names to be considered.
    data_dict : dict
        An optional `items_additional_data_dict` type dictionary to update.
    keys_list : list
        An optional `items_additional_data_keys` type list to update.

    Returns
    =======
    data_dict : dict
        An `items_additional_data_dict` type dictionary.
    keys_list : list
        An `items_additional_data_keys` type list.
    """
```

Additional sections follow the same pattern:

```python
    Notes
    =====
    - If calling many times, consider `self.get_gene_info_for_each_position`
    - This function will not work unless gene functions are initialized first

    Examples
    ========
    >>> import tracemalloc
    >>> import anvio.utils as utils
    >>> tracemalloc.start()
    >>> snap = tracemalloc.take_snapshot()
    >>> utils.display_top_memory_usage(snap)
```

Short one-liner docstrings (for simple methods) are also fine:

```python
def get_all_genome_names_in_gene_clusters_dict(self, gene_clusters_dict):
    """Returns all genome names found in a `gene_clusters_dict`"""
```

Never add docstrings to code you didn't write. Only add or update them when making substantive changes to the function.

### Error Messages

Error messages are **conversational and specific**. They explain what went wrong and, where possible, what the user or developer can do about it. Self-deprecating humor is welcome since anvi'o users know that there are people behind the codebase and they are as frustrated with these issues as the users themselves:

```python
raise ConfigError(f"Oh no :( There is a contigs database here at '{self.db_path}', but it seems to be broken :/ "
                  f"It is very likely that the process that was trying to create this database failed, "
                  f"and left behind this unfinished thingy. Well, anvi'o believes it is best if you make "
                  f"it go away with fire, and try again.")
```

When an error comes from deep inside a class method, it is common to prefix the message with the function name so it is traceable without `--debug`:

```python
raise ConfigError("get_corresponding_codon_order_in_gene :: pos_in_contig must be of type 'int'")
```

Both `%` formatting and f-strings are used. Either is acceptable; match the style of the surrounding code.

### Progress and Error Interaction

**Always call `progress.end()` before raising an exception** when a progress bar is active. Failing to do so leaves a broken progress line on the terminal:

```python
self.progress.new('Processing splits')
self.progress.update('...')

for split_name in split_names:
    if split_name not in expected_splits:
        self.progress.end()   # ← must come before raise
        raise ConfigError("Split '%s' was not found." % split_name)
    ...

self.progress.end()
```

### Inline Comments

Comments above lines explain intent, not mechanics. They are written as sentences (capitalized, ending with a period or not — both are fine):

```python
# learn the number of categories for the function source
function_source_categories = set([])

# we can't use this strategy if there are many categories
if len(function_source_categories) > 10:
    raise ConfigError(...)
```

For longer functions, ALL-CAPS section labels mark major phases:

```python
# BUSINESS TIME
for split in splits:
    ...

# ACTION
self.progress.new(...)
```

`####...####` banners (as seen in `anvio/tables/__init__.py`) are used to separate logical sections in large files:

```python
####################################################################################################
#
#     TABLE DESCRIPTIONS FOR THE PROFILE DATABASE
#
####################################################################################################
```

### Default Mutable Arguments

The codebase uses `set([])` (not `set()`) and `[]` / `{}` as default argument values in many function signatures — this is a known Python anti-pattern but is pervasive in the existing code. Do not change existing signatures; follow the pattern when adding new parameters to existing functions. For new standalone functions, prefer `None` as default and initialize inside:

```python
# existing pattern (leave as is in old code)
def init_split_sequences(self, min_contig_length=0, split_names_of_interest=set([])):

# preferred in new code
def my_new_function(self, items=None):
    if items is None:
        items = set()
```

### Triggering Lazy Properties

To force a lazy-loaded property to load without using its value, assign to `_`:

```python
# this triggers the LazyProperty descriptor to load the data and cache it
_ = self.split_name_to_genes_in_splits_entry_ids
```

This pattern appears frequently before code that depends on a lazy property being populated as a side effect.

### String Formatting

The codebase has a long history of `%`-style formatting, and it remains common in older code. **Do not go back and change old `%`-style strings to f-strings** — that kind of sweep creates noisy diffs with no functional benefit. Simply leave existing code as-is and use f-strings in any new code you write:

```python
# older style — leave it alone when you see it
self.run.info('Pan DB', 'Initialized: %s (v. %s)' % (self.pan_db_path, version))

# f-strings for new code
self.run.info('Pan DB', f'Initialized: {self.pan_db_path} (v. {version})')
```

Never mix both styles in a single expression.

**Multi-line f-strings:** when a string is long enough to wrap across multiple lines, put the `f` prefix on *every* line — including lines that contain no interpolated variables. This keeps the lines visually uniform and makes it easy to add a variable to any segment later without also having to add the `f`:

```python
# correct — f on every line
raise ConfigError(f"Something went wrong while processing '{sample_name}'. "
                  f"This is likely due to a mismatch between the profile database "
                  f"and the contigs database. Please make sure both were generated "
                  f"from the same FASTA file and try again.")

# wrong — f only where there happens to be a variable
raise ConfigError(f"Something went wrong while processing '{sample_name}'. "
                  "This is likely due to a mismatch between the profile database "
                  "and the contigs database. Please make sure both were generated "
                  "from the same FASTA file and try again.")
```

Multi-line strings use implicit concatenation inside parentheses (not `\` continuation). This applies to both `%`-style and f-string multi-line literals.

---

## Important Shared Modules

| Module | Role |
|---|---|
| `anvio/terminal.py` | `Run`, `Progress`, `Timer`, `pluralize`, `pretty_print` — all user-facing output |
| `anvio/errors.py` | Exception hierarchy — all error raising goes through here |
| `anvio/db.py` | `DB` class — low-level SQLite3 with retry logic, numpy adapters, versioning |
| `anvio/dbops.py` | `ContigsDatabase`, `ProfileDatabase`, `PanDatabase`, `ContigsSuperclass`, `ProfileSuperclass`, `LazyProperty`, `DBClassFactory` |
| `anvio/tables/__init__.py` | All table schema definitions |
| `anvio/utils.py` | Hundreds of utilities — file I/O, format conversion, `run_command`, `rev_comp`, etc. |
| `anvio/filesnpaths.py` | `is_file_exists`, `is_output_dir_writable`, `get_temp_directory_path`, etc. |
| `anvio/constants.py` | `essential_data_fields_for_anvio_profiles`, `gene_call_types`, BAM fetch filters, `max_depth_for_coverage` |
| `anvio/version.py` | All version numbers (single source of truth) |
| `anvio/argparse.py` | Wrapper around `ArgumentParser` — use this, not stdlib directly |
| `anvio/__init__.py` | Global flags, argument dict `D` with all 200+ CLI argument definitions |

---

## Database Key Concepts

### Contigs Database

The foundation. Created by `anvi-gen-contigs-database` from a FASTA file.

Key tables:

- `contig_sequences` — raw DNA per contig
- `contigs_basic_info` — length, gc_content, num_splits per contig
- `splits_basic_info` — contigs are split into ~20kb chunks ("splits") for display; this table stores split coordinates
- `genes_in_contigs` — gene calls (gene_callers_id, start, stop, direction, call_type)
- `gene_amino_acid_sequences` — translated sequences
- `gene_functions` — functional annotations (COGs, KOfams, Pfams, etc.)
- `hmm_hits` — HMMER search results
- `scg_taxonomy` / `trna_taxonomy` — taxonomic assignments
- `self` — key-value metadata (db_type, version, creation_date, etc.)

### Profile Database

Created by `anvi-profile` (one per BAM file) and merged by `anvi-merge`. Contains coverage and variability per split per sample.

Key tables:

- `mean_coverage`, `detection`, `variability`, etc. — view tables (one per coverage metric per sample)
- `variable_nucleotides` — SNVs
- `variable_codons` — codon-level variability
- `indels` — insertion/deletion events
- `collections_of_splits`, `collections_bins_info` — bin/collection data
- `item_additional_data`, `layer_additional_data` — miscellaneous data layers

The `self` table stores `merged`, `samples`, `contigs_db_hash` (for linkage verification).

### Splits vs. Contigs

A fundamental design decision: long contigs are split into ~20 kb chunks called "splits" for the visualization interface. Most display and analysis operates on splits, not contigs. Gene-to-split mapping is stored in `genes_in_splits`.

### Collections and Bins

Users group splits into "bins" and name groups of bins "collections". Stored in `collections_of_splits` and `collections_bins_info` tables (present in both contigs and profile databases).

---

## Non-Obvious Design Decisions

1. **Splits are the unit of display, not contigs.** Even short contigs get one split. The split_length default (~20000 bp) is configurable at contigs DB creation time and cannot be changed after.

2. **The `self` table is a key-value store.** Every database has a `self` table. Access it with `db.get_meta_value('key')` and `db.set_meta_value('key', value)`. It holds DB type, version, creation date, and analysis-specific metadata.

3. **Database version bumps require migration scripts.** Version numbers are in `anvio/version.py`. Any schema change that increments a version needs a migration script in `anvio/migrations/` and must be wired into `anvi-migrate`. Do not increment versions casually.

4. **`ROWID` prepending.** For tables where the first column is not unique, `db.py` automatically prepends `ROWID as "entry_id"` when reading. This is controlled by `tables.is_table_requires_unique_entry_id(table_name)`. This means some tables get an implicit `entry_id` column on read that doesn't exist in the schema definition.

5. **The argument dictionary `D` in `__init__.py`.** All ~200+ CLI arguments across all programs are defined once in `anvio/__init__.py` in the dict `D`. Each entry has `metavar`, `help`, and sometimes `required`, `type`, etc. CLI modules call `A('--arg-name')` to pull from this shared definition. This ensures consistent argument naming and help text across all programs.

6. **`fastalib.py` is a custom FASTA parser.** Anvi'o does not use BioPython or other libraries for FASTA I/O in hot paths. `anvio/fastalib.py` is a purpose-built, performance-oriented FASTA reader/writer.

7. **numpy types need SQLite adapters.** `db.py` registers `sqlite3.register_adapter(numpy.int64, int)` etc. at module import time. If you add new numpy types being stored, add adapters here.

8. **`multiprocess` not `multiprocessing`.** See "Coding Patterns" above. Always `import multiprocess as multiprocessing` if you need multiprocessing.

9. **Auxiliary data is separate.** Per-nucleotide coverage arrays are stored in a separate auxiliary `.db` file (HDF5-based, handled by `auxiliarydataops.py`), not in the main profile database. This avoids bloating the SQLite file.

---

## Things That Should Not Be Changed Without Extreme Care

- **Database version numbers** (`anvio/version.py`) — changing without a migration script breaks all existing databases
- **Table schemas** (`anvio/tables/__init__.py`) — breaking changes need version bumps + migrations
- **The `self` table key conventions** — many programs read specific keys; changing names silently breaks them
- **`essential_data_fields_for_anvio_profiles`** in `constants.py` — these are the view tables created by profiling; removing or renaming breaks profile databases
- **The `gene_callers_id` numbering** — these IDs appear in multiple tables and are the primary key linking genes across all analyses
- **The split naming convention** (`contig_name_split_00001`) — splits are referenced by name throughout databases and display states
- **The `--debug` / `--quiet` / `--no-progress` global flag mechanism** — checked at import time via `sys.argv`; cannot be set after import

---

## Developer Workflow

### Installation (in development mode)

See https://anvio.org/install/

The development mode installation will allow editing the code and immediately testing it without re-installing anything.

### Running Tests

Tests are shell scripts that run actual anvi'o commands end-to-end for component testing. They are under,

```bash
cd anvio/tests
bash run_component_tests_for_metagenomics.sh
bash run_component_tests_for_pangenomics.sh
```

There are no pytest unit tests in the traditional sense — testing is integration-based via shell scripts.

### Documentation (`anvio/docs/`)

Anvi'o maintains hand-written documentation for every program and artifact (data type) that is rendered and published at **https://anvio.org/help/**. The docs live inside the package itself under three subdirectories:

```
anvio/docs/
  programs/    # one .md file per anvi-* program
  artifacts/   # one .md file per anvi'o artifact/data type
  workflows/   # one .md file per anvi'o Snakemake workflow
```

**Programs** correspond 1-to-1 with `anvi-*` CLI commands. Each file describes what the program does, how to use it, its key flags, and worked examples. **Artifacts** are anvi'o's named data types (e.g. `contigs-db`, `profile-db`, `collection`, `fasta`) — the same names used in `__requires__` and `__provides__` in every CLI module. Artifact docs explain what that data type is and how it is created and consumed.

Docs use two special conventions:

- `%(artifact-name)s` — cross-links to another artifact or program's help page (e.g. `%(contigs-db)s`, `%(anvi-run-hmms)s`)
- `{{ codestart }}` / `{{ codestop }}` — wraps shell command examples (rendered as styled code blocks on the website; use these instead of fenced code blocks for commands the user is meant to run)

**The rule: docs must stay in sync with the code.** Specifically:

- **Adding a new program** → create `anvio/docs/programs/anvi-my-program.md`. The file should describe what the program does, its inputs and outputs, and include at least one usage example. Run `anvi-script-gen-help-pages` to see if the new markdown document compiles well.
- **Changing a program's behavior** — new flags, changed inputs/outputs, altered logic that affects what users see — → update the corresponding `anvio/docs/programs/` file to match.
- **Adding or renaming an artifact** → create or update the corresponding `anvio/docs/artifacts/` file.

Stale docs that no longer match the code are actively harmful because they are what users read first. Treat docs as part of the same commit as the code change, not an afterthought.

### Adding a New CLI Program

1. Create `anvio/cli/my_program.py` with the standard module-level metadata (`__requires__`, `__provides__`, etc.) and a `main(args)` function.
2. Add the entry point to `pyproject.toml` under `[project.scripts]`.
3. Re-run `pip install -e .` to register the entry point.
4. Create `anvio/docs/programs/anvi-my-program.md` (see Documentation section above).

### Debugging

- Pass `--debug` to any program to get full tracebacks from `AnvioError` exceptions.
- Pass `--display-db-calls` to print all SQL statements.
- Use `anvio.P(some_dict)` (defined in `__init__.py`) for pretty-printed JSON debugging.
- Use `anvio.SHOW_CALLER()` to print the call stack.

---

## Dependencies Worth Knowing

These may change over time:

- `numpy==1.24.1` — pinned; many array operations throughout
- `pandas==1.4.4` — pinned; used for tabular data
- `scikit-learn==1.2.2` — pinned; clustering, ordination
- `matplotlib==3.5.1` — pinned; static figure generation
- `bottle` — lightweight web framework for the interactive interface
- `snakemake` — workflow engine
- `pysam` — BAM file reading
- `pyrodigal_gv` — gene calling (default caller, replaces prodigal)
- `multiprocess` — fork of multiprocessing using dill
- `colored` — terminal color codes (`Fore`, `Back`, `Style` in terminal.py)
- `networkx==3.1` — graph operations (reaction networks, programs network)
- `ete3` — phylogenetic tree handling

Python requirement is strict: `==3.10.*` (checked at import, enforced in `pyproject.toml`).