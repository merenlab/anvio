# pyANI-plus Python 3.13 prototype validation

This note records compatibility evidence for the `pyani-plus` backend,
which is now a required dependency in the later ANI feature branch. The Python 3.13 port supports ANIb and ANIm; the evidence below does
not cover every workflow or production-scale input.

## Current installation recipe

From the ANI feature checkout, use `python -m pip install -e .`.
The required dependencies include `pyani-plus>=1.0.1,<2` and `click>=8,<9`;
no optional extra is needed. `anvi-run-workflow` explicitly registers the
bundled logger for ordinary editable installs. Standalone
`snakemake --logger anvio` discovery remains unsupported in that mode;
wheel installs retain automatic discovery. ANIb needs BLAST, and ANIm
needs MUMmer on PATH. Backend execution currently supports Linux only.

## Historical installation evidence

The following commands and versions describe earlier validation snapshots,
including the original optional-backend prototype. They are retained as
evidence, not as the current required-dependency recipe. The tested
combined recipe used a fresh CPython 3.13.13 venv and resolved both applications
in one pip invocation:

```bash
python3.13 -m venv /tmp/anvio-pyani-joint-313
/tmp/anvio-pyani-joint-313/bin/python -m pip install \
  -e /tmp/anvio-pyani-plus-prototype 'pyani-plus==1.0.1' \
  'typer==0.27.2' 'click==8.3.1'
/tmp/anvio-pyani-joint-313/bin/python -m pip check
```

Standalone Click is explicit because pyani-plus imports it directly, while
current Typer no longer brings that package into this installation. An initial
unconstrained pyani-plus installation failed with `ModuleNotFoundError: No
module named 'click'`. Typer 0.15.4 with Click 8.1.8 was an initial workaround;
subsequent actual runs established that Typer 0.27.2 with standalone Click 8.3.1
works, including the private worker CLI. No Typer upper bound is required by
the tested combination. Exact versions above describe the reproducible test,
not permanent constraints in anvi'o dependency metadata.

The original runtime required `networkx==3.1`, conflicting with pyani-plus
1.0.1's `networkx>=3.4.2`. Runtime commit `39307922b` changes this to
`networkx>=3.4.2,<4` after graph checks at 3.4.2 and 3.7. The feature is rebased
onto that runtime commit. The original isolated `/tmp/pyani-plus-313` runs below
used Typer 0.15.4/Click 8.1.8 and remain historical evidence; they no longer imply
that isolation or the old Typer workaround is required. An isolated environment
can instead install `pyani-plus==1.0.1`, `typer==0.27.2`, and `click==8.3.1` and
be selected through `--pyani-plus-program`.

ANIm also needs MUMmer on `PATH`. Validation installed it in a separate conda
prefix, without changing the anvi'o or pyani-plus Python environments:

```bash
/home/bcoltman/miniforge3/bin/conda create -p /tmp/pyani-mummer -y \
  --override-channels -c conda-forge -c bioconda 'mummer=3.23'
export PATH="/tmp/pyani-mummer/bin:$PATH"
```

The installed conda package was MUMmer 3.23 (`mummer 3.23
pl5321h503566f_21`); `nucmer --version` reports the NUCmer engine version as
3.1.

## Synthetic fixtures

Both smoke runs used these deterministic FASTA inputs, generated with seed
713:

```bash
mkdir -p /tmp/pyani-real/inputs
python3.13 - <<'PY'
import random
from pathlib import Path

r = random.Random(713)
base = ''.join(r.choice('ACGT') for _ in range(6120))
mutant = list(base)
next_base = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
for i in range(0, len(mutant), 50):
    mutant[i] = next_base[mutant[i]]
sequences = {
    'base': base,
    'duplicate': base,
    'mutant': ''.join(mutant),
    'short': base[:3060],
    'unrelated': ''.join(r.choice('ACGT') for _ in range(6120)),
}
out = Path('/tmp/pyani-real/inputs')
for name, sequence in sequences.items():
    (out / f'{name}.fa').write_text(f'>{name}\n{sequence}\n')
PY
```

## ANIb smoke evidence

The ANIb run used BLASTN+ 2.17.0 from the existing anvi'o environment at
`/home/bcoltman/miniforge3/envs/anvio-dev-313/bin/blastn`. Direct pyANI-plus runs
prepended `/tmp/pyani-plus-313/bin` for Snakemake/helper discovery and included
the anvi'o environment bin directory for BLAST. For an independent environment,
install BLAST+ separately (for example `conda install -p /tmp/pyani-mummer
-c conda-forge -c bioconda blast=2.17.0`) and include its bin directory on `PATH`.
Five local FASTA files were used:

| Input | Length | Construction |
| --- | ---: | --- |
| `base.fa` | 6,120 nt | Base synthetic sequence |
| `duplicate.fa` | 6,120 nt | Identical to base |
| `mutant.fa` | 6,120 nt | Base with approximately 2% substitutions |
| `short.fa` | 3,060 nt | Half-length sequence |
| `unrelated.fa` | 6,120 nt | Unrelated sequence |

The initial command created the database, but failed inside the sandbox.
The subsequent resume completed run 1 and marked it `Done`:

```bash
/tmp/pyani-plus-313/bin/pyani-plus anib /tmp/pyani-real/inputs \
  --database /tmp/pyani-real/anib.sqlite --create-db
# Resume outside the socket-restricted sandbox after the first failure:
/tmp/pyani-plus-313/bin/pyani-plus resume --database /tmp/pyani-real/anib.sqlite
```

The output contained 25 directional comparisons for five genomes. Spot checks
in the database matched the fixture construction: base versus its duplicate
reported 100% ANI; base versus mutant reported approximately 98.022% ANI with
6,118 aligned nucleotides; base versus the half-length exact sequence reported
100% identity with 50% query coverage. Mutant query versus the half-length
sequence reported approximately 98.006% ANI with 3,059 aligned nucleotides and
approximately 50% query coverage. The unrelated genome had no alignment
result. The successful resume exited without an error and logged completion.

The first invocation in this restricted WSL environment failed when the
Snakemake progress helper attempted to start Python's multiprocessing
`forkserver`: binding its Unix socket raised `PermissionError: [Errno 1]
Operation not permitted`. The database was resumable; invoking `pyani-plus
resume --database /tmp/pyani-real/anib.sqlite` completed the run in 2.47 seconds.
This initial failure reflects the environment's socket restrictions and is
recorded so the successful resume is not mistaken for a clean first attempt.

## ANIm smoke evidence

The ANIm run used the same five fixtures generated above, pyani-plus 1.0.1,
Python 3.13.13, Snakemake 9.27.0, and the MUMmer installation above:

```bash
export PATH="/tmp/pyani-mummer/bin:$PATH"
/tmp/pyani-plus-313/bin/pyani-plus anim /tmp/pyani-real/inputs \
  --database /tmp/pyani-real/anim.sqlite --create-db
```

The run completed and logged run 1 as `Done`. Exported matrices in
`/tmp/pyani-real/anim-export` retained all five genome labels, including the
duplicate sequence and unrelated genome. They reported base versus mutant ANI
of approximately 0.9800621, with 6,119 aligned nucleotides and 122 similarity
errors. Base versus the half-length exact sequence had identity 1.0 and query
coverage 0.5; in the reverse direction coverage was 1.0. Comparisons involving
the unrelated genome and any other genome had blank/NA identity, coverage,
alignment-length, and error fields. In the mutant-query to short-subject row,
ANIm reported approximately 0.9800588 identity, 3,059 aligned nucleotides,
61 errors, and approximately 50% query coverage.

For this ANIm result, pyani-plus computes weighted identity as
`sum((reference_aligned_length + query_aligned_length) - 2 * errors) /
sum(reference_aligned_length + query_aligned_length)`. This weights both sides
of each alignment and differs from a reference-only indel weighting. The
observed scores therefore do not establish strict numeric equivalence with
legacy PyANI. The original prototype was optional. The later ANI feature now requires
pyANI-plus; these bounded tests still do not establish broader compatibility.

## Validation boundary and adoption notes

The combined environment and small fixtures establish installation, CLI
startup, `pip check`, and ANIb/ANIm execution on Python 3.13.13. A five-genome
repository fixture with 250,000 bases per genome also completed both methods;
it is a reduced test fixture, not a production genome set. A complete ANIb run
updated a copied Pan Database with all 25 pairwise identity values and layer
orders. Dereplication imported that complete matrix and placed all five genomes
in one cluster at the configured 0.99 threshold.

Legacy PyANI 0.2.12 under Python 3.10.15 was run on the same fixture and on the
small substitution/half-length fixtures. ANIb identity, query coverage and
alignment lengths matched within serialization precision. ANIm had small
differences on the 6,120-base fixture and on the repository fixture (up to
529 bp in an alignment length); the implementations must not be treated as
numerically interchangeable. Identity values are fractions from 0 to 1,
alignment coverage is directional query coverage, and alignment lengths are
base pairs. For example, repository-fixture ANIb coverage was 0.959376 from
E_faecalis_6240 to E_faecalis_6563 and 0.954780 in the reverse direction.

pyANI-plus 1.0.1 exposes no public thread limit and its local workflow invokes
Snakemake with `--cores all`. On Linux, anvi'o launches the pyANI-plus commands
through `taskset` with a selected subset of the current process's allowed CPU
IDs. `--num-threads` defaults to 1, keeps the global `ANVIO_THREADS` default
behavior, and is capped to the number of allowed CPUs when it is larger. The
workflow can use all CPUs visible within that allocation.
Anvi'o logs the requested count, allowed CPU count, and actual selected IDs.
This controls CPU affinity inherited by child processes, not every OS thread.
The pyANI-plus backend currently supports Linux only and is unsupported on
macOS and Windows. The allocation requires Linux, `os.sched_getaffinity`, and
`taskset`; it fails before starting pyANI-plus when any requirement is
unavailable. Legacy pyANI
and fastANI retain their existing thread behavior.

## Feature CLI checks

The combined CPython 3.13.13 environment ran both ANIb and ANIm through
`anvi-compute-genome-similarity --program pyANI` on five packaged E. faecalis
test genomes. Each FASTA exported from its ContigsDB contains 250,000 bases
(1.25 Mb total); this is a reduced test fixture rather than a full production
genome set. Both methods wrote six complete matrices for all five labels.

A complete ANIb run with `--pan-db` updated a copied Pan Database. Readback
assertions found all 25 identity values and matched each exported matrix value
to the stored numeric value within 1e-12; all six ANI layer orders were present.
Dereplication imported the complete ANIb results at a configured 0.99 threshold
and clustered all five genomes together. After the importer was changed to
initialize its ANI driver only when calculating new results, this same
`--ani-dir` operation succeeded with a PATH containing neither ANI executable.

Legacy PyANI 0.2.12 under Python 3.10.15 was run on the same test genomes and
on the small substitution/half-length FASTAs. On these inputs, ANIb identity,
directional query coverage and alignment lengths matched to serialization
precision. ANIm outputs differed: on the test genomes, identity varied by up
to 7.51e-7, query coverage by 0.002116, and alignment length by 529 bp. The
small fixture's base and half-length sequence show directional coverage 0.5
versus 1.0 and lengths 3060 versus 6120 bp. These comparisons do not establish
equivalence between backend implementations.

Incomplete ANIb/ANIm matrices retain blank undefined comparisons and omit
Newick trees. PanDB updates reject incomplete matrices before mutation, and
dereplication rejects them rather than inventing ANI values.

## NetworkX compatibility before the combined install

Runtime commit `39307922b` changes only the NetworkX dependency range and its
architecture inventory. Current graph source from runtime parent `0852402b0`
was tested with isolated overlays at NetworkX 3.4.2 and 3.7 using CPython 3.13.13.
No graph API changes were necessary; the old DirectedForce/Edmonds code is no
longer present in this source. Existing environments were not modified.

For each version, the current `anvi-pan-genome-graph` CLI processed the prepared
five-genome PAN/storage fixture, including graph construction, topological
layout, region summaries, distances, and graph database output. Both outputs
contained 383 nodes, 405 edges, 32 regions, and 20 genome-distance rows.
Independent review confirmed the source selection, versions, and emitted
artifacts. The recorded command was:

```bash
# v is nx342 or nx37; the isolated venv contains the corresponding version.
cd /tmp/nx-graph-$v
MPLCONFIGDIR=/tmp/nx-mpl-$v PYTHONPATH=/tmp/anvio-pyani-runtime \
  /tmp/anvio-$v/bin/python \
  /home/bcoltman/miniforge3/envs/anvio-dev-313/bin/anvi-pan-genome-graph \
  -p TEST-PAN.db -g TEST-GENOMES.db -e external-genomes.txt -n NX_$v
```

The migration check loaded the current `pan-graph/v5_to_v6.py` module and
called `_migrate_edges` on an in-memory v5-shaped SQLite table. It verified one
reversed edge was retained/flipped and one was cut because `nx.has_path` found
a cycle, including expected genome memberships. Both versions passed. This
validates the migration's NetworkX path operations; it is not a complete
database-version migration. Full evidence and the reproduction fixture are
`/tmp/anvio-networkx-validation.txt` and `/tmp/nx-migration-smoke.py`.

## Combined environment results

The fresh `/tmp/anvio-pyani-joint-313` environment has
`include-system-site-packages = false`. A single joint pip installation from
the feature source completed, and `pip check` reported no broken requirements.
From working directory `/tmp`, anvi'o imported from the feature worktree while
pyani-plus, Typer, Click and Numba imported from this environment. Numba was
installed through anvi'o's declared dependency, without manual preloading.

| Component | Tested version |
| --- | --- |
| Python | 3.13.13 |
| pyani-plus | 1.0.1 |
| NetworkX | 3.7 |
| Typer | 0.27.2 |
| standalone Click | 8.3.1 |
| Numba | 0.68.0 |
| NumPy | 2.5.3 |
| pandas | 3.0.6 |
| Snakemake | 9.27.0 |

Public and private pyani-plus help and worker options passed independent
review. Additional isolated tests routed private workers through Typer 0.27.2
and completed ANIb, ANIm, and matrix export with explicit database, run ID,
label, executor and log options. Typer 0.25.1 with Click 8.3.1 also completed
ANIb/export, but no fallback constraint is needed by the current tested recipe.
The former Typer 0.15.4 recipe is superseded by explicitly installing Click.

Using the combined environment's Python and `pyani-plus` from the same bin
directory, both related-sequence FASTA CLI runs completed with five matrices
and five Newick trees. Assertions confirmed both duplicate names, all four
leaf labels, diagonal values, identity fractions, and query coverage 0.5 versus
1.0 for the half-length sequence. Outputs are under
`/tmp/pyani-real/joint-matching-anib` and `joint-matching-anim`.

Both five-genome genome-similarity CLI runs also completed in the combined
environment. Each wrote six matrices retaining all five labels and blank
undefined pairs, including the derived full identity matrix; all six Newick
trees were omitted with an explicit warning. Assertions checked every output
matrix for scientific null preservation. Outputs are under
`/tmp/pyani-real/joint-genomes-anib` and `joint-genomes-anim`. Ten targeted tests
also passed using the combined interpreter. Existing shared file helpers
emitted ResourceWarnings about unclosed handles; no test failed.

### Linux CPU-affinity smoke

On Linux, the current adapter was exercised with pyANI-plus 1.0.1 and the same
five small genomes. ANIb with `--num-threads 1` completed and logged
`taskset -c 0`; ANIm with `--num-threads 2` completed and logged
`taskset -c 0,1`. Both runs wrote all six matrices with the five input labels
in matching row and column order, and retained undefined comparisons (48 blank
cells across the six matrices per run). The logs record requested, detected,
and allocated CPU counts and explicitly describe affinity rather than an OS
thread limit. The adapter outputs for both methods were compared with
unrestricted raw pyANI-plus exports: all five scientific matrices matched in
labels, nulls, and values within 1e-10, and the derived full-identity matrices
were checked including nulls. The comparison record is
`/tmp/pyani-plus-affinity-smoke/adapter-scientific-comparison.txt`.

`anvi-dereplicate-genomes` also completed an ANIb calculation from four related
fixture genomes with `ANVIO_THREADS=2` and no explicit thread flag. It logged
allocation to CPUs 0 and 1, wrote six complete four-genome matrices, clustered
the three redundant genomes, and reported one representative. The focused
adapter/parser/CLI suite passed 32 tests, including default and `None` API
budgets, invalid budgets, non-contiguous allowed CPU IDs, capping, missing
Linux support checks, and all three CLI help descriptions. Artifacts are under
`/tmp/pyani-plus-affinity-smoke`; logs include `genomes-anib-escalated.log`,
`genomes-anim-escalated.log`, `derep-anib-pyani.log`, `derep-anib` outputs and
`unit-tests.log`.

The initial non-escalated Snakemake attempt failed at forkserver Unix-socket
creation under the sandbox and produced no ANI comparisons; it is excluded from
the successful-run claims above.

Standalone `click` is a required additional package for the tested pyani-plus
CLI installation. `click==8.3.1` and `typer==0.27.2` in the recipe pin the exact
validated versions for reproduction, not a claim that other versions are
incompatible. At the time of this historical run, anvi'o's default dependencies did not
include pyani-plus or Click. The later ANI feature now requires both.

During the earlier optional-backend prototype, the extra was installed directly with `python -m pip install --no-build-isolation -e ".[pyani-plus]"` in the joint environment (no manual pyANI-plus or Click arguments); `pip check` remained clean. Logs: `/tmp/pyani-plus-adoption/extra-install.log` and `extra-pipcheck.log`. That historical extra used API-major dependency bounds; exact versions above are a validation snapshot.
