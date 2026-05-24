# The Anvi'o Tool Registry

The registry is a small system that centralizes metadata about every third-party program and database anvi'o depends on, and provides a single execution layer that can dispatch those programs to different runtimes (i.e., a local binary in `PATH`, a Docker container, a Singularity image, or an isolated conda environment) without any driver needing to know which one.

The problem it solves is scattered, ad hoc tool invocation which we have been relying upon until now. Before the registry, every driver in `anvio/drivers/` hard-coded its binary name, called `utils.run_command()` or `subprocess.Popen()` directly, and duplicated tool-availability checks. Switching a single analysis to run inside a container required touching every driver that analysis uses. With the registry, that switch is one environment variable:

```bash
export ANVIO_RUNNER=docker
anvi-run-hmms ...
```

or

```bash
anvi-run-hmms --runner docker ...
```

etc.

---

## Key files and directories

Here is where things related to the registry and its downstream use are:

```
anvio/registry/
  third-party-programs.yaml
  anvio-programs.yaml
  databases.yaml

anvio/exec.py (the draft module in whihc the Registry and AnvioExecutor classes live)
```

---

## The YAML files

### `third-party-programs.yaml`

One entry per external tool. The entry key is the **registry name**, the string you pass to `AnvioExecutor(...)`. It does not have to match the binary name (e.g., the entry for the MUMmer suite is `mummer`, not `nucmer`).

```yaml
diamond:
  authors: Buchfink et al.
  citation: doi:10.1038/s41592-021-01101-x
  url: https://github.com/bbuchfink/diamond
  version: 2.1.9
  runners:
    conda: bioconda::diamond={version}
    docker: buchfink/diamond:v{version}
    singularity: https://depot.galaxyproject.org/singularity/diamond:{version}--HASH_PLACEHOLDER
```

**Fields**

| Field | Required | Notes |
|---|---|---|
| `authors` | yes | Last name(s), e.g. `Buchfink et al.` |
| `citation` | yes | DOI string |
| `url` | yes | Project home page |
| `version` | yes | The version pinned in the standard conda environment |
| `runners.conda` | recommended | `channel::package={version}` (the `{version}` token is resolved at load time) |
| `runners.docker` | optional | Full image tag; `{version}` is resolved; use `HASH_PLACEHOLDER` for the biocontainers build hash |
| `runners.singularity` | optional | Full image URL; same rules as docker |

**`{version}` resolution.** The `Registry` class replaces every `{version}` token in runner strings at load time. Callers always receive fully-resolved strings from `registry.get_runner()`.

**`HASH_PLACEHOLDER`.** Biocontainers image tags include a build hash (e.g. `3.4--hdbdd923_2`). Fill these in from [https://biocontainers.pro](https://biocontainers.pro) or by browsing `https://depot.galaxyproject.org/singularity/`. `AnvioExecutor` raises a `ConfigError` with a clear message if you try to use a backend whose image tag still contains `HASH_PLACEHOLDER`.

**Tool suites.** When a package provides multiple binaries (e.g. `hmmer` provides `hmmbuild`, `hmmpress`, `hmmscan`, `hmmsearch`), use a single entry for the package. The `AnvioExecutor` key is the package name; the argv list you pass to `run_command()` specifies the actual binary (this is a bit finnicky as is, and we can register all binaries that come with a package in the YAML file at some point).

---

### `anvio-programs.yaml`

Metadata for anvi'o programs, used for citation generation. The `anvi-*` wildcard entry provides defaults merged into every specific entry.

```yaml
anvi-*:
  authors: Anvi'o Authors
  citation: doi:10.1038/s41564-020-00834-3

anvi-estimate-metabolism:
  authors: Veseli et al.
  citation: doi:10.1038/s41564-023-01439-2
  uses:
    - kegg-data
  example_text: >-
    Metabolic pathway completeness was estimated using {self}
    ({self.citation}) against the KEGG database ({kegg-data.citation}).
```

**Fields**

| Field | Notes |
|---|---|
| `authors` | Last name(s) or GitHub handles for anvi'o developers |
| `citation` | DOI of the primary paper describing this program's method |
| `uses` | List of `databases.yaml` keys this program relies on |
| `example_text` | Optional methods-section template; `{self}`, `{self.<attr>}`, and `{<key>.<attr>}` tokens are resolved by `ProgramEntry.get_methods_text()` |

Fields in a specific entry always override the wildcard defaults. A program not listed by name still gets the wildcard defaults (main anvi'o citation) when queried through `registry.get_program()`.

---

### `databases.yaml`

Reference databases that anvi'o downloads and manages. These are distinct from tools (as in, a database is a data artifact, not an executable).

```yaml
kegg-data:
  authors: Kanehisa et al.
  provider: third-party
  description: Kyoto Encyclopedia of Genes and Genomes (KEGG)
  citation: doi:10.1093/nar/gkw1092
  url: https://www.genome.jp/kegg/
```

**Fields**

| Field | Notes |
|---|---|
| `provider` | `third-party` (external group) or `anvio` (generated/distributed by anvi'o) |
| `description` | The name of the database (as official as possible) |
| `citation` | DOI of the paper that should be cited when the database is used |

---

## The Python API

### `Registry`

`anvio.exec.registry` is a module-level singleton. It is **not** populated at import time; the three YAML files are parsed on the first call to any lookup method, then cached for the life of the process. Thread-safe via a double-checked lock.

Here is how it looks:

```python
from anvio.exec import registry

# Third-party tools returns a ToolEntry or None
tool   = registry.get_tool('diamond')
runner = registry.get_runner('diamond', 'docker')   # → resolved image string or None
runner = registry.get_runner('diamond', 'conda')    # → 'bioconda::diamond=2.1.9'
print(tool.version)            # '2.1.9'
print(tool.runners['conda'])   # 'bioconda::diamond=2.1.9'

# Anvi'o programs returns a ProgramEntry or None
prog = registry.get_program('anvi-estimate-metabolism')
print(prog.citation)           # doi:10.1038/s41564-023-01439-2
print(prog.authors)            # Veseli et al.
print(prog.uses)               # ['kegg-data']
print(prog.get_methods_text()) # Which currently takes this text:
#
#   Metabolic pathway completeness was estimated using {self}
#   ({self.citation}) against the {kegg-data.description} ({kegg-data.citation}).
#
# and renders it as the following
#
#   Metabolic pathway completeness was estimated using
#   anvi-estimate-metabolism (doi:10.1038/s41564-023-01439-2)
#   against the Kyoto Encyclopedia of Genes and Genomes (KEGG)
#   (doi:10.1093/nar/gkw1092).

# Reference databases returns a DatabaseEntry or None
db = registry.get_database('kegg-data')
print(db.citation)             # doi:10.1093/nar/gkw1092
print(db.description)          # Kyoto Encyclopedia of Genes and Genomes
print(db.provider)             # 'third-party'
```

All getters return `None` for unknown keys. They never raise.

---

### `AnvioExecutor`

The dispatcher. Construct one per tool invocation (cheap; no subprocess is started on construction) and call `run_command()`.

```python
from anvio.exec import AnvioExecutor

executor = AnvioExecutor(
    tool_name='diamond',   # registry key
    backend=None,          # None → read $ANVIO_RUNNER → default 'native'
    work_dir=None,         # None → os.getcwd(); used as bind-mount root for containers
    run=self.run,          # pass your instance's run/progress objects
    progress=self.progress,
)

result = executor.run_command(
    cmd=['diamond', 'blastp', '-q', query, '-d', db, '-o', out],
    raise_on_error=True,   # raises CommandError on non-zero exit (default)
)

# result is a RunResult namedtuple:
result.stdout      # str for captured stdout
result.stderr      # str for captured stderr
result.returncode  # int
result.cmd         # list for the final argv executed (includes container wrappers)
```

**Backend selection precedence** (highest to lowest):

1. Explicit `backend=` constructor argument
2. `$ANVIO_RUNNER` environment variable
3. `'native'` (run directly from `PATH`)

**`get_conda_env_name()`** returns the expected conda environment name so users can create it before switching to the conda backend:

```python
executor = AnvioExecutor('diamond')
print(executor.get_conda_env_name())   # 'anvio-tool-diamond-2.1.9'
```

We also will have a global parameter (like `--debug` which sets `anvio.DEBUG`) called `--runner` to set `anvio.RUNNER`.

---

## Backend reference

### `native` (default)

Runs the binary directly from `PATH` via `subprocess.run`. Zero overhead. This is what every existing driver already does, just routed through a uniform return type.

```python
executor = AnvioExecutor('diamond')
result = executor.run_command(['diamond', 'blastp', ...])
```

---

### `docker`

Wraps the command in `docker run --rm` with the `work_dir` bind-mounted at the same absolute path inside the container:

```
docker run --rm \
  -v /data/run1:/data/run1 \
  -w /data/run1 \
  -u 1000:1000 \
  buchfink/diamond:v2.1.9 \
  diamond blastp ...
```

File paths in `cmd` resolve identically inside and outside the container because the directory is mounted at the same path (the same strategy Nextflow uses for process work directories).

**Requirements:**
- Docker daemon must be running (`docker info` must succeed)
- The image tag must be filled in (no `HASH_PLACEHOLDER`)
- All input and output files must live under `work_dir`

Set `work_dir` explicitly when your files are not under `os.getcwd()`:

```python
executor = AnvioExecutor('diamond', backend='docker', work_dir='/tmp/xxx')
```

---

### `singularity`

Wraps the command in `singularity exec` with `work_dir` bound:

```
singularity exec \
  --bind /data/run1 \
  --pwd  /data/run1 \
  https://depot.galaxyproject.org/singularity/diamond:2.1.9--h43eeafb_0 \
  diamond blastp ...
```

Same file-path resolution assumptions as docker.

**Requirements:**
- `singularity` must be in `PATH`
- The image URL must have its `HASH_PLACEHOLDER` filled in

---

### `conda`

Looks up the conda package spec in the registry, derives a canonical environment name (e.g. `anvio-tool-diamond-2.1.9`), checks that the environment exists, then wraps with `conda run`:

```
mamba run --no-capture-output --name anvio-tool-diamond-2.1.9 -- diamond blastp ...
```

Prefers `mamba` over `micromamba` over `conda` in that order. Does **not** create the environment automatically since environment creation is slow and should be an explicit user action. If the environment is missing, the error message includes the exact `conda create` command needed:

```
conda create -n anvio-tool-diamond-2.1.9 -c bioconda -c conda-forge diamond=2.1.9
```

---

## Migration guide

The goal is to replace direct `utils.run_command()` calls and raw `subprocess.Popen()` calls in driver classes with `AnvioExecutor`, one driver at a time. Existing drivers keep working until migrated. Nothing about `utils.run_command()` changes.

There are multiple patterns in the codebase, and how we run commands clearly changed over time. Here is a proposed guideline for every pattern:

### Pattern 1: `utils.run_command()` (most drivers use this)

This is the dominant pattern: build an argv list, call `utils.run_command`, check whether an expected output file was produced. The log file is written as a side effect.

**Before (from `anvio/drivers/mcl.py`):**

```python
cmd_line = ['mcl',
            self.mcl_input_file_path,
            '--abc',
            '-I', self.inflation,
            '-o', self.clusters_file_path,
            '-te', self.num_threads]

self.run.info('mcl cmd', ' '.join([str(x) for x in cmd_line]), quiet=True)

utils.run_command(cmd_line, log_file_path)

self.check_output(self.clusters_file_path, 'makedb')
```

**After:**

```python
from anvio.exec import AnvioExecutor

cmd_line = ['mcl',
            self.mcl_input_file_path,
            '--abc',
            '-I', self.inflation,
            '-o', self.clusters_file_path,
            '-te', self.num_threads]

executor = AnvioExecutor('mcl', run=self.run, progress=self.progress, log=log_file_path)
result = executor.run_command(cmd_line)

self.check_output(self.clusters_file_path, 'makedb')
```

The key difference: `utils.run_command` tees stdout and stderr to a log file and discards the captured bytes. But `AnvioExecutor` captures them in `result.stdout` and `result.stderr` and, when `log=` is provided, writes a structured log entry automatically (including the full command, backend, return code, and both output streams). The log is always written, even when the command fails, so it is available for debugging before the `CommandError` propagates. For drivers that do not need a persistent log file, omit the `log=` argument entirely.

---

### Pattern 2: `subprocess.Popen` with stdout captured (FastTree)

FastTree writes the Newick tree to stdout, which the existing driver reads via `Popen.communicate()`.

**Before (from `anvio/drivers/fasttree.py`):**

```python
from subprocess import Popen, PIPE

fasttree = Popen(['FastTree'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
output = fasttree.communicate(input=input_file.read())
output_stdout = output[0].decode().rstrip()
output_stderr = output[1].decode().splitlines()

run.info("Version", output_stderr[0])
# ... parse stderr for warnings ...
```

**After:**

```python
from anvio.exec import AnvioExecutor

executor = AnvioExecutor('fasttree', run=self.run)
result = executor.run_command(['FastTree', '-nt', '-gtr', input_file_path])

tree_string = result.stdout.rstrip()
stderr_lines = result.stderr.splitlines()

self.run.info("Version", stderr_lines[0])
# ... parse stderr_lines as before ...
```

**Note on stdin piping.** The existing FastTree driver passes the alignment via stdin. The current `AnvioExecutor.run_command()` does not support stdin injection. Pass the input as a file argument instead (FastTree accepts a filename directly), or use the `native` backend with a direct `subprocess.run(input=data, ...)` call for now. Stdin support can be added to `AnvioExecutor` when needed, or better yet, we can make things work with files and stop passing stuff through stdin (I had implemented this to overcome the I/O overhead of dealing with files due to a slow NFS at the MBL; maybe it is a problem of a distant past now).

---

### Pattern 3: `utils.run_command_and_get_output()` (more recent drivers use this)

This pattern already captures stdout as a string. The migration is straightforward.

**Before:**

```python
output = utils.run_command_and_get_output(
    ['samtools', 'view', '-c', bam_path],
    raise_on_error=True,
)
count = int(output.strip())
```

**After:**

```python
from anvio.exec import AnvioExecutor

executor = AnvioExecutor('samtools', run=self.run)
result = executor.run_command(['samtools', 'view', '-c', bam_path])
count = int(result.stdout.strip())
```

---

### Pattern 4: `utils.is_program_exists()` at driver construction

Drivers currently call `utils.is_program_exists('diamond')` in `__init__` to fail early with a helpful message when the binary is missing. Keep this call for the `native` backend. For other backends the binary is not in `PATH` by design, so the check would always fail. Skip it when `backend != 'native'`:

```python
from anvio.exec import AnvioExecutor, DEFAULT_BACKEND

class Diamond:
    def __init__(self, ...):
        if DEFAULT_BACKEND == 'native':
            utils.is_program_exists('diamond')

        self.executor = AnvioExecutor('diamond', run=self.run)
```

---

### What not to migrate at this stage

These cases are not a good fit for `AnvioExecutor` in the current form of the class (which can change, but I think shouldn't):

- **STDIN piping** (`utils.run_command_STDIN`, `Popen(stdin=PIPE)`): drivers that stream data into a tool via stdin (FAMSA, muscle) need either a stdin parameter on `run_command()` or a different approach. These must be taken care of, or left alone until that feature is added.
- **Streaming / long-running processes** (`utils.start_command`, returning a `Popen` object): some drivers start a process and monitor it while it runs. `AnvioExecutor.run_command` waits for completion. We either will have to address this and change driver behavior, or leave them as is also.

---

## Adding a new tool

Adding a new tool is easy.

1. **Add an entry to `third-party-programs.yaml`**:

   ```yaml
   mynewtool:
     authors: Smith et al.
     citation: doi:10.1000/xyz123
     url: https://github.com/example/mynewtool
     version: 1.2.3
     runners:
       conda: bioconda::mynewtool={version}
       docker: quay.io/biocontainers/mynewtool:{version}--HASH
       singularity: https://depot.galaxyproject.org/singularity/mynewtool:{version}--HASH
   ```

2. **Use `AnvioExecutor` in the driver**:

   ```python
   from anvio.exec import AnvioExecutor

   class MyNewTool:
       def __init__(self, run=run, progress=progress):
           self.run = run
           self.progress = progress

       def run_analysis(self, input_file, output_file):
           executor = AnvioExecutor('mynewtool', run=self.run)
           result = executor.run_command(
               ['mynewtool', '--input', input_file, '--output', output_file]
           )
           self.run.info('MyNewTool output', output_file)
   ```

3. **Fill in the `HASH` and confirm the URL for the program** by looking up the image on [https://biocontainers.pro](https://biocontainers.pro) or `https://depot.galaxyproject.org/singularity/mynewtool/`.

4. **Add a program entry to `anvio-programs.yaml`** if this tool is exposed through an `anvi-*` CLI program (see that file for examples).

---

## Environment variable reference

| Variable | Values | Effect |
|---|---|---|
| `ANVIO_RUNNER` | `native`, `docker`, `singularity`, `conda` | Sets the default backend for all `AnvioExecutor` instances in the process. Overridable per-instance via `backend=` argument. |

---

## Design notes

**YAML and not Python**, because the metadata (citations, versions, image tags) is consumed by documentation generators, future CLI tools (`anvi-check-programs`, methods-section generators), and potentially external scripts. YAML is human-editable, diff-friendly, and parseable without importing the full anvi'o package.

**Lazy singleton for `Registry`** because most anvi'o programs never need to query the registry at all, so loading three YAML files on every import would add latency for nothing. The lazy load is transparent to callers.

**No auto-pulling missing container images** because silently downloading a 500 MB image mid-analysis would make users hate us. The executor should fail loudly if the image is absent, and the error message should tell the user exactly which image to pull if they requested a docker executor.

**No auto-creating conda envs** because the same thing applies: creating conda envs can take a lot of time, and it should be the user who makes that call.

**Relationship to `utils.run_command`.** `AnvioExecutor` is additive: `utils.run_command` and all existing drivers are untouched. Migration is incremental: drivers can be updated one at a time. The only shared state is the `registry` singleton, which is read-only after load.
