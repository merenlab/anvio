# Per-rule conda environments for anvi'o workflows

Anvi'o workflows call a number of third-party programs (mappers, assemblers, QC tools).
The long-term goal is that **none** of these need to live in the anvi'o environment itself:
each rule that wraps a third-party tool can instead be pointed at its own standalone conda
environment, which Snakemake builds and activates only for that rule.

This directory ships one environment file per tool, with a pinned, tested version.

## Tools covered

| Rule config key | Tool (executable)        | Env file            |
|-----------------|--------------------------|---------------------|
| `bowtie`        | Bowtie2 (`bowtie2`)      | `bowtie.yaml`       |
| `minimap2`      | minimap2 (`minimap2`)    | `minimap2.yaml`     |
| `megahit`       | MEGAHIT (`megahit`)      | `megahit.yaml`      |
| `metaspades`    | SPAdes (`metaspades.py`) | `metaspades.yaml`   |
| `idba_ud`       | IDBA-UD (`idba_ud`)      | `idba_ud.yaml`      |
| `flye`          | Flye (`flye`)            | `flye.yaml`         |
| `filtlong`      | Filtlong (`filtlong`)    | `filtlong.yaml`     |
| `nanoplot`      | NanoPlot (`NanoPlot`)    | `nanoplot.yaml`     |
| `fastqc`        | FastQC (`fastqc`)        | `fastqc.yaml`       |
| `multiqc`       | MultiQC (`multiqc`)      | `multiqc.yaml`      |

## How to use

**By default a rule runs its tool from your `$PATH`**: set nothing, and anvi'o expects the program
on your `$PATH` and checks for it before starting. To have a rule use a conda environment instead,
set exactly one of the three mutually exclusive options below. Anvi'o raises a ConfigError if you
set more than one.

1. **`use_anvio_conda_yaml`** (boolean; default `false`) — set to `true` to use the env file anvi'o
   ships for this rule (the file in *this* directory). The path is resolved at runtime from the
   installed anvi'o location, so the config stays reproducible across machines — no hard-coded repo
   path. This is the easiest way to get a pinned, tested version without installing anything:

   ```json
   "flye":     { "run": true, "use_anvio_conda_yaml": true },
   "nanoplot": { "run": true, "run_on_raw": true, "use_anvio_conda_yaml": true }
   ```

   Anvi'o auto-adds `--use-conda` to the Snakemake command when a rule uses a conda YAML, so you
   normally don't pass it yourself. The first run builds the env; later runs reuse the cache.

2. **`conda_yaml`** — path to *your own* env YAML (e.g. a customized copy of one of these).
   Snakemake builds the environment from it. Because this is an explicit path it is machine-specific;
   prefer `use_anvio_conda_yaml` when you just want the shipped env:

   ```json
   "flye": { "run": true, "conda_yaml": "/abs/path/to/my-flye.yaml" }
   ```

3. **`conda_env`** — the name of an existing conda environment you have already created
   (e.g. `conda env create -f flye.yaml -n my-flye`). The rule then runs the tool via
   `conda run -n <name> ...` (no `--use-conda` needed):

   ```json
   "flye": { "run": true, "conda_env": "my-flye" }
   ```

Set at most one of these per rule; if more than one is in effect anvi'o errors rather than guessing.

To run a tool from your **`$PATH`** — the default — just leave all three off
(`use_anvio_conda_yaml: false`, and `conda_yaml` / `conda_env` empty).

## Notes

- `nanoplot.yaml` pins only the NanoPlot version and leaves the rest of the dependency tree to
  the solver, which resolves an importable environment; pinning `python`/`python-kaleido` breaks
  the solve. See that file's header for details.
- `flye`/`minimap2` versions track the presets validated in
  `anvio/workflows/lr_technology_presets.yaml`.
- Still on `$PATH` (no conda hook yet): `samtools`, `krakenuniq`, `centrifuge`, `emapper`.
  These are the remaining rules to migrate as this effort continues.
