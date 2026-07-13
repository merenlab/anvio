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
| `fastqc_sr`     | FastQC (`fastqc`)        | `fastqc_sr.yaml`    |
| `multiqc`       | MultiQC (`multiqc`)      | `multiqc.yaml`      |

## How to use

Three mutually exclusive options per rule, set under the rule's block in your workflow config.
Anvi'o raises a ConfigError if you set more than one.

1. **`use_anvio_conda_yaml`** (boolean, **the default — `true`**) — use the env file anvi'o ships
   for this rule (the file in *this* directory). The path is resolved at runtime from the installed
   anvi'o location, so the config stays reproducible across machines — no hard-coded repo path.
   This is the recommended option, and because it is the default you usually don't set anything:

   ```json
   "flye":     { "run": true },
   "nanoplot": { "run": true, "run_on_raw": true }
   ```

   Anvi'o auto-adds `--use-conda` to the Snakemake command when a rule uses a conda YAML, so you
   normally don't pass it yourself. The first run builds the env; later runs reuse the cache.

2. **`conda_yaml`** — path to *your own* env YAML (e.g. a customized copy of one of these).
   Snakemake builds the environment from it. Because this is an explicit path it is machine-specific;
   prefer `use_anvio_conda_yaml` when you just want the shipped env. Requires
   `use_anvio_conda_yaml: false`:

   ```json
   "flye": { "run": true, "use_anvio_conda_yaml": false, "conda_yaml": "/abs/path/to/my-flye.yaml" }
   ```

3. **`conda_env`** — the name of an existing conda environment you have already created
   (e.g. `conda env create -f flye.yaml -n my-flye`). The rule then runs the tool via
   `conda run -n <name> ...` (no `--use-conda` needed). Requires `use_anvio_conda_yaml: false`:

   ```json
   "flye": { "run": true, "use_anvio_conda_yaml": false, "conda_env": "my-flye" }
   ```

Because `use_anvio_conda_yaml` defaults to `true`, options 2 and 3 require you to ALSO set
`use_anvio_conda_yaml: false` on the same rule (otherwise anvi'o errors — it won't guess).

To run a tool from your **`$PATH`** instead (the classic behavior), set `use_anvio_conda_yaml: false`
and leave `conda_yaml` / `conda_env` empty.

## Notes

- `nanoplot.yaml` pins `python-kaleido=0.2.1` and `python <3.13`: NanoPlot 1.42 imports the
  legacy `kaleido.scopes` API that `python-kaleido >= 1.0` removed, and an unpinned solve
  otherwise pulls python 3.14 + kaleido 1.x and crashes on import.
- `flye`/`minimap2` versions track the presets validated in
  `anvio/workflows/lr_technology_presets.yaml`.
- Still on `$PATH` (no conda hook yet): `samtools`, `krakenuniq`, `centrifuge`, `emapper`.
  These are the remaining rules to migrate as this effort continues.
