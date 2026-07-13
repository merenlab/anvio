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

Two options per rule, set under the rule's block in your workflow config:

1. **`conda_yaml`** — path to one of these YAML files. Snakemake builds the environment
   from it. Run the workflow with `--use-conda`:

   ```json
   "flye":     { "run": true, "conda_yaml": "/path/to/anvio/workflows/conda_envs/flye.yaml" },
   "nanoplot": { "run": true, "run_on_raw": true, "conda_yaml": ".../conda_envs/nanoplot.yaml" }
   ```

   ```
   anvi-run-workflow -w metagenomics -c config.json \
       --additional-params --use-conda --conda-frontend mamba
   ```

2. **`conda_env`** — the name of an existing conda environment you have already created
   (e.g. `conda env create -f flye.yaml -n my-flye`). The rule then runs the tool via
   `conda run -n <name> ...`; no `--use-conda` needed.

Set only one of `conda_yaml` / `conda_env` per rule. If neither is set, the tool must be
on `$PATH` (the legacy behavior).

## Notes

- `nanoplot.yaml` pins `python-kaleido=0.2.1` and `python <3.13`: NanoPlot 1.42 imports the
  legacy `kaleido.scopes` API that `python-kaleido >= 1.0` removed, and an unpinned solve
  otherwise pulls python 3.14 + kaleido 1.x and crashes on import.
- `flye`/`minimap2` versions track the presets validated in
  `anvio/workflows/lr_technology_presets.yaml`.
- Still on `$PATH` (no conda hook yet): `samtools`, `krakenuniq`, `centrifuge`, `emapper`.
  These are the remaining rules to migrate as this effort continues.
