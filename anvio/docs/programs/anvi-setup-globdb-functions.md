This program **downloads and sets up a local copy of the GlobAA gene family database** for use in functional annotation with %(anvi-run-globdb-functions)s. It produces a %(globdb-data)s artifact.

### Basic usage

{{ codestart }}
anvi-setup-globdb-functions
{{ codestop }}

{:.warning}
We recommend using `--num-threads` to speed up the DIAMOND database build step.

If you already have a %(globdb-data)s artifact and want to re-download and rebuild everything from scratch:

{{ codestart }}
anvi-setup-globdb-functions --reset
{{ codestop }}

### Custom data directory

By default, anvi'o stores the GlobAA data in a location inside the anvi'o package directory. If you do not have write access to that location, or if you want to keep the data elsewhere, use:

{{ codestart }}
anvi-setup-globdb-functions --globdb-data-dir /path/to/your/directory
{{ codestop }}

You can also set the environment variable `ANVIO_GLOBAA_DATA_DIR` to your preferred path so anvi'o will use it automatically without requiring the `--globdb-data-dir` flag each time:

{{ codestart }}
export ANVIO_GLOBAA_DATA_DIR=/path/to/your/directory
anvi-setup-globdb-functions
{{ codestop }}

### What happens during setup

1. The GlobAA data package (that is maintained by GlobDB folk, including Daan Speth et al) is downloaded and extracted.
2. Every gene family YAML file is validated for required fields (`gene_family`, `description`, `version`, and `cutoffs` including `lasr`, `selfmax`, `selfmin`, and `matrix`).
3. All per-family FASTA files are concatenated into a single `GlobAA.faa` (with GAA identifiers prepended to sequence headers).
4. All per-family YAML files are merged into a single `GlobAA.yaml`.
5. A DIAMOND search database is built from `GlobAA.faa`.
