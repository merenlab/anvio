This artifact represents **a local copy of the GlobDBFunctions gene family database** set up on your computer for use with %(anvi-run-globdb-functions)s.

The source of %(globdb-data)s is [GlobDB](https://globdb.org/home), a curated database of a dereplicated set of 'species' representative genomes of Bacteria and Archaea.

The GlobDB Functions database, which underlies %(globdb-data)s artifact in the anvi'o ecosystem, stores functional information on protein families at a fine-grained level, with each gene family corresponding to a protein or subunit of a protein complex with a specific physiological function. Its increasing coverage aims to focus on gene families relevant for microbiology and microbial ecology where the physiological function has not yet been experimentally determined, and thus will be complementary to the KEGG database. The GlobDB Functions database is also complementary to other general databases for functional annotation such as PFAM and COG, which often include broader annotations. The gene family metadata includes the KEGG, COG, and PFAM annotation where available, and whether a gene family has members that have been experimentally characterized.

GlobDBFunctions gene families are defined using the amino acid sequence toolkit [(AASTK)](https://github.com/dspeth/aastk), and include synteny information derived from consensus genomic context analyses. Per sequence metadata including other annotations, genome taxonomy, availability of the organism in culture, and environmental detection is included. AASTK allows for automatic updates of the gene family datasets with each annual release of the GlobDB, ensuring that the GlobDB Functions gene family databases remain up to date over time.

{:.warning}
Please note that GlobDB Functions database is in development, and the tools associated with this database may not be in their final form if you are accessing them through `anvio-dev`. We do appreciate your patience, and look forward to receiving your input to improve the anvi'o tools around this resource.

{:.notice}
GlobDB developers welcome user created/curated submissions for functions to maximize the utility of the GlobDB Functions. Please reach out to [Daan Speth](https://dome.univie.ac.at/speth) for more information.

The GlobDB Functions database contains amino acid sequences and empirically-derived quality cutoffs for a curated set of microbial gene families. Each gene family is identified by a GAA accession number (e.g., `GAA00000001`) and carries gene-family-level Local Alignment Score Ratio (LASR) cutoffs that allow anvi'o to assess the quality of each annotation hit.

## What it contains

After running %(anvi-setup-globdb-functions)s, the data directory will contain:

- `GlobDB.faa` — all gene family sequences concatenated into a single FASTA file (headers prefixed with the GAA identifier for fast lookup)
- `GlobDB-gene-family-data.yaml` — a master YAML file containing gene family metadata and LASR cutoffs for every gene family
- `GlobDB-synteny-data.yaml` — a master YAML file containing synteny information for each gene family where this is available
- `GlobDB-validation-data.yaml` — a master YAML file containing experimental validation data for each gene family where this is available
- `DB_DIAMOND/GlobDB.dmnd` — a DIAMOND search database built from `GlobDB.faa`
- `.VERSION` — an internal version file used by anvi'o to detect incompatible setups

## Set it up

{{ codestart }}
anvi-setup-globdb-functions
{{ codestop }}

To force a complete re-download and rebuild:

{{ codestart }}
anvi-setup-globdb-functions --reset
{{ codestop }}

To store the data in a custom location instead of the anvi'o default:

{{ codestart }}
anvi-setup-globdb-functions --globdb-data-dir /path/to/your/directory
{{ codestop }}

You can also set the environment variable `ANVIO_GLOBDB_DATA_DIR` to your custom path so you do not need to pass it every time.
