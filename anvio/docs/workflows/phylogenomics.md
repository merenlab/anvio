The phylogenomics workflow starts with one or more %(contigs-db)s files and extracts a set of genes defined by HMMs. It then concatenates the recovered sequences, aligns them, trims the alignment, and infers a phylogenomic tree.

The workflow is meant for cases where you want to build a tree from homologous genes already annotated in anvi'o contigs databases. It can use internal genomes, external genomes, or a mix of both, as long as they are provided through the workflow config.

## Required input

The phylogenomics workflow requires a %(workflow-config)s file. You can generate a default config like this:

{{ codestart }}
anvi-run-workflow -w phylogenomics \
                  --get-default-config config.json
{{ codestop }}

The workflow config typically includes:

1. `project_name`, which is used as the prefix for workflow outputs.
2. `internal_genomes` and/or `external_genomes`, which point to the genomes or metagenomes that should be used.
3. Parameters for `anvi_get_sequences_for_hmm_hits`, `trimal`, and `iqtree`.

An example minimal config looks like this:

```json
{
    "workflow_name": "phylogenomics",
    "config_version": "3",
    "project_name": "phylo_project",
    "internal_genomes": "internal-genomes.txt",
    "external_genomes": "external-genomes.txt",
    "anvi_get_sequences_for_hmm_hits": {
        "--return-best-hit": true,
        "--align-with": "famsa",
        "--concatenate-genes": true,
        "--get-aa-sequences": true,
        "--hmm-sources": "Bacteria_71"
    },
    "trimal": {
        "-gt": 0.5
    },
    "iqtree": {
        "threads": 8,
        "-m": "WAG",
        "-bb": 1000
    },
    "output_dirs": {
        "PHYLO_DIR": "01_PHYLOGENOMICS",
        "LOGS_DIR": "00_LOGS"
    }
}
```

The `project_name` is mandatory. The workflow uses it to name the output FASTA, alignment, and tree files.

## Run it

Create a workflow graph first if you want to inspect the plan:

{{ codestart }}
anvi-run-workflow -w phylogenomics \
                  -c config.json \
                  --save-workflow-graph
{{ codestop }}

Then run the workflow:

{{ codestart }}
anvi-run-workflow -w phylogenomics \
                  -c config.json
{{ codestop }}

If everything completes successfully, you should end up with a concatenated amino acid FASTA, a trimmed alignment, and a final tree in the phylogenomics output directory.

## Output structure

The workflow writes its main outputs under `01_PHYLOGENOMICS/` by default.

Typical files include:

```text
01_PHYLOGENOMICS/
├── PROJECT-proteins.fa
├── PROJECT-proteins_GAPS_REMOVED.fa
└── PROJECT-proteins_GAPS_REMOVED.fa.contree
```

The intermediate files represent the main stages:

1. `anvi-get-sequences-for-hmm-hits` extracts the target proteins.
2. `trimal` trims the multiple sequence alignment.
3. `iqtree` infers the phylogenomic tree.

Workflow logs are written under `00_LOGS/phylogenomics` by default. Logs are organized by rule name, and the workflow also writes a tab-delimited manifest named `00_LOGS/phylogenomics/phylogenomics-workflow-manifest.tsv` that records whether each job succeeded or failed and points to the relevant rule log.

## Notes

This workflow inherits the contigs workflow, so the same contigs database setup and log organization conventions apply. If you are building your phylogeny from anvi'o HMM hits, the `--return-best-hit` and `--concatenate-genes` settings are usually the important ones to review first.
