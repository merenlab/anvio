
This workflow takes a set of genomes -- either as %(fasta)s files, as already-existing %(contigs-db)s files, or a mix of both -- and creates a pangenome (%(pan-db)s) you can visualize with %(anvi-display-pan)s. Along the way it can also compute a phylogenetic tree for your genomes and order the layers of your pangenome with it, and estimate genome similarity (i.e. ANI) across your genomes.

The pangenomics workflow inherits the contigs workflow, so everything the contigs workflow can do for you (reformatting FASTA files, gene calling, functional annotation, SCG taxonomy, and so on) is also available here through the very same %(workflow-config)s parameters. On top of that, it adds the following steps:

1. %(anvi-gen-genomes-storage)s to build a %(genomes-storage-db)s from your genomes.
2. %(anvi-pan-genome)s to compute gene clusters and produce a %(pan-db)s.
3. Optionally, a phylogenetic tree (from gene clusters or from HMM hits) that is imported into the %(pan-db)s as %(misc-data-layer-orders)s.
4. Optionally, %(anvi-compute-genome-similarity)s to add %(genome-similarity)s estimates into the %(pan-db)s.

{:.warning}
Since this workflow inherits the contigs workflow, the default %(workflow-config)s will try to annotate your contigs databases with NCBI COGs, KEGG KOfams, and SCG taxonomy. If you have not yet run %(anvi-setup-ncbi-cogs)s, %(anvi-setup-kegg-data)s, and %(anvi-setup-scg-taxonomy)s on your system, you will get cryptic errors from this workflow. You can avoid this by first running those setup programs (which is done only once for every anvi'o installation), **or** by explicitly setting `run=false` for the `anvi_run_ncbi_cogs`, `anvi_run_kegg_kofams`, and `anvi_run_scg_taxonomy` rules in your config file.

## Required input

The workflow needs a %(workflow-config)s file. To get a default one:

{{ codestart }}
anvi-run-workflow -w pangenomics \
                  --get-default-config config-pangenomics-default.json
{{ codestop }}

The resulting file is quite long since it includes every parameter of every rule inherited from the contigs workflow. It is worth investigating its contents to see everything you can tweak, but only two things in it are mandatory:

1. `project_name` -- the workflow will refuse to run without it, and uses it to name most of its outputs.
2. `internal_genomes` and/or `external_genomes` -- you must provide at least one of them (see below).

### Describing your genomes

There are two ways to get your genomes into this workflow, and you can combine them freely.

**If you are starting from FASTA files**, you need to (1) manually create a %(fasta-txt)s file that describes the name and location of each FASTA file, (2) point the `fasta_txt` variable in your config to it, and (3) *also* give the `external_genomes` variable a file name. That last step is easy to miss:

{:.notice}
The %(external-genomes)s file **does not have to exist** when you start. Anvi'o will generate it for you from the information in your %(fasta-txt)s, but you still have to tell anvi'o what to call it. If you set `fasta_txt` without setting `external_genomes`, the workflow will stop and complain.

**If you are starting from contigs databases you already have**, simply point `external_genomes` to an existing %(external-genomes)s file and/or `internal_genomes` to an existing %(internal-genomes)s file, and set `"fasta_txt": ""` so the workflow does not look for FASTA files it does not need. If `fasta_txt` points to a path that does not exist on your disk, the workflow will stop with an error.

Here is a minimal config that starts from FASTA files, skips the annotation steps that require setting up external databases, and computes genome similarity:

```json
{
    "workflow_name": "pangenomics",
    "config_version": "6",
    "project_name": "TEST",
    "fasta_txt": "fasta.txt",
    "external_genomes": "external-genomes.txt",
    "anvi_run_hmms": {
        "run": false
    },
    "anvi_run_ncbi_cogs": {
        "run": false
    },
    "anvi_compute_genome_similarity": {
        "run": true
    },
    "output_dirs": {
        "FASTA_DIR": "01_FASTA",
        "CONTIGS_DIR": "02_CONTIGS",
        "PAN_DIR": "03_PAN",
        "LOGS_DIR": "00_LOGS"
    }
}
```

{:.notice}
By default, the name of your %(pan-db)s comes from `project_name`. If you set `--project-name` under the `anvi_pan_genome` rule, that value will be used for the pangenome instead (and anvi'o will warn you about it), while the %(genomes-storage-db)s will still be named after `project_name`.

## Run it

To see if everything looks alright, you can generate a 'workflow graph' for your config file and input files:

{{ codestart }}
anvi-run-workflow -w pangenomics \
                  -c config-pangenomics.json \
                  --save-workflow-graph
{{ codestop }}

{:.notice}
Please note that the generation of this workflow graph requires the usage of a program called [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)). If you are using MAC OSX, you can use [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) by installing [graphviz](http://www.graphviz.org/) through `brew` or `conda`, and if you are on WSL, you can install `graphviz` through `sudo apt install`.

You can also ask the workflow which programs it is going to need, before it needs them:

{{ codestart }}
anvi-run-workflow -w pangenomics \
                  -c config-pangenomics.json \
                  --list-dependencies
{{ codestop }}

If everything looks alright, run the workflow:

{{ codestart }}
anvi-run-workflow -w pangenomics \
                  -c config-pangenomics.json
{{ codestop }}

At the end of it all you should have a %(genomes-storage-db)s and a %(pan-db)s in your `03_PAN` directory, which you can visualize right away:

{{ codestart }}
anvi-display-pan -g 03_PAN/TEST-GENOMES.db \
                 -p 03_PAN/TEST-PAN.db
{{ codestop }}

## Adding a phylogeny to your pangenome

If you set the `sequence_source_for_phylogeny` variable in your config, the workflow will additionally compute a phylogenetic tree for your genomes and import it into your %(pan-db)s as a layer order, so you can order the genomes in the pangenome display by their phylogenetic relationships rather than by gene cluster presence/absence. There are two valid sources:

* `gene_clusters` -- the tree is computed from gene clusters in your pangenome, exported with %(anvi-get-sequences-for-gene-clusters)s.
* `hmm` -- the tree is computed from HMM hits in your contigs databases, exported with %(anvi-get-sequences-for-hmm-hits)s.

In both cases the sequences are then trimmed with `trimal` and the tree is inferred with `iqtree`, exactly as in the phylogenomics workflow. The resulting tree is imported into the %(pan-db)s under the name you give with the `tree_name` parameter of the `import_phylogenetic_tree_to_pangenome` rule (which defaults to `phylogeny`).

Here are the relevant parts of a config that builds the phylogeny from single-copy core gene clusters:

```json
{
    "sequence_source_for_phylogeny": "gene_clusters",
    "anvi_get_sequences_for_gene_clusters": {
        "--min-num-genomes-gene-cluster-occurs": 5,
        "--max-num-genes-from-each-genome": 1,
        "--add-into-items-additional-data-table": "GCs_for_phylogeny",
        "--concatenate-gene-clusters": true,
        "--align-with": "famsa"
    },
    "import_phylogenetic_tree_to_pangenome": {
        "tree_name": "phylogeny_gene_clusters"
    }
}
```

{:.notice}
The parameters above are the ones that determine which gene clusters make it into your phylogeny, so they deserve a careful look. Asking for gene clusters that occur in all of your genomes (`--min-num-genomes-gene-cluster-occurs`, here set to `5` because this example has 5 genomes) with a single gene from each (`--max-num-genes-from-each-genome`) is a common way to approximate single-copy core genes.

And here is the equivalent for a phylogeny from HMM hits:

```json
{
    "sequence_source_for_phylogeny": "hmm",
    "anvi_get_sequences_for_hmm_hits": {
        "--hmm-sources": "Bacteria_71",
        "--return-best-hit": true,
        "--concatenate-genes": true,
        "--get-aa-sequences": true,
        "--min-num-bins-gene-occurs": 5,
        "--align-with": "muscle"
    },
    "import_phylogenetic_tree_to_pangenome": {
        "tree_name": "phylogeny_hmms"
    }
}
```

If you use `hmm` as your sequence source, make sure the `anvi_run_hmms` rule is set to `run=true` so the HMM hits your tree needs actually exist in your contigs databases.

## Computing genome similarity

The `anvi_compute_genome_similarity` rule is **not** run by default. Setting `run=true` for it will run %(anvi-compute-genome-similarity)s on your genomes and store the resulting %(genome-similarity)s data in your %(pan-db)s, so that similarity metrics show up as additional layers in the pangenome display. Any extra parameters for the program go into `additional_params`:

```json
{
    "anvi_compute_genome_similarity": {
        "run": true,
        "threads": 5,
        "additional_params": "--program pyANI --method ANIm"
    }
}
```

You can choose between `pyANI` and `fastANI` for ANI (the former is more accurate, the latter is much faster), or `sourmash` to compute mash distances instead.

{:.notice}
Even though `sourmash` is one of the options here, we don't recommend using it for genome comparisons -- it excels at other tasks -- and it remains available only as a legacy option. Please see %(anvi-compute-genome-similarity)s for the parameters each of these programs accepts.

## Output structure

With the minimal config above and a `project_name` of `TEST`, the workflow produces something like this:

```text
01_FASTA/                          # reformatted FASTA files (if anvi_script_reformat_fasta ran)
02_CONTIGS/                        # one contigs-db per genome
03_PAN/
├── TEST-GENOMES.db                # the genomes storage
├── TEST-PAN.db                    # the pangenome
└── TEST-ANI-OUTPUT/               # genome similarity results (if anvi_compute_genome_similarity ran)
00_LOGS/
└── pangenomics/
```

If you also asked for a phylogeny, the intermediate files of the tree computation land in the `PHYLO_DIR` directory (`01_PHYLOGENOMICS` by default):

```text
01_PHYLOGENOMICS/
├── TEST-GC-sequences.fa               # or TEST-proteins.fa, when the source is `hmm`
├── TEST-proteins_GAPS_REMOVED.fa      # the trimal output
└── TEST-proteins_GAPS_REMOVED.fa.contree   # the iqtree output, imported into the pan-db
```

You will also see a few empty `.done` files in `03_PAN`. These are flag files the workflow uses to keep track of steps that do not produce a single predictable output file (such as importing the tree into the %(pan-db)s), and you can safely ignore them.

All of these directory names are set in the `output_dirs` section of your config file, so feel free to change them.

Workflow logs will be under `00_LOGS/pangenomics` by default. Rule logs are organized by rule name, and the workflow writes a tab-delimited manifest at `00_LOGS/pangenomics/pangenomics-workflow-manifest.tsv` that points to each job's log and records whether it succeeded or failed.
