
This workflow is extremely useful if you have one or more %(fasta)s files that describe one or more contig sequences for your genomes or assembled metagenomes, and all you want to turn them into %(contigs-db)s files.

{:.warning}
If you have not yet run anvi'o programs %(anvi-setup-ncbi-cogs)s and %(anvi-setup-scg-taxonomy)s on your system yet, you will get a cryptic error from this workflow if you run it with the default %(workflow-config)s. You can avoid this by first running these two anvi'o programs to setup the necessary databases (which is done only once for every anvi'o installation), **or** set the rules for COG functions and/or SCG taxonomy to `run=false` explicitly.

To start things going with this workflow, first ask anvi'o to give you a default %(workflow-config)s file for the contigs workflow:

```bash
anvi-run-workflow -w contigs \
                  --get-default-config config-contigs-default.json
```

This will generate a file in your work directory called `config-contigs-default.json`. You should investigate its contents, and familiarize youself with it. It should look something like this, but much longer:
and you could examine its content to find out all possible options to tweak. We included a much simpler config file, `config-contigs.json`, in the mock data package for the sake of demonstrating how the contigs workflow works:

```json
{
    "workflow_name": "contigs",
    "config_version": "2",
    "fasta_txt": "fasta.txt",
    "output_dirs": {
        "FASTA_DIR": "01_FASTA",
        "CONTIGS_DIR": "02_CONTIGS",
        "LOGS_DIR": "00_LOGS"
    }
}
```

The only mandatory thing you need to do is to (1) manually create a %(fasta-txt)s file to describe the name and location of each FASTA file you wish to work with, and (2) make sure the `fasta_txt` variable in your %(workflow-config)s point to the location of your %(fasta-txt)s.

To see if everything looks alright, you can simply run the following command, which should generate a 'workflow graph' for you, given your config file parameters and input files:

```bash
anvi-run-workflow -w contigs \
                  -c config-contigs.json \
                  --save-workflow-graph
```

For the example config file shown above, this command will generate something similar to this:

[![DAG-contigs](../../images/workflows/contigs/DAG-contigs.png)]( ../../images/workflows/contigs/DAG-contigs.png){:.center-img .width-50}

{:.notice}
Please note that the generation of this workflow graph requires the usage of a program called [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)). If you are using MAC OSX, you can use [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) by installing [graphviz](http://www.graphviz.org/) through `brew` or `conda`.

If everything looks alright, you can run this workflow the following way:

```bash
anvi-run-workflow -w contigs \
                  -c config-contigs.json
```

If everything goes smoothly, you should see happy messages flowing on your screen, and at the end of it all you should see your contigs databases are generated and annotated properly. At the end of this process, you will have all your %(contigs-db)s files in the `02_CONTIGS` directory (as per the instructions in the config file, which you can change). You can use the program %(anvi-display-contigs-stats)s on one of them to see if everything makes sense.

Workflow logs will be under `00_LOGS/contigs` by default. Rule logs are organized by rule name, and the contigs workflow writes a tab-delimited manifest at `00_LOGS/contigs/contigs-workflow-manifest.tsv` that points to each job's log and records whether it succeeded or failed.

## Predicting protein structures

The contigs workflow can also predict the structures of the genes in each %(contigs-db)s, and store them in a %(structure-db)s using %(anvi-gen-structure-database)s. This step is off by default. To turn it on, set `run` to `true` for the rule `anvi_gen_structure_database` in your %(workflow-config)s. Every parameter of the prediction (including the `--engine`, which can be `modeller` or `colabfold`) lives under this rule:

```json
"anvi_gen_structure_database": {
    "run": true,
    "threads": 4,
    "--engine": "colabfold",
    "--colabfold-conda-env": "colabfold",
    "--colabfold-db": "/path/to/colabfold_db"
}
```

The workflow writes one %(structure-db)s per FASTA file into `03_STRUCTURE`, e.g., `03_STRUCTURE/SAMPLE_01-STRUCTURE.db`. The prediction does not hold up the rest of the workflow: it only needs the %(contigs-db)s, and annotating the %(contigs-db)s afterwards will not trigger a new prediction.

### Choosing genes

By default, anvi'o predicts a structure for every gene of every %(contigs-db)s, and warns you when it is about to do so. A single genome has thousands of genes, so this can take a very long time. To limit the prediction to some genes, add a `structure_genes_of_interest` column to your %(fasta-txt)s that points to a %(genes-of-interest-txt)s for each FASTA file:

|name|path|structure_genes_of_interest|
|:--|:--|:--|
|SAMPLE_01|path/to/sample_01.fa|genes_01.txt|
|SAMPLE_02|path/to/sample_02.fa||

Here, anvi'o predicts structures for the genes listed in `genes_01.txt` for `SAMPLE_01`, and for every gene of `SAMPLE_02`.

Since gene caller ids only exist once a %(contigs-db)s is built, you can first run the workflow without structure prediction, pick your genes, add the column to your %(fasta-txt)s, turn on `anvi_gen_structure_database`, and then run the workflow again: anvi'o will only run the missing structure steps. Changing a %(genes-of-interest-txt)s later will re-predict the structures of that %(contigs-db)s.

### Splitting ColabFold into CPU and GPU steps

ColabFold generates multiple sequence alignments (MSA) on the CPU, then predicts structures on the GPU. With a local ColabFold database (`--colabfold-db`), you can set `split_msa_and_predict` to `true` for `anvi_gen_structure_database` to run these as two separate jobs: `colabfold_msa` (with `--only-msa`) and `colabfold_predict` (with `--only-predict`). This way, a cluster can run each one on the appropriate node. Each of these two rules has its own `threads`, and takes every other parameter from `anvi_gen_structure_database`. The MSAs go into a temporary `03_STRUCTURE/SAMPLE_01-COLABFOLD-MSA` directory that the workflow removes once the %(structure-db)s is built.

### GPU resources

The rules that predict structures with ColabFold (`colabfold_predict`, and `anvi_gen_structure_database` when `--engine` is `colabfold`) declare a GPU through the snakemake resource `gpu=1`. MODELLER and `colabfold_msa` declare `gpu=0`. You can change this number with the `gpu` parameter of `anvi_gen_structure_database` and `colabfold_predict`.

It is up to you to map this resource onto your scheduler, e.g., through your snakemake cluster profile. On a single machine with one GPU, you can make sure that only one prediction runs at a time with `--additional-params --resources gpu=1 nodes=N`, where `N` is the total number of threads for your workflow (anvi'o sets `nodes` itself only when you do not pass your own `--resources`).
