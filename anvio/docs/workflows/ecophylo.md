The ecophylo workflow starts with a user-defined target gene ([HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms)) and a list of assembled genomes and/or metagenomes and results in an %(interactive)s interface that includes (1) a phylogenetic analysis of all genes found in genomes and metagenomes that match to the user-defined target gene, and (2) the distribution pattern of each of these genes across metagenomes if the user provided metagenomic short reads to survey.

The user-defined target genes can be described by an %(hmm-list)s. Furthermore, the assemblies of genomes and/or metagenomes to search these genes can be passed to the workflow via the artifacts %(external-genomes)s and %(metagenomes)s, respectively. Finally, the user can also provide a set of metagenomic short reads via the artifact %(samples-txt)s to recover the distribution patterns of genes.

In a standard run, ecophylo first identifies matching genes based on their [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms)s, then clusters them based on sequence similarity at a threshold defined by the user, and finally selects a representative sequence from each cluster that contains more than two genes. Next, ecophylo calculates a phylogenetic tree to infer evolutionary associations between these sequences to produce a NEWICK-formatted %(dendrogram)s. If the user provided a %(samples-txt)s for metagenomic [read recruitment](https://anvio.org/vocabulary/#read-recruitment), the workflow will also perform a [read recruitment](https://anvio.org/vocabulary/#read-recruitment) step to recover and store coverage statistics of the final set of genes for downstream analyses in the form of a %(profile-db)s. The completion of the workflow will yield all files necessary to explore the results through an anvi'o %(interactive)s interface and investigate associations between ecological and evolutionary relationships between target genes. The workflow can use any [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) that models amino acid sequences. Using [single-copy core genes](https://anvio.org/vocabulary/#single-copy-core-gene-scg) such as Ribosomal Proteins will yield taxonomic profiles of metagenomes *de facto*.

The ecophylo workflow has 2 modes which can be designated in the %(workflow-config)s by changing the input files that are provided: [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) and [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments). In [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes), the sequences will be used to calculate a phylogenetic tree. In [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments), the sequences will be used to calculate a phylogenetic tree and be additionally profiled via [read recruitment](https://anvio.org/vocabulary/#read-recruitment) across user-provided metagenomes.

## Required input

The ecophylo workflow requires the following files:

- %(workflow-config)s: This allows you to customize the workflow step by step. Here is how you can generate the default version:

{{ codestart }}
anvi-run-workflow -w ecophylo \
                  --get-default-config config.json
{{ codestop }}


{:.notice}
Here is a tutorial walking through more details regarding the ecophylo %(workflow-config)s file: coming soon!

- %(hmm-list)s: This file designates which HMM should be used to extract the target gene from your %(contigs-db)s.
- %(metagenomes)s and/or %(external-genomes)s: These files hold the assemblies where you are looking for the target gene. Genomes in %(external-genomes)s can be reference genomes, [SAGs](https://anvio.org/vocabulary/#single-amplified-genome-sag), and/or [MAGs](https://anvio.org/vocabulary/#metagenome-assembled-genome-mag).

## tree-mode: Insights into the evolutionary patterns of target genes 

This is the simplest implementation of ecophylo where only an amino acid based phylogenetic tree is calculated. The workflow will extract the target gene from input assemblies, cluster and pick representatives, then calculate a phylogenetic tree based on the amino acid representative sequences. There are two sub-modes of [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) which depend on how you pick representative sequences, [NT-mode](#nt-mode) or [AA-mode](#aa-mode) where extracted genes associated nucleotide version (NT) or the amino acid (AA) can be used to cluster the dataset and pick representatives, respectively.

### NT-mode

**Cluster and select representative genes based on NT sequences.**

This is the default version of [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) where the extracted gene sequences are clustered based on their associated NT sequences. This is done to prepare for [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments),  where adequate sequence distance is needed between gene NT sequences to prevent [non-specific-read-recruitment](https://anvio.org/vocabulary/#non-specific-read-recruitment). The translated amino acid versions of the NT sequence clusters are then used to calculate an AA based phylogenetic tree. This mode is specifically useful to see what the gene phylogenetic tree will look like before the [read recruitment](https://anvio.org/vocabulary/#read-recruitment) step in [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments),  (for gene phylogenetic applications of ecophylo please see [AA-mode](#Cluster based on AA sequences - AA-mode)). If everything looks good you can add in your %(samples-txt)s and continue with [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments) to add metagenomic [read recruitment](https://anvio.org/vocabulary/#read-recruitment) results.

Here is what the start of the ecophylo %(workflow-config)s should look like if you want to run [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes):

```bash
{
    "metagenomes": "metagenomes.txt",
    "external_genomes": "external-genomes.txt",
    "hmm_list": "hmm_list.txt",
    "samples_txt": ""
}
```

### AA-mode

**Cluster and select representative genes based on AA sequences. If you are interested specifically in gene phylogenetics, this is the mode for you!**

This is another sub-version of [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) where representative sequences are chosen via AA sequence clustering.

To initialize [AA-mode](#aa-mode), go to the rule `cluster_X_percent_sim_mmseqs` in the ecophylo %(workflow-config)s and turn "AA_mode" to true:

```bash
{
    "metagenomes": "metagenomes.txt",
    "external_genomes": "external-genomes.txt",
    "hmm_list": "hmm_list.txt",
    "samples_txt": ""
    "cluster_X_percent_sim_mmseqs": {
        "AA_mode": true,
    }
}
```

{:.notice}
Be sure to change the `--min-seq-id` of the `cluster_X_percent_sim_mmseqs` rule to the appropriate clustering threshold depending if you are in [NT-mode](#nt-mode) or [AA-mode](#aa-mode).

## profile-mode: Insights into the ecological and evolutionary patterns of target genes and environments

[profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments),  is an extension of default [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) ([NT-mode](#nt-mode)) where NT sequences representatives are profiled with metagenomic reads from user provided metagenomic samples. This allows for the simultaneous visualization of phylogenetic and ecological relationships of genes across metagenomic datasets.

Additional required files:
- %(samples-txt)s

To initialize [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments), , add the path to your %(samples-txt)s to your ecophylo %(workflow-config)s:

```bash
{
    "metagenomes": "metagenomes.txt",
    "external_genomes": "external-genomes.txt",
    "hmm_list": "hmm_list.txt",
    "samples_txt": "samples.txt"
}
```

## Config file options

Ecophylo will sanity check all input files that contain %(contigs-db)ss before the workflow starts. This can take a while especially if you are working with 1000's of genomes. If you want to skip sanity checks for %(contigs-db)ss in your %(external-genomes)s and/or %(metagenomes)s then adjust your config to the following:

```bash
{
    "run_genomes_sanity_check": false
}
```