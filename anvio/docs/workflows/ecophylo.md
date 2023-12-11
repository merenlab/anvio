The ecophylo workflow starts with a user-defined target gene family defined by an [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) and a list of assembled genomes and/or metagenomes. The final output is an %(interactive)s interface that includes (1) a phylogenetic analysis of all genes detected by the HMM in genomes and/or metagenomes, and (2) the distribution pattern of each of these genes across metagenomes if the user provided metagenomic short reads to survey.

While the 'user-defined [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms)' is passed to ecophylo via the %(hmm-list)s artifact, the input assemblies of genomes and/or metagenomes to query using the [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) are passed to the workflow via the artifacts %(external-genomes)s and %(metagenomes)s, respectively. Finally, the user can also provide a set of metagenomic short reads via the artifact %(samples-txt)s to recover the distribution patterns of genes across samples.

Ecophylo first identifies homologous genes based on the input [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms), clusters matching sequences based on a user-defined sequence similarity threshold, and finally selects a representative sequence from each cluster that contains more than two genes. The final set of representative genes are filtered for QC at multiple steps of the workflow which is discussed later in this document in the section "[Quality control and processing of hmm-hits](#Quality control and processing of hmm-hits)". After this step, the ecophylo workflow can continue with one of two modes that the user defines in the %(workflow-config)s: The so-called [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) or the so-called [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments).

In the [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes), the user must provide an %(hmm-list)s and %(metagenomes)s and/or %(external-genomes)s, and the workflow will stop after extracting representative sequences and calculating a phylogenetic tree (without any insights into the ecology of sequences through a subsequent step of metagenomic [read recruitment](https://anvio.org/vocabulary/#read-recruitment)). In contrast, the [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments) will require an additional file: %(samples-txt)s. In this mode the workflow will continue with the profiling of representative sequences via read recruitment across user-provided metagenomes to recover and store coverage statistics. The completion of the workflow will yield all files necessary to explore the results in downstream analyses to investigate associations between ecological and evolutionary relationships between target genes.

The ecophylo workflow can leverage any [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) that models amino acid sequences. If the user chooses an [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) for a [single-copy core gene](https://anvio.org/vocabulary/#single-copy-core-gene-scg), such as ribosomal protein, the workflow will yield multi-domain taxonomic profiles of metagenomes *de facto*.

## Required input

The minimum requirements of the ecophylo workflow are the following:

- %(workflow-config)s: This allows you to customize the workflow step by step. Here is how you can generate the default version:

{{ codestart }}
anvi-run-workflow -w ecophylo \
                  --get-default-config config.json
{{ codestop }}

- %(hmm-list)s: This file designates which HMM should be used to extract the target gene from your %(contigs-db)s. Please note that the ecophylo workflow can only process one gene family at a time i.e. %(hmm-list)s can only contain one HMM. If you would like to process multiple gene families from the same input assemblies then you will need to re-run the workflow with a separate %(hmm-list)s.
- %(metagenomes)s and/or %(external-genomes)s: These files hold the assemblies where you are looking for the target gene. Genomes in %(external-genomes)s can be reference genomes, [SAGs](https://anvio.org/vocabulary/#single-amplified-genome-sag), and/or [MAGs](https://anvio.org/vocabulary/#metagenome-assembled-genome-mag).

## Quality control and processing of hmm-hits

[Hidden Markov Models](https://anvio.org/vocabulary/#hidden-markov-models-hmms) are the crux of the ecophylo workflow and will determine the sensitivity and specificity of the gene family hmm-hits you seek to investigate. However, not all %(hmm-hits)s are created equal. Just how BLAST can detect spurious hits with [high-scoring segment pairs](https://www.ncbi.nlm.nih.gov/books/NBK62051/), an HMM search can yield non-homologous hits as well. To address this, we have a series of parameters you can adjust in the %(workflow-config)s to fine tune the input set of %(hmm-hits)s that ecophylo will process.

### HMM alignment coverage filtering

The first step to removing bad %(hmm-hits)s is to filter out hits with low quality alignment coverage. This is done with the rule `filter_hmm_hits_by_model_coverage` which leverages %(anvi-script-filter-hmm-hits-table)s. We recommend 80%% model coverage filter for most cases. However, it is always recommended to explore the distribution of model coverage with any new HMM which will help you determine a proper cutoff (citation). To adjust this parameter, go to the `filter_hmm_hits_by_model_coverage` rule and change the parameter `--min-model-coverage`.

{:.notice}
Some full gene length HMM models align to a single hmm-hit independently at different coordinates when there should only be one annotation. To merge these independent alignment into one HMM alignment coverage stat, set `--merge-partial-hits-within-X-nts` to any distance between the hits for which you would like to merge and add it to the rule `filter_hmm_hits_by_model_coverage` under `additional_params`.

### conservative-mode: complete open-reading frames only

Genes predicted from genomes and metagenomes can be partial or complete depending on whether a stop and stop codon is detected. Even if you filter out %(hmm-hits)s with bad alignment coverage as discussed above, HMMs can still detect low quality hits with good alignment coverage and homology statistics due to partial genes. Unfortunately, partial genes can lead to spurious phylogenetic branches and/or inflate the number of observed populations or functions in a given set of genomes/metagenomes.

To remove partial genes from the ecophylo analysis, the user can assign `true` for `--filter-out-partial-gene-calls` parameter so that only complete open-reading frames are processed.

{:.notice}
What is below is the default settings in the ecophylo %(workflow-config)s file.

```bash
{
    "filter_hmm_hits_by_model_coverage": {
        "threads": 5,
        "--min-model-coverage": 0.8,
        "--filter-out-partial-gene-calls": true,
        "additional_params": ""
    },
}
```

### discovery-mode: ALL open-reading frames

However, maybe you're a risk taker, a maverick explorer of metagenomes. Complete or partial you accept all genes and their potential tree bending shortcomings! In this case, set `--filter-out-partial-gene-calls false` in the %(workflow-config)s.

{:.notice}
Simultaneously exploring complete and partial ORFs will increase the distribution of sequence lengths and thus impact sequence clustering. We recommend adjusting `cluster_X_percent_sim_mmseqs` to `"--cov-mode": 1` to help insure ORFs of all length properly cluster together. Please refer to the [MMseqs2 user guide description of --cov-mode](https://mmseqs.com/latest/userguide.pdf) for more details.

```bash
{
    "filter_hmm_hits_by_model_coverage": {
        "threads": 5,
        "--min-model-coverage": 0.8,
        "--filter-out-partial-gene-calls": false,
        "additional_params": ""
    },
      "cluster_X_percent_sim_mmseqs": {
      "threads": 5,
      "--min-seq-id": 0.94,
      "--cov-mode": 1,
      "clustering_threshold_for_OTUs": [
          0.99,
          0.98,
          0.97
      ],
      "AA_mode": false
    },
}
```

Now that you have fine tuned the gene family input into the ecophylo workflow, it's time to decide what output best fits your science question at hand.

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

## Miscellaneous config file options

Ecophylo will sanity check all input files that contain %(contigs-db)ss before the workflow starts. This can take a while especially if you are working with 1000's of genomes. If you want to skip sanity checks for %(contigs-db)ss in your %(external-genomes)s and/or %(metagenomes)s then adjust your %(workflow-config)s to the following:

```bash
{
    "run_genomes_sanity_check": false
}
```

The ecophylo workflow by default uses [FastTree](http://www.microbesonline.org/fasttree/) to calculate the output phylogenetic tree. This is because the workflow was designed to be run on large genomic datasets that could yield thousands of input sequences. However, if you like to run [IQ-TREE](https://github.com/Cibiv/IQ-TREE) adjust your %(workflow-config)s to the following:

```bash
{
    "fasttree": {
        "run": "",
        "threads": 5
    },
    "iqtree": {
        "threads": 5,
        "-m": "MFP",
        "run": true,
        "additional_params": ""
    },
}
```