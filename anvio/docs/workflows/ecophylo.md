The EcoPhylo workflow is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow run by the anvi'o script %(anvi-run-workflow)s. Briefly, the workflow extracts a target protein from any number of assemblies in %(metagenomes)s and/or %(external-genomes)s using a user-designated [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) from %(hmm-list)s. Next, the workflow clusters the sequences and selects representatives from each cluster. After that, the workflow calculates a phylogenetics tree and performs metagenomic read recruitment (if %(samples-txt)s is provided) to produce an anvi'o %(interactive)s interface to explore the proteins phylogenetics and co-occurrence across metagenomes (via simultaneous visualization of a phylogenetic tree and read recruitment results). The workflow can use any protein-based [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) including [single-copy core genes](https://anvio.org/vocabulary/#single-copy-core-gene-scg) to taxonomically profile metagenomes or any functional protein to explore variants across samples. 

The EcoPhylo workflow has 2 modes which can be designated in the %(workflow-config)s by changing the input files that are provided: `tree-mode` and `profile-mode`. In `tree-mode`, the sequences will be used to calculate a phylogenetic tree. In `profile-mode`, the sequences will be used to calculate a phylogenetic tree and be additionally profiled via read recruitment across user-provided metagenomes. 

## Required input

The EcoPhylo workflow requires the following files:

- %(workflow-config)s: This allows you to customize the workflow step by step. Here is how you can generate the default version:

{{ codestart }}
anvi-run-workflow -w ecophylo --get-default-config config.json
{{ codestop }}

{:.notice}
Here is a tutorial walking through more details regarding the EcoPhylo %(workflow-config)s file: coming soon!

- %(hmm-list)s: This file designates which HMM should be used to extract the target protein from your %(contigs-db)s.  
- %(metagenomes)s and/or %(external-genomes)s: These files hold the assemblies where you are looking for the target protein. Genomes in %(external-genomes)s can be reference genomes, [SAGs](https://anvio.org/vocabulary/#single-amplified-genome-sag), and/or [MAGs](https://anvio.org/vocabulary/#metagenome-assembled-genome-mag). 

## Want to explore phylogenetic relationships of proteins across assemblies? Tree-mode

This is the simplest implementation of EcoPhylo where only an amino acid based phylognetic tree is calculated. The workflow will extract the target protein from input assemblies, cluster and pick representatives, then calculate a phylogenetic tree based on the amino acid representative sequences. There are two sub-modes of `tree-mode` which depend on how you pick representative sequences, `NT-mode` or `AA-mode` where extracted proteins associated nucleotide version (NT) or the amino acid (AA) can be used to cluster the dataset and pick representatives, respectively.

### NT-mode

**Cluster and select representative proteins based on NT sequences.**

This is the default version of `tree-mode` where the extracted protein sequences are clustered based on their associated NT sequences. This is done to prepare for `profile-mode` where adequate NT sequence distance is needed between protein NT sequences to prevent [non-specific-read-recruitment](https://anvio.org/vocabulary/#non-specific-read-recruitment). The translated amino acid versions of the NT sequence clusters are then used to calculate an AA based phylogenetic tree. This mode is specifically useful to see what the protein phylogenetic tree will look like before the read recruitment step in `profile-mode` (for protein phylogenetic applications of EcoPhylo please see [AA-mode](#Cluster based on AA sequences - AA-mode)). If everything looks good, you can add in your %(samples-txt) and continue with `profile-mode` to add metagenomic read recruitment results. 

Here is what the start of the EcoPhylo %(workflow-config)s should look like if you want to run `tree-mode`:

```bash
{
    "metagenomes": "metagenomes.txt",
    "external_genomes": "external-genomes.txt",
    "hmm_list": "hmm_list.txt",
    "samples_txt": ""
}
```

### AA-mode

**Cluster and select representative proteins based on AA sequences. If you are interested specifically in protein phylogenetics, this is the mode for you!**

This is another sub-version of `tree-mode` where representative sequences are chosen via AA sequence clustering.

To initialize `AA-mode`, go to the rule `cluster_X_percent_sim_mmseqs` in the EcoPhylo %(workflow-config)s and turn "AA_mode" to true:

```bash
{
    "metagenomes": "metagenomes.txt",
    "external_genomes": "external-genomes.txt",
    "hmm_list": "hmm_list.txt",
    "samples_txt": ""
    "cluster_X_percent_sim_mmseqs": {
        "AA_mode": True,
    }
}
```

{:.notice}
Be sure to change the `--min-seq-id` of the `cluster_X_percent_sim_mmseqs` rule to the appropriate clustering threshold depending if you are in `NT-mode` or `AA-mode`.

## Want to track proteins across metagenomic samples via read recruitment? Profile-mode

`profile-mode` is an extension of default `tree-mode` where NT sequences representatives are profiled with metagenomic reads from user provided metagenomic samples. This allows for the simultaneous visualization of phylogenetic and ecological relationships of proteins across metagenomic datasets. 

Additional required files:
- %(samples-txt)s

To initialize `profile-mode`, add the path to your %(samples-txt)s to your EcoPhylo %(workflow-config)s:

```bash
{
    "metagenomes": "metagenomes.txt",
    "external_genomes": "external-genomes.txt",
    "hmm_list": "hmm_list.txt",
    "samples_txt": "samples.txt"
}
```