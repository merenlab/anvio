The ecophylo workflow starts with a user-provided target gene family defined by an [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) and a list of assembled genomes and/or metagenomes. The final output is an %(interactive)s interface that includes (1) a phylogenetic analysis of all genes detected by the HMM in genomes and/or metagenomes, and (2) the distribution pattern of each of these genes across metagenomes if the user provided metagenomic short reads to survey.

The 'user-provided [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms)' is passed to ecophylo via the %(hmm-list)s file, and the input assemblies of genomes and/or metagenomes to query using the [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) are passed to the workflow via the files %(external-genomes)s and %(metagenomes)s, respectively. Finally, the user can also provide a set of metagenomic short reads via a %(samples-txt)s to recover the distribution patterns of genes across samples.

Ecophylo first identifies homologous genes based on the input [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms), clusters matching sequences based on a user-defined sequence similarity threshold, and finally selects a representative sequence from each cluster that contains more than two genes. The final set of representative genes are filtered for QC at multiple steps of the workflow which is discussed later in this document in the section "[Quality control and processing of hmm-hits](#Quality control and processing of hmm-hits)". After this step, the ecophylo workflow can continue with one of two modes that the user defines in the %(workflow-config)s: The so-called [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes) or the so-called [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments).

In the [tree-mode](#tree-mode-insights-into-the-evolutionary-patterns-of-target-genes), the user must provide an %(hmm-list)s and %(metagenomes)s and/or %(external-genomes)s, and the workflow will stop after extracting representative sequences and calculating a phylogenetic tree (without any insights into the ecology of sequences through a subsequent step of metagenomic [read recruitment](https://anvio.org/vocabulary/#read-recruitment)). In contrast, the [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments) will require an additional file: %(samples-txt)s. In this mode the workflow will continue with the profiling of representative sequences via read recruitment across user-provided metagenomes to recover and store coverage statistics. The completion of the workflow will yield all files necessary to explore the results in downstream analyses to investigate associations between ecological and evolutionary relationships between target genes.

The ecophylo workflow can leverage any [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) that models amino acid sequences. If the user chooses an [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms) for a [single-copy core gene](https://anvio.org/vocabulary/#single-copy-core-gene-scg), such as ribosomal protein, the workflow will yield multi-domain taxonomic profiles of metagenomes *de facto*.

{:.notice}
If you have never run an anvi'o snakemake workflow, please checkout the [anvi'o snakemake workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/). This is where you can learn the basics about how anvi'o leverages Snakemake to process data. In fact, the EcoPhylo workflow uses the anvi'o metagenomics workflow to profile protein families across metagenomes. Two birds, two workflows?

## Required input

The minimum requirements of the ecophylo workflow are the following:

- %(workflow-config)s: This allows you to customize the workflow step by step. Here is how you can generate the default version:

{{ codestart }}
anvi-run-workflow -w ecophylo \
                  --get-default-config config.json
{{ codestop }}

- %(hmm-list)s: This file designates which HMM should be used to extract the target gene from your %(contigs-db)s. Please note that the ecophylo workflow can only process one gene family at a time i.e. %(hmm-list)s can only contain one HMM. If you would like to process multiple gene families from the same input assemblies then you will need to re-run the workflow with a separate %(hmm-list)s.
- %(metagenomes)s and/or %(external-genomes)s: These files hold the assemblies where you are looking for the target gene. Genomes in %(external-genomes)s can be reference genomes, [SAGs](https://anvio.org/vocabulary/#single-amplified-genome-sag), and/or [MAGs](https://anvio.org/vocabulary/#metagenome-assembled-genome-mag).

## A quick tour of the output directory structure

The ecophylo workflow produces a ton of intermediate files that can be useful for you to explore your data! Here is a basic look at the directory structure after successfully running the workflow:

```bash
$ tree ECOPHYLO_WORKFLOW -L 1
├── 00_LOGS
├── 01_REFERENCE_PROTEIN_DATA
├── 02_NR_FASTAS
├── 03_MSA
├── 04_SEQUENCE_STATS
├── 05_TREES
├── 06_MISC_DATA
├── METAGENOMICS_WORKFLOW
```

Let's dive into some key intermediate files!

`01_REFERENCE_PROTEIN_DATA/`

This directory contains data extracted from each individual %(contigs-db)s the user provides via the %(external-genomes)s and/or %(metagenomes)s files. (The wide part of the funnel)

- `ASSEMBLY-PROTEIN-external_gene_calls.tsv`: %(external-gene-calls)s for each open reading frame in the analyzed %(contigs-db)s
- `ASSEMBLY-PROTEIN-external_gene_calls_renamed.tsv`: %(external-gene-calls)s renamed and subsetted for target protein
- `ASSEMBLY-PROTEIN-hmm_hits.faa` and `ASSEMBLY-PROTEIN-hmm_hits.fna`: target protein sequences extracted with %(anvi-get-sequences-for-hmm-hits)s
- `ASSEMBLY-PROTEIN-hmm_hits_renamed.faa` and `ASSEMBLY-PROTEIN-hmm_hits_renamed.fna`: fasta files with renamed headers
- `ASSEMBLY-PROTEIN-orfs.fna`: output fasta from %(anvi-get-sequences-for-gene-calls)s
- `ASSEMBLY-PROTEIN-reformat_report_AA.txt` and `ASSEMBLY-PROTEIN-reformat_report_nt.txt`
- `ASSEMBLY-dom-hmmsearch/`: this directory contains all the homologs extracted from your input %(contigs-db)s (genomes and/or metagenomic assembles) with the user-provided [HMM](https://anvio.org/vocabulary/#hidden-markov-models-hmms). Here are some key files:
  - `hmm.domtable` contains the raw output the domain hits table from `hmmsearch` run by %(anvi-run-hmms)s
  - `hmm_hits.txt` is the hmm-hits table stored in the associated %(contigs-db)s
  - `hmm_hits_filtered.txt` is the filtered version `hmm_hits.txt` based on user-defined parameters


`02_NR_FASTAS/`

This directory is where all the data from individual %(contigs-db)s is combined and contains all the protein family clustering information from the workflow. (The narrow part of the funnel)

- `PROTEIN-all.faa` and `PROTEIN-all.fna` contains ALL the amino acid and nucleotide sequences that made it past the initial filtering steps and will be clustered
- `PROTEIN-mmseqs_NR_cluster.tsv`: mmseqs cluster output file. VERY helpful for looking inside clusters (column 1: representative sequence; column 2: cluster member)
- `PROTEIN-references_for_mapping_NT.fa`: Input ORFs used for the metagenomics workflow
- `PROTEIN-AA_subset.fa`: translated sequences from `PROTEIN-references_for_mapping_NT.fa`
- `PROTEIN-external_gene_calls_all.tsv`: external-gene calls file for `PROTEIN-references_for_mapping_NT.fa`

`03_MSA/`

This directory contains all the intermediate files from multiple sequences alignment steps.

- `PROTEIN-aligned.fa`: raw output from MSA
- `PROTEIN_aligned_trimmed.fa`: trimmed MSA from `trimal`
- `PROTEIN_aligned_trimmed_filtered.fa`: subsetted MSA removing sequences with x > 50%% gaps
- `PROTEIN_gaps_counts.tsv`: Table counting number of gaps per sequence in MSA

`04_SEQUENCE_STATS/`

This directory contains information regarding the number of sequences filtered at various steps of the workflow.

`PROTEIN_stats.tsv`: tracks number of sequences filtered at different steps of the workflow. Here are definitions for each rule:
    - `combine_sequence_data`: total number of sequences before clustering
    - `cluster_X_percent_sim_mmseqs`: number of cluster representative sequences
    - `remove_sequences_with_X_percent_gaps`: number of sequences left after filtering during MSA
    - `99_percent`: number of clusters after clustering at 99%% nucleotide similarity
    - `98_percent`: number of clusters after clustering at 98%% nucleotide similarity

`05_TREES/`

Here we have all things phylogenetics in the workflow.

- `PROTEIN-PROFILE.db`:
- `PROTEIN.nwk`: PROTEIN phylogenetic tree
- `PROTEIN_renamed.faa`:renamed PROTEIN phylogenetic tree fasta file
- `PROTEIN_renamed.nwk`: renamed PROTEIN phylogenetic tree to be imported in Metagenomics workflow merged profile
- `PROTEIN_renamed_all.faa`: renamed fasta file with ALL proteins

`06_MISC_DATA/`

This directory contains miscellaneous data created from the flow to help you interpret the phylogeography of your target PROTEIN.

- `PROTEIN_misc.tsv`: This file contains basic information about each representative sequence in the workflow including: 
  - contigs_db_type: genome or metagenomic assembly
  - genomic_seq_in_cluster: YES/NO a sequence that originate from an input genome is in the cluster
  - cluster_size: number of sequences in mmseqs cluster

These files are the output of %(anvi-estimate-scg-taxonomy)s and will only be there if you are explore the phylogeography of a compatible single-copy core gene
- `PROTEIN_scg_taxonomy_data.tsv`: 
- `PROTEIN_estimate_scg_taxonomy_results-RAW-LONG-FORMAT.txt`: 

`METAGENOMICS_WORKFLOW/`

This directory contains the output of the EcoPhylo sequences profiled with metagenomes with the [anvi'o metagenomics workflow](https://anvio.org/help/main/workflows/metagenomics/).


## Quality control and processing of hmm-hits

[Hidden Markov Models](https://anvio.org/vocabulary/#hidden-markov-models-hmms) are the crux of the ecophylo workflow and will determine the sensitivity and specificity of the gene family hmm-hits you seek to investigate. However, not all %(hmm-hits)s are created equal. Just how BLAST can detect spurious hits with [high-scoring segment pairs](https://www.ncbi.nlm.nih.gov/books/NBK62051/), an HMM search can yield non-homologous hits as well. To address this, we have a series of parameters you can adjust in the %(workflow-config)s to fine tune the input set of %(hmm-hits)s that ecophylo will process.

### HMM alignment coverage filtering

The first step to removing bad %(hmm-hits)s is to filter out hits with low quality alignment coverage. This is done with the rule `filter_hmm_hits_by_model_coverage` which leverages %(anvi-script-filter-hmm-hits-table)s. This tool uses the output of hmmsearch to filter out hits basedon the model and/or gene coverage. We recommend 80%% model coverage filter for most cases. However, it is always recommended to explore the distribution of model coverage with any new HMM which will help you determine a proper cutoff (citation). To adjust this parameter, go to the `filter_hmm_hits_by_model_coverage` rule and change the parameter `--min-model-coverage`. You can also adjust the gene coverage by change the parameter `--min-gene-coverage`. This can help remove ORFs with outlier lengths but completely depends on the HMM you are using.

{:.notice}
Please consider exploring the distribution of alignment coverages before choosing HMM alignment coverage filtering values. [Interproscan](https://www.ebi.ac.uk/interpro/) is a great way to visualize how publicly available HMMs align to proteins. Additionally, you can parse the domtblout files from hmmsearch to explore these values in high throughput. 

```bash
{
    "filter_hmm_hits_by_model_coverage": {
        "threads": 5,
        "--min-model-coverage": 0.8,
        "--min-gene-coverage": 0.5,
        "additional_params": ""
    },
}
```

{:.notice}
Some full gene length HMM models align to a single hmm-hit independently at different coordinates when there should only be one annotation. To merge these independent alignment into one HMM alignment coverage stat, set `--merge-partial-hits-within-X-nts` to any distance between the hits for which you would like to merge and add it to the rule `filter_hmm_hits_by_model_coverage` under `additional_params`.

```bash
{
    "filter_hmm_hits_by_model_coverage": {
        "additional_params": "--merge-partial-hits-within-X-nts"
    },
}
```

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

### Multiple sequence alignment step with MUSCLE

One step of ecophylo is to perform a multiple sequence alignment of the recruited homologs and depending on your application, this could recruit thousands of ORFs which make the MSA a challenging feat. By default, the ecophylo is designed for quick insights, and thus the %(workflow-config)s file uses MUSCLE parameters to perform a large MSA, swiftly: 

```bash
"align_sequences": {
    "threads": 5,
    "additional_params": "-maxiters 1 -diags -sv -distance1 kbit20_3"
},
```

However, these parameters may not be optimal for your use case. For example, maybe you are trying to explore branches patterns of a specific protein family and would prefer to have mulitple interations of the MSA. Please explore the MUSCLE documentation to [documentation](https://www.drive5.com/muscle/muscle.html) customize the MSA step for your needs. You can replace the `additional_params` with whatever MUSCLE parameters that are best for you. 

### discovery-mode: ALL open-reading frames

However, maybe you're a risk taker, a maverick explorer of metagenomes. Complete or partial you accept all genes and their potential tree bending shortcomings! In this case, set `--filter-out-partial-gene-calls false` in the %(workflow-config)s.

{:.notice}
Simultaneously exploring complete and partial ORFs will increase the distribution of sequence lengths and thus impact sequence clustering. We recommend adjusting `cluster_X_percent_sim_mmseqs` to `"--cov-mode": 1` to help insure ORFs of all length properly cluster together. Please refer to the [MMseqs2 user guide description of --cov-mode](https://mmseqs.com/latest/userguide.pdf) for more details.

#FIXME: we ALWAYS recommend --cov-mode 1 to group protein fragments as well as overextended ORFs caused by early or late stop codons.

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

{:.notice}
It's common that not all genomes or metagenomes will have the gene family of interest either due to it not being detect by the input HMM or filtered out during the QC steps. Please check this log file for %(contigs-db)s that did not contain your gene family of interest: `00_LOGS/contigDBs_with_no_hmm_hit_*.log`

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

### Visualize the output

To visualize the output of [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments), run %(anvi-interactive)s on the %(contigs-db)s and %(profile-db)s located in the `METAGENOMICS_WORKFLOW` directory.

```bash
PROTEIN=""
anvi-interactive -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db \
                 -c METAGENOMICS_WORKFLOW/03_CONTIGS/"${PROTEIN}"-contigs.db \
                 --manual
```

Just want a quick look at the tree without read recruitment results?

```bash
PROTEIN=""
anvi-interactive -t 05_TREES/"${PROTEIN}"/"${PROTEIN}"_renamed.nwk \
                 -p 05_TREES/"${PROTEIN}"/"${PROTEIN}"-PROFILE.db \
                 --fasta 05_TREES/"${PROTEIN}"/"${PROTEIN}"_renamed.faa \
                 --manual 
```

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

### Visualize the output

To visualize the output of [profile-mode](#profile-mode-insights-into-the-ecological-and-evolutionary-patterns-of-target-genes-and-environments), run %(anvi-interactive)s on the %(contigs-db)s and %(profile-db)s located in the `METAGENOMICS_WORKFLOW` directory.

{{ codestart }}
PROTEIN=""
anvi-interactive -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db \
                 -c METAGENOMICS_WORKFLOW/03_CONTIGS/"${PROTEIN}"-contigs.db
{{ codestop }}

## Manual curation of the ecophylo phylogeny

Calculating a phylogeny from ORFs recruited from a metagenomic assembly can result in some unnatural long branches. This can be caused by a variety of reasons including a misassembled sequences that ecophylo couldn't removed automatically :(
    
To remove branches from the phylogenetic tree in the ecophylo interface, you can manually curate the tree and reimport it into anvi'o. At the moment, anvi'o does not have an automated program to do this but here is a workflow:

{:.notice}
This is just a code outline. Please adjust parameters for the various steps to match your specific needs.

**Step 1.** Make collection of bad branches

Make a collection of branches you would like to remove and safe it!

**Step 2.** Export collection and remove those sequences from the protein fasta file

```bash
HOME_DIR="ECOPHYLO"
PROTEIN=""
cd $HOME_DIR

mkdir SUBSET_TREE

# Export default collection and collection of bad branches
anvi-export-collection -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db -C DEFAULT --output-file-prefix SUBSET_TREE/DEFAULT
anvi-export-collection -C bad_branches -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db --output-file-prefix SUBSET_TREE/bad_branches

# Clean fasta headers
cut -f 1 SUBSET_TREE/bad_branches | sed 's|_split_00001||' > SUBSET_TREE/bad_branches_headers.txt

anvi-script-reformat-fasta 02_NR_FASTAS/"${PROTEIN}"/"${PROTEIN}"-AA_subset.fa --exclude-ids bad_branches_headers.txt \ 
                                                                               -o SUBSET_TREE/"${PROTEIN}"-AA_subset_remove_bad_branches.fa

ALIGNMENT_PREFIX=""${PROTEIN}"-AA_subset_remove_bad_branches"

# Align
clusterize "muscle -in SUBSET_TREE/"${PROTEIN}"-AA_subset_remove_bad_branches.fa -out SUBSET_TREE/"${ALIGNMENT_PREFIX}".faa -maxiters 1" -n 15 -o 00_LOGS/align.log

# Trim
clusterize "trimal -in SUBSET_TREE/"${ALIGNMENT_PREFIX}".faa -out SUBSET_TREE/"${ALIGNMENT_PREFIX}"_trimmed.faa -gappyout" -n 5 -o 00_LOGS/trim.log

# Calculate tree
clusterize "FastTree SUBSET_TREE/"${ALIGNMENT_PREFIX}"_trimmed_filtered.faa > SUBSET_TREE/"${ALIGNMENT_PREFIX}"_trimmed_filtered_FastTree.nwk" -n 15 -o 00_LOGS/FastTree.log
```

**Step 3.** Use `anvi-split` to remove bad branches from the ecophylo interface

```bash
PROTEIN=""

grep -v -f SUBSET_TREE/bad_branches_headers.txt SUBSET_TREE/collection-DEFAULT.txt | sed 's|EVERYTHING|EVERYTHING_curated|' > SUBSET_TREE/my_bins.txt

anvi-import-collection SUBSET_TREE/my_bins.txt -C curated \
                                               -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db \
                                               -c METAGENOMICS_WORKFLOW/03_CONTIGS/"${PROTEIN}"-contigs.db
                        
anvi-split -C curated \
           --bin-id EVERYTHING_curated \
           -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db \
           -c METAGENOMICS_WORKFLOW/03_CONTIGS/"${PROTEIN}"-contigs.db \
           --output-dir SUBSET_TREE/"${PROTEIN}"_curated
```

**Step 4.** Add the string "_split_00001" to each tree leaf to import it back into the interface

{{ codestart }}
packages <- c("tidyverse", "ape", "phytools", "glue")
suppressMessages(lapply(packages, library, character.only = TRUE))

add_split_string_to_tree <- function(IN_PATH, OUT_PATH) {
  
  # Import tree
  tree <- read.tree(IN_PATH)
  
  # Create DF with tree tip metadata
  tree_tip_metadata <- tree$tip.label %%>%% 
    as_tibble() %%>%%
    rename(tip_label = value) %%>%%
    mutate(tip_label = str_c(tip_label, "_split_00001"))
  
  
  tree$tip.label <- tree_tip_metadata$tip_label
  
  # Write DF 
  print(OUT_PATH)
  write.tree(tree, file = OUT_PATH)
}

PROTEIN=""
add_split_string_to_tree(IN_PATH = glue("{PROTEIN}_trimmed_filtered_FastTree.nwk"),
                         OUT_PATH = glue("{PROTEIN}_trimmed_filtered_FastTree_ed.nwk"))
{{ codestop }}

**Step 5.** Revisualize the subsetted tree

et voila!


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


## Common questions

### The ecophylo workflow died at the run_metagenomics_workflow rule and printed this message in the log file, what should I do?

```bash
Error in rule run_metagenomics_workflow:
    jobid: 0
    input: ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/metagenomics_config.json
    output: ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/metagenomics_workflow.done

RuleException:
CalledProcessError in file /Users/mschechter/github/anvio/anvio/workflows/ecophylo/rules/profile_mode.smk, line 87:

Command 'set -euo pipefail;  cd ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW && anvi-run-workflow -w metagenomics -c metagenomics_config.json --additional-params  --rerun-incomplete --latency-wait 100 --keep-going &> 00_LOGS/run_metagenomics_workflow.log && cd -' returned non-zero exit status 1.
  File "/Users/mschechter/github/anvio/anvio/workflows/ecophylo/rules/profile_mode.smk", line 87, in __rule_run_metagenomics_workflow
  File "/Users/mschechter/miniconda3/envs/anvio-dev/lib/python3.10/concurrent/futures/thread.py", line 58, in run
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```

One explanation for this error is none of the metagenomes in the %(samples-txt)s. you provided mapped any reads to extracted target proteins. To test this run the following command and see if you get this error:

```bash
$ grep -B 10 "Nothing to merge" ECOPHYLO_WORKFLOW/METAGENOMICS_WORKFLOW/00_LOGS/run_metagenomics_workflow.log
Command 'set -euo pipefail;  echo Nothing to merge for Ribosomal_S11. This should only happen if all profiles were empty (you can check the log file: 00_LOGS/Ribosomal_S11-anvi_merge.log to see if that is indeed the case). This file was created just so that your workflow would continue with no error (snakemake expects to find these output files and if we don't create them, then it will be upset). As we see it, there is no reason to throw an error here, since you mapped your metagenome to some fasta files and you got your answer: whatever you have in your fasta file is not represented in your  metagenomes. Feel free to contact us if you think that this is our fault. sincerely, Meren Lab >> 00_LOGS/Ribosomal_S11-anvi_merge.log' returned non-zero exit status 1.
```

### Can I run multiple hmms on the same data?

Yes! But sadly, not at the same time, and anvi'o feels really bad about that :(

To be clear, you can run one complete workflow, then change the path in the `"hmm_list"` parameter in the config file to a different %(hmm-list)s then rerun the workflow on the same data in the same directory. For example, this ecophylo directory contains the outputs of Ribosomal_L16 and Ribosomal_S11 over the same data:

```bash
$ tree ECOPHYLO_WORKFLOW/ -L 2

ECOPHYLO_WORKFLOW/
├── 01_REFERENCE_PROTEIN_DATA
│   ├── E_facealis_MAG
│   ├── Enterococcus_faecalis_6240
│   ├── Enterococcus_faecium_6589
│   ├── S_aureus_MAG
│   └── co_assembly
├── 02_NR_FASTAS
│   ├── Ribosomal_L16
│   └── Ribosomal_S11
├── 03_MSA
│   ├── Ribosomal_L16
│   └── Ribosomal_S11
├── 04_SEQUENCE_STATS
│   ├── Ribosomal_L16
│   └── Ribosomal_S11
├── 05_TREES
│   ├── Ribosomal_L16
│   ├── Ribosomal_L16_combined.done
│   ├── Ribosomal_S11
│   └── Ribosomal_S11_combined.done
├── 06_MISC_DATA
│   ├── Ribosomal_L16_estimate_scg_taxonomy_results-RAW-LONG-FORMAT.txt
│   ├── Ribosomal_L16_misc.tsv
│   ├── Ribosomal_L16_scg_taxonomy_data.tsv
│   ├── Ribosomal_S11_estimate_scg_taxonomy_results-RAW-LONG-FORMAT.txt
│   ├── Ribosomal_S11_misc.tsv
│   └── Ribosomal_S11_scg_taxonomy_data.tsv
├── METAGENOMICS_WORKFLOW
│   ├── 00_LOGS
│   ├── 03_CONTIGS
│   ├── 04_MAPPING
│   ├── 05_ANVIO_PROFILE
│   ├── 06_MERGED
│   ├── 07_SUMMARY
│   ├── Ribosomal_L16_ECOPHYLO_WORKFLOW_state.json
│   ├── Ribosomal_L16_add_default_collection.done
│   ├── Ribosomal_L16_state_imported_profile.done
│   ├── Ribosomal_S11_ECOPHYLO_WORKFLOW_state.json
│   ├── Ribosomal_S11_add_default_collection.done
│   ├── Ribosomal_S11_state_imported_profile.done
│   ├── fasta.txt
│   ├── metagenomics_config.json
│   ├── metagenomics_workflow.done
│   └── samples.txt
├── Ribosomal_L16_anvi_estimate_scg_taxonomy_for_SCGs.done
├── Ribosomal_S11_anvi_estimate_scg_taxonomy_for_SCGs.done
```

To visualize the results of the different ecophylo runs, just change the paths to include the different proteins:

```bash
PROTEIN=""
anvi-interactive -c METAGENOMICS_WORKFLOW/03_CONTIGS/"${PROTEIN}"-contigs.db -p METAGENOMICS_WORKFLOW/06_MERGED/"${PROTEIN}"/PROFILE.db
```

However, if you are interested in comparing the outputs of different parameters on the same protein, make a new %(workflow-config)s file, otherwise ecophylo will try and fail to overwrite your original run.

To do this, first make your new %(workflow-config)s file:
```bash
cp config_RP_L16.json config_RP_L16_new_parameters.json
```

Next, adjust any parameters you want!

Finally, edit the `"HOME"` string to a new path to ensure you make a new directory structure, like this:
```bash
# edit
"output_dirs": {
    "HOME": "ECOPHYLO_WORKFLOW_new_parameters",
    "EXTRACTED_RIBO_PROTEINS_DIR": "ECOPHYLO_WORKFLOW/01_REFERENCE_PROTEIN_DATA",
    "RIBOSOMAL_PROTEIN_FASTAS": "ECOPHYLO_WORKFLOW/02_NR_FASTAS",
    "MSA": "ECOPHYLO_WORKFLOW/03_MSA",
    "RIBOSOMAL_PROTEIN_MSA_STATS": "ECOPHYLO_WORKFLOW/04_SEQUENCE_STATS",
    "TREES": "ECOPHYLO_WORKFLOW/05_TREES",
    "MISC_DATA": "ECOPHYLO_WORKFLOW/06_MISC_DATA",
    "SCG_NT_FASTAS": "ECOPHYLO_WORKFLOW/07_SCG_NT_FASTAS",
    "RIBOSOMAL_PROTEIN_FASTAS_RENAMED": "ECOPHYLO_WORKFLOW/08_RIBOSOMAL_PROTEIN_FASTAS_RENAMED",
    "LOGS_DIR": "00_LOGS"
},
```

### Can I add more genomes and metagenomes to my analysis?

Yes you can add more genomes and metagenomes in your %(metagenomes)s, %(external-genomes)s, and %(samples-txt)s

BUT, you need to do a couple of steps first so that Snakemake can restart all the processes and maintain as much data as possible:

```bash
HOME_DIR="ECOPHYLO_WORKFLOW_asdf"
PROTEIN="Ribosomal_S11"
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/03_CONTIGS/
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/05_ANVIO_PROFILE/
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/06_MERGED/
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/07_SUMMARY/
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/"${PROTEIN}"_ECOPHYLO_WORKFLOW_state.json
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/"${PROTEIN}"_add_default_collection.done
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/"${PROTEIN}"_state_imported_profile.done
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/fasta.txt
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/metagenomics_config.json
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/metagenomics_workflow.done
rm -rf "${HOME_DIR}"/METAGENOMICS_WORKFLOW/samples.txt
```