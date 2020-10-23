# Snakemake workflow for tRNA-seq datasets

## Contents

- [Snakemake workflow for tRNA-seq datasets](#snakemake-workflow-for-trna-seq-datasets)
  - [Contents](#contents)
  - [Introduction](#introduction)
  - [Standard Usage](#standard-usage)
    - [Pre-Quality Control Input](#pre-quality-control-input)
    - [Post-Quality Control Input](#post-quality-control-input)

## Introduction

The pipeline includes the following steps:

1. Quality control tRNA-seq datasets using [Illumina-utils](https://github.com/merenlab/illumina-utils/) (*optional step*).
2. Construct a database of annotated tRNA sequences from each dataset using anvi-find-tRNA.
3. Cluster tRNA sequences from all datasets to choose seed sequences using [MMseqs2](https://github.com/soedinglab/MMseqs2).
4. Construct an Anvi'o tRNASEEDS database using anvi-gen-tRNAseeds-database.
5. Map tRNA sequences from each dataset to tRNA seed sequences using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
6. Construct a database of nucleotide modifications in the seed sequences using anvi-find-tRNA-mods.
7. Refine seeds by splitting seed clusters into clusters representing true sequence variants rather than modification-induced variants using anvi-refine-tRNA-seeds, map tRNA sequences from each dataset to the new seeds using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and construct a database of nucleotide modifications in the new seed sequences using anvi-find-tRNA-mods (*optional step*).
8. Construct Anvi'o tRNAPROFILE databases from each dataset's BAM file using anvi-tRNA-profile.
9. Merge tRNAPROFILE databases from each pair of demethylase-treated/untreated datasets using anvi-merge-tRNA-treatments (*optional step*).
10. Merge resulting tRNAPROFILE databases using anvi-merge-tRNA.

## Standard Usage

tRNA-seq samples can be split for the analysis of nucleotide methylation -- one split is left untreated, while the other is treated with demethylase. It is *optional* for the input to this pipeline to include paired demethylase-treated/untreated datasets.

### *Pre-Quality Control Input*

All input datasets should be described in a tab-delimited file formatted like the pre-QC example file available [here](mock_files_for_merenlab_tRNAseq_workflow/preQC-samples.txt). When a sample doesn't have paired demethylase-treated/untreated datasets, enter `NA` in the appropriate places.

The file describes the names of each sample, which datasets from each sample are demethylase-treated/untreated, and the relative or absolute paths to each dataset's paired-end FASTQ files (for now, we do not support single-end reads or interleaved FASTQ files).

### *Post-Quality Control Input*

All input datasets should be described in a tab-delimited file formatted like the post-QC example file available [here](mock_files_for_merenlab_tRNAseq_workflow/postQC-samples.txt). When a sample doesn't have paired demethylase-treated/untreated datasets, enter `NA` in the appropriate place.

The file describes the names of each sample, which datasets from each sample are demethylase-treated/untreated, and the relative or absolute path to each dataset's FASTA file.

[Back to Table of Contents](#contents)
