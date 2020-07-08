---
layout: post
title: pan-db [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-pan-genome](../../programs/anvi-pan-genome)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-analyze-synteny](../../programs/anvi-analyze-synteny)</span> <span class="artifact-r">[anvi-compute-genome-similarity](../../programs/anvi-compute-genome-similarity)</span> <span class="artifact-r">[anvi-display-pan](../../programs/anvi-display-pan)</span> <span class="artifact-r">[anvi-export-items-order](../../programs/anvi-export-items-order)</span> <span class="artifact-r">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span> <span class="artifact-r">[anvi-export-state](../../programs/anvi-export-state)</span> <span class="artifact-r">[anvi-get-enriched-functions-per-pan-group](../../programs/anvi-get-enriched-functions-per-pan-group)</span> <span class="artifact-r">[anvi-get-sequences-for-gene-clusters](../../programs/anvi-get-sequences-for-gene-clusters)</span> <span class="artifact-r">[anvi-import-collection](../../programs/anvi-import-collection)</span> <span class="artifact-r">[anvi-import-items-order](../../programs/anvi-import-items-order)</span> <span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-import-state](../../programs/anvi-import-state)</span> <span class="artifact-r">[anvi-meta-pan-genome](../../programs/anvi-meta-pan-genome)</span> <span class="artifact-r">[anvi-show-collections-and-bins](../../programs/anvi-show-collections-and-bins)</span> <span class="artifact-r">[anvi-summarize](../../programs/anvi-summarize)</span></p>

## Description

An anvi’o database that contains **key information associated with your gene clusters**.

This is vital for its pangenomic analsis, hence the name. 

This is the output of the program <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/programs/anvi-pan-genome)</span>, which can be run after you've created a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span> with the genomes you want to analyze. That script will identify your gene clusters and put them in a pan-db. 

Then, you can use the pan database to run a variety of other functions. You can think of it as a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> focused on the gene clusters of multiple genomes. Some of these functions include <span class="artifact-n">[anvi-compute-genome-similarity](/software/anvio/help/programs/anvi-compute-genome-similarity)</span>, <span class="artifact-n">[anvi-analyze-synteny](/software/anvio/help/programs/anvi-analyze-synteny)</span>, and <span class="artifact-n">[anvi-get-enriched-functions-per-pan-group](/software/anvio/help/programs/anvi-get-enriched-functions-per-pan-group)</span>. You can also view the data in a pan-db using <span class="artifact-n">[anvi-display-pan](/software/anvio/help/programs/anvi-display-pan)</span>. 

## Advanced information for programmers

While it is possible to read and write a given anvi'o pan database through SQLite functions directly, one can also use anvi'o libraries to initiate a pan database to read from. 

### Initiate a pan database instance

``` python
import argparse

from anvio.dbops import PanSuperclass

args = argparse.Namespace(pan_db="PAN.db", genomes_storage="GENOMES.db")

pan_db = PanSuperclass(args)

```

### Gene clusters dictionary

Once an instance from `PanSuperclass` is initiated, the following member function will give access to gene clusters:

``` pyton
pan_db.init_gene_clusters()
print(pan_db.gene_clusters)
```

```
{
  "GC_00000001": {
    "Genome_A": [19, 21],
    "Genome_B": [30, 32],
    "Genome_C": [122, 125],
    "Genome_D": [44, 42]
  },
  "GC_00000002": {
    "Genome_A": [123],
    "Genome_B": [176],
    "Genome_C": [175],
    "Genome_D": []
  },
  (...)
  "GC_00000036": {
    "Genome_A": [],
    "Genome_B": [24],
    "Genome_C": [],
    "Genome_D": []
  }
  (...)
```

Each item in this dictionary is a gene cluster describes anvi'o gene caller ids of each gene from each genome that contributes to this cluster.

### Sequences in gene clusters

```
gene_clusters_of_interest = set(["GC_00000006", "GC_00000036"])
gene_cluster_sequences = pan_db.get_sequences_for_gene_clusters(gene_cluster_names= gene_clusters_of_interest)

print(gene_cluster_sequences)
```

```
{
  "GC_00000006": {
    "Genome_A": {
      23: "MDVKKGWSGNNLND--NNNGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    },
    "Genome_B": {
      34: "MDVKKGWSGNNLND--NNNGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    },
    "Genome_C": {
      23: "MDVKKGWSGNNLNDWVNNNGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    },
    "Genome_D": {
      23: "MDVKKGWSGNNLNDWVNNAGSFTLFNAYLPQAKLANEAMHQKIMEMSAKAPNATMSITGHSLGTMISIQAVANLPQAD"
    }
  },
  "GC_00000036": {
    "Genome_A": {},
    "Genome_B": {
      24: "MSKRHKFKQFMKKKNLNPMNNRKKVGIILFATSIGLFFLFAFRTTYIVATGKVAGVSLKEKTA"
    },
    "Genome_C": {},
    "Genome_D": {}
  }
}
```


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/pan-db.md) to update this information.

