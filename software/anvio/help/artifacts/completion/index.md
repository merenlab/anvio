---
layout: post
title: completion [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-estimate-genome-completeness](../../programs/anvi-estimate-genome-completeness)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"></p>

## Description

An estimate of the completeness of purity of a genome based on single-copy core genes.

{:.notice}
See [this blog post](http://merenlab.org/2016/06/09/assessing-completion-and-contamination-of-MAGs/) for more information, and [this paper](https://doi.org/10.1038/nbt.3893) for the community standards for metagenome-assembled and single-amplified genomes.

There are two essential features to this metric: **completion**, and **redundancy**.

### Completion

A rough estimate of how completely a set of contigs represents a full genome based on the presence or absence of single-copy core genes (SCGs) they contain. 

SCGs are a set of special genes that occur in every single genome once and once only. So theoretically, the higher the percentage of SCGs found in a genome bin, the more likely that the bin represents a complete genome. Of course, SCGs are typically determined by analyzing isolate genomes that are available to find out which genes match to this criterion, hence, the accuracy of their predictions may be limited when this approach is applied to genome bins that represent populations from poorly studies clades of life. Even for genomes of well-studied organisms, our methods to identify these genes in genomes may prevent us from getting to 100% completeness.

### Redundancy

A measure of how many copies of each single-copy core gene (SCG) are found in a genome or a genome bin.

Usually, we expect to have only one copy of each of these genes (that’s why they’re called ‘single-copy’), and for this reason, redundancy of SCGs is commonly used as an estimate the level of potential ‘contamination’ within a bin (i.e., higher values of redundancy may indicate that more than one population may be contributing to a given genome bin).

However, interpretations of ‘contamination’ as a function of redundant occurrence of SCGs may not be straightforward as some genomes may have multiple copies of generally single-copy core genes, hence we prefer not to draw conclusions about contamination right away. In addition, lack of redundancy does not necessarily mean the lack of contamination, since contaminant contigs that do not include SCGs will not be in the radar of these estimates.

### Attention

Regardless of their utility to gain quick insights, single-copy core genes are mere approximations to understanding the quality of a genome and [SCGs cannot ensure the absence of contamination or level of true completion](https://doi.org/10.1101/gr.258640.119).

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/completion.md) to update this information.

