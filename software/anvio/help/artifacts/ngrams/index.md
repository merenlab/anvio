---
layout: post
title: ngrams [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-analyze-synteny](../../programs/anvi-analyze-synteny)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"></p>

## Description

An <span class="artifact-n">[ngrams](/software/anvio/help/artifacts/ngrams)</span> object is a DataFrame that contains count data of synteny patterns collected from a group of similar loci or genomes. It is produced by running <span class="artifact-n">[anvi-analyze-synteny](/software/anvio/help/programs/anvi-analyze-synteny)</span> when given a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span> and an annotation source.

An `ngram` is a group of neighboring genes that include precisely `n` genes, inspired by the term ngram in [linguistics and natural language processing](https://en.wikipedia.org/wiki/N-gram). This object was inspired by kmer count tables but is inherently different because it is counting adjacent genes and not nucleotides.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/ngrams.md) to update this information.

