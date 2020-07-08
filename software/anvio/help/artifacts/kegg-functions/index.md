---
layout: post
title: kegg-functions [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-kegg-kofams](../../programs/anvi-run-kegg-kofams)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span></p>

## Description

[Kegg Orthology](https://www.genome.jp/kegg/ko.html) (KO) functional annotations, produced by finding HMM hits to the KEGG KOfam database.

You can annotate a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> with these KEGG functions by running <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/programs/anvi-run-kegg-kofams)</span>.

Another program that relies on these annotations is <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>, which uses them to determine presence and completeness of metabolic pathways that are defined by KOs.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/kegg-functions.md) to update this information.

