---
layout: post
title: kegg-db [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-setup-kegg-kofams](../../programs/anvi-setup-kegg-kofams)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-run-kegg-kofams](../../programs/anvi-run-kegg-kofams)</span></p>

## Description

A **directory of data** downloaded from the [KEGG database resource](https://www.kegg.jp/) for use in function annotation and metabolism estimation.

It is created by running the program <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/programs/anvi-setup-kegg-kofams)</span>. Not everything from KEGG is included in this directory, only the information relevant to downstream programs. The most critical components of this directory are KOfam HMM profiles and the `MODULES.db` database which contains information on metabolic pathways as described in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html).

Programs that rely on this data directory include <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/programs/anvi-run-kegg-kofams)</span> and <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/kegg-db.md) to update this information.

