---
layout: post
title: hmm-hits [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-hmms](../../programs/anvi-run-hmms)</span> <span class="artifact-p">[anvi-scan-trnas](../../programs/anvi-scan-trnas)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"></p>

## Description

The search results for an <span class="artifact-n">[hmm-source](/software/anvio/help/artifacts/hmm-source)</span> in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. Essentially, this is the part of an individual <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> that handles the HMM data (In anvi'o, this is usually information about gene annotation). 

Upon creation, a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> will not contain any HMM results. In order to populate it, users can run <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/programs/anvi-run-hmms)</span> to use a custom or default anvio <span class="artifact-n">[hmm-source](/software/anvio/help/artifacts/hmm-source)</span>. The program <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/programs/anvi-scan-trnas)</span> also populates a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>'s hmm-hits.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/hmm-hits.md) to update this information.

