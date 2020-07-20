---
layout: post
title: coverages-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-splits-and-coverages](../../programs/anvi-export-splits-and-coverages)</span> <span class="artifact-p">[anvi-script-get-coverage-from-bam](../../programs/anvi-script-get-coverage-from-bam)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"></p>

## Description

This is a text file containing **the average coverage for each contig in each sample** that was in the <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> that you used when you ran <span class="artifact-n">[anvi-export-splits-and-coverages](/software/anvio/help/programs/anvi-export-splits-and-coverages)</span>. 

This is a tab-delimited file where each row describes a specific contig and each column describes one of your samples. Each cell contains the average coverage of that contig in that sample. 

This artifact is really only used when taking information out of anvi'o, so enjoy your coverage information :) 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/coverages-txt.md) to update this information.

