---
layout: post
title: raw-bam-file [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/BAM.png" alt="BAM" style="width:100px; border:none" />

A BAM-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


Most likely provided by the user.


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-init-bam](../../programs/anvi-init-bam)</span></p>

## Description

This is a **<span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span> (which contains aligned sequence data) that has not yet been indexed and sorted**. 

### What does being "indexed" mean? 

Think of your BAM file as a long, complex book. In order to get the most out of it when trying to perform analysis, it will be super helpful to have a table of contents. Indexing your BAM file basically creates a second file that serves as an external table of contents, so that anvi'o doesn't have to keep looking through the entire BAM file during analysis. 

You can tell whether or not your BAM file is indexed based on the presence of this second file, which will have the same title as your BAM file, but end with the extension `.bai`. For example, if your directory contained these files:

<div class="codeblock" markdown="1">
Lake_Michigan_Sample_1.bam
Lake_Michigan_Sample_1.bam.bai
Lake_Michigan_Sample_2.bam 
</div>

then you would still need to index `Lake_Michigan_Sample_2.bam`. 

### How do you index a BAM file?

You can either do this directly using samtools, or you can just run the anvi'o program <span class="artifact-n">[anvi-init-bam](/software/anvio/help/programs/anvi-init-bam)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/raw-bam-file.md) to update this information.

