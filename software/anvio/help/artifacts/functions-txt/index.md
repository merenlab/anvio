---
layout: post
title: functions-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-functions](../../programs/anvi-export-functions)</span> <span class="artifact-p">[anvi-search-functions](../../programs/anvi-search-functions)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-functions](../../programs/anvi-import-functions)</span></p>

## Description

This artifact **contains information about the functions of its contigs in a tab delimated text file.**

This file is formatted with a single gene per row with the following columns: 
1. The gene caller ID
2. The source (the database that you got this function data from)
3. The assession number of the gene (optional)
4. The function information
1. The e-value (optional, but helpful if you plan to use the interactive interface and have more than one function for a single gene)

For an example, check out [this lovely page](http://merenlab.org/2016/06/18/importing-functions/#simple-matrix). 

This is primarily used to import and export function information from an anvi'o project. See [this page](http://merenlab.org/2016/06/18/importing-functions/) for more information on importing function information. 

It is also the output of <span class="artifact-n">[anvi-search-functions](/software/anvio/help/programs/anvi-search-functions)</span> which searches for specific terms in your functional annotations. For example, you could search for all genes related to "ribosomes" and get a functions-txt of all of those genes. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/functions-txt.md) to update this information.

