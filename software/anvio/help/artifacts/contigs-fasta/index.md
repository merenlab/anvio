---
layout: post
title: contigs-fasta [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-contigs](../../programs/anvi-export-contigs)</span> <span class="artifact-p">[anvi-export-splits-and-coverages](../../programs/anvi-export-splits-and-coverages)</span> <span class="artifact-p">[anvi-script-reformat-fasta](../../programs/anvi-script-reformat-fasta)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-gen-contigs-database](../../programs/anvi-gen-contigs-database)</span></p>

## Description

A <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span> is a <span class="artifact-n">[fasta](/software/anvio/help/artifacts/fasta)</span> file that is suitable to be used by <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/programs/anvi-gen-contigs-database)</span> to create a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>.

The most critical requirement for this file is that **it must have simple deflines**. If your <span class="artifact-n">[fasta](/software/anvio/help/artifacts/fasta)</span> file doesn't have simple deflines, it is not a proper <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span>. If you intend to use this file with anvi'o, **you must fix your FASTA file prior to mapping**.

Take a look at your deflines prior to mapping, and remove anything that is not a digit, an ASCII letter, an underscore, or a dash character. Here are some example deflines that are not suitable for a <span class="artifact-n">[fasta](/software/anvio/help/artifacts/fasta)</span> to be considered a <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span>

``` bash
>Contig-123 length:4567 
>Another defline 42
>gi|478446819|gb|JN117275.2|
```

And here are some OK ones:

``` bash
>Contig-123
>Another_defline_42
>gi_478446819_gb_JN117275_2
```

The program <span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/programs/anvi-script-reformat-fasta)</span> can do this automatically for you.

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/contigs-fasta.md) to update this information.

