---
layout: post
title: metagenomes [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


Most likely provided by the user.


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span></p>

## Description

A metagenome is any set of sequences that collectively describes multiple different populations (rather than just one genome) and has been converted into a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>.

Metagenomes file format enables anvi'o to work with one or more metagenomes. A TAB-delimited external genomes file will be composed of at least the following two columns:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

In some cases, (<span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/programs/anvi-estimate-scg-taxonomy)</span>, for example), you may also want to provide the <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> that is associated with the <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. Then the metagenomes file will be composed of three columns:

|name|contigs_db_path|profile_db_path|
|:--|:--|:--|
|Name_01|/path/to/contigs-01.db|/path/to/profile.db|
|Name_02|/path/to/contigs-02.db|/path/to/profile.db|
|Name_03|/path/to/contigs-03.db|/path/to/profile.db|
|(...)|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see, **<span class="artifact-n">[internal-genomes](/software/anvio/help/artifacts/internal-genomes)</span>** and **<span class="artifact-n">[external-genomes](/software/anvio/help/artifacts/external-genomes)</span>**.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/metagenomes.md) to update this information.

