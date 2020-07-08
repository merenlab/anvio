---
layout: post
title: single-profile-db [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-profile](../../programs/anvi-profile)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-import-taxonomy-for-layers](../../programs/anvi-import-taxonomy-for-layers)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-merge](../../programs/anvi-merge)</span></p>

## Description

An anvi'o database that contains the same information as a merged <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span>, namely **key information about the mapping of short reads *in a single sample* to your contigs.** 

You can think of this as a extension of a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> that contains information about how your contigs align with a single one of your individual samples. If you have more than one sample, you'll probably want to use <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span> to merge your databases into a merged <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span>. The vast majority of programs that use the profile database will also ask for the contigs database associated with it. 

A single profile database contains information about how the short reads in a single BAM-file map to the contigs in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. Specificially, a profile database contains 
* the coverage and abundance per nucleotide position for each contig 
* varience of various kinds (single-nucleotide, single-codon, and single-amino acid)
* structural variance (ex insertions and deletions)

Once created, a single profile database is almost interchangable with a <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> (even though the names can be a little confusing. Think of a single-profile-db as a type of profile-db, since it has only a few differences). The main differences between the two are as follows: 
* You cannot run <span class="artifact-n">[anvi-cluster-contigs](/software/anvio/help/programs/anvi-cluster-contigs)</span> or <span class="artifact-n">[anvi-mcg-classifier](/software/anvio/help/programs/anvi-mcg-classifier)</span> on a single profile db, since these two programs look at the alignment data in many samples. 
* You can run <span class="artifact-n">[anvi-import-taxonomy-for-layers](/software/anvio/help/programs/anvi-import-taxonomy-for-layers)</span> on a single profile database but not a merged one. 
* You can only run <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span> on a single profile database.

If you want to look inside a single profile database, you can do so using <span class="artifact-n">[anvi-interactive](/software/anvio/help/programs/anvi-interactive)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/single-profile-db.md) to update this information.

