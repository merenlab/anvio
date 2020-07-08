---
layout: post
title: profile-db [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-merge](../../programs/anvi-merge)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-3dev](../../programs/anvi-3dev)</span> <span class="artifact-r">[anvi-cluster-contigs](../../programs/anvi-cluster-contigs)</span> <span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-estimate-genome-completeness](../../programs/anvi-estimate-genome-completeness)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span> <span class="artifact-r">[anvi-export-collection](../../programs/anvi-export-collection)</span> <span class="artifact-r">[anvi-export-items-order](../../programs/anvi-export-items-order)</span> <span class="artifact-r">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span> <span class="artifact-r">[anvi-export-splits-and-coverages](../../programs/anvi-export-splits-and-coverages)</span> <span class="artifact-r">[anvi-export-state](../../programs/anvi-export-state)</span> <span class="artifact-r">[anvi-gen-fixation-index-matrix](../../programs/anvi-gen-fixation-index-matrix)</span> <span class="artifact-r">[anvi-gen-gene-level-stats-databases](../../programs/anvi-gen-gene-level-stats-databases)</span> <span class="artifact-r">[anvi-gen-variability-profile](../../programs/anvi-gen-variability-profile)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-get-short-reads-from-bam](../../programs/anvi-get-short-reads-from-bam)</span> <span class="artifact-r">[anvi-import-collection](../../programs/anvi-import-collection)</span> <span class="artifact-r">[anvi-import-items-order](../../programs/anvi-import-items-order)</span> <span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-import-state](../../programs/anvi-import-state)</span> <span class="artifact-r">[anvi-inspect](../../programs/anvi-inspect)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-refine](../../programs/anvi-refine)</span> <span class="artifact-r">[anvi-show-collections-and-bins](../../programs/anvi-show-collections-and-bins)</span> <span class="artifact-r">[anvi-split](../../programs/anvi-split)</span> <span class="artifact-r">[anvi-summarize](../../programs/anvi-summarize)</span></p>

## Description

An anvi'o database that **contains key information about the mapping of short reads *from multiple samples* to your contigs.** 

You can think of this as a extension of a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> that contains information about how your contigs align with each of your samples. The vast majority of programs that use the profile database will also ask for the contigs database associated with it. 

A profile database contains information about how short reads map to the contigs in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. Specificially, a profile database contains 
* the coverage and abundance per nucleotide position for each contig 
* varience of various kinds (single-nucleotide, single-codon, and single-amino acid)
* structural variance (ex insertions and deletions)

This information is neccessary to run anvi'o programs like <span class="artifact-n">[anvi-cluster-contigs](/software/anvio/help/programs/anvi-cluster-contigs)</span>, <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>, and <span class="artifact-n">[anvi-gen-gene-level-stats-databases](/software/anvio/help/programs/anvi-gen-gene-level-stats-databases)</span>. You can also interact with a profile database using programs like <span class="artifact-n">[anvi-interactive](/software/anvio/help/programs/anvi-interactive)</span>.

Specifically, a profile-db refers to a profile database that contains the data from several samples -- in other words, the result of running <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span> on several <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>. However, since a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span> has a lot of the functionality of a profile-db, it might be easier to think of a profile database as a header referring to both single-profile-dbs and profile-dbs (which can also be called merged-profile-dbs). For simplicity sake, since most users are dealing with multiple samples, the name was shortened to just profile-db. The following are a list of differences in functionality between a single profile database and a merged profile database:
* You can only run <span class="artifact-n">[anvi-cluster-contigs](/software/anvio/help/programs/anvi-cluster-contigs)</span> or <span class="artifact-n">[anvi-mcg-classifier](/software/anvio/help/programs/anvi-mcg-classifier)</span> on a merged profile database (or profile-db), since they look at the allignment data in many samples 
* You cannot run <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span> or <span class="artifact-n">[anvi-import-taxonomy-for-layers](/software/anvio/help/programs/anvi-import-taxonomy-for-layers)</span> on a merged profile database, only on a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>.

## How to make a profile database

### If you have multiple samples 
1. Prepare your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>
2. Run $(anvi-profile)s with an appropriate <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span>. The output of this will give you a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>. You will need to do this for each of your samples, converted into a <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span> with your short reads.
3. Run <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span> on your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> (from step 1) and your <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>s (from step 2). The output of this is a profile-db 

### If you have a single sample
1. Prepare your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>
2. Run $(anvi-profile)s with an appropriate <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span>. The output of this will give you a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>. You can see that page for more information, but essentially you can use a single-profile-db instead of a profile database to run most anvi'o functions. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/profile-db.md) to update this information.

