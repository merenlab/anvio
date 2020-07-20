---
layout: page
title: anvi-profile [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Creates a single anvi&#x27;o profile database.                    When it is run on a BAM file, depending on the user parameters, the program quantifies                   coverage per nucleotide position (and averages them  per contig), calculates                   single-nucleotide, single-codon, and single-amino acid variants, as well as structural variants                   such as insertion and deletions and stores these data into appropriate tables.

See **[program help menu](../../../vignette#anvi-profile)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[single-profile-db](../../artifacts/single-profile-db)</span> <span class="artifact-p">[misc-data-item-orders](../../artifacts/misc-data-item-orders)</span> <span class="artifact-p">[variability-profile](../../artifacts/variability-profile)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[bam-file](../../artifacts/bam-file)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


This program **creates a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span> from a <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>**. 

Once you have a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>, you can run programs like <span class="artifact-n">[anvi-cluster-contigs](/software/anvio/help/programs/anvi-cluster-contigs)</span>, <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>, and <span class="artifact-n">[anvi-gen-gene-level-stats-databases](/software/anvio/help/programs/anvi-gen-gene-level-stats-databases)</span>, as well as use the interactive interface with <span class="artifact-n">[anvi-interactive](/software/anvio/help/programs/anvi-interactive)</span>. If you want to run these same contigs against multiple BAM files (because you have multiple samples), you'll combine your <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>s after you've created them all using <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span>. See the pages for <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> for more you can do with these artifiacts. 

In short, this program runs various analyses on the contigs in your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> and how they relate to the sample information stored in the <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span> you provided. It then stores this information into a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>. Specifically, this program calculates 
* coverage per nucleotide position (if you're unsure what coverage refers to, check out [this page](http://merenlab.org/vocabulary/#coverage))
* single-nucleotide, single-codon, and single-amino acid variants (You can find all of those terms on the vocab page linked above, as well as a more detailed explaination [here](http://merenlab.org/2015/07/20/analyzing-variability/#an-intro-to-single-nucleotidecodonamino-acid-variation))
* structural variance such as insertions or deletions 

## Basic Usage

### Inputs 

This program takes in an [indexed](https://merenlab.org/software/anvio/help/programs/anvi-init-bam) <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span> and a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. The BAM file contains the short reads from a single sample that will be used to create the profile database. Thus, here is a standard run with default parameters: 

<div class="codeblock" markdown="1">
anvi&#45;profile &#45;i <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/artifacts/bam&#45;file)</span> \
            &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> 
</div>

Alternatively, if you lack mapping data, you can add the flag `--blank-profile` so that you can still get the functionality of a profile database. 

<div class="codeblock" markdown="1">
anvi&#45;profile &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span>  \ 
            &#45;&#45;blank&#45;profile
</div>

### Checking your BAM file: always a good idea 

If you want to first check your BAM file to see what contigs it contains, just use the flag `--list-contigs` to see a comprehensive list. 

### Profiling a subset of contigs

*Note: To profile a subset of contigs based on their characterists (for example, only contigs of a certain length or that have a certain coverage), see the section on "contig specifications"*

By default, anvi'o will use every contig in your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. However, if you wish to focus specifically on a subset of these contigs, you can do so. Just provide a file that contains only the names of the contigs you want to analyze, one per line, using the tag `--contigs-of-interest`.

For example, you could run

<div class="codeblock" markdown="1">
anvi&#45;profile &#45;c Ross_sea_contigs.db  \ 
            &#45;&#45;blank&#45;profile
            &#45;&#45;contigs&#45;of&#45;interest contigs_i_like.txt
</div>

Where `contigs_i_like.txt` looks like this: 

    SF15-RossSeacontig4922
    SF15-RossSeacontig702

## Analysis Parameters
 
Changing these will affect the way that your sequences are analyzed. 

Keep in mind that if you plan to merge your resulting <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span> with others later in the project, you'll want to keep these parameters consistent. 

### Contig Specification 

To profile only contigs within a specific length, you can use the flags `--min-contig-length` and `-max-contig-length`. By default, the minimum length for analysis is 1000 and there is no maximum length. You can also profile only contigs that have a certain average coverage. 

#### Specifications for your BAM file

You can also ignore reads in your BAM file with a percent identity to the reference less than some threshold using the flag `--min-percent-identity`.  By default, all reads are used. 

For example, the following code will only look at contigs longer than 2000 nts and will ignore BAM file reads with less than 95 percent identity to the reference:

<div class="codeblock" markdown="1">
anvi&#45;profile &#45;c Ross_sea_contigs.db  \ 
            &#45;i bam_file.bam \
            &#45;&#45;min&#45;contig&#45;length 2000 \
            &#45;&#45;min&#45;percent&#45;identity 95 
</div>

### Hierarchical Clustering 

#### To cluster or not to cluster? 

By default, anvi'o will not try to cluster your splits (since it takes quite a bit of runtime) unless you are using the tag `--blank-profile`. If you don't want to run this, use the tag `--skip-hierarchical-clustering`. 

If you're planning to later merge this sample with others, it is better to perform clustering while running <span class="artifact-n">[anvi-merge](/software/anvio/help/programs/anvi-merge)</span> than at this stage. 

However, if you want to bin this single sample or otherwise want clustering to happen, just use the tag `--cluster-contigs`. 

If you do plan to cluster, you can set a custom distance metric or a custom linkage method. 

### Variability 

Anvi-profile will throw away variability data below certain thresholds to reduce noise. After all, if you have a single C read at a position with a 1000X coverage where all other reads are T, this is probably not a variance position that you want to investigate further. By default, it will not analyze positions will coverage less than 10X, and it will not report every position. 

However, you can change the coverage threshold using the  `--min-coverage-for-variability` flag. You can also report every variability position using the flag `--report-variability-full`. 

For example, if your data quality was poor and you wanted to view every single position for variance, you could call the following: 

<div class="codeblock" markdown="1">
anvi&#45;profile &#45;c Ross_sea_contigs.db  \ 
            &#45;i bam_file.bam \
            &#45;&#45;min&#45;coverage&#45;for&#45;variability 1 \
            &#45;&#45;report&#45;variability&#45;full
</div>

## Other Parameters 

You should provide the sample name with the flag `-S` and can provide a description of your project using the `--description` tag followed by a text file. These will help anvi'o name output files and will show up in the anvi'o interfaces down the line. 

You can characterize the codon frequencies of genes in your sample at the cost of some runtime. Despite time being money, codon frequency analysis can be helpful downstream. Simply add the tag `--profile-SCVs` and watch the magic happen. 

Alternatively, you can choose not to store insertion and deletion data or single nucleotide variance data.

If you know the limits of your system, you can also multithread this program. See the program help menu for more information. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-profile.md) to update this information.


## Additional Resources


* [The usage of the profiler in metagenomics workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-profile) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
