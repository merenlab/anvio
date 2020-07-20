---
layout: page
title: anvi-init-bam [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Sort/Index BAM files.

See **[program help menu](../../../vignette#anvi-init-bam)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[bam-file](../../artifacts/bam-file)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[raw-bam-file](../../artifacts/raw-bam-file)</span></p>

## Usage


This program sorts and indexes your BAM files: essentially converting <span class="artifact-n">[raw-bam-file](/software/anvio/help/artifacts/raw-bam-file)</span> into <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span>, which are ready to be used in anvi'o. 

If you're unsure what a BAM file is, check out the <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span> page or [this file](https://samtools.github.io/hts-specs/SAMv1.pdf), written by the developers of samtools. For a descritption of what indexing a BAM file does, check out the page for <span class="artifact-n">[raw-bam-file](/software/anvio/help/artifacts/raw-bam-file)</span>. 

To run this program, just provide a path to the bam files that you want to index. For example, 

<div class="codeblock" markdown="1">
anvi&#45;init&#45;bam <span class="artifact&#45;n">[raw&#45;bam&#45;file](/software/anvio/help/artifacts/raw&#45;bam&#45;file)</span> 
</div>

You can also multithread this to shorten runtime with the flag `-T` followed by the desired number of threads if your system is capable of this. 

To see it in action (plus a description on how to run it on an entire folder), check out [this page](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-init-bam). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-init-bam.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-init-bam) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
