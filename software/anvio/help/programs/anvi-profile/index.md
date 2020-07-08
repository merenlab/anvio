---
layout: page
title: anvi-profile [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Creates a single anvi&#x27;o profile database. The default input to this program is a BAM file.                   When it is run on a BAM file, depending on the user parameters, the program quantifies                   coverage per nucleotide position (and averages them out per contig), calculates                   single-nucleotide, single-codon, and single-amino acid variants, as well as structural variants                   such as insertion and deletions and stores these data into appropriate tables. Anvi&#x27;o single                   profiles can be merged by the program `anvi-merge`.

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





{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-profile.md) to update this information.


## Additional Resources


* [The usage of the profiler in metagenomics workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-profile) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
