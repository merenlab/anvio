---
layout: page
title: anvi-gen-fixation-index-matrix [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate a pairwise matrix of a fixation indices between samples.

See **[program help menu](../../../vignette#anvi-gen-fixation-index-matrix)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[fixation-index-matrix](../../artifacts/fixation-index-matrix)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[structure-db](../../artifacts/structure-db)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span> <span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-gen-fixation-index-matrix.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [Utilizing fixation index to study SAR11 population structure](http://merenlab.org/data/sar11-saavs/#generating-distance-matrices-from-fixation-index-for-saavs-and-snvs-data)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-fixation-index-matrix) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
