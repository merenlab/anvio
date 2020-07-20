---
layout: page
title: anvi-script-gen-hmm-hits-matrix-across-genomes [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

A simple script to generate a TAB-delimited file that reports the frequency of HMM hits for a given HMM source across contigs databases..

See **[program help menu](../../../vignette#anvi-script-gen-hmm-hits-matrix-across-genomes)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[hmm-source](../../artifacts/hmm-source)</span> <span class="artifact-p">[hmm-hits](../../artifacts/hmm-hits)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-script-gen-hmm-hits-matrix-across-genomes.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-hmm-hits-matrix-across-genomes) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
