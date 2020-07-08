---
layout: page
title: anvi-refine [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start an anvi&#x27;o interactive interactive to manually curate or refine a genome, whether it is a metagenome-assembled, single-cell, or an isolate genome.

See **[program help menu](../../../vignette#anvi-refine)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[bin](../../artifacts/bin)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-refine.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [Refining a bin](http://merenlab.org/2015/05/11/anvi-refine/)

* [Notes on genome refinement](http://merenlab.org/2017/05/11/anvi-refine-by-veronika/)

* [A case study: Inspecting the genomic link between Archaea and Eukaryota](http://merenlab.org/2017/01/03/loki-the-link-archaea-eukaryota/)

* [A demo](https://www.youtube.com/watch?v=vXPKP5vKiBM)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-refine) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
