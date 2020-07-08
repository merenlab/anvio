---
layout: page
title: anvi-import-misc-data [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Populate additional data or order tables in pan or profile databases for items and layers, OR additional data in contigs databases for nucleotides and amino acids (the Swiss army knife-level serious stuff).

See **[program help menu](../../../vignette#anvi-import-misc-data)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[misc-data-items](../../artifacts/misc-data-items)</span> <span class="artifact-p">[misc-data-layers](../../artifacts/misc-data-layers)</span> <span class="artifact-p">[misc-data-layer-orders](../../artifacts/misc-data-layer-orders)</span> <span class="artifact-p">[misc-data-nucleotides](../../artifacts/misc-data-nucleotides)</span> <span class="artifact-p">[misc-data-amino-acids](../../artifacts/misc-data-amino-acids)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[genes-db](../../artifacts/genes-db)</span> <span class="artifact-r">[pan-db](../../artifacts/pan-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[misc-data-items-txt](../../artifacts/misc-data-items-txt)</span> <span class="artifact-r">[dendrogram](../../artifacts/dendrogram)</span> <span class="artifact-r">[phylogeny](../../artifacts/phylogeny)</span> <span class="artifact-r">[misc-data-layers-txt](../../artifacts/misc-data-layers-txt)</span> <span class="artifact-r">[misc-data-layer-orders-txt](../../artifacts/misc-data-layer-orders-txt)</span> <span class="artifact-r">[misc-data-nucleotides-txt](../../artifacts/misc-data-nucleotides-txt)</span> <span class="artifact-r">[misc-data-amino-acids-txt](../../artifacts/misc-data-amino-acids-txt)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-import-misc-data.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [A primer on anvi&#x27;o misc data tables](http://merenlab.org/2017/12/11/additional-data-tables/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-misc-data) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
