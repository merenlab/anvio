---
layout: page
title: anvi-interactive [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start an anvi&#x27;o server for the interactive interface.

See **[program help menu](../../../vignette#anvi-interactive)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection](../../artifacts/collection)</span> <span class="artifact-p">[bin](../../artifacts/bin)</span> <span class="artifact-p">[interactive](../../artifacts/interactive)</span> <span class="artifact-p">[svg](../../artifacts/svg)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[single-profile-db](../../artifacts/single-profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[genes-db](../../artifacts/genes-db)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span> <span class="artifact-r">[view-data](../../artifacts/view-data)</span> <span class="artifact-r">[dendrogram](../../artifacts/dendrogram)</span> <span class="artifact-r">[phylogeny](../../artifacts/phylogeny)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-interactive.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [A beginners tutorial on anvi&#x27;o interactive interface](http://merenlab.org/tutorials/interactive-interface/)

* [How to add more data to a display for layers and items](http://merenlab.org/2017/12/11/additional-data-tables/)

* [An overview of interactive data types](http://merenlab.org/2016/02/27/the-anvio-interactive-interface/)

* [Anvi&#x27;o &#x27;views&#x27; demystified](http://merenlab.org/2017/05/08/anvio-views/)

* [Working with SVG files from the interactive interface](http://merenlab.org/2016/10/27/high-resolution-figures/)

* [Running remote anvi&#x27;o interactive interfaces on your local computer](http://merenlab.org/2018/03/07/working-with-remote-interative/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-interactive) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
