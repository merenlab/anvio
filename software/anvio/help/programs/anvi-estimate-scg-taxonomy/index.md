---
layout: page
title: anvi-estimate-scg-taxonomy [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Estimates taxonomy at genome and metagenome level. This program is the entry point to estimate taxonomy for a given set of contigs (i.e., all contigs in a contigs database, or contigs described in collections as bins). For this, it uses single-copy core gene sequences and the GTDB database.

See **[program help menu](../../../vignette#anvi-estimate-scg-taxonomy)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[genome-taxonomy](../../artifacts/genome-taxonomy)</span> <span class="artifact-p">[genome-taxonomy-txt](../../artifacts/genome-taxonomy-txt)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[scgs-taxonomy](../../artifacts/scgs-taxonomy)</span> <span class="artifact-r">[collection](../../artifacts/collection)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span> <span class="artifact-r">[metagenomes](../../artifacts/metagenomes)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-estimate-scg-taxonomy.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [Usage examples and warnings](http://merenlab.org/scg-taxonomy)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-estimate-scg-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
