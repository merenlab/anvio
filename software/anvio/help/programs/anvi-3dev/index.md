---
layout: page
title: anvi-3dev [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Interactively visualize sequence variants on protein structures.

See **[program help menu](../../../vignette#anvi-3dev)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive)</span> <span class="artifact-p">[summary](../../artifacts/summary)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[structure-db](../../artifacts/structure-db)</span> <span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-3dev.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-3dev) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
