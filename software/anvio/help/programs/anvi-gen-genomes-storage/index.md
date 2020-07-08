---
layout: page
title: anvi-gen-genomes-storage [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Create a genome storage from internal and/or external genomes for a pangenome analysis.

See **[program help menu](../../../vignette#anvi-gen-genomes-storage)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[genomes-storage-db](../../artifacts/genomes-storage-db)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[external-genomes](../../artifacts/external-genomes)</span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes)</span></p>

## Usage


### From external genomes

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;genomes&#45;storage &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/artifacts/external&#45;genomes)</span> \
                         &#45;o <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/artifacts/genomes&#45;storage&#45;db)</span>
</div>

### From internal genomes

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;genomes&#45;storage &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/artifacts/internal&#45;genomes)</span> \
                         &#45;o <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/artifacts/genomes&#45;storage&#45;db)</span>
</div>

### From internal and external genomes

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;genomes&#45;storage &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/artifacts/internal&#45;genomes)</span> \
                         &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/artifacts/external&#45;genomes)</span> \
                         &#45;o <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/artifacts/genomes&#45;storage&#45;db)</span>
</div>

See also <span class="artifact-n">[anvi-pan-genome](/software/anvio/help/programs/anvi-pan-genome)</span>, which computes a pangenome from a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-genomes-storage.md) to update this information.


## Additional Resources


* [A tutorial on pangenomics](http://merenlab.org/2016/11/08/pangenomics-v2/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-genomes-storage) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
