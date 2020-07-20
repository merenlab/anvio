---
layout: page
title: anvi-get-codon-frequencies [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Get amino acid or codon frequencies of genes in a contigs database.

See **[program help menu](../../../vignette#anvi-get-codon-frequencies)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[codon-frequencies-txt](../../artifacts/codon-frequencies-txt)</span> <span class="artifact-p">[aa-frequencies-txt](../../artifacts/aa-frequencies-txt)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


This program **calculates the frequency of each codon and amino acid in your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>**. 




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-codon-frequencies.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-codon-frequencies) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
