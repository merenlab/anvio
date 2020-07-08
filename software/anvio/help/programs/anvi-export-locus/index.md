---
layout: page
title: anvi-export-locus [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program helps you cut a &#x27;locus&#x27; from a larger genetic context (e.g., contigs, genomes). By default, anvi&#x27;o will locate a user-defined anchor gene, extend its selection upstream and downstream based on the --num-genes argument, then extract the locus to create a new contigs database. The anchor gene must be provided as --search-term, --gene-caller-ids, or --hmm-sources. If --flank-mode is designated, you MUST provide TWO flanking genes that define the locus region (Please see --flank-mode help for more information). If everything goes as plan, anvi&#x27;o will give you individual locus contigs databases for every matching anchor gene found in the original contigs database provided. Enjoy your mini contigs databases!.

See **[program help menu](../../../vignette#anvi-export-locus)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[locus-fasta](../../artifacts/locus-fasta)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-export-locus.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-locus) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
