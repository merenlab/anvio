---
layout: page
title: anvi-search-functions [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Search functions in an anvi&#x27;o contigs database or genomes storage. Basically, this program searches for one or more search terms you define in functional annotations of genes in an anvi&#x27;o contigs database, and generates multiple reports. The default report simply tells you which contigs contain genes with functions matching to serach terms you used, useful for viewing in the interface. You can also request a much more comprehensive report, which gives you anything you might need to know for each hit and serach term.

See **[program help menu](../../../vignette#anvi-search-functions)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[functions-txt](../../artifacts/functions-txt)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db)</span></p>

## Usage


This program **searches for keywords in the function annotations of your database.** 

You can use this program to look for specific functon keywords in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>, <span class="artifact-n">[genomes-storage-db](/software/anvio/help/artifacts/genomes-storage-db)</span> or <span class="artifact-n">[pan-db](/software/anvio/help/artifacts/pan-db)</span>. For example, say you wanted your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> to search for genes that encoded some type of kinase. You could call 

<div class="codeblock" markdown="1">
anvi&#45;search&#45;functions &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;search&#45;terms kinase
</div>

By default, the output will be a fairly barren <span class="artifact-n">[functions-txt](/software/anvio/help/artifacts/functions-txt)</span>, only telling you which contigs contain genes that matched your search. This will be most helpful as an additional layer in the anvi'o interactive interface, so you can quickly see where the kinase-encoding genes are in the genome. To do this, run anvi-interactive with the `--aditional-layer` parameter with the <span class="artifact-n">[functions-txt](/software/anvio/help/artifacts/functions-txt)</span>. 

However, you can also request a much more comprehensive output that contains much more information, including the matching gene's caller id, functional annotation source and full function name. 

For example, to run the same search as above, but with a more comprehensive output, you could call 

<div class="codeblock" markdown="1">
anvi&#45;search&#45;functions &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;search&#45;terms kinase
            &#45;&#45;full&#45;report kinase_information.txt
            &#45;&#45;include&#45;sequences
            &#45;&#45;verbose
</div>

Following this run, the file `kinase_information.txt` will contain comprehensive information about the matching genes, including their sequences. 

You can also search for multiple terms at the same time, or for terms from only specific annotation sources. For example, if you only wanted Pfam hits with functions related to kinases or phosphatases, you could call 

<div class="codeblock" markdown="1">
anvi&#45;search&#45;functions &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;search&#45;terms kinase,phosphatase
            &#45;&#45;annotation&#45;sources Pfam
            &#45;&#45;full&#45;report kinase_phosphatase_information.txt
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-search-functions.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-search-functions) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
