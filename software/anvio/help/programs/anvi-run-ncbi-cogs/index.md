---
layout: page
title: anvi-run-ncbi-cogs [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run NCBI&#x27;s COGs to associate genes in an anvi&#x27;o contigs database with functions. COGs database was been designed as an attempt to classify proteins from completely sequenced genomes on the basis of the orthology concept. It is no longer actively developed, however, it is still very effective for daily needs. You may want to consider Pfams or the eggNOG database for more comprehensive functional insights.

See **[program help menu](../../../vignette#anvi-run-ncbi-cogs)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[functions](../../artifacts/functions)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[cogs-data](../../artifacts/cogs-data)</span></p>

## Usage


This program **associates genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> with functions using NCBI's [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/).**

Before you run this program, you'll have to set up the COGs database on your computer with the program <span class="artifact-n">[anvi-setup-ncbi-cogs](/software/anvio/help/programs/anvi-setup-ncbi-cogs)</span>.  

As mentioned above, the COGs database is no longer actively added to, so might also want to consider using a separate database. As of yet, anvi'o does not have a program to accesss the eggNOG database (instructions to use this database to get function information are [here](http://merenlab.org/2016/06/18/importing-functions/#eggnog-database--emapper)), but does have the functionality to use the Pfams database (check out <span class="artifact-n">[anvi-run-pfams](/software/anvio/help/programs/anvi-run-pfams)</span> for more information). 

To run, you'll need to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. If you stored the <span class="artifact-n">[cogs-data](/software/anvio/help/artifacts/cogs-data)</span> that you got from running <span class="artifact-n">[anvi-setup-ncbi-cogs](/software/anvio/help/programs/anvi-setup-ncbi-cogs)</span> in a custom location, you'll need to provide that path as well. The output is a <span class="artifact-n">[functions](/software/anvio/help/artifacts/functions)</span> artifact. 

By default, this program uses DIAMOND in the "fast" setting for database searching. To instead run in "sensitive" mode, just call: 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;cog&#45;data&#45;dir <span class="artifact&#45;n">[cogs&#45;data](/software/anvio/help/artifacts/cogs&#45;data)</span> \
            &#45;&#45;sensitive
</div>

You can also use blastp to search, by running: 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;search&#45;with blastp
</div>

*Note: without the flag `--cog-data-dir`, anvi'o will just search in the default location.*




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-ncbi-cogs.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-ncbi-cogs) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
