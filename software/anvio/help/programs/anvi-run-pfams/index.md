---
layout: page
title: anvi-run-pfams [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run Pfam on Contigs Database.

See **[program help menu](../../../vignette#anvi-run-pfams)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[functions](../../artifacts/functions)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[pfams-data](../../artifacts/pfams-data)</span></p>

## Usage


This program **associates genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> with functions using the EBI's [Pfam database](https://pfam.xfam.org/).** 

Before you run this program, you'll have to set up the Pfam database on your computer with the program <span class="artifact-n">[anvi-setup-pfams](/software/anvio/help/programs/anvi-setup-pfams)</span>.  

The Pfam database is based on protein sequences, so anvi'o will convert your genetic information into protein sequences and then use HMMs to compare them to the database. (*Unsure what an HMM is? Check out [our vocab page](http://merenlab.org/vocabulary/#hmm)*)

To run, you'll need to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. If you stored the <span class="artifact-n">[pfams-data](/software/anvio/help/artifacts/pfams-data)</span> that you got from running <span class="artifact-n">[anvi-setup-pfams](/software/anvio/help/programs/anvi-setup-pfams)</span> in a custom location, you'll need to provide that path as well. The output is a <span class="artifact-n">[functions](/software/anvio/help/artifacts/functions)</span> artifact. 

Here is a default run: 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;pfam&#45;data&#45;dir <span class="artifact&#45;n">[pfams&#45;data](/software/anvio/help/artifacts/pfams&#45;data)</span> \
</div>

By default, this uses hmmsearch to run HMMs. You can choose to use hmmscan instead by running

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
            &#45;&#45;pfam&#45;data&#45;dir <span class="artifact&#45;n">[pfams&#45;data](/software/anvio/help/artifacts/pfams&#45;data)</span> \
            &#45;&#45;hmmer&#45;program hmmscan
</div>

See [this article](https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/) for a discussion on the performance of the two HMMER programs. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-pfams.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-pfams) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
