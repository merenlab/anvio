---
layout: page
title: anvi-scan-trnas [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Identify and store tRNA genes in a contigs database.

See **[program help menu](../../../vignette#anvi-scan-trnas)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[hmm-hits](../../artifacts/hmm-hits)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


This program identifies the tRNA genes in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> and stores them in an <span class="artifact-n">[hmm-hits](/software/anvio/help/artifacts/hmm-hits)</span>. 

To run, just provide a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> that you want to look through. 

<div class="codeblock" markdown="1">
anvi&#45;scan&#45;trnas &#45;c CONTIGS_DB
</div>

### Customizing the cut off score

What counts as a tRNA gene? That could be up to you. 

The default minimum score for a gene to be counted is 20.  However, you can set this cutoff to anywhere between 0-100. This value is actually used by the module tRNAScan-SE, so view their documentation for details (here is a [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768409/) to their paper). For example, to find more non-cononical tRNA genes, a user could lower the cutoff score to 10 as follows:

<div class="codeblock" markdown="1">
anvi&#45;scan&#45;trnas &#45;c CONTIGS_DB &#45;&#45;trna&#45;cut&#45;off&#45;score 10
</div>

### Other options 

- It is easy to modify where the outputs will go:

    - Use the parameter `--log-file` to provide a path for the output messages to go.
    - Use the parameter `--trna-hits-file` to provide a path for the raw tRNA scan data to go. 
- Like many other anvi'o programs, you can use the tag `--just-do-it` to not have to look at questions or warnings
- You can also try to multithread whenever possible by setting the `--num-threads` parameter (it is 1 by default). This can be used to speed up runtime, but please be aware of your system and its limitations before trying this. 

### Understanding the output 

Essentially, the output of this program states the probability that each gene is a tRNA gene. See <span class="artifact-n">[hmm-hits](/software/anvio/help/artifacts/hmm-hits)</span> for more information. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-scan-trnas.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-scan-trnas) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
