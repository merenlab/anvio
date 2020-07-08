---
layout: page
title: anvi-run-kegg-kofams [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run KOfam HMMs on an anvi&#x27;o contigs database.

See **[program help menu](../../../vignette#anvi-run-kegg-kofams)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[kegg-functions](../../artifacts/kegg-functions)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[kegg-db](../../artifacts/kegg-db)</span></p>

## Usage


Essentially, this program uses the KEGG database to annotate functions and metabolic pathways in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. More specifically, <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/programs/anvi-run-kegg-kofams)</span> annotates a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> with HMM hits from KOfam, a database of KEGG Orthologs (KOs). You must set up these HMMs on your computer using <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/programs/anvi-setup-kegg-kofams)</span> before you can use this program.

Briefly, what this program does is extract all the gene calls from the <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> and checks each one for hits to the KOfam HMM profiles in your <span class="artifact-n">[kegg-db](/software/anvio/help/artifacts/kegg-db)</span>. This can be time-consuming given that the number of HMM profiles is quite large, and especially time-consuming if the number of genes in the <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> is also large. Multi-threading is a good idea if you have the computational capability to do so.

Many HMM hits will be found, most of them weak. The weak hits will by default be eliminated according to the score thresholds provided by KEGG; that is, only hits with scores above the threshold for a given KO profile will be annotated in the <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. It is perfectly normal to notice that the number of raw hits found is many, many times larger than the number of annotated KO hits in your database.

In the <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> functions table, annotated KO hits (<span class="artifact-n">[kegg-functions](/software/anvio/help/artifacts/kegg-functions)</span>) will have the source `KOfam`.

Running this program is a pre-requisite for metabolism estimation with <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>. Note that if you are planning to run metabolism estimation, it must be run with the same <span class="artifact-n">[kegg-db](/software/anvio/help/artifacts/kegg-db)</span> that is used in this program to annotate KOfam hits.

### Standard usage

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db
</div>

### Use a specific non-default KEGG data directory

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db &#45;&#45;kegg&#45;data&#45;dir /path/to/directory/KEGG
</div>

### Run with multiple threads

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db &#45;T 4
</div>

### Use a different HMMER program
By default, <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/programs/anvi-run-kegg-kofams)</span> uses `hmmsearch` to find KO hits. If for some reason you would rather use a different program (`hmmscan` is also currently supported), you can do so.

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db &#45;&#45;hmmer&#45;program hmmscan
</div>

### Keep all HMM hits
Usually, this program parses out weak HMM hits and keeps only those that are above the score threshold for a given KO. If you would like to turn off this behavior and keep all hits (there will be _a lot_ of weak ones), you can follow the example below:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;kegg&#45;kofams &#45;c CONTIGS.db &#45;&#45;keep&#45;all&#45;hits
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-kegg-kofams.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-kegg-kofams) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
