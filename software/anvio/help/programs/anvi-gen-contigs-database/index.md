---
layout: page
title: anvi-gen-contigs-database [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate a new anvi&#x27;o contigs database.

See **[program help menu](../../../vignette#anvi-gen-contigs-database)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-db](../../artifacts/contigs-db)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-fasta](../../artifacts/contigs-fasta)</span> <span class="artifact-r">[external-gene-calls](../../artifacts/external-gene-calls)</span></p>

## Usage


The input for this program is a <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span>, which should contain one or more sequences. These sequences may belong to a single genome, or could be contigs obtained from an assembly.

Make sure the input file file matches to requirements of <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span>. If you are planning to use the resulting contigs-db with <span class="artifact-n">[anvi-profile](/software/anvio/help/programs/anvi-profile)</span>, it is essential that you convert your <span class="artifact-n">[fasta](/software/anvio/help/artifacts/fasta)</span> file to a properly formatted <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span> *before* you perform the read recruitment.

An anvi'o contigs database will keep all the information related to your sequences: positions of open reading frames, k-mer frequencies for each contig, functional and taxonomic annotation of genes, etc. The contigs database is one of the most essential component of anvi'o.

When run on a <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span> this program will,

* **Compute k-mer frequencies** for each contig (the default is `4`, but you can change it using `--kmer-size` parameter if you feel adventurous).

* **Soft-split contigs** longer than 20,000 bp into smaller ones (you can change the split size using the `--split-length`). When the gene calling step is not skipped, the process of splitting contigs will consider where genes are and avoid cutting genes in the middle. For very very large assemblies this process can take a while, and you can skip it with `--skip-mindful-splitting` flag.

* **Identify open reading frames** using [Prodigal](http://prodigal.ornl.gov/), UNLESS, (1) you have used the flag `--skip-gene-calling` (no gene calls will be made) or (2) you have provided <span class="artifact-n">[external-gene-calls](/software/anvio/help/artifacts/external-gene-calls)</span>.


### Create a contigs database from a FASTA file

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;contigs&#45;database &#45;f <span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/artifacts/contigs&#45;fasta)</span> \
                          &#45;o <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span>
</div>

### Create a contigs database with external gene calls

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;contigs&#45;database &#45;f <span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/artifacts/contigs&#45;fasta)</span> \
                          &#45;o <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                          &#45;&#45;external&#45;gene&#45;calls <span class="artifact&#45;n">[external&#45;gene&#45;calls](/software/anvio/help/artifacts/external&#45;gene&#45;calls)</span>
</div>

See <span class="artifact-n">[external-gene-calls](/software/anvio/help/artifacts/external-gene-calls)</span> for the description and formatting requirements of this file.

If user-provided or anvi'o-calculated amino acid sequences contain internal stop codons, anvi'o will yield an error. The following command will persist through this error:

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;contigs&#45;database &#45;f <span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/artifacts/contigs&#45;fasta)</span> \
                          &#45;o <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                          &#45;&#45;external&#45;gene&#45;calls <span class="artifact&#45;n">[external&#45;gene&#45;calls](/software/anvio/help/artifacts/external&#45;gene&#45;calls)</span> \
                          &#45;&#45;ignore&#45;internal&#45;stop&#45;codons
</div>

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-contigs-database.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-contigs-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
