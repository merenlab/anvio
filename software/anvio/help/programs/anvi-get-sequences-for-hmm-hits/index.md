---
layout: page
title: anvi-get-sequences-for-hmm-hits [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Get sequences for HMM hits from many inputs.

See **[program help menu](../../../vignette#anvi-get-sequences-for-hmm-hits)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[genes-fasta](../../artifacts/genes-fasta)</span> <span class="artifact-p">[concatenated-gene-alignment-fasta](../../artifacts/concatenated-gene-alignment-fasta)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes)</span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes)</span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source)</span></p>

## Usage


This program can work with anvi'o <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>, <span class="artifact-n">[external-genomes](/software/anvio/help/artifacts/external-genomes)</span>, or <span class="artifact-n">[internal-genomes](/software/anvio/help/artifacts/internal-genomes)</span> files to return sequences for HMM hits identified through the default anvi'o <span class="artifact-n">[hmm-source](/software/anvio/help/artifacts/hmm-source)</span>s (such as the domain-specific single-copy core genes) or user-defined <span class="artifact-n">[hmm-source](/software/anvio/help/artifacts/hmm-source)</span>s (such as HMMs for specific antibiotic resistance gene families or any other targets).

Using it with single-copy core genes in default anvi'o HMMs make it a very versatile tool for phylogenomics as the user can define specific sets of genes to be aligned and concatenated.


### Learn available HMM sources

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;&#45;list&#45;hmm&#45;sources

AVAILABLE HMM SOURCES
&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;&#61;=
&#42; 'Bacteria_71' (type: singlecopy; num genes: 71)
&#42; 'Archaea_76' (type: singlecopy; num genes: 76)
&#42; 'Protista_83' (type: singlecopy; num genes: 83)
&#42; 'Ribosomal_RNAs' (type: Ribosomal_RNAs; num genes: 12)
</div>

### Get all sequences in a given HMM source

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;&#45;hmm&#45;source Bacteria_71 \
                                &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/artifacts/genes&#45;fasta)</span>
</div>

### Learn available genes in a given HMM source

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;&#45;hmm&#45;source Bacteria_71 \
                                &#45;&#45;list&#45;available&#45;gene&#45;names

&#42; Bacteria_71 [type: singlecopy]: ADK, AICARFT_IMPCHas, ATP&#45;synt, ATP&#45;synt_A,
Chorismate_synt, EF_TS, Exonuc_VII_L, GrpE, Ham1p_like, IPPT, OSCP, PGK,
Pept_tRNA_hydro, RBFA, RNA_pol_L, RNA_pol_Rpb6, RRF, RecO_C, Ribonuclease_P,
Ribosom_S12_S23, Ribosomal_L1, Ribosomal_L13, Ribosomal_L14, Ribosomal_L16,
Ribosomal_L17, Ribosomal_L18p, Ribosomal_L19, Ribosomal_L2, Ribosomal_L20,
Ribosomal_L21p, Ribosomal_L22, Ribosomal_L23, Ribosomal_L27, Ribosomal_L27A,
Ribosomal_L28, Ribosomal_L29, Ribosomal_L3, Ribosomal_L32p, Ribosomal_L35p,
Ribosomal_L4, Ribosomal_L5, Ribosomal_L6, Ribosomal_L9_C, Ribosomal_S10,
Ribosomal_S11, Ribosomal_S13, Ribosomal_S15, Ribosomal_S16, Ribosomal_S17,
Ribosomal_S19, Ribosomal_S2, Ribosomal_S20p, Ribosomal_S3_C, Ribosomal_S6,
Ribosomal_S7, Ribosomal_S8, Ribosomal_S9, RsfS, RuvX, SecE, SecG, SecY, SmpB,
TsaE, UPF0054, YajC, eIF&#45;1a, ribosomal_L24, tRNA&#45;synt_1d, tRNA_m1G_MT,
Adenylsucc_synt
</div>

### Get sequences for some sequences in a given HMM source

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;&#45;hmm&#45;source Bacteria_71 \
                                &#45;&#45;gene&#45;names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/artifacts/genes&#45;fasta)</span>
</div>

### Get HMM hits in bins of a collection

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/artifacts/profile&#45;db)</span> \
                                &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/artifacts/collection)</span>
                                &#45;&#45;hmm&#45;source Bacteria_71 \
                                &#45;&#45;gene&#45;names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/artifacts/genes&#45;fasta)</span>
</div>

### Get amino acid sequences for HMM hits

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/artifacts/profile&#45;db)</span> \
                                &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/artifacts/collection)</span>
                                &#45;&#45;hmm&#45;source Bacteria_71 \
                                &#45;&#45;gene&#45;names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                &#45;&#45;get&#45;aa&#45;sequences \
                                &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/artifacts/genes&#45;fasta)</span>
</div>

### Get HMM hits independently aligned and concatenated

The resulting file can be used for phylogenomics analyses via <span class="artifact-n">[anvi-gen-phylogenomic-tree](/software/anvio/help/programs/anvi-gen-phylogenomic-tree)</span> or through more sophisticated tools for curating alignments and computing trees.

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                                &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/artifacts/profile&#45;db)</span> \
                                &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/artifacts/collection)</span>
                                &#45;&#45;hmm&#45;source Bacteria_71 \
                                &#45;&#45;gene&#45;names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                &#45;&#45;get&#45;aa&#45;sequences \
                                &#45;&#45;concatenate&#45;genes \
                                &#45;&#45;return&#45;best&#45;hit
                                &#45;o <span class="artifact&#45;n">[genes&#45;fasta](/software/anvio/help/artifacts/genes&#45;fasta)</span>
</div>


### Want to play?

You can play with this program using the anvi'o data pack for the [infant gut data](/tutorials/infant-gut) and by replacing the parameters above with appropriate ones in the following commands.

Download the latest version of the data from here:

[doi:10.6084/m9.figshare.3502445](https://doi.org/10.6084/m9.figshare.3502445)

Unpack it:

<div class="codeblock" markdown="1">
tar &#45;zxvf INFANTGUTTUTORIAL.tar.gz && cd INFANT&#45;GUT&#45;TUTORIAL
</div>

Import the collection `merens`:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;import&#45;collection](/software/anvio/help/programs/anvi&#45;import&#45;collection)</span> additional&#45;files/collections/merens.txt \
                       &#45;p PROFILE.db \
                       &#45;c CONTIGS.db \
                       &#45;C merens
</div>

Then run the program to,

### Learn available HMM sources

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits &#45;p PROFILE.db \
                                &#45;c CONTIGS.db \
                                &#45;C merens \
                                &#45;o OUTPUT.fa \
                                &#45;&#45;hmm&#45;source Campbell_et_al \
                                &#45;&#45;gene&#45;names Ribosomal_L27,Ribosomal_L28,Ribosomal_L3 \
                                &#45;&#45;return&#45;best&#45;hit \
                                &#45;&#45;get&#45;aa&#45;sequences \
                                &#45;&#45;concatenate
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-sequences-for-hmm-hits.md) to update this information.


## Additional Resources


* [A tutorial on anvi&#x27;o phylogenomics workflow](http://merenlab.org/2017/06/07/phylogenomics/)

* [A detailed application of phylogenomics to place a new genome on a tree](http://merenlab.org/data/parcubacterium-in-hbcfdna/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-sequences-for-hmm-hits) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
