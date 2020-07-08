---
layout: page
title: anvi-estimate-metabolism [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Reconstructs metabolic pathways and estimates pathway completeness for a given set of contigs.

See **[program help menu](../../../vignette#anvi-estimate-metabolism)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[kegg-metabolism](../../artifacts/kegg-metabolism)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[kegg-db](../../artifacts/kegg-db)</span> <span class="artifact-r">[kegg-functions](../../artifacts/kegg-functions)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[collection](../../artifacts/collection)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes)</span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes)</span> <span class="artifact-r">[metagenomes](../../artifacts/metagenomes)</span></p>

## Usage


<span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span> predicts the metabolic capabilities of organisms based on their genetic content. It relies upon <span class="artifact-n">[kegg-functions](/software/anvio/help/artifacts/kegg-functions)</span> and metabolism information from the KEGG resource, which is stored in <span class="artifact-n">[kegg-db](/software/anvio/help/artifacts/kegg-db)</span>.

The metabolic pathways that this program currently considers are those defined by KOs in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). Each KO represents a gene function, and a KEGG module is a set of KOs that collectively carry out the steps in a metabolic pathway. Therefore, for this to work, you need to have annotated your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> with hits to the KEGG KOfam database by running <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/programs/anvi-run-kegg-kofams)</span> prior to using this program.

Given a properly annotated <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>, this program determines which KOs are present and from those determines the completeness of each KEGG module. The results are described in a set of output text files, collectively referred to as <span class="artifact-n">[kegg-metabolism](/software/anvio/help/artifacts/kegg-metabolism)</span>.

## Running metabolism estimation on a single contigs database

There are several possible inputs to this program. For single genomes - isolate genomes or MAGs, for example - you can provide a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. If your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> describes a metagenome rather than a single genome, you can provide the flag `--metagenome-mode`. In metagenome mode, KOfam hits in the <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> are analyzed as though they belong to one collective genome, despite the fact that the sequences represent multiple different populations. Alternatively, if you have binned your metagenome sequences into separate populations and would like metabolism estimation to be run separately on each bin, you can provide a <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> and a <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span>.

### Estimation for a single genome

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db
</div>

### Estimation for a metagenome

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;metagenome&#45;mode
</div>

### Estimation for bins in a metagenome

You can estimate metabolism for each bin in a <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span>:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;p PROFILE.db &#45;C COLLECTION_NAME
</div>

You can also provide a specific <span class="artifact-n">[bin](/software/anvio/help/artifacts/bin)</span> in that <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span> to run on:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;p PROFILE.db &#45;C COLLECTION_NAME &#45;b BIN_NAME
</div>

Or, you can provide a specific list of bins in a text file:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;p PROFILE.db &#45;C COLLECTION_NAME &#45;B bin_ids.txt
</div>

Each line in the `bin_ids.txt` file should be a bin name from the <span class="artifact-n">[collection](/software/anvio/help/artifacts/collection)</span>.

## Running metabolism estimation on multiple contigs databases

If you have a set of contigs databases of the same type (ie, all of them are single genomes, or all are binned metagenomes), you can analyze them all at once. What you need to do is put the relevant information for each <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> into a text file and pass that text file to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>. The program will then run estimation individually on each contigs database in the file. The estimation results for each database will be aggregated and printed to the same output file(s).

### Estimation for multiple single genomes

Multiple single genomes (also known as <span class="artifact-n">[external-genomes](/software/anvio/help/artifacts/external-genomes)</span>) can be analyzed with the same command by providing an external genomes file to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>. To see the required format for the external genomes file, see <span class="artifact-n">[external-genomes](/software/anvio/help/artifacts/external-genomes)</span>.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;e external&#45;genomes.txt
</div>

### Estimation for multiple metagenomes

Multiple metagenomes can be analyzed with the same command by providing a metagenomes input file. Metagenome mode will be used to analyze each contigs database in the file. To see the required format for the external genomes file, see <span class="artifact-n">[metagenomes](/software/anvio/help/artifacts/metagenomes)</span>.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;M metagenomes.txt
</div>

### Estimation for multiple bins in different metagenomes

If you have multiple bins (also known as <span class="artifact-n">[internal-genomes](/software/anvio/help/artifacts/internal-genomes)</span>) across different collections or even different metagenomes, they can be analyzed with the same command by providing an internal genomes file to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span>. To see the required format for the external genomes file, see <span class="artifact-n">[internal-genomes](/software/anvio/help/artifacts/internal-genomes)</span>.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt
</div>

## Adjusting module completion threshold

KEGG module completeness is computed as the percentage of steps in the metabolic pathway that are 'present' based on the KOs found in the contigs database. If this completeness is larger than a certain percentage, then the entire module is considered to be 'present' in the genome or metagenome. By default, this module completion threshold is 0.75; that is, 75 percent of the KOs in a module must have a KOfam hit in the contigs database in order for the module to be considered 'complete' as a whole. This threshold can be adjusted.

### Changing the module completion threshold

In this example, we change the threshold to 50 percent.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;module&#45;completion&#45;threshold 0.5
</div>

## Controlling output

<span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/programs/anvi-estimate-metabolism)</span> can produce a variety of output files. All will be prefixed with the same string, which by default is "kegg-metabolism".

### Changing the output file prefix

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;O my&#45;cool&#45;prefix
</div>


This program has two major output options - long format (tab-delimited) output files and matrices.

Long format output has several preset "modes" as well as a "custom" mode in which the user can define the contents of the output file. Multiple modes can be used at once, and each requested "mode" will result in a separate output file. The default output mode is "modules" mode.

You can find more details on the output format by looking at <span class="artifact-n">[kegg-metabolism](/software/anvio/help/artifacts/kegg-metabolism)</span>.

### Viewing available output modes

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;list&#45;available&#45;modes
</div>

### Using a non-default output mode

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes kofam_hits
</div>

### Using multiple output modes

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes kofam_hits,modules
</div>

### Viewing available output headers for 'custom' mode

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;list&#45;available&#45;output&#45;headers
</div>

### Using custom output mode

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes custom &#45;&#45;custom&#45;output&#45;headers kegg_module,module_name,module_is_complete
</div>


Matrix format is only available when working with multiple contigs databases. Several output matrices will be generated, each of which describes one KEGG module statistic such as completion score or presence/absence.  

### Obtaining matrix-formatted output

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt &#45;&#45;matrix&#45;format
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-estimate-metabolism.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-estimate-metabolism) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
