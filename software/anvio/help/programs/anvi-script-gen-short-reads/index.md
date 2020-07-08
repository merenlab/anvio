---
layout: page
title: anvi-script-gen-short-reads [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate short reads from contigs. Useful to reconstruct mock data sets from already assembled contigs.

See **[program help menu](../../../vignette#anvi-script-gen-short-reads)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[configuration-ini](../../artifacts/configuration-ini)</span></p>

## Usage


This program uses already assembled contigs to create a mock list of short reads. You can then use these short reads to reassemble your data in order to test alternative reassembly programs or analysis methods as a positive control. 

Basically, this attempts to undo the assembly and produce a data set that could have been directly received from laboratory sequencing. While the computer's mock short reads won't be perfect, they can be used to make sure your analysis pipeline is working from step 1. 

## Example Usage

This program takes an INI file - a form of text file containing various information. For this program, the example provided in the anvi'o test suite looks like this: 

```ini
[general]
short_read_length = 10
error_rate = 0.05
coverage = 100
contig = CTGTGGTTACGCCACCTTGAGAGATATTAGTCGCGTATTGCATCCGTGCCGACAAATTGCCCAACGCATCGTTCCTTCTCCTAAGTAATTTAACATGCGT
```

Note that this file contains both the contig that you want to break down, and various information about the short reads that you want to create. To run this program, just call 

    anvi-script-gen-short-reads <span class="artifact-n">[configuration-ini](/software/anvio/help/artifacts/configuration-ini)</span>
    --output-file-path <span class="artifact-n">[short-reads-fasta](/software/anvio/help/artifacts/short-reads-fasta)</span>
    
The resulting FASTA file with short reads will cover the `contig` with short reads that are 10 nts long at 100X coverage. There will also be an error-rateof 0.05, to mimic the sequencing errors you would get from sequencing in the wet lab. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-short-reads.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-short-reads) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
