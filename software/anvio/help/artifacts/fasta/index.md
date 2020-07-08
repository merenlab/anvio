---
layout: post
title: fasta [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span></p>


## Required or used by

<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span> <span class="artifact-r">[anvi-script-compute-ani-for-fasta](../../programs/anvi-script-compute-ani-for-fasta)</span> <span class="artifact-r">[anvi-script-reformat-fasta](../../programs/anvi-script-reformat-fasta)</span></p>

## Description

A FASTA-formatted file that does not necessarily meet the standards of a <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span>.

<span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/programs/anvi-script-reformat-fasta)</span> can turn a regular fasta into a <span class="artifact-n">[contigs-fasta](/software/anvio/help/artifacts/contigs-fasta)</span>, which anvi'o will be able to utilize better.

### But what is a FASTA file? 

A FASTA file contains sequences (in this case, nucleotide sequences, though they can also describe peptide sequences) that are formatted as follows: 

    >SEQUENCE_ID VARIOUS_SEQUENCE_DATA
    SEQUENCE
    
The `VARIOUS_SEQUENCE_DATA` region can contain data such as the NCBI taxon ID, gi assession number, a text description of the sequence, or the start and end positions if the sequence is a portion of a larger sample. All of this information is optional. 

The sequence itself is written in standard IUPAC format (though it can be written in lower-case letters).  

For a concrete example, you can download sequences from the NCBI database in FASTA format. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/fasta.md) to update this information.

