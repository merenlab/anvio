The goal of this script is to identify whether a given genome has high metabolic independence (HMI) or not. We developed it in response to findings from [our study on fecal microbiota transplant (which is currently a pre-print and in review)](https://doi.org/10.1101/2021.03.02.433653), because we needed a quantitative way of characterizing the concept of HMI and generalizing it to new datasets.

We expect that users interested in this script will mostly be those seeking to reproduce our work or to see how HMI applies to their collection of genomes. But we know that anvi'o users are quite creative in applying our tools in ways beyond the developers' limited imaginations, so feel free to use this script in the context of other research questions if you think it will be appropriate. :)

## HMI? What? Please explain.

To summarize, a microbe has high metabolic independence (HMI) when its genome encodes, with high completeness, many metabolic pathways for biosynthesis of key molecules such as amino acids, cofactors, nucleotides, lipids, etc -- and this makes it fairly robust to stresses in its environment that can disrupt microbial communities.  In contrast, low metabolic independence (LMI) is characterized by the complete absence (or low completeness) of these pathways, meaning that an LMI microbe can only access these critical metabolites (which are often necessary for survival) if they are produced by other microbes in the community.

But which specific metabolic pathways are needed to have high metabolic independence? And how complete do these pathways need to be? Well, nobody _really_ knows yet as this phenomenon hasn't yet been studied in enough detail. So we cannot pretend to have exactly the right answers to these questions (if there is even just one "right" answer). But we do have some insights into these questions from our recent work, and you can see the relevant sections below to read about how we select metabolic pathways, the set that we use, and how we choose a scoring threshold for classification. But it is important to remember that you are not limited to using the pathways and threshold(s) that we use - this script accepts this information as parameters, so you can choose what to use based on your own insights and/or hypotheses.

For a bit of additional (and possibly helpful) context, the observation leading to this concept of HMI/LMI occurred when we were studying colonization success after fecal microbiota transplant (FMT) in people with recurrent _C. dificile_ infections. The microbial populations that successfully colonized (and survived long-term) in the FMT recipients in our study all had relatively larger genomes, and were enriched in a set of 33 mostly-biosynthetic metabolic pathways, when we compared them to microbial populations that failed to colonize. We then looked at the completeness of these pathways in metagenome-assembled genomes (MAGs) from a variety of publicly-available gut metagenomes, and saw that the pathways were more complete in the MAGs obtained from people with gastrointestinal stresses (like Crohn's disease, ulcerative colitis, and pouchitis). From this information, the concept of HMI emerged as an ecological phenomenon that supports microbial survival in stressful environments. We invite interested readers to consider [our pre-print](https://doi.org/10.1101/2021.03.02.433653) for more details.

## What does the script do?

It is very basic (at least in its current iteration). You give the script a genome fasta file, a set of metabolic pathways (as [KEGG MODULE](https://www.genome.jp/kegg/module.html) accessions), and a threshold. The script will tell you 1) whether it is an 'HMI' genome or an 'LMI' genome, and 2) the genome's 'HMI score'. To do this, it will go though the following steps:

1. Create a %(contigs-db)s for the genome, using the programs %(anvi-script-reformat-fasta)s followed by %(anvi-gen-contigs-database)s.
2. Annotate the gene calls in this genome with KEGG KOfams using the program %(anvi-run-kegg-kofams)s and the KEGG database (%(kegg-data)s) that you provided.
3. Use %(anvi-estimate-metabolism)s to estimate the completeness of each metabolic pathway in the provided list. This program provides a fractional completeness score for each pathway, which will have a value between 0 and 1 (inclusive). Note that this step also depends on the provided KEGG database, and that we currently use the pathwise completeness score for this step (you can find [an explanation of pathwise completeness here](https://anvio.org/help/main/programs/anvi-estimate-metabolism/#technical-details)).
4. Add up all the completeness scores to get one sum (this is the 'HMI score').
5. Classify the genome: if the HMI score is greater than or equal to the provided threshold, the genome is labeled as 'HMI'. Otherwise, it is classified as 'LMI'.

The most important things to remember is that 'HMI score' is just a sum of completeness scores, and the threshold is just a number representing the lower boundary of 'high' completeness. The last two steps are the only things that are unique to this particular script, which simply calls other anvi'o programs to do the rest of the work. If you want to know more details about the other steps, you can check the individual program pages that are linked above.

## How do I use the script?

Like this:
```
anvi-script-classify-hmi-genomes -f %(fasta)s \
                                 --module-list modules.txt \
                                 --threshold 20 \
                                 --kegg-data-dir %(kegg-data)s \
```

Above, the `modules.txt` is a file in which each line contains a different module accession number, as defined by KEGG. Note that these modules need to be present in the KEGG database that you pass to the `--kegg-data-dir` parameter, so make sure you have a good version of this database (see %(anvi-setup-kegg-kofams)s for help with setting up different versions of KEGG, including the option to [download the most up-to-date data directly from the KEGG website](https://anvio.org/help/main/programs/anvi-setup-kegg-kofams/#getting-the-most-up-to-date-kegg-data-downloading-directly-from-kegg)).  

The threshold is a number that you set to differentiate between HMI and LMI genomes. It will be applied to the sum of completeness scores of the metabolic pathways in the `modules.txt` file, so you should set its value appropriately according to how many pathways are in your list, and how complete you expect each of those pathways to be in HMI genomes. More advice on this below.

Note that the script can take a while to run. Gene annotation with KOfams is the bottleneck step. To make it run as fast as possible, we recommend increasing the number of threads using the `-T` parameter (within the reasonable limits of your computer system).

### How do I select which metabolic modules to use?

What a great question! 

### How do I select a threshold score?
