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

### Example output

When the script is finished running, you should see some text like the following on your terminal screen:
```
CLASSIFICATION RESULT
===============================================
Your genome has now been classified as an LMI genome. It had a score of
6.5214285714285705, where the score is simply the sum of all pathwise
completeness scores of the modules in the list you provided.
```
Hopefully it is clear that the classification (LMI or HMI) and the score will be different from genome to genome.

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

What a great question! The short but not very helpful answer is that you are free to define HMI and choose your metabolic pathways in a way that fits your research question. But the longer and more practical answer is: here we will describe how we selected our set of pathways, and you are welcome to use our set, or to use our method for selecting a set of modules using your own data.

It's important to remember that you will need some justification for how you select modules, which will depend on your knowledge of your dataset. So maybe spend some time thinking about what high metabolic independence means in your research context before you start.

**Here is what we did.** In our study of colonization after FMT, we had two groups of microbial genomes: one group (the good colonizers) which we realized had much higher metabolic capacity than the other group (the poor colonizers). To figure out exactly which metabolic pathways were differentiating between the two groups, we ran %(anvi-compute-metabolic-enrichment)s on the metabolism data from these genomes (obtained via %(anvi-estimate-metabolism)s) and selected the KEGG modules that were enriched in the first group (the good colonizers). Specifically, in our case this translates to filtering the enrichment output for modules with 1) a q-value less than 0.05, 2) the good colonizer group as its 'associated group', and 3) at least 75%% completeness in at least 50%% of the group members. This resulted in the following list of 33 KEGG modules:

```
M00049
M00050
M00007
M00140
M00005
M00083
M00120
M00854
M00527
M00096
M00048
M00855
M00022
M00844
M00051
M00082
M00157
M00026
M00526
M00015
M00019
M00432
M00018
M00570
M00126
M00115
M00028
M00924
M00122
M00125
M00023
M00631
M00061
```

You are welcome to copy and paste this list directly into a file and use it as your `--module-list` input, if you want to use the same set of modules that we did. Note that the definitions for these modules are coming from the %(kegg-data)s version with a hash of `45b7cc2e4fdc`, which is a static snapshot of KEGG from December 2020 (`v2020-12-23`), so you will need to provide this same version to the `--kegg-data-dir` parameter for things to work correctly.

**But I want to do it myself!** If you instead want to use a similar method for obtaining a _de novo_ set of modules from your own data, here are the key components:

1. You need two groups of genomes, one that is made up of HMI genomes and one that is made up of LMI genomes.
2. You should estimate the completeness of all metabolic pathways of potential interest in these genomes.
3. Finally, you need some way of determining which modules are (mostly) present in high completeness in the HMI group, but largely not complete in the LMI group. This is the list you would provide to the script for classifying additional genomes as HMI or LMI.

The exact methods of estimating metabolism and computing enrichment are up to you, but note that this strategy only works if you already have genomes that you have manually identified as HMI or LMI.

### How do I select a threshold score?

This question has a similar answer to the previous one, in that we can tell you how we did it for our study, but it doesn't mean everyone has to do it that way. Science is fun like that!

The key, again, is to already have some idea of what HMI means, in a given genome. This usually means that you have some test genomes that you've manually identified as HMI and LMI, and you know which metabolic pathways represent high metabolic independence.

Remember that our very basic script will add up all of the completeness scores for each metabolic module in your list. Then all you need to do is find a number for this sum that distinguishes between your HMI genomes and your LMI ones.

**Here is what we did.** We had our two groups of genomes, good and poor colonizers. We knew that all of the good colonizers were HMI genomes, and that most (not all) of the poor colonizers were LMI genomes. So what we did was to add up the completeness scores of our 33 metabolic pathways (which, you may remember from the previous section, were enriched in the good colonizer group) and plot these sums for each genome, separated by group. This is the plot we got:

![A plot of HMI scores for the good colonizer genomes and the poor colonizer genomes](../../images/FMT_HMI_score_plot.png)

Here, the 'HIG' label corresponds to the group with all HMI genomes, and the 'LOW' label corresponds to the other group (with mostly LMI genomes). We used this plot to choose a threshold score of 20, which does a good job of distinguishing between the 'HIG' group (all of these genomes have scores above this value) and the 'LOW' group (most of these genomes have scores below this value). That's all!

**A very general rule of thumb.** Since individual completeness scores have a maximum of 1 (for 100%% complete), the maximum value of the sum will be _n_, where _n_ is the number of metabolic pathways in your list. So you should be selecting a threshold between 0 and _n_, but most likely on the higher end of that range, since high metabolic independence is generally defined as _high_ completeness scores across the set of pathways). It will depend on your data, of course.
