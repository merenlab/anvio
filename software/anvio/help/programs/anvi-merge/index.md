---
layout: page
title: anvi-merge [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Merge multiple anvio profiles.

See **[program help menu](../../../vignette#anvi-merge)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-p">[misc-data-item-orders](../../artifacts/misc-data-item-orders)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[single-profile-db](../../artifacts/single-profile-db)</span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span></p>

## Usage


The main function of anvi-merge is to convert multiple <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>s into a single <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span> (also called a merged profile database). Basically, this takes the alignment data from each sample (each contained in its own <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>) and combines them into a single database that anvi'o can look through more easily. 

### Overview: How to run anvi-merge
1. Set up your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>. See that page for more information
1. Use <span class="artifact-n">[anvi-profile](/software/anvio/help/programs/anvi-profile)</span> to create a <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span> using each of your samples (formatted into <span class="artifact-n">[bam-file](/software/anvio/help/artifacts/bam-file)</span>) *(Note: for each of these you'll need to use the same <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>)*
1. Use anvi-merge to merge those <span class="artifact-n">[single-profile-db](/software/anvio/help/artifacts/single-profile-db)</span>s into a single database, called a <span class="artifact-n">[profile-db](/software/anvio/help/artifacts/profile-db)</span>. This will look something like the following:

        anvi-merge -c cool_contigs.db \
                    Single_profile_db_1 Single_profile_db_2 \
                    -o cool_contigs_merge
                    
    This will put all of the output files (the final profile-db as well as a <span class="artifact-n">[misc-data-item-orders](/software/anvio/help/artifacts/misc-data-item-orders)</span> containing various other information) into the folder `cool_contigs_merge `.
        

## Other Parameters

You must give anvi-merge your contigs database and single profile databases. However, you can also provide more information or give addtional instructions.

### Hierarchical Clustering 
#### To run or not to run? 
* Use the flag `--skip-hierarchical-clustering` to turn hierarchical clustering off. This will save on computation time, but will skip out on creating the tree of contigs at the center of the interactive interface. If you have more than 25,000 splits in the final profile, this will be set automatically. 
* Use the flag `--enforce-hierarchical-clustering` to turn hierarchical clustering back on. This will take a long time, but will produce a lovely contigs tree for the interactive interface. 
#### Additional parameters
* Provide a custom distance metric for clustering using the flag `--distance.` (The default is euclidean)
* Provide a custom linkage method for clustering using the flag `--linkage.` (The default is ward)
### Providing additional information
* Provide the sample name using the flag `-S`. If you don't, anvi'o will come up with one, but it probably won't be any good. 
* Provide a text file in markdown to describe the project using the flag `--description`. This will show up when you later use the interactive interface to analyze your profile-db. 
### Output Information
* Provide an output destination with the flag `-o`.
* Add the flag `-W` to overwrite existing files in that directory. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-merge.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-merge) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.
