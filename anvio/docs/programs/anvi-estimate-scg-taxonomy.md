This program **quickly makes taxonomy estimates for genomes, metagenomes, or bins stored in your %(contigs-db)s.**

This is the final step in the scg-taxonomy workflow (described in its entirety [here](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/)). Before running this program, you'll need to have run both %(anvi-setup-scg-taxonomy)s and %(anvi-run-scg-taxonomy)s on the %(contigs-db)s you're working with for this project. 

[This tutorial](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal) also includes a comprehensive overview of what this program can do. See that page for more information on all of the features described below. 

Keep in mind that this uses single-copy core genes and their hits in  [GTDB](https://gtdb.ecogenomic.org/), so it will not work well in bins with low completion or for Eukaryotic organisms. 

This program is implicitly run in the interactive interface, when you turn on "Realtime taxonomy estimation for bins (whenever possible)." So, if you've ever wondered where those estimates were coming from, now you know. 

So, what can this program do?

### 1. Estimate the taxonomy of a single genome

By default, this program wll assume your %(contigs-db)s contains only one genome, and will try to use the single-copy core genes (that were associated with taxonomy when you ran %(anvi-run-scg-taxonomy)s) to try to identify the taxonomy of your genome. 

When you run 

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s
{{ codestop }}

It will give you the best taxonomy hit for your genome. If you would like to see how it got there (by looking at the hits for each of the single-copy core genes), just use the `--debug` flag to see more information, as so: 

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --debug 
{{ codestop }}

### 2. Estimate the taxa within a metagenome 

By running this program in metagenome mode, it will assume that your %(contigs-db)s contains multiple genomes and will try to give you an overview of the taxa within it. To do this, it will determine which single-copy core gene has the most hits in your contigs (for example `Ribosomal_S6`), and then will look at the taxnomy hits for that gene across your contigs. The output will be this list of taxonomy results. 

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode 
{{ codestop }}

If you want to look at a specific gene (instead of the one with the most hits), you can also tell it to do that. For example, to tell it to look at Ribosomal_S9, run

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode  \
                           --scg-name Ribosomal_S9
{{ codestop }}

### 3. Look at relative abundance across multiple samples 

If you provide a %(profile-db)s or %(single-profile-db)s, then you'll be able to look at the relative abundance of your taxonomy hits (through a single-copy core gene) across your samples. Essentially, this adds additional columns to your output (one per sample) that descrbe the relative abundance of each hit in each sample. 

Running this will look something like this, 
{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode  \
                           --p %(profile-db)s \
                           --compute-scg-coverages
{{ codestop }}

For an example output, take a look at [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#contigs-db--profile-db).

### 4. Estimate the taxonomy of your bins 

This program basically looks at each of the %(bin)ss in your %(collection)s as a single genome and tries to assign it taxonomy information. To do this, simply provide a collection, like this:

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --C %(collection)s
{{ codestop }}

You can also look at the relative abundances across your samples at the same time, by running something like this: 

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --C %(collection)s  \
                           --p %(profile-db)s \
                           --compute-scg-coverages
{{ codestop }}

### 5. Look at multiple metagenomes at the same time

You can even use this program to look at multiple metagenomes by providing a %(metagenomes)s artifact. This is useful to get an overview of what kinds of taxa might be in your metagenomes, and what kinds of taxa they share. 

Running this

{{ codestart }}
anvi-estimate-scg-taxonomy --metagenomes %(metagenomes)s \
                           --output-file-prefix EXAMPLE
{{ codestop }}

will give you an output file containing all taxonomic levels found and their coverages in each of your metagenomes. 

For a concrete example, check out [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#many-contigs-dbs-for-many-metagenomes). 
