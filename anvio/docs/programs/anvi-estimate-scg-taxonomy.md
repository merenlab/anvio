This program makes **quick taxonomy estimates for genomes, metagenomes, or bins stored in your %(contigs-db)s** using single-copy core genes.

You can run this program on an anvi'o contigs database only if you already have setup the necessary databases to assign taxonomy on your computer by running %(anvi-setup-scg-taxonomy)s and annotated the %(contigs-db)s you are working with using %(anvi-run-scg-taxonomy)s, which are described in greater detail in [this document](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/)), which also offers a [comprehensive overview](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal) of what %(anvi-estimate-scg-taxonomy)s can do.

Keep in mind that the scg-taxonomy framework currently uses single-copy core genes found in [GTDB](https://gtdb.ecogenomic.org/) genomes, thus it will not work well for low-completion, viral, or eukaryotic genomes.

This same functionality %(anvi-estimate-scg-taxonomy)s is implicitly accessed thorugh the anvi'o %(interactive)s interface, when you turn on real-time taxonomy estimation for bins. So, if you've ever wondered where those estimates were coming from, now you know.

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
                           --metagenome-mode \
                           --scg-name Ribosomal_S9
{{ codestop }}

### 3. Look at relative abundance of taxa across samples

If you provide a merged %(profile-db)s or %(single-profile-db)s, then you'll be able to look at the relative abundance of your taxonomy hits (through a single-copy core gene) across your samples. Essentially, this adds additional columns to your output (one per sample) that descrbe the relative abundance of each hit in each sample.

Running this will look something like this,
{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode \
                           -p %(profile-db)s \
                           --compute-scg-coverages
{{ codestop }}

For an example output, take a look at [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#contigs-db--profile-db).

### 4. Estimate the taxonomy of your bins

This program basically looks at each of the %(bin)ss in your %(collection)s as a single genome and tries to assign it taxonomy information. To do this, simply provide a collection, like this:

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           -C %(collection)s
{{ codestop }}

You can also look at the relative abundances across your samples at the same time, by running something like this:

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           -C %(collection)s \
                           -p %(profile-db)s \
                           --compute-scg-coverages
{{ codestop }}

Pro tip: you can use the output that emerges from the following output,

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           -p %(profile-db)s \
                           -C %(collection)s \
                           -o TAXONOMY.txt
{{ codestop }}

to display the taxonomy of your bins in the anvi'o interactive interface in **collection mode**:

{{ codestart }}
%(anvi-interactive)s -c %(contigs-db)s \
                     -p %(profile-db)s \
                     -C %(collection)s \
                     --additional-layers TAXONOMY.txt
{{ codestop }}

That simple.

### 5. Look at multiple metagenomes at the same time

You can even use this program to look at multiple metagenomes by providing a %(metagenomes)s artifact. This is useful to get an overview of what kinds of taxa might be in your metagenomes, and what kinds of taxa they share.

Running this

{{ codestart }}
anvi-estimate-scg-taxonomy --metagenomes %(metagenomes)s \
                           --output-file-prefix EXAMPLE
{{ codestop }}

will give you an output file containing all taxonomic levels found and their coverages in each of your metagenomes.

For a concrete example, check out [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#many-contigs-dbs-for-many-metagenomes).
