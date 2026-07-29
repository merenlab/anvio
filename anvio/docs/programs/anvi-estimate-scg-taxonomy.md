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

Without a %(profile-db)s, the output will report the **frequency** of each taxon — i.e., how many times that taxon was detected across the SCG hits in your contigs database. If you only care whether a taxon is present or absent rather than how many times it was detected, you can use the `--presence-absence-only` flag to get a binary report instead:

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode \
                           --presence-absence-only
{{ codestop }}

### 3. Look at the taxonomic composition as a tree

The table above answers the question "what is the taxonomy of each of these hits?", which is not the same question as "what is in here?". Since each row spells out an entire lineage, a taxon that occurs 30 times occupies 30 nearly identical rows, and the shape of the community is nowhere to be seen. If that is what you are after, add the flag `--tree-output`:

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode \
                           --tree-output
{{ codestop }}

which will show you the very same results as a hierarchical tree instead (please remember that the tree here has no phylogenetic order or meaning in its current implementation):

```
All Ribosomal_S11 copies (96)
├── Bacteria (94)
│   ├── Campylobacterota (36)
│   │   └── Campylobacteria (36)
│   │       └── Campylobacterales (36)
│   ├── Desulfobacterota (16)
│   │   ├── Desulfobulbia (9)
│   │   │   └── Desulfobulbales (9)
│   │   ├── Desulfobacteria (3)
│   │   │   └── Desulfobacterales (3)
(...)
└── Archaea (2)
    └── Halobacteriota (2)
        └── Methanosarcinia (2)
            └── Methanosarcinales (2)
```

The number next to each taxon is the number of things that were assigned to it **or to anything under it**, which is why the numbers of a node's children always add up to the number of the node itself. What is being counted depends on how you ran the program: copies of the single-copy core gene anvi'o surveyed in metagenome mode (`Ribosomal_S11` in the example above), bins if you gave anvi'o a %(collection)s, and genomes otherwise. The root of the tree spells out which one it is, so you never have to guess.

Hits that could not be resolved all the way down show up in an explicit `Unknown_*` node at the level where they ran out of names (as in `Unknown_genera`) rather than being quietly dropped, so nothing goes missing from the tree.

By default the tree does not go deeper than genus names, since species-level assignments from single-copy core genes are often absent or not particularly trustworthy. If you want a different cutoff, use `--tree-output-level`:

{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode \
                           --tree-output \
                           --tree-output-level t_family
{{ codestop }}

{:.notice}
Please note that `--tree-output-level` has nothing to do with the `--taxonomic-level` parameter, and only influences the tree that is displayed in your terminal.

Two small things to know about this flag. First, it only changes what is displayed: if you also ask for an output file with `-o`, that file will contain the usual TAB-delimited table, exactly as it would have without `--tree-output`. Second, since a tree is the only thing this flag produces, anvi'o will refuse to run it together with `--quiet` (which would show you nothing) or `--as-markdown` (which would mangle the characters anvi'o uses to draw the tree).

If you are working with a %(metagenomes)s or %(external-genomes)s file, `--tree-output` will give you a single tree in which the results from every one of your inputs are merged together. This is not really the appropriate display for many (meta)genomes at once, since a tree gives you no way to tell which input a taxon came from, and anvi'o will warn you about that -- but if a bird's eye view of everything at once is what you want, there it is.

### 4. Look at relative abundance of taxa across samples

If you provide a merged %(profile-db)s or %(single-profile-db)s, then you'll be able to look at the relative abundance of your taxonomy hits (through a single-copy core gene) across your samples. Essentially, this adds additional columns to your output (one per sample) that descrbe the relative abundance of each hit in each sample.

Running this will look something like this,
{{ codestart }}
anvi-estimate-scg-taxonomy -c %(contigs-db)s \
                           --metagenome-mode \
                           -p %(profile-db)s \
                           --compute-scg-coverages
{{ codestop }}

For an example output, take a look at [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#contigs-db--profile-db).

If you add `--tree-output` to a command like the one above, each node of the tree will report a coverage value in addition to its count. That value is the total coverage of everything under that node, summed across all of your samples. A tree has no room for one column per sample, so if you need per-sample numbers, ask for an output file with `-o`.

### 5. Estimate the taxonomy of your bins

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

### 6. Look at multiple metagenomes at the same time

You can even use this program to look at multiple metagenomes by providing a %(metagenomes)s artifact. This is useful to get an overview of what kinds of taxa might be in your metagenomes, and what kinds of taxa they share.

Running this

{{ codestart }}
anvi-estimate-scg-taxonomy --metagenomes %(metagenomes)s \
                           --output-file-prefix EXAMPLE
{{ codestop }}

will give you an output file containing all taxonomic levels found and their coverages in each of your metagenomes.

For a concrete example, check out [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#many-contigs-dbs-for-many-metagenomes).

### 7. Estimate taxonomy across multiple genomes using an external genomes file

You can also run this program on a set of genomes described in an %(external-genomes)s file:

{{ codestart }}
anvi-estimate-scg-taxonomy --external-genomes %(external-genomes)s \
                           --output-file-prefix EXAMPLE
{{ codestop }}

If you want to treat each genome as a metagenome (i.e., report the taxonomic composition within each genome rather than a single consensus taxonomy for it), add the `--metagenome-mode` flag:

{{ codestart }}
anvi-estimate-scg-taxonomy --external-genomes %(external-genomes)s \
                           --metagenome-mode \
                           --output-file-prefix EXAMPLE
{{ codestop }}

### A note on SCG selection when working with multiple contigs databases

When you use `--external-genomes` or `--metagenomes` together with `--metagenome-mode`, anvi'o must pick a **single SCG to use consistently across all contigs databases**. This is critical for result comparability: using different SCGs for different databases would make the outputs impossible to compare.

If you do not explicitly name an SCG with `--scg-name-for-metagenome-mode`, anvi'o will automatically select the SCG that is most frequent across all your contigs databases combined, and will report which one it chose. You can inspect per-SCG frequencies beforehand with `--report-scg-frequencies` to make an informed choice.

If instead you want to use the most frequent SCG independently for each contigs database (which would make results incomparable across databases), you should run %(anvi-estimate-scg-taxonomy)s on each contigs database separately in `--metagenome-mode`.
