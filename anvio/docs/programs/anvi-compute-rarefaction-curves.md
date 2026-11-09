The program %(anvi-compute-rarefaction-curves)s goes through all genomes in a given %(pan-db)s and calculates rarefaction curves for all gene clusters and core gene clusters. It also computes the [Heaps' Law](https://en.wikipedia.org/wiki/Heaps'_law) fit to model the relationship between genome sampling and the number of new gene clusters discovered for you to have a more comprehensive reporting of your pangenome.

It also works on a %(pan-graph-db)s, in which case the thing that accumulates as genomes are added is the *synteny-aware gene cluster* (SynGC) rather than the gene cluster. See [Rarefaction curves for a pangenome graph](#rarefaction-curves-for-a-pangenome-graph) below.

### On the utility of rarefaction curves and Heaps' Law fit

Rarefaction curves are helpful in the analysis of pangenome as they help visualize the *discovery rate of new gene clusters* as a function of increasing number of genomes. While a steep curve suggests that many new gene clusters are still being discovered, indicating incomplete coverage of the potential gene cluster space, a curve that reaches a plateau suggests sufficient sampling of gene cluster diversity.

However, rarefaction curves have inherent limitations. Because genome sampling is often biased and unlikely to fully capture the true genetic diversity of any taxon, rarefaction analysis provides only dataset-specific insights. Despite these limitations, rarefaction curves remain a popular tool for characterizing whether a pangenome is relatively 'open' (with continuous gene discovery) or 'closed' (where new genome additions contribute few or no new gene clusters). As long as you take such numerical summaries with a huge grain of salt, it is all fine.

Fitting Heaps' Law to the rarefaction curve provides a quantitative measure of pangenome openness. The *alpha* value derived from Heaps' Law (sometimes referred to as *gamma* in the literature) reflects how the number of new gene clusters scales with increasing genome sampling. There is no science to define an absolute threshold for an open or a closed pangenome. However, pangenomes with alpha values below 0.3 tend to be relatively closed, and those above 0.3 tend to be relatively open. Higher alpha values will indicate increasingly open pangenomes and lower values will identify progressively closed ones.

### Running the program

The simplest for of the command will look like this,

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s
{{ codestop }}

Running this program on a %(pan-db)s will,

* Report the Heaps' Law alpha value in its output,
* Generate an SVG file for the visualization of the rarefaction curves for the whole pangenome and core gene clusters with all the necessary information embedded in the figure,
* Generate four additional text files that represent the exact data used to visualize rarefaction curves (both averages and iterations for GC gain per genome calculations for the whole pangenome and for the core genome).

The program will use the 'project name' information stored in the %(pan-db)s as a 'prefix' to all resulting files, and the output will look like this:

```
Num genomes found ............................: 5
Num iterations to run ........................: 100
Output file prefix ...........................: TEST
Heaps' Law parameters estimated ..............: K=245.3049, alpha=0.2484

OUTPUT FILES
===============================================
Rarefaction curves ...........................: TEST-rarefaction-curves.svg
GC gain per genome for core (averages) .......: TEST-rarefaction-core-averages.txt
GC gain per genome for core (each iteration) .: TEST-rarefaction-core-iterations.txt
GC gain per genome for all (averages) ........: TEST-rarefaction-pangenome-averages.txt
GC gain per genome for all (each iteration) ..: TEST-rarefaction-pangenome-iterations.txt

```

You can change the prefix to something else,

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --output-file-prefix MY_PREFIX
{{ codestop }}

or you can ask the program to *not* generate any output files and simply report the Heaps' Law parameters:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --skip-output-files
{{ codestop }}

You can also determine the number of random sampling to be conducted through the `--iteration` parameter. The default is 100. Going above this value will unlikely refine the results, but going below 10 will have a negative influence since the fit will be affected by small amount of sampling:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --iterations 50
{{ codestop }}

### How about Pangeome Graphs?

Everything above works the same way for a %(pan-graph-db)s as `-p` takes either kind of database:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-graph-db)s
{{ codestop }}

The only difference is *what is being counted*. A node in a pangenome graph is a **synteny-aware gene cluster** (SynGC): a cluster that contains at most one gene per genome, and whose members were found in the same place relative to their neighbors. So the curves describe the discovery rate of new SynGCs rather than of new gene clusters, and the terminal output and figure say `SynGC` where they would otherwise say gene cluster:

```
Input database type ..........................: pan-graph
Number of genomes found ......................: 5
Number of SynGCs found .......................: 30
Number of iterations to run ..................: 100
Heaps' Law parameters estimated ..............: K=18.6778, alpha=0.2958

OUTPUT FILES
===================================================
Rarefaction curves ...............................: TEST-rarefaction-curves.svg
SynGC gain per genome for core (averages) ........: TEST-rarefaction-core-averages.txt
SynGC gain per genome for core (each iteration) ..: TEST-rarefaction-core-iterations.txt
SynGC gain per genome for all (averages) .........: TEST-rarefaction-pangenome-averages.txt
SynGC gain per genome for all (each iteration) ...: TEST-rarefaction-pangenome-iterations.txt
```

{:.warning}
How the calcualtions of rarefaction curves for GC sand SynGCs differ, and how these differences affect the interpretations of Heaps' fits is an important intellectual question for which we offer no insights here (but ask you to think about for your own dataset). A single gene cluster of the source pangenome can be split into several SynGCs by %(anvi-pan-genome-graph)s, since genes that cluster together by sequence similarity but occur in different syntenic contexts end up in different nodes. A %(pan-graph-db)s curve therefore will almost always sit *above* the %(pan-db)s curve it was derived from, and its alpha is a different number describing a different quantity. What does it mean, then? Well, we have our own opinions, but we would like to not contaminate your thinking with them.

The %(rarefaction-curves)s text files use the same column names (`avg_num_gene_clusters`, `GeneClusters`) whichever database they came from. This is deliberate, so that anything parsing these files does not have to care which kind of pangenome produced them, but it does mean the column header will say `gene_clusters` while holding SynGC counts.
