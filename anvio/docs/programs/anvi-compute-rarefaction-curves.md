The program %(anvi-compute-rarefaction-curves)s processes all genomes in a given %(pan-db)s and calculates rarefaction curves for all gene clusters and core gene clusters. It also computes the [Heaps' Law](https://en.wikipedia.org/wiki/Heaps'_law) fit to model the relationship between genome sampling and the discovery of new gene clusters, providing comprehensive reporting of your pangenome characteristics.

### On the utility of rarefaction curves and Heaps' Law fit

Rarefaction curves are valuable tools in pangenome analysis as they help visualize the *discovery rate of new gene clusters* as a function of increasing numbers of genomes. A steep curve suggests that many new gene clusters are still being discovered, indicating incomplete coverage of the potential gene cluster space, while a curve that reaches a plateau suggests sufficient sampling of gene cluster diversity.

However, rarefaction curves have inherent limitations. Because genome sampling is often biased and unlikely to fully capture the true genetic diversity of any taxon, rarefaction analysis provides only dataset-specific insights. Despite these limitations, rarefaction curves remain a popular tool for characterizing whether a pangenome is relatively 'open' (with continuous gene discovery) or 'closed' (where new genome additions contribute few or no novel gene clusters). Such numerical summaries should be interpreted with appropriate caution.

Fitting Heaps' Law to the rarefaction curve provides a quantitative measure of pangenome openness. The *alpha* value derived from Heaps' Law (sometimes referred to as *gamma* in the literature) reflects how the number of new gene clusters scales with increasing genome sampling. There is no scientific consensus defining an absolute threshold for distinguishing open versus closed pangenomes. However, pangenomes with alpha values below 0.3 tend to be relatively closed, while those above 0.3 tend to be relatively open. Higher alpha values indicate increasingly open pangenomes, and lower values identify progressively closed ones.

### Running the program

The simplest form of the command will look like this:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s
{{ codestop }}

Executing this program on a %(pan-db)s will:

* Report the Heaps' Law alpha value in its output,
* Generate an SVG file for visualizing the rarefaction curves for the complete pangenome and core gene clusters with all necessary information embedded in the figure,
* Generate four additional text files that represent the exact data used to visualize rarefaction curves (both averages and iterations for gene cluster gain per genome calculations for the whole pangenome and for the core genome).

The program will use the 'project name' information stored in the %(pan-db)s as a 'prefix' for all resulting files, and the output will appear as follows:

```
Number of genomes found ......................: 5
Number of iterations to run ..................: 100
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

You can modify the prefix to a custom value:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --output-file-prefix MY_PREFIX
{{ codestop }}

or you can instruct the program to *not* generate any output files and simply report the Heaps' Law parameters:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --skip-output-files
{{ codestop }}

You can also specify the number of random sampling iterations to be conducted using the `--iteration` parameter. The default is 100. Increasing this value beyond 100 is unlikely to significantly refine the results, but reducing it below 10 will negatively impact the analysis since the fit will be affected by insufficient sampling:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --iterations 50
{{ codestop }}