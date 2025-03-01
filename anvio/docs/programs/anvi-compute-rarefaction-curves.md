The program %(anvi-compute-rarefaction-curves)s goes through all genomes in a given %(pan-db)s and calculates rarefaction curves for all gene clusters and core gene clusters. It also computes the Heaps' Law fit to model the relationship between genome sampling and the number of new gene clusters discovered for you to have a more comprehensive reporting of your pangenome.

### On the utility of rarefaction curves and Heaps' Law fit

You can use the program {% include PROGRAM name="anvi-compute-rarefaction-curves" %} to calculate and report rarefaction curves and Heaps' Law fit for your pangenome. As described in the program help page, rarefaction curves are helpful in the analysis of pangenome as they help visualize the *discovery rate of new gene clusters* as a function of increasing number of genomes. While a steep curve suggests that many new gene clusters are still being discovered, indicating incomplete coverage of the potential gene cluster space, a curve that reaches a plateau suggests sufficient sampling of gene cluster diversity.

However, rarefaction curves have inherent limitations. Because genome sampling is often biased and unlikely to fully capture the true genetic diversity of any taxon, rarefaction analysis provides only dataset-specific insights. Despite these limitations, rarefaction curves remain a popular tool for characterizing whether a pangenome is relatively 'open' (with continuous gene discovery) or 'closed' (where new genome additions contribute few or no new gene clusters). As long as you take such numerical summaries with a huge grain of salt, it is all fine.

Fitting [Heaps' Law](https://en.wikipedia.org/wiki/Heaps'_law) to the rarefaction curve provides a quantitative measure of pangenome openness. The *alpha* value derived from Heaps' Law (sometimes referred to as *gamma* in the literature) reflects how the number of new gene clusters scales with increasing genome sampling. There is no science to define an absolute threshold for an open or a closed pangenome. However, pangenomes with alpha values below 0.3 tend to be relatively closed, and those above 0.3 tend to be relatively open. Higher alpha values will indicate increasingly open pangenomes and lower values will identify progressively closed ones.

### Running the program

The simplest for of the command will look like this,

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s
{{ codestop }}

Which will only report the Heaps' Law alpha value for downstream reporting.

But you can also determine the number of random sampling to be conducted through the `--iteration` parameter. The default is 100. Going above this value will unlikely refine the results, but going below 10 will have a negative influence since the fit will be affected by small amount of sampling:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --iterations 50
{{ codestop }}

When an output file is provided, the program will store the rarefaction curve visualizations in a file:

{{ codestart }}
anvi-compute-rarefaction-curves -p %(pan-db)s \
                                --iterations 50 \
                                --output-file rarefactions.svg
{{ codestop }}

Please note, the file extension (e.g., `.pdf`, `.svg`, `.png`, etc.) will determine the resulting file format.
