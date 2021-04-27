This artifact is the output tables that are displayed when you run %(anvi-estimate-scg-taxonomy)s or %(anvi-estimate-trna-taxonomy)s. 

By default, they won't be outputed anywhere, just displayed in the terminal for your viewing pleasure. If you want them in a tab-delimited file (as a %(genome-taxonomy-txt)s), just provide the `-o` or the `-O` prefix and anvi'o will do that for you.

The content of these tables will depend on how you ran %(anvi-estimate-trna-taxonomy)s or %(anvi-estimate-scg-taxonomy)s. [This blog post](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal) gives you examples of what this looks like for each of the input scenarios for anvi-estimate-scg-taxonomy. Anvi-estimate-scg-taxonomy's output is very similar, just with the results coming from different gene types. They will also be briefly described below. 

When you run %(anvi-estimate-scg-taxonomy)s or %(anvi-estimate-scg-taxonomy)s on 

- a single genome, this table will contain a single line telling you the taxonomy estimate for your genome. It will also show the number of single-copy core genes or tRNA genes that support this estimate. If you run the `--debug` flag, it will also display the hits for all of the single-copy core genes.  
- a single metagenome, this table will list all of the hits for the chosen single-copy core gene or anticodon (by default, the one with the most hits) and their taxonomy information.   
- a %(contigs-db)s and %(profile-db)s with the flag `--compute-scg-coverages`, additional columns will be added that describe the coverage values for your single-copy core gene or tRNA gene hits across your samples.   
- a %(collection)s, this table will show you each of your bins, and the best taxonomy estimate for each one, similarly to how it's displayed for a run on a single genome. 
- a %(metagenomes)s artifact, this table will give a gene entry ID, its taxonomy, and its corresponding coverage in your metagenomes. This format is essentially identical to the output for a single metagenome. If you provide the flag `--matrix-format`, then it will list taxonomy information in each row, and tell you the coverage of each in each of your metagenomes.   

This may sound confusing, but it is easier to understand when looking at the functionality of %(anvi-estimate-scg-taxonomy)s and the comprehensive examples given on [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal).
