This artifact **contains information about the functions of the genes in your %(contigs-db)s.**

It is a table within your %(contigs-db)s that contains functional annotations for your genes. 

To get one of these for your %(contigs-db)s, you can either import it (using %(anvi-import-functions)s) or make one yourself by running your contigs against one of two databases available in anvi'o:
* NCBI [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) -- see %(anvi-run-ncbi-cogs)s for instructions
* EBI's [Pfam database](https://pfam.xfam.org/) -- see %(anvi-run-pfams)s for instructions
* The [KOfamKOALA database](https://www.genome.jp/tools/kofamkoala/) -- see %(anvi-run-kegg-kofams)s for instructions. Functions specifically from the KOfam database are used to run %(anvi-estimate-metabolism)s and are also called %(kegg-functions)s.

You can also use EggNOG or InterProScan and then import the results into anvi'o, as described in [this blog post](http://merenlab.org/2016/06/18/importing-functions/).

This is used to run %(anvi-analyze-synteny)s. 

You can also use %(anvi-export-functions)s to obtain a file containing these functional annotations through a %(functions-txt)s artifact. 
