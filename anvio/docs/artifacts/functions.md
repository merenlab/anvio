This is an artifact that describes **annotation of genes in your %(contigs-db)s with functions**.

Broadly used across anvi'o, functions are one of the most essential pieces of information stored in any %(contigs-db)s. To see what annotation sources for functions are available in a given %(contigs-db)s or %(genomes-storage-db)s, you can use the program %(anvi-db-info)s.

To populate a given %(contigs-db)s with functions, anvi'o includes multiple programs that can annotate genes using various sources of annotation. These programs include,

* %(anvi-run-ncbi-cogs)s, which uses NCBI's [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/),
* %(anvi-run-pfams)s, which uses EBI's [Pfam database](https://pfam.xfam.org/),
* %(anvi-run-cazymes)s, which uses the dbCAN [CAZyme HMMs](https://bcb.unl.edu/dbCAN2/download/Databases/)
* %(anvi-run-kegg-kofams)s, which uses the [Kyoto Encyclopedia of Genes and Genomes](https://www.genome.jp/kegg/) (KEGG) database and produces %(kegg-functions)s, which is the necessary annotation information that can be used by the program %(anvi-estimate-metabolism)s.

In addition, you can use the program %(anvi-import-functions)s with a simple %(functions-txt)s to import functions from any other annotation source, or to import any ad hoc, user-defined function to later access through anvi'o interfaces or programs.

{:.notice}
You can use %(anvi-import-functions)s also to import functions from EggNOG or InterProScan as described in [this blog post](http://merenlab.org/2016/06/18/importing-functions/).

You can also use %(anvi-export-functions)s to obtain a file containing these functional annotations through a %(functions-txt)s artifact, and use %(anvi-display-functions)s to show the distribution of functions across multiple %(contigs-db)ss.
