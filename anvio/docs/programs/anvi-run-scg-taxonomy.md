This program **associates the single-copy core genes in your %(contigs-db)s with taxnomy information.**  

Once this information is stored in your %(contigs-db)s (in the form of a %(scgs-taxonomy)s artifact), you can run %(anvi-estimate-scg-taxonomy)s or use the %(anvi-interactive)s and enable "Realtime taxonomy estimate for bins." Check out [this tutorial](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/) for more information. 

In order to run this program, you'll need a %(scgs-taxonomy-db)s, which you can set up by running %(anvi-setup-scg-taxonomy)s. 

### What does this program do? 

In short, this program searches all of the single-copy core genes that it uses for this workflow (which are the 22 listed on [this page](https://github.com/merenlab/anvio/tree/master/anvio/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES)) against the [GTDB](https://gtdb.ecogenomic.org/) databases that you downloaded, and stores hits in your %(contigs-db)s. In other words, it finds your single-copy core genes and assigns them taxonomy. This way, it can use these single-copy core genes later to estimate the taxnomy of larger groups of contigs that include these single-copy core genes when you run %(anvi-estimate-scg-taxonomy)s. 

### Sweet. How do I run it? 

{{ codestart }}
anvi-run-scg-taxonomy -c %(contigs-db)s
{{ codestop }}

In case you're running this on a genome and not getting any hits, you have the option to try lowering the percent identity required for a hit (as long as you're careful with it). The default value is 90 percent. 

{{ codestart }}
anvi-run-scg-taxonomy -c %(contigs-db)s \
                      --min-percent-identity 70
{{ codestop }}
