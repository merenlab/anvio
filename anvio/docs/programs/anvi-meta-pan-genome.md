This program integrates the information from an %(internal-genomes)s artifact into a %(pan-db)s, creating a metapangenome. 

A metapangenome contains both the information in a metagenome (i.e. their abundances in different samples as described in your %(profile-db)s) and the information in a pangenome (i.e. the gene clusters in your dataset). This is useful because you are able to observe which gene cluster patterns are present in certain environments. For an example of a metapangenomic workflow, take a look [here](http://merenlab.org/data/prochlorococcus-metapangenome/) (this tutorial was written before this program, but the insights persist). 

To use this program, provide a %(pan-db)s and %(genomes-storage-db)s pair, as well as an %(internal-genomes)s.

{{ codestart }}
anvi-meta-pan-genome -p %(pan-db)s \
                     -g %(genomes-storage-db)s \
                     -i %(internal-genomes)s 
{{ codestop }}

However, when integrating metagenomic and pangenomic data together, you'll get a lot of data. You can set two additional parameters to help you filter out data that doesn't mean certain standards:

- `--fraction-of-median-coverage`: this threshold removes genes with less than this fraction of the median coverage. The default is 0.25. So, for example, if the median coverage in your data was 100X, this would remove all genes with coverage less than 25X. 
- `--min-detection`: this threshold removes genomes with detection less than this value in all samples. The default is 0.5.
