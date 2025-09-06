This program **computes both the geometric homogeneity and functional homogeneity for the gene clusters stored in a %(pan-db)s.** 

*Geometric homogeneity* and *functional homogeneity* are anvi'o-specific metrics that describe how similar genes within a gene cluster are to each other in different ways. Geometric homogeneity analyzes the positions of gaps in aligned residues without considering specific amino acids, while functional homogeneity examines point mutations to amino acids and compares the chemical similarity of the resulting amino acids. See [this page](http://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters) for more detailed information. 

You can execute this program as follows: 

{{ codestart }}
anvi-compute-gene-cluster-homogeneity -p %(pan-db)s \
                                      -g %(genomes-storage-db)s \
                                      -o path/to/output.txt \
                                      --store-in-db
{{ codestop }}

This execution will store the output directly in the database and provide it as a separate file at the specified output path. 

The analysis can also be restricted to specific gene clusters by providing a gene cluster ID, list of gene cluster IDs, %(collection)s, or %(bin)s. 

To reduce runtime, you can enable the `--quick-homogeneity` option, which skips the horizontal geometric homogeneity analysis (i.e., it will not examine alignments within individual genes). This approach provides faster execution but with reduced accuracy for detailed analyses. 

The following example demonstrates usage of this flag while restricting analysis to a specific collection: 

{{ codestart }}
anvi-compute-gene-cluster-homogeneity -p %(pan-db)s \
                                      -g %(genomes-storage-db)s \
                                      -o path/to/output.txt \
                                      --store-in-db \ 
                                      -C %(collection)s \
                                      --quick-homogeneity 
{{ codestop }}

Multithreading capabilities are also available for users who require parallel processing.
