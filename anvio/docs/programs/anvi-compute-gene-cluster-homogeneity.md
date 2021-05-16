This program **computes both the geometric homogeneity and functional homogeneity for the gene clusters in a %(pan-db)s.** 

*Geometric homogeneity* and *functional homogeneity* are anvi'o specific terms that describe how similar genes within a gene cluster are to each other in different ways. Briefly, geometric homogeneity compares the positions of gaps in the aligned residues without considering specific amino acids, and functional homogeneity examines point mutations to amino acids and compares how similar the resulting amino acids are chemically. See [this page](http://merenlab.org/2016/11/08/pangenomics-v2/#inferring-the-homogeneity-of-gene-clusters) for more details. 

You can run this program as so: 

{{ codestart }}
anvi-compute-gene-cluster-homogeneity -p %(pan-db)s \
                                      -g %(genomes-storage-db)s \
                                      -o path/to/output.txt \
                                      --store-in-db
{{ codestop }}

This run will put the output directly in the database, as well as provide it as a separate file as the specified output path. 

You also have the option to calculate this information about only specific gene clusters, either by providing a gene cluster ID, list of gene cluster IDs, %(collection)s or %(bin)s. 

To save on runtime, you can also enable `--quick-homogeneity`, which will not check for horizontal geometric homogenity (i.e. it will not look at alignments within a single gene). This will be less accurate for detailed analyses, but it will run faster. 

Here is an example run that uses this flag and only looks at a specific collection: 

{{ codestart }}
anvi-compute-gene-cluster-homogeneity -p %(pan-db)s \
                                      -g %(genomes-storage-db)s \
                                      -o path/to/output.txt \
                                      --store-in-db \ 
                                      -C %(collection)s \
                                      --quick-homogeneity 
{{ codestop }}

You can also use multithreading if you're familiar with that. 
