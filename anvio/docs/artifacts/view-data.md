View data refers to a matrx where each column represents a specific sample and each row describes some attribute of that sample (most often a sequence's abundance per sample). 

For example, in the [pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions), the `PROCHLORO-functions-occurrence-frequency.txt` is a view-data. 

You can use this to compute a distance matrix to generate a dendrogram (using %(anvi-matrix-to-newick)s) or direclty input it to %(anvi-interactive)s to visualize the distribution of your items across samples. 
