This program converts a distance matrix (computed from a %(view-data)s artifact) into a %(dendrogram)s. 

It uses the numerical data in a %(view-data)s to compute a distance matrix behind the scenes, and then runs some hierarchical clustering to create a %(dendrogram)s for all of your items. 

With all default parameters, a run would look like this:

{{ codestart }}
anvi-matrix-to-newick -o path/for/%(dendrogram)s \ 
                      %(view-data)s 
{{ codestop }}

If your input file has your samples as rows instead of columns, just add the flag `--transpose`. 

You can also ask for an additional output file: the order of the items in the resulting dendrogram as a %(misc-data-items-order)s in LIST format. To get this, simply provide a path to its desired location  with `--items-order-file`. 

Additionally, for hierarchical clustering, you can change the distance metric (a full list of the available metrics can be found [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)) or the linkage method (though this is not recommended, the list of options can be found [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)).
