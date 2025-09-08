You can send any matrix file to this program to get a %(dendrogram)s from it.

An example run would look like this:

{{ codestart }}
anvi-matrix-to-newick TAB_DELIMITED_DATA.txt \
                      %(dendrogram)s
{{ codestop }}

By default, %(anvi-matrix-to-newick)s will cluster rows. With the flag `--transpose`, it will cluster columns.

See [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html) a list of distance metrics you can use, and [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html) a list of linkage methods you can use.

%(anvi-matrix-to-newick)s can handle missing data, but in that case the program will not normalize your data and will assume that it is already normalized.
