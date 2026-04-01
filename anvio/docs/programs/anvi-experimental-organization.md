This program uses an anvi'o %(clustering-configuration)s file to access various data sources in anvi'o databases to produce a hierarchical clustering dendrogram for items.

It is especially powerful when the user wishes to create a hierarchical clustering of contigs or gene clusters using only a specific set of samples. If you would like to see an example usage of this program see the article on [combining metagenomics with metatranscriptomics](https://merenlab.org/2015/06/10/combining-omics-data/).

### How does it work?

A %(clustering-configuration)s file tells the program which data matrices to use and how to process them. The program reads all matrices described in the config, scales and normalizes them as instructed, and then merges them into a single combined matrix by concatenating the columns from each matrix. This final merged matrix is then used to perform hierarchical clustering, producing a %(dendrogram)s.

A simple run of this program looks like this:

{{ codestart }}
anvi-experimental-organization %(clustering-configuration)s \
                               -c %(contigs-db)s \
                               -p %(profile-db)s \
                               -N my_organization \
                               -o %(dendrogram)s
{{ codestop }}

If you don't want to store the result in your %(profile-db)s, use the `--skip-store-in-db` flag:

{{ codestart }}
anvi-experimental-organization %(clustering-configuration)s \
                               -c %(contigs-db)s \
                               --skip-store-in-db \
                               -o %(dendrogram)s
{{ codestop }}

You can use the `--dry-run` flag to check whether the program can parse the config file and find the relevant data sources without actually performing the clustering:

{{ codestart }}
anvi-experimental-organization %(clustering-configuration)s \
                               -c %(contigs-db)s \
                               --skip-store-in-db \
                               --dry-run
{{ codestop }}

### Exporting the merged data matrix

In some cases, you may want to see the actual data that goes into the clustering. Since the program combines multiple data matrices into one before clustering, the final form of this merged matrix may not be immediately obvious to the user. But it can be recovered using the `--export-merged-matrix` flag with any of the clustering configurations.

For instance, running the program this way will export the combined and scaled matrix as a TAB-delimited file while still performing the clustering and storing the result in the database:

{{ codestart }}
anvi-experimental-organization %(clustering-configuration)s \
                               -c %(contigs-db)s \
                               -o %(dendrogram)s \
                               --export-merged-matrix merged_matrix.txt
{{ codestop }}

The resulting file will contain one row per item and columns from all input matrices after they have been normalized, log-transformed, and scaled according to the config. This can be useful for debugging, or for generating dendrograms using other software.

You can also recover the matrix *without* running the clustering step, which may be costly for large datasets, and without storing anything in the database -- as in "just scale the data, merge it, give it back to me as a TAB-delimited file, and stop there":

{{ codestart }}
anvi-experimental-organization %(clustering-configuration)s \
                               -c %(contigs-db)s \
                               --dry-run \
                               --export-merged-matrix merged_matrix.txt
{{ codestop }}

This will scale and merge the matrices but skip the hierarchical clustering entirely.

Let's assume you have a merged %(profile-db)s, and you wish to get the scaled and merged matrix for tetranucleotide frequency and coverage data. This is exactly what you could do to recover the primary data that would go into the clustering step for that dataset::

```bash
# locate the tnf-cov config file in the anvi'o source directory
tnf_cov=$(python -c "from pathlib import Path; import anvio; print(Path(anvio.__file__).parent / 'data/clusterconfigs/merged/tnf-cov')")

# run the experimental organization to get the merged matrix for tnf-cov without doing the clustering
anvi-experimental-organization $tnf_cov \
                               -c CONTIGS.db \
                               -p PROFILE.db \
                               --dry-run \
                               --export-merged-matrix tnf_cov_merged_matrix.txt
```

Et voilà! You will have a file called `tnf_cov_merged_matrix.txt` that contains the merged and scaled data for tetranucleotide frequency and coverage for all contigs in your dataset.
