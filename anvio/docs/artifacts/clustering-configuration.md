A clustering configuration file is **a special file that tells anvi'o how to combine multiple data sources into a single merged matrix prior to hierarchical clustering**.

These configuration files are primarily used internally by anvi'o when it merges profile databases to calculate hierarchical clustering of contigs, when it generates pangenomes to calculate hierarchical clustering of gene clusters, and so on. But they can also be utilized by %(anvi-experimental-organization)s to re-run such clustering tasks with the default clustering configuration files or user-defined ones.

Anvi'o ships with a set of default clustering configurations under `anvio/data/clusterconfigs/`, organized into subdirectories for different use cases (`merged/`, `single/`, `pan/`, etc.). You can browse these configurations online at [github.com/merenlab/anvio/tree/master/anvio/data/clusterconfigs](https://github.com/merenlab/anvio/tree/master/anvio/data/clusterconfigs), or find them on your computer wherever anvi'o is installed. For instance, you can run the following command in your terminal and it will show you the exact path for the clustering configuration files on your system:

```
python -c "from pathlib import Path; import anvio; print(Path(anvio.__file__).parent / 'data/clusterconfigs')"
```

These built-in configurations are good starting points for writing your own.

### File format

A clustering configuration file uses the [INI format](https://en.wikipedia.org/wiki/INI_file) with a `[general]` section and one or more matrix sections. Here is an example that combines tetranucleotide frequency and coverage data:

{{ codestart }}
[general]
distance = euclidean
linkage = ward

[TNF !CONTIGS.db::kmer_contigs]

[Coverage PROFILE.db::mean_coverage_contigs]
table_form = view
normalize = False
log = True
{{ codestop }}

### The general section

The `[general]` section can contain the following optional parameters:

- **`distance`**: The distance metric to use for hierarchical clustering (e.g., `euclidean`, `cityblock`, `cosine`). If not specified, anvi'o will use its default.
- **`linkage`**: The linkage method to use for hierarchical clustering (e.g., `ward`, `complete`, `average`). If not specified, anvi'o will use its default.
- **`output_file`**: A file path to write the resulting Newick tree to.
- **`name`**: A single-word name for the resulting clustering.

### Matrix sections

Each additional section describes one data matrix. The section header has two parts separated by a space: the **alias** (a short name for reporting) and the **matrix source** (a file path or a database reference).

{{ codestart }}
[alias matrix_source]
{{ codestop }}

The matrix source can be:

- A plain file name (resolved relative to the input directory),
- A database table reference like `PROFILE.db::mean_coverage_contigs`, or
- A database table reference with a `!` prefix like `!CONTIGS.db::kmer_contigs`, which tells anvi'o to resolve the database path from the calling program rather than looking for a literal file.

Each matrix section can contain the following optional parameters:

- **`normalize`** (default: `True`): Whether to normalize the vectors in this matrix prior to merging. Set to `False` to skip normalization.
- **`log`** (default: `False`): Whether to apply a log transformation (`log10(x + 1)`) to the values. Set to `True` to enable.
- **`columns_to_use`**: A comma-separated list of column names to include from this matrix. If not specified, all columns that pass anvi'o's essential field check are used.
- **`table_form`**: Set to `view` if the matrix source refers to a view table in a profile database rather than a regular table.
- **`ratio`**: An integer that describes the weight of this matrix relative to others (used in some experimental clustering workflows).

### How the merging works

When multiple matrices are listed in a configuration file, anvi'o processes them as follows:

1. Each matrix is loaded and filtered to a common set of rows (items) shared across all matrices.
2. Each matrix is independently scaled: normalization and/or log transformation is applied as specified.
3. The scaled matrices are concatenated column-wise into a single merged matrix, so each row contains the combined features from all input matrices.
4. Hierarchical clustering is performed on this merged matrix.

This means that a matrix with more columns will inherently have more influence on the final clustering. Understanding this is important when designing your configurations. You can use the `--export-merged-matrix` flag of %(anvi-experimental-organization)s to inspect the final merged matrix that goes into clustering.
