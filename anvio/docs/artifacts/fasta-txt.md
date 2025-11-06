This is a file used by %(anvi-run-workflow)s that lists the name and path of all of the input %(fasta)s files.

In its simplest form, a %(fasta-txt)s is a TAB-delmited file with two columns for `name` and `path`. Here is an example:

|name|path|
|:--|:--|
|SAMPLE_01|path/to/sample_01.fa|
|SAMPLE_02|path/to/sample_02.fa|

Paths can be absolute or relative, and FASTA files can be compressed or not. That's all up to you.

One of the primary users of the %(fasta-txt)s is the [anvi'o snakemake workflows](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/), and to make it more compatible to complex workflow needs, %(fasta-txt)s supports the following additional columns to provide more information for each FASTA file when available, such as %(external-gene-calls)s file and/or a %(functions-txt)s.

Here is an example with those additional columns:

|name|path|external_gene_calls|gene_functional_annotation|
|:--|:--|:--|:--|
|SAMPLE_01|path/to/sample_01.fa|%(external-gene-calls)s_01.txt|%(functions-txt)s_01.txt|
|SAMPLE_02|path/to/sample_02.fa|%(external-gene-calls)s_02.txt|%(functions-txt)s_02.txt|

For more information, check out the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)
