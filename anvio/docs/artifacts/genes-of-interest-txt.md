Some anvi'o programs that analyze genes allow the user to specify gene caller IDs to work on (i.e., instead of all genes in a %(contigs-db)s). This can sometimes be done via the command line parameter `--gene-caller-ids`, which accepts a list (often comma-delimited) of gene caller IDs, and sometimes via the command line parameter `--genes-of-interest`, which accepts a file in which the gene caller IDs are stored.

The input file for `--genes-of-interest` looks like the following example, with every line of the file corresponding to a single, integer gene caller ID:

```
5
13
79
148
206
```

Ideally, each gene caller ID would match to those stored in the %(contigs-db)s currently being analyzed. But if this is not the case, some anvi'o programs will throw you an error to let you know of any mismatches.