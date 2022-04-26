Output text files produced by %(anvi-estimate-metabolism)s that describe the presence of **user-defined metabolic pathways** in a %(contigs-db)s.

These files are exactly the same format as those described by %(kegg-metabolism)s, but in addition to (or instead of) information on KEGG modules and KEGG Orthologs, they contain information on user-defined metabolic pathways" (and their component enzymes), as described in %(user-modules-data)s.

## How to get to this output?

You should first read the page on %(user-modules-data)s to learn how to define and set up your own metabolic pathways for use in anvi'o. The program that generates this output is %(anvi-estimate-metabolism)s, and you should run that program with the `--user-modules` parameter to make sure the resulting text files contains the information on your user-defined metabolic pathways. There are two main ways to do it (which are also described on the %(anvi-estimate-metabolism)s help page):

1. To get files describing user-defined metabolic modules _in addition to_ KEGG modules, just use the `--user-modules` parameter.

2. To get files describing _only_ user-defined metabolic modules (instead of KEGG stuff), use both `--user-modules` and `--only-user-modules` parameters.

## What do these files look like?

Check out the %(kegg-metabolism)s page for a comprehensive description of the file formats and various options to customize them. The examples on that page show KEGG data, but the format is the same for user-defined data.
