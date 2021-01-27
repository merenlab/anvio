This program **associates genes in your %(contigs-db)s with functions using NCBI's [Clusters of Orthologus Groups (COGs) database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/).**

Before you run this program, you'll have to set up the COGs database on your computer with the program %(anvi-setup-ncbi-cogs)s.

To run, you'll need to provide a %(contigs-db)s. If you stored the %(cogs-data)s that you got from running %(anvi-setup-ncbi-cogs)s in a custom location, you'll need to provide that path as well. The output is a %(functions)s artifact.

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --cog-data-dir path/to/%(cogs-data)s
{{ codestop }}

Without the flag `--cog-data-dir`, anvi'o will just search in the default location.

By default, this program uses DIAMOND in the "fast" setting for database searching. To instead run in "sensitive" mode, just call:

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --sensitive
{{ codestop }}

You can also use blastp to search, by running:

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --search-with blastp
{{ codestop }}


