This program **annotates genes in your %(contigs-db)s with functions using NCBI's [Clusters of Orthologus Groups (COGs) database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/).**

This program assumes that the user has successfully set up the COGs database on their computer using the anvi'o program %(anvi-setup-ncbi-cogs)s.

The only critical parameter to %(anvi-run-ncbi-cogs)s is a %(contigs-db)s. The program will store its output in the %(contigs-db) sas a %(functions)s artifact.

If the %(cogs-data)s was stored at a specific path when %(anvi-setup-ncbi-cogs)s was run, then providing that path using the `--cog-data-dir` parameter is also necessary.

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --cog-data-dir path/to/%(cogs-data)s
{{ codestop }}

Without the flag `--cog-data-dir`, anvi'o will just search in the default location.

### Choosing a different database version
If you want to annotate your genes with a non-default version of %(cogs-data)s, provide that version after the `--cog-version` parameter:

{{ codestart }}
anvi-setup-ncbi-cogs --cog-version COG14
{{ codestop }}

### Choosing a different search program
By default, this program uses `diamond` to search for hits to the database. You can also use `blastp`` to search, by running:

{{ codestart }}
anvi-run-ncbi-cogs -c %(contigs-db)s \
            --search-with blastp
{{ codestop }}
