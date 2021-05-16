This program exports your gene calls given a %(contigs-db)s and a gene caller. The output of this is a %(gene-calls-txt)s. 

To see the gene callers available in your contigs database, run 

{{ codestart }}
anvi-export-gene-calls -c %(contigs-db)s \
                       --list-gene-callers
{{ codestop }}

By default, this list will probably include [Prodigal](https://github.com/hyattpd/Prodigal), which identifies genes when creating a %(contigs-db)s. For in this example, we'll use export Prodigal-identified genes. Note that you can also get genes from more than one source by providing several gene-callers in a comma-delimited list.  

Then, you can export all of your gene callers in an orderly fashion by running 

{{ codestart }}
anvi-export-gene-calls -c %(contigs-db)s \
                       --gene-caller Prodigal \
                       -o path/to/output
{{ codestop }}

This will put a %(gene-calls-txt)s in the path you specified containing all of your Prodigal genes. 

If you don't want to display the amino acid sequences of each gene (they can crowd the file very quickly if you don't want to see them), just add the flag `--skip-sequence-reporting`
