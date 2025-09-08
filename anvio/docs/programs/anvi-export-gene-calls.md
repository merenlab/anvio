The purpose of this program is to exports your gene calls in a given %(contigs-db)s and a gene caller, in the form of a %(gene-calls-txt)s. 

To see the gene callers available in your contigs database, you can use %(anvi-db-info)s or use this program with the following flag: 

{{ codestart }}
anvi-export-gene-calls -c %(contigs-db)s \
                       --list-gene-callers
{{ codestop }}

Running this will export all of your gene calls identified by the gene caller [pyrodigal-gv](https://github.com/althonos/pyrodigal-gv) (assuming it is in your %(contigs-db)s:

{{ codestart }}
anvi-export-gene-calls -c %(contigs-db)s \
                       --gene-caller Prodigal \
                       -o %(gene-calls-txt)s
{{ codestop }}

{:.notice}
You can export genes from more gene callers by providing a comma-separated list of gene caller names.

If you don't want to display the amino acid sequences of each gene (they can crowd the file very quickly if you don't want to see them), you can add the following flag:

{{ codestart }}
anvi-export-gene-calls -c %(contigs-db)s \
                       --gene-caller Prodigal \
                       --skip-sequence-reporting \
                       -o %(gene-calls-txt)s
{{ codestop }}

## Advanced uses

This program can take a lot of time and memory when working with very large %(contigs-db)s files (such as those that are more than 10 Gb in file size or more than 10 million contigs).

In that case you can export your gene calls the following way within minutes and a small memory space.

First open your %(contigs-db)s:

{{ codestart }}
sqlite3 %(contigs-db)s
{{ codestop }}

Then run these lines,

{{ codestart }}
.mode csv 
.headers on 
.out %(gene-calls-txt)s
select gene_callers_id, contig, start, stop, direction, partial from genes_in_contigs;
{{ codestop }}

You can also continue with these lines to get the amino acid sequences for them:

{{ codestart }}
.mode csv 
.headers on 
.out AMINO-ACID-SEQUENCES.txt
select * from genes_in_contigs;
{{ codestop }}
