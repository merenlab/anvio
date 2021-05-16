This program allows you to **export the sequences of your gene calls** from a %(contigs-db)s or %(genomes-storage-db)s in the form of a %(genes-fasta)s. 

If you want other information about your gene calls from a %(contigs-db)s, you can run %(anvi-export-gene-calls)s (which outputs a %(gene-calls-txt)s) or get the coverage and detection information with %(anvi-export-gene-coverage-and-detection)s.

### Running on a contigs database

You can run this program on a %(contigs-db)s like so:

{{ codestart }}
anvi-get-sequences-for-gene-calls -c %(contigs-db)s \
                                  -o path/to/output
{{ codestop }}

This is create a %(genes-fasta)s that contains every gene in your contigs database. If you only want a specific subset of genes, you can run the following: 

{{ codestart }}
anvi-get-sequences-for-gene-calls -c %(contigs-db)s \
                                  -o path/to/output \
                                  --gene-caller-ids 897,898,1312 \
                                  --delimiter ,
{{ codestop }}

Now the resulting %(genes-fasta)s will contain only those three genes. 

You also have the option to report the output in [gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), report extended deflines for each gene, or report amino acid sequences instead of nucleotide sequences.

### Running on a genomes storage database

You can also get the sequences from gene calls in a %(genomes-storage-db)s, like so:

{{ codestart }}
anvi-get-sequences-for-gene-calls -g %(genomes-storage-db)s \
                                  -o path/to/output
{{ codestop }}

This will create a %(genes-fasta)s that contains every gene in your genomes storage database. To focus on only a subset of the genomes contained in your database, use the flag `--genome-names`. You can provide a comma-delimited list of genome names or a flat text file that contains one genome per line. Alternatively, you could provide a list of gene-caller-ids as specified above. 

You also have the option to report the output in [gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), report extended deflines for each gene, or report amino acid sequences instead of nucleotide sequences.
