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

Please note that this program allows you to format the deflines of the resulting FASTA file to a great extent. For this, it uses a set of previously-defined variables you can use to define a template of your liking. You can learn about the available variables, you can include the following flag in your command:

{{ codestart }}
anvi-get-sequences-for-gene-calls -c %(contigs-db)s \
                                  -o path/to/output \
                                  --list-defline-variables
{{ codestop }}

which will give you an output similar to this:

```
WARNING
===============================================
Here are the variables you can use to provide a user-defined defline template:

* {gene_caller_id}
* {contig_name}
* {start}
* {stop}
* {direction}
* {length}
* {contigs_db_project_name}

Remember, by default, anvi'o will only use '{gene_caller_id}' to format the
deflines of FASTA files it produces.
```

Now you can use the `--defline-format` parameter with any combination of these variables to control the output FASTA file deflines. Here are some examples using the %(contigs-db)s file included in the Infant Gut Dataset, and the contents of the resulting FASTA file:

```
anvi-get-sequences-for-gene-calls -c INFANT-GUT-TUTORIAL/CONTIGS.db \
                                  --gene-caller-ids 77 \
                                  -o 77.fa \
                                  --defline-format "{gene_caller_id}"
```

```
>77
ATGAGTAATTTATTGCGAGCAGAAGGCGTGTCTTATCAAGTAAATGGTCGTAGCATTCTTTCTGATATTGATTTGTCATTTGAAACAGGCAGCAATACAACAATTGTTGGTCCTTCAGGT
AGCGGGAAAAGTACATTTTTAAAAATTTTATCTTCATTATTAAGTCCTACAGAAGGCGAAATTTTTTATCAAGAAGCGCCAATTACTACAATGCCAATCGAAACATACCGCCAAAAGGTT
TCTTATTGTTTTCAGCAGCCAACTTTATTTGGTGAAACCGTGTATGATAATTTGTTATTTCCATTTACCGTCAGACAAGAAGCGTTTAATCAGGAAAAAGTCGTGGCATTACTCCAACAA
GTGAAATTGCCCGCTGCTTATCTTGAAAAGAAAATAGCCGAACTCTCTGGTGGTGAGCGACAACGGGTTGCTTTGCTACGAAACATTATTTTTGTACCAGATGTTTTATTATTAGACGAA
GTTACAACGGGATTAGATGAAGAAAGCAAACAGATTGTCAATCAATTGTTAAACCAATTAAACAAAGAGCAAGGAGTCACGCTGGTTCGTGTCACGCATGATACCGAAGAAATTCAGCAA
GCACAGCAAGTGATTCGTATTGTAGCAGGAAAGGTGGCGCCGACAGATGGATTTAGCAGTTAA
```

The simple option shown above also is the default defline format anvi'o uses. Here is a more sophisticated example:

```
anvi-get-sequences-for-gene-calls -c INFANT-GUT-TUTORIAL/CONTIGS.db \
                                  --gene-caller-ids 77 \
                                  -o 77.fa \
                                  --defline-format "{contigs_db_project_name}_{contig_name}_{gene_caller_id}"
```

```
>Infant_Gut_Contigs_from_Sharon_et_al_Day17a_QCcontig1_77
ATGAGTAATTTATTGCGAGCAGAAGGCGTGTCTTATCAAGTAAATGGTCGTAGCATTCTTTCTGATATTGATTTGTCATTTGAAACAGGCAGCAATACAACAATTGTTGGTCCTTCAGGT
AGCGGGAAAAGTACATTTTTAAAAATTTTATCTTCATTATTAAGTCCTACAGAAGGCGAAATTTTTTATCAAGAAGCGCCAATTACTACAATGCCAATCGAAACATACCGCCAAAAGGTT
TCTTATTGTTTTCAGCAGCCAACTTTATTTGGTGAAACCGTGTATGATAATTTGTTATTTCCATTTACCGTCAGACAAGAAGCGTTTAATCAGGAAAAAGTCGTGGCATTACTCCAACAA
GTGAAATTGCCCGCTGCTTATCTTGAAAAGAAAATAGCCGAACTCTCTGGTGGTGAGCGACAACGGGTTGCTTTGCTACGAAACATTATTTTTGTACCAGATGTTTTATTATTAGACGAA
GTTACAACGGGATTAGATGAAGAAAGCAAACAGATTGTCAATCAATTGTTAAACCAATTAAACAAAGAGCAAGGAGTCACGCTGGTTCGTGTCACGCATGATACCGAAGAAATTCAGCAA
GCACAGCAAGTGATTCGTATTGTAGCAGGAAAGGTGGCGCCGACAGATGGATTTAGCAGTTAA
```

And finally, here is an example that is quite comprehensive in its request from anvi'o:

```
anvi-get-sequences-for-gene-calls -c INFANT-GUT-TUTORIAL/CONTIGS.db \
                                  --gene-caller-ids 77 \
                                  -o 77.fa \
                                  --defline-format "{gene_caller_id} contig_name:{contig_name}|gene_start:{start}|gene_stop:{stop}|gene_length:{length}|gene_direction:{direction}"
```

```
>77 contig_name:Day17a_QCcontig1|gene_start:67276|gene_stop:67939|gene_length:663|gene_direction:f
ATGAGTAATTTATTGCGAGCAGAAGGCGTGTCTTATCAAGTAAATGGTCGTAGCATTCTTTCTGATATTGATTTGTCATTTGAAACAGGCAGCAATACAACAATTGTTGGTCCTTCAGGT
AGCGGGAAAAGTACATTTTTAAAAATTTTATCTTCATTATTAAGTCCTACAGAAGGCGAAATTTTTTATCAAGAAGCGCCAATTACTACAATGCCAATCGAAACATACCGCCAAAAGGTT
TCTTATTGTTTTCAGCAGCCAACTTTATTTGGTGAAACCGTGTATGATAATTTGTTATTTCCATTTACCGTCAGACAAGAAGCGTTTAATCAGGAAAAAGTCGTGGCATTACTCCAACAA
GTGAAATTGCCCGCTGCTTATCTTGAAAAGAAAATAGCCGAACTCTCTGGTGGTGAGCGACAACGGGTTGCTTTGCTACGAAACATTATTTTTGTACCAGATGTTTTATTATTAGACGAA
GTTACAACGGGATTAGATGAAGAAAGCAAACAGATTGTCAATCAATTGTTAAACCAATTAAACAAAGAGCAAGGAGTCACGCTGGTTCGTGTCACGCATGATACCGAAGAAATTCAGCAA
GCACAGCAAGTGATTCGTATTGTAGCAGGAAAGGTGGCGCCGACAGATGGATTTAGCAGTTAA
```

You also have the option to report the output in [gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), report extended deflines for each gene, or report amino acid sequences instead of nucleotide sequences.

### Running on a genomes storage database

You can also get the sequences from gene calls in a %(genomes-storage-db)s, like so:

{{ codestart }}
anvi-get-sequences-for-gene-calls -g %(genomes-storage-db)s \
                                  -o path/to/output
{{ codestop }}

This will create a %(genes-fasta)s that contains every gene in your genomes storage database. To focus on only a subset of the genomes contained in your database, use the flag `--genome-names`. You can provide a comma-delimited list of genome names or a flat text file that contains one genome per line. Alternatively, you could provide a list of gene-caller-ids as specified above.

You also have the option to report the output in [gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), report extended deflines for each gene, or report amino acid sequences instead of nucleotide sequences.
