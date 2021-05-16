This program **exports the contig sequences from a %(contigs-db)s**, outputting them as a %(contigs-fasta)s. It also has the ability to output the sequences of your splits instead. 

You can run this program as follows: 

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s
{{ codestop }}

To run it on only a named subset of your contigs, you can provide a list of contigs as a separate file (in the same format as a %(splits-txt)s). For example: 

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s \
                    --contigs-of-interest my_favorite_contigs.txt 
{{ codestop }}

where `my_favorite_contigs.txt` looks like this:

    contig_0001
    contig_0005
    contig_0035
    
### Splits mode

Want to look at your splits instead of your contigs? Just run with the flag `splits-mode` attached. 

{{ codestart }}
anvi-export-contigs -c %(contigs-db)s \
                    -o path/to/%(contigs-fasta)s \
                    --splits-mode
{{ codestop }}
