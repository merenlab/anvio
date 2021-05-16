This program **calculates the frequency of each codon or amino acid of every gene in your %(contigs-db)s**. 

To run with all standard parameters, simply provide a %(contigs-db)s and path for the output file as follows: 

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \ 
                -o name/of/output_file.txt 
{{ codestop }}

The output of this is a %(codon-frequencies-txt)s that counts the number of times each codon appears in all of your genes.

If instead you want to calculate the data for the amino acids, run 

{{ codestart }}
anvi-get-codon-frequencies -c %(contigs-db)s \ 
                -o name/of/output_file.txt  \
                --return-AA-frequencies-instead \
                --gene-caller-id MY_FAVORITE_GENE
{{ codestop }}

In this example, the flag `gene-caller-id` means that it will only count the amino acid frequencies of a single gene, namely `MY_FAVORITE_GENE`.

You can also return the data as a percent of the total number of codons or amino acids in the gene (with the flag `--percent-normalize`) or calculate the percent that each codon encoding the same amino acid appears in the gene (for example, 0.4 GCT and 0.6 GCC for alanine) (with the flag `--merens-codon-normalization`). 
