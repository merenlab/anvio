This program **provides consensus sequences for the genes within a %(contigs-db)s and %(profile-db)s pair**.

In other words, this collapses variability by assigning the most abundant nucleotide in your sample at each position, giving single consensus sequences for each gene for each sample. 

A basic run of this program will resemble the following: 

{{ codestart }}
anvi-gen-gene-consensus-seuqences -p %(profile-db)s \
                                  -c %(contigs-db)s \
                                  -o %(genes-fasta)s 
{{ codestop }}

The default output is a %(genes-fasta)s, but you can also get a tab-delimited output matrix by adding the flag  `--tab-delimited`.

You also have the option to focus on a subset of the data in your %(contigs-db)s and %(profile-db)s by providing either: 

- A list of gene caller IDs (either as a parameter or through a file with one gene caller ID put line)
- A list of samples to focus on (as a file with a single sample name per line) 

### Additional Parameters 

- You have the option to change the variability engine (i.e. to codons), where variability at this level will be resolved. 
- To compress all variability profiles for each of your samples for a single gene, use the flag `--conpress samples`. This way, the program will only report one consensus sequence for each gene instead of reporting one for each sample. 
- You can get consensus sequences for each contig instead of for each gene with `--contigs-mode`
- To report all consensus sequences (even when there are no variable positions), activate `--quince-mode`
