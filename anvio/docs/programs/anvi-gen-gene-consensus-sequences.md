This program **provides consensus sequences for genes within a %(contigs-db)s and %(profile-db)s pair**.

In other words, this program collapses variability by assigning the most abundant nucleotide in your sample at each position, generating single consensus sequences for each gene in each sample. 

A basic execution of this program will resemble the following: 

{{ codestart }}
anvi-gen-gene-consensus-sequences -p %(profile-db)s \
                                  -c %(contigs-db)s \
                                  -o %(genes-fasta)s 
{{ codestop }}

The default output is a %(genes-fasta)s, but you can also obtain tab-delimited output by adding the flag `--tab-delimited`.

You also have the option to focus on a subset of the data in your %(contigs-db)s and %(profile-db)s by providing either: 

- A list of gene caller IDs (either as a parameter or through a file with one gene caller ID per line)
- A list of samples to focus on (as a file with a single sample name per line) 

### Additional Parameters 

- You have the option to change the variability engine (i.e., to codons), where variability at this level will be resolved. 
- To compress all variability profiles for each of your samples for a single gene, use the flag `--compress-samples`. This way, the program will report only one consensus sequence for each gene instead of reporting one for each sample. 
- You can obtain consensus sequences for each contig instead of for each gene using `--contigs-mode`
- To report all consensus sequences (even when there are no variable positions), activate `--quince-mode`
