This file contains **the frequency of each amino acid for some reference context in your %(contigs-db)s**.  

This is a tab-delimited table where each column represents an amino acid and each row represents a specific reference context (most often this will be a gene after running %(anvi-get-codon-frequencies)s). The numbers will either refer to counts of each amino acid or precent normalizations depending on the parameters with which you ran %(anvi-get-codon-frequencies)s. 

You can also use %(anvi-get-aa-counts)s to get this information for a %(bin)s, %(collection)s, or %(splits-txt)s. 

### Example

    gene_caller_id  Ala Arg Thr Asp ...
        1           0   0   1   2
        2           1   0   0   2
        .
        .
        .
