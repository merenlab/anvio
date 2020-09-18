This is a NEWICK-formatted tree that describes the phylogenic relationships of your data. 

{:.notice}
Wondering what the NEWICK format is? Then you're in luck! It has its own [Wikipedia page](https://en.wikipedia.org/wiki/Newick_format).

### How to get one of these? 

You can use %(anvi-gen-phylogenomic-tree)s to create a phylogeny based on a series of genes. 

As discussed on the page for %(anvi-gen-phylogenomic-tree)s, you can also use an external program to get a NEWICK-formatted tree and use that. 

### What can I do with it? 

Firstly, you can use it to reorder elements of the interactive interface. To import this to rearrange the orders that your items appear (in other words, as the central phylogenetic tree when you open the interface), import it using %(anvi-import-items-order)s. To import this as a tree describing your layers (the concentric circles in the anvi'o interface), convert this to a %(misc-data-layer-orders-txt)s and use the program %(anvi-import-misc-data)s.

Secondly, as done in the [Phylogenetics tutorial](http://merenlab.org/2017/06/07/phylogenomics/#working-with-fasta-files), you can open it in the interactive interface without an associated %(contigs-db)s. To do this, run %(anvi-interactive)s as so:

{{ codestart }}
anvi-interactive -t %(phylogeny)s \
                 --title "Phylogenomics Tutorial" \
                 --manual
{{ codestop }}

This will create an empty %(profile-db)s to store any %(bin)ss you create and other such data. You can also add various information, such as taxonomy hits, as done in that same [Phylogenetics tutorial](http://merenlab.org/2017/06/07/phylogenomics/#working-with-fasta-files). 
