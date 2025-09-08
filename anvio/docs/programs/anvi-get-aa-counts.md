Similarly to %(anvi-get-codon-frequencies)s, this program counts the number of times each amino acid occurs in a given sequence, whether that's a %(collection)s, %(bin)s, set of contigs (listed in a %(splits-txt)s), or a set of genes. The output of this is a %(aa-frequencies-txt)s. 

There are four possible things you can count the amino acid frequencies in: 
* All of the contigs in a %(contigs-db)s
* A series of %(bin)ss
* A list of contigs
* A list of genes

Examples for each are below.

### Option 1: all contigs in a contigs-db

To count the amino acids in all of the contigs in a %(contigs-db)s, you can just provide the %(contigs-db)s of interest, as so:

{{ codestart }}
anvi-get-aa-counts -c %(contigs-db)s \
                   -o path/to/%(aa-frequencies-txt)s
{{ codestop }}

### Option 2: a series of bins in a collection 

To count the amino acid frequencies for a series of %(bin)ss, you'll need to provide three additional parameters: the %(profile-db)s that you used for binning, the %(collection)s that your bins are contained in, and a text file that describes which bins you are interested in. This text file should have only one bin ID per line. 

So, your run would look something like this: 

{{ codestart }}
anvi-get-aa-counts -c %(contigs-db)s \
                   -o path/to/%(aa-frequencies-txt)s \
                   -p %(profile-db)s \
                   -C %(collection)s \
                   -B my_favorite_bins.txt
{{ codestop }}

`my_favorite_bins.txt` would look something like this:

    bin_00001
    bin_00004
    
### Option 3: a list of contigs

Just provide a %(splits-txt)s file that lists the contigs you want to look at. 

{{ codestart }}
anvi-get-aa-counts -c %(contigs-db)s \
                   -o path/to/%(aa-frequencies-txt)s \
                   --contigs-of-interest %(splits-txt)s
{{ codestop }}

### Option 4: a list of genes 

Just provide a list of gene caller ids, straight into the terminal, like so:

{{ codestart }}
anvi-get-aa-counts -c %(contigs-db)s \
                   -o path/to/%(aa-frequencies-txt)s \
                   --gene-caller-ids gene_1,gene_2,gene_3
{{ codestop }}
