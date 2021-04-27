This aptly-named program **gets the sequences for the gene clusters stored in a %(pan-db)s and returns them as either a %(genes-fasta)s or a %(concatenated-gene-alignment-fasta)s** (which you can use to run %(anvi-gen-phylogenomic-tree)s). This gives you advanced access to your gene clusters, which you can take out of anvi'o, use for phylogenomic analyses, or do whatever you please with. 

You also have the option to output the sequences of your choice as a %(misc-data-items)s (with `add-into-items-additional-data-table`), which can be added to the %(interactive)s interface as additional layers. 

While the number of parameters may seem daunting, many of the options just help you specify exactly which gene clusters you want to get the sequences  from. 

### Running on all gene clusters

Here is a basic run, that will  export alignments for every single gene cluster found in the %(pan-db)s as amino acid sequences :

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s
{{ codestop }}

To get the DNA sequences instead, just add `--report-DNA-sequences`. 

### Exporting only specific gene clusters

#### Part 1: Choosing gene clusters by collection, bin, or name

You can export only the sequences for a specific %(collection)s or %(bin)s with the parameters `-C` or `-b` respectively. You also have the option to display the collections and bins available in your %(pan-db)s with `--list-collections` or `--list-bins`

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s \
                                     -C %(collection)s 
{{ codestop }}

Alternatively, you can export the specific gene clusters by name, either by providing a single gene cluster ID or a file with one gene cluster ID per line. For example: 

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s \
                                     --gene-cluster-ids-file gene_clusters.txt
{{ codestop }}

where `gene_clusters.txt` contains the following:

    GC_00000618
    GC_00000643
    GC_00000729

#### Part 2: Choosing gene clusters by their attributes

These parameters are used to exclude gene clusters that don't reach certain thresholds and are applies on top of filters already applied (for example, you can use these to exclude clusters within a specific bin). 

Here is a list of the different filters that you can use to exclude some subsection of your gene clusters:

- min/max number of genomes that the gene cluster occurs in. 
- min/max number of genes from each genome. For example, you could exclude clusters that don't appear in every genome 3 times, or get single-copy genes by setting `max-num-genes-from-each-genome` to 1. 
- min/max [geometric homogenity index](http://merenlab.org/2016/11/08/pangenomics-v2/#geometric-homogeneity-index) 
- min/max [functional homogenity index](http://merenlab.org/2016/11/08/pangenomics-v2/#functional-homogeneity-index)
- min/max combined homogenity index 

For example, the following run on a %(genomes-storage-db)s that contains 50 genomes will report only the single-copy core genes with a functional homogenity index above 0.25:

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s \
                                     --max-num-genes-from-each-genome 1 \
                                     --min-num-genomes-gene-cluster-occurs 50 \
                                     --min-functional-homogenity-index 0.25 
{{ codestop }}

You can also exclude genomes that are missing some number of the gene clusters that you're working with by using the paramter `--max-num-gene-clusters-missing-from-genome`. 

For each of these parameters, see the program's help menu for more information. 

### Fun with phylogenomics! 

To get a %(concatenated-gene-alignment-fasta)s (which you can use to run %(anvi-gen-phylogenomic-tree)s), use the parameter `--concatenate-gene-clusters`

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s \
                                     --concatenate-gene-clusters
{{ codestop }}

Here, you also have the option to specify a specific aligner (or list the available aligners), as well as provide a NEXUS formatted partition file, if you so choose. 
