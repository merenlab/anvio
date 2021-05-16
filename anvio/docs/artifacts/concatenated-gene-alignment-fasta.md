This file **contains the alignment information for multiple genes across different organisms**.

Basically, a single gene alignment compares a single gene's sequence across multiple organisms. For example, you could align some specific rRNA sequence across all of the organisms in your sample. This alignment highlights both mutations and insertions and deletions (indicated with dashes). 

Clustal programs do a great job of visualizing this data, by color coding it. Here is an example from Anvi'o's pangenome display: 

![A lovely clustal-like alignment from the anvi'o pangenome display](../../images/example_alignment.png)

A concatenated gene alignment fasta contains multiple of these gene alignments, in order to generate a tree based off of multiple genes. 

This information can then be used to generate a phylogenomic tree using %(anvi-gen-phylogenomic-tree)s or through programs like [FastTree](http://www.microbesonline.org/fasttree/). 

In Anvi'o, this is an output of %(anvi-get-sequences-for-gene-clusters)s (for generating a tree based off of gene clusters in your workflow) as well as %(anvi-get-sequences-for-hmm-hits)s (for generating a tree based off of the genes that got HMM hits). 

