This program computes the detection of genes (inputted as a %(bin)s) across your samples, so that you can visualize them in the %(interactive)s interface. 

This program is used in [the metapangenomic workflow](https://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes) on genes with metagenomes as samples to visually identify the environmental core genes and accessory genes. 

### Inputs  

Essentially, you provide a %(contigs-db)s and %(profile-db)s pair, as well as the %(bin)s you want to look at, and this program will  search each gene in your bin against the samples denoted in your %(profile-db)s: 

{{ codestart }}
anvi-script-gen-distribution-of-genes-in-a-bin -c %(contigs-db)s \ 
                                               -p %(profile-db)s \
                                               -C %(collection)s \
                                               -b %(bin)s 
{{ codestop }}

There are two other parameters that you can set to focus the genes that you're looking at: 
- The minimum detection required for a gene to be included (by default, a gene must have a detection value of `0.5` in at least one of your samples)
-The minimum coverage required for a gene to be included (by default, a gene must have a total coverage of `0.25` times the mean total coverage in your data) 

### Outputs

This program will produce two outputs: 

1. `[your bin name]-GENE-COVs.txt`, which is a %(view-data)s artifact. This is a matrix where each row represents a gene, each column represents one of your samples, and the cells each contain a coverage value. 
2. `[your bin name]-ENV-DETECTION.txt`, which is a %(misc-data-layers)s. It is a two-column file, where each row is a gene and and the second column describes whether or not that gene is systematically detected in your samples. Thus, this can be added as an additional layer in the interface that describes describes which genes are detected in your samples. (as an example, see the outermost layer [here](https://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes))

Thus, after running this program on a bin with name `BIN_NAME`, you can run 

{{ codestart }}
%(anvi-interactive)s -d BIN_NAME-GENE-COVs.txt \
                 -A BIN_NAME-ENV-DETECTION.txt \
                 --manual \
                 -p %(profile-db)s
{{ codestop }}                                                   

This will visually show you the coverage and detection of your genes across your samples in the %(interactive)s interface (simlarly to [this figure](https://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes)). 
