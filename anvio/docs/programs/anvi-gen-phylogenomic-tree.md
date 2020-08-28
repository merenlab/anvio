This program generates a NEWICK-formatted phylogenomic tree (see %(phylogeny)s) based on a given %(concatenated-gene-alignment-fasta)s. 

As mentioned in the [phylogenetics tutorial](http://merenlab.org/2017/06/07/phylogenomics/), it currently only has the option to use [FastTree](http://microbesonline.org/fasttree/) to do so, but be aware that there are many other programs that you can do this with. Some of the options we are familiar with (and are not yet represented in `anvi-gen-phylogenomic-tree`) include [MrBayes](http://mrbayes.sourceforge.net/), [MEGA](http://www.megasoftware.net/), and PHYLIP, [among many others](http://evolution.genetics.washington.edu/phylip/software.html#methods), most of which will happily take a %(concatenated-gene-alignment-fasta)s. 

Anyway, running this program is simple. Just provide the %(concatenated-gene-alignment-fasta)s with all of the genes that you want to use and the output file path for your %(phylogeny)s:

{{ codestart }}
anvi-gen-phylogenomic-tree -f %(concatenated-gene-alignment-fasta)s \
                           -o PATH/TO/%(phylogeny)s
{{ codestop }}
