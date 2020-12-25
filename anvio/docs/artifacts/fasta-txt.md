This is a file used by %(anvi-run-workflow)s that lists the name and path of all of the input %(fasta)s files. 

As of now, this file is used in the %(contigs-workflow)s, %(pangenomics-workflow)s, and [the reference mode](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#references-mode) of the %(metagenomics-workflow)s.

This file will look something like this: 

    name        path
    SAMPLE_01   path/to/sample_01.fa
    SAMPLE_02   path/to/sample_02.fa
    
Note that the paths can be either absolute or relative, and the fasta files can be either compressed or not. That's all up to you. 

To input more information, for each file you can also specify an %(external-gene-calls)s file and/or a %(functions-txt)s. Just provide those with additional columns, as so: 

    name        path                    external_gene_calls             gene_functional_annotation
    SAMPLE_01   path/to/sample_01.fa    %(external-gene-calls)s_01.txt  %(functions-txt)s_01.txt
    SAMPLE_02   path/to/sample_02.fa    %(external-gene-calls)s_02.txt  %(functions-txt)s_02.txt

For more information, check out the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)
