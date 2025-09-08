This program lets you export selections of your %(contigs-db)s around all occurances of a user-defined anchor gene. 

The output of this is a folder that contains a separate %(contigs-db)s for the region around each hit of the anchor gene. (In fact, you'll get a FASTA file, %(contigs-db)s, %(profile-db)s, and a copy of the runlog).

For example, you could specify the recognition site for a specific enzyme and use this program to pull out all potential sites where that enzyme could bind. 

### Required Parameters

You'll need to provide a %(contigs-db)s (of course), as well as the name of the output directory and a prefix to use when naming all of the output databases. 

You can define the region of interest either by defining the two flanking genes or by searching for an anchor gene and defining a number of genes around this gene that you want to look at. For example, if you set `num-genes` as 1, then each locus will contain the gene of interest, a gene upstream of it, and a gene downstream of it, for a total of three genes. 

### Defining the region of interest

There are four ways to indicate the desired anchor gene:

1. Provide a search term in the functional annotations of all of your genes. (If you're trying to find a gene with a vague function, you might want to use %(anvi-search-functions)s to find out which genes will show up first. Alternatively, you can you %(anvi-export-functions)s to look at a full list of the functional annotaitons in this database). 

    {{ codestart }}
    anvi-export-locus -c %(contigs-db)s \
                      --num-genes 2 \
                      -o GLYCO_DIRECTORY \
                      -O Glyco \
                      --search-term "Glycosyltransferase involved in cell wall bisynthesis" \ 
    {{ codestop }}
    
    You also have the option to specify an annotation source with the flag `--annotation source`

2.  Provide a specific gene caller ID. 

    {{ codestart }}
    anvi-export-locus -c %(contigs-db)s \
                      --num-genes 2 \
                      -o output_directory \
                      -O GENE_1 \
                      --gene-caller-ids 1
    {{ codestop }}

3. Provide a search term for the HMM source annotations. To do this, you must also specify an hmm-source. (You can use the flag `--list-hmm-sources` to list the available sources). 

    {{ codestart }}
    anvi-export-locus -c %(contigs-db)s \
                      --num-genes 2 \
                      -o Ribosomal_S20p \
                      -O Ribosomal_S20p \
                      --use-hmm \
                      --hmm-source Bacteria_71 \
                      --search-term Ribosomal_S20p
    {{ codestop }}
    
    4. Run in `flank-mode` and provide two flanking genes that define the locus region.
    
    {{ codestart }}
    anvi-export-locus -c %(contigs-db)s \
                      --flank-mode \
                      -o locus_output \
                      -O gyclo_to_acyl \
                      --search-term "Glycosyltransferase involved in cell wall bisynthesis","Acyl carrier protein" \ 
    {{ codestop }}

### Additional Options 

You can also remove partial hits, ignore reverse complement hits, or overwrite all files in a pre-existing output. 
