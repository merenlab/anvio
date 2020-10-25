This program uses the user's similarity metric of choice to calculate the distance between the input genomes. 

The currently available distance metrics include: 
- [PyANI](https://github.com/widdowquinn/pyani)) to calculate the average nucleotide identity (ANI) (i.e. what portion of orthologous gene pairs align)
-[fastANI](https://github.com/ParBLiSS/FastANI) also to calcualte the ANI but at a faster speed (at the drawback of a slight reduction in accuracy)
- [sourmash](https://sourmash.readthedocs.io/en/latest/) to calculate the mash distance between genomes

### Input/Output 

The expected input is any combination of %(external-genomes)s, %(internal-genomse)s, and text files that contains paths to %(fasta)s files that describe each of your genomes. This is a tab-delimited file with two columns (`name` and `path` to the fasta files, each of which is assumed to be a single genome),

You also have the option to provide a %(pan-db)s, in which case the output data will be written to the database as %(misc-data-layers)s and %(misc-data-layer-orders)s data. This was done in the [pangenomic tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#computing-the-average-nucleotide-identity-for-genomes-and-other-genome-similarity-metrics-too). Else, you'll need to provide an output path for the %(genome-similarity)s files. 

Here is an example run with pyANI from an %(external-genomes)s without any parameter changes:

    anvi-compute-genome-similarity -e %(external-genomes)s \
                                   -o path/for/%(genome-simliarty)s \ 
                                   --program pyANI

### Genome similarity metrics: parameters

#### pyANI

You have the option to change any of the follow parameters:

- The method used for alignment. The options are:
    -`ANIb` (default): uses [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)+ to align 1020 nt fragments of the inputs
    -`ANIm`: uses [MUMmer](http://mummer.sourceforge.net/) to align
    -`ANIblastall`: Uses legacy [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) to align 1020 nt fragments 
    -`TETRA`: caclulates tetranucleotide frequencies for each input 
    
- The minimum alignment fraction (all percent identity scores lower than this will be automatically set to 0). The default is 0. When using this, you can also use `--significant-alignment-length` to overwrite this discard for sequences with alignmens longer than a specific absolute length 
-Similarly, you can discard all results less than some full percent identity (which accounts for what total percent of the two sequences of interest aligned). 

#### fastANI

You can change any of the following fastANI parameters:

-The kmer size. The default is 16.     

-The fragement length. The default is 30. 

-The minimum number of fragments for a result to count. The default is 50. 

#### sourmash

You have the option to change the `kmer-size`. This value should depend on the relationship between your samples. The default is 31 ([as recommended by sourmash for genus-level distances](https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html), but we found that 13 most closely parallels that results from an ANI alignment. 

You can also set the compression ratio for your fasta files. Decreasing this from the default (1000) will decrease sensitivity. 

### Other Parameters

If requested (i.e. a %(pan-db)s is not provided), this program outputs similarity matrix files, which can be clustered into a %(dendrogram)s. You can choose to change the distance metric or linkage algorithm for hierarchical clustering. 

If your getting a lot of debug/output messages, you can turn them off with `--just-do-it` or helpfully store them into a file with `--log-file`. 

