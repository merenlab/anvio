This program calculates similarity between input genomes using the specified similarity metric.

The currently available programs for calculating similarity metrics can be selected using the `--program` parameter:
- [PyANI](https://github.com/widdowquinn/pyani)) to calculate the average nucleotide identity (ANI) (i.e., what proportion of orthologous gene pairs align)
- [fastANI](https://github.com/ParBLiSS/FastANI) also to calculate the ANI but at faster speeds (with a slight reduction in accuracy)
- [sourmash](https://sourmash.readthedocs.io/en/latest/) to calculate the mash distance between genomes. Though this option is provided, we do not recommend using sourmash for genome comparisons as it excels at other tasks; it remains available as a legacy option.

### Input/Output

The expected input is any combination of %(external-genomes)s, %(internal-genomes)s, and text files containing paths to %(fasta)s files that describe each of your genomes. This is a tab-delimited file with two columns (`name` and `path` to the fasta files, each of which is assumed to represent a single genome).


The program outputs a directory containing %(genome-similarity)s data. The specific contents depend on how similarity scores are computed (specified with `--program`), but generally include tab-separated files of similarity scores between genomes and related metrics.


You also have the option to provide a %(pan-db)s, in which case the output data will additionally be stored in the database as %(misc-data-layers)s and %(misc-data-layer-orders)s data. This approach was demonstrated in the [pangenomic tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#computing-the-average-nucleotide-identity-for-genomes-and-other-genome-similarity-metrics-too).  

Here is an example run with pyANI using %(external-genomes)s with default parameters: 

{{ codestart }}
anvi-compute-genome-similarity -e %(external-genomes)s \
                               -o path/for/%(genome-similarity)s \
                               --program pyANI
{{ codestop }}

### Genome similarity metrics: parameters

Parameters have been organized according to the `--program` you select.

#### pyANI

You can modify any of the following parameters:

- The method used for alignment. The available options are:
    - `ANIb` (default): uses [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)+ to align 1020 nt fragments of the inputs
    - `ANIm`: uses [MUMmer](http://mummer.sourceforge.net/) to perform alignment
    - `ANIblastall`: Uses legacy [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) to align 1020 nt fragments
    - `TETRA`: Alignment-free approach. This calculates similarity scores by comparing tetranucleotide frequencies for each input

- The minimum alignment fraction (all percent identity scores lower than this threshold will be set to 0). The default is 0.


- If you want to retain alignments that are long despite not passing the minimum alignment fraction filter, you can specify a `--significant-alignment-length` to override `--min-alignment-fraction`.


- Similarly, you can discard all results below a specified full percent identity threshold (percent identity of aligned segments × aligned fraction).


#### fastANI

You can modify any of the following fastANI parameters:

* The k-mer size. The default is 16.

* The fragment length. The default is 30.

* The minimum number of fragments required for a result to be considered valid. The default is 50.

#### sourmash

You can modify the `kmer-size` parameter. This value should depend on the evolutionary relationships between your samples. The default is 31 ([as recommended by sourmash for genus-level distances](https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html)), but we found that 13 most closely parallels the results from ANI alignment.  

You can also set the compression ratio for your fasta files. Decreasing this from the default (1000) will decrease sensitivity.  

### Other Parameters 

Once calculated, the similarity matrix is used to create dendrograms via hierarchical clustering, which are stored in the output directory (and in the %(pan-db)s, if provided). You can modify the distance metric or linkage algorithm used for this clustering process.


If you are receiving excessive debug or output messages, you can disable them with `--just-do-it` or redirect them to a file with `--log-file`.



