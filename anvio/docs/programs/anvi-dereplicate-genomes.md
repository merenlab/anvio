This program uses the user's similarity metric of choice to identify genomes that are highly similar to each other. 

This is done by analyzing the results of %(anvi-compute-genome-similairty)s to cluster your input genomes based on the provided threshold. By default, the program will report one fasta file with a representative sequence for each cluster. 

#### Input Options 

You have two options for the input to this program: 
- the results of %(anvi-compute-genome-similarity)s (a %(genome-similarity)s directory). If you used `fastANI` or `pyANI` when you ran %(anvi-compute-genome-similarity)s, provide this using the parameter `--ani-dir`; if you used sourmash, use the parameter `--mash-dir`. 
- an %(internal-genomes)s, %(external-genomes)s or a series of %(fasta)s files (each of which represents a genome), in which case anvi'o will run %(anvi-compute-genome-similarity)s for you.  When providing these inputs, you can also provide any of the parameters that %(anvi-compute-genome-similarity)s can take, including the `--program` you want to use (out of  [PyANI](https://github.com/widdowquinn/pyani), [fastANI](https://github.com/ParBLiSS/FastANI),  [sourmash](https://sourmash.readthedocs.io/en/latest/)) and their parameters. Details about all of this can be found in the help menu for %(anvi-compute-genome-similiarty)s.

#### Output Options 

The output of this program is a directory containing {stuff}. You can also choose to report all genome fasta files (including from redundant genomes) (with `--report-all`) or report no fasta files (with `--skip-fasta-report`). 

#### Required Parameters and Example Runs

You are required to set the threshold for two genomes to be considered redundant and put in the same cluster. 

For example, if you had the results from an %(anvi-compute-genome-similarity)s run where you had used `PyANI` and wanted the threshold to be 90 percent, you would run: 

{{ codestart }}
anvi-dereplictate-genomes --ani-dir %(genome-similiarty)s \ 
                          -o path/to/output \
                          --similiarity-threshold 0.90
{{ codestop }}

If instead you hadn't yet run %(anvi-compute-genome-similarity)s and instead wanted to cluster the genomes in your %(external-genomes)s file with similarity 85 percent or more (no fasta files necessary) using sourmash, you could run: 

{{ codestart }}
anvi-dereplictate-genomes -e %(external-genomes)s \ 
                          --skip-fasta-report \
                          --program sourmash \
                          -o path/to/output \
                          --similiarity-threshold 0.85 
{{ codestop }}

#### Other parameters

You can change how anvi'o picks the representative sequence from each cluster with the parameter `--representative-method`. For this you have three options:

- `Qscore`: picks the genome with highest completition and lowest redundency 
-`length`: picks the longest genome in the cluster
-`Centrality` (default): picks the genome with highest average similiarty to every other genome in the cluster

You can also choose to skip checking genome hashes (which will warn you if you have identical sequences in separate genomes with different names), provide a log path for debug messages or use multithreading. 

