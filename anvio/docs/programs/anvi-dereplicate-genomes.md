This program uses the specified similarity metric to identify genomes that are highly similar to each other and groups them into redundant clusters. The program identifies representative sequences for each cluster and outputs them as %(fasta)s files.


#### Input Options 

You have two options for providing input to this program: 

- the results of %(anvi-compute-genome-similarity)s (a %(genome-similarity)s directory). If you used `fastANI` or `pyANI` when you ran %(anvi-compute-genome-similarity)s, provide this using the parameter `--ani-dir`; if you used sourmash, use the parameter `--mash-dir`. 

- an %(internal-genomes)s, %(external-genomes)s or a series of %(fasta)s files (each representing a genome), in which case anvi'o will run %(anvi-compute-genome-similarity)s for you. When providing these inputs, you can also specify any of the parameters that %(anvi-compute-genome-similarity)s accepts, including the `--program` you want to use (from [PyANI](https://github.com/widdowquinn/pyani), [fastANI](https://github.com/ParBLiSS/FastANI), or [sourmash](https://sourmash.readthedocs.io/en/latest/)) and their associated parameters. Details about all available options can be found in the help menu for %(anvi-compute-genome-similarity)s.

#### Output Format 

By default, the output of this program is a directory containing two descriptive text files (the cluster report and fasta report) and a subdirectory called `GENOMES`:

- The cluster report is a tab-delimited text file where each row describes a cluster. This file contains four columns: the cluster name, the number of genomes in the cluster, the representative genome of the cluster, and a list of the genomes that belong to the cluster. Here is an example describing 11 genomes organized into three clusters:

|**cluster**|**size**|**representative**|**genomes**|
|:--|:--|:--|:--|
|cluster_000001|1|G11_IGD_MAG_00001|G11_IGD_MAG_00001|
|cluster_000002|8|G11_IGD_MAG_00012|G08_IGD_MAG_00008,G33_IGD_MAG_00011,G01_IGD_MAG_00013,G06_IGD_MAG_00023,G03_IGD_MAG_00021,G05_IGD_MAG_00014,G11_IGD_MAG_00012,G10_IGD_MAG_00010|
|cluster_000003|2|G03_IGD_MAG_00011|G11_IGD_MAG_00013,G03_IGD_MAG_00011|

- The subdirectory `GENOMES` contains fasta files describing the representative genome from each cluster. For example, if your original set of genomes contained two identical genomes, this program would cluster them together, and the `GENOMES` folder would only include one of their sequences. 

- The fasta report describes the fasta files contained in the subdirectory `GENOMES`. By default, this describes the representative sequence of each final cluster. It specifies the genome name, its source, its cluster (and the representative sequence of that cluster), and the path to its fasta file in `GENOMES`. For the example above, the fasta report would appear as follows:

|**name**|**source**|**cluster**|**cluster_rep**|**path**|
|:--|:--|:--|:--|:--|
|G11_IGD_MAG_00001|fasta|cluster_000001|G11_IGD_MAG_00001|GENOMES/G11_IGD_MAG_00001.fa|
|G11_IGD_MAG_00012|fasta|cluster_000002|G11_IGD_MAG_00012|GENOMES/G11_IGD_MAG_00012.fa|
|G03_IGD_MAG_00011|fasta|cluster_000003|G03_IGD_MAG_00011|GENOMES/G03_IGD_MAG_00011.fa|

You can also choose to report all genome fasta files (including redundant genomes) (with `--report-all`) or report no fasta files (with `--skip-fasta-report`). This would modify the fasta files included in `GENOMES` and the genomes mentioned in the fasta report. The cluster report would remain unchanged.

#### Required Parameters and Example Runs

You are required to set the threshold for two genomes to be considered redundant and placed in the same cluster. 

For example, if you had the results from an %(anvi-compute-genome-similarity)s run where you used `pyANI` and wanted the threshold to be 90 percent, you would execute: 

{{ codestart }}
anvi-dereplicate-genomes --ani-dir %(genome-similarity)s \ 
                         -o path/to/output \
                         --program pyANI \
                         --similarity-threshold 0.90
{{ codestop }}

If instead you had not yet run %(anvi-compute-genome-similarity)s and wanted to cluster the genomes in your %(external-genomes)s file with 85 percent or greater similarity (without generating fasta files) using sourmash, you could execute: 

{{ codestart }}
anvi-dereplicate-genomes -e %(external-genomes)s \ 
                         --skip-fasta-report \
                         --program sourmash \
                         -o path/to/output \
                         --similarity-threshold 0.85 
{{ codestop }}

#### Other parameters

You can modify how anvi'o selects the representative sequence from each cluster using the `--representative-method` parameter. Three options are available:

- `Qscore`: selects the genome with the highest completion and lowest redundancy
- `length`: selects the longest genome in the cluster
- `centrality` (default): selects the genome with the highest average similarity to every other genome in the cluster

Additional options include skipping genome hash verification (which will warn you if you have identical sequences in separate genomes with different names), providing a log path for debug messages, or using multithreading (relevant only if not providing `--ani-dir` or `--mash-dir`).


