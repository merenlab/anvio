This aptly-named program **gets the sequences for the gene clusters stored in a %(pan-db)s and returns them as either a %(genes-fasta)s or a %(concatenated-gene-alignment-fasta)s**, which can directly go into the program %(anvi-gen-phylogenomic-tree)s for phylogenomics. This gives you advanced access to your gene clusters, so you can take them out of anvi'o and do whatever you please with them.

The program parameters also include a large collection of advanced filtering options. Using these options you can scrutinize your gene clusters in creative and precise ways. Using the combination of these filters you can focus on single-copy core gene clusters in a pangenome, or those occur only as singletons, or paralogs that contain more than a given number of sequences, and so on. Once you are satisfied with the output a given set of filters generate, you can add the matching gene clusters a %(misc-data-items)s with the flag `--add-into-items-additional-data-table`, which can be added to the %(interactive)s interface as additional layers when you visualize your %(pan-db)s using the program %(anvi-display-pan)s 

By default, %(anvi-get-sequences-for-gene-clusters)s will generate a single output file. But you can ask the program to report every gene cluster that match to your filters as a separate FASTA file depending on your downstream analyses.

While the number of parameters this powerful program can utilize may seem daunting, many of the options just help you specify exactly from which gene clusters you want to get sequences. 

### Running on all gene clusters

Here is an example that shows the simplest possible run, which will export sequences for every single gene cluster found in the %(pan-db)s as amino acid sequences:

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s
{{ codestop }}

{:.notice}
The program will report the DNA sequences if the flag `--report-DNA-sequences` is used.

### Splitting gene clusters into their own files

The command above will put all gene cluster sequences in a single output %(fasta)s file. If you would like to report each gene cluster in a separate FASTA file, it is also an option thanks to the flag `--split-output-per-gene-cluster`. This optional reporting throught this flag applies to all commands shown on this page. For instance, the following command will report every gene cluster as a separate FASTA file in your directory,

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     --split-output-per-gene-cluster
{{ codestop }}

where the output files and paths will look like this in your work directory:

```
GC_00000001.fa
GC_00000002.fa
GC_00000003.fa
GC_00000004.fa
GC_00000005.fa
GC_00000006.fa
GC_00000007.fa
GC_00000008.fa
GC_00000009.fa
GC_00000010.fa
(...)
```

You can use the parameters `--output-file-prefix` to add file name prefixes to your output files. For instance, the following command,

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     --split-output-per-gene-cluster \
                                     --output-file-prefix MY_PROJECT
{{ codestop }}

will result in the following files in your work directory:

```
MY_PROJECT_GC_00000001.fa
MY_PROJECT_GC_00000002.fa
MY_PROJECT_GC_00000003.fa
MY_PROJECT_GC_00000004.fa
MY_PROJECT_GC_00000005.fa
MY_PROJECT_GC_00000006.fa
MY_PROJECT_GC_00000007.fa
MY_PROJECT_GC_00000008.fa
MY_PROJECT_GC_00000009.fa
MY_PROJECT_GC_00000010.fa
(...)
```

You can also use the parameter `--output-file-prefix` to store files in different directories. For instance, the following command (note the trailing `/` in the `--output-file-prefix`),

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     --split-output-per-gene-cluster \
                                     --output-file-prefix A_TEST_DIRECTORY/
{{ codestop }}

will result in the following files:

```
A_TEST_DIRECTORY/GC_00000001.fa
A_TEST_DIRECTORY/GC_00000002.fa
A_TEST_DIRECTORY/GC_00000003.fa
A_TEST_DIRECTORY/GC_00000004.fa
A_TEST_DIRECTORY/GC_00000005.fa
A_TEST_DIRECTORY/GC_00000006.fa
A_TEST_DIRECTORY/GC_00000007.fa
A_TEST_DIRECTORY/GC_00000008.fa
A_TEST_DIRECTORY/GC_00000009.fa
A_TEST_DIRECTORY/GC_00000010.fa
(...)
```

In contrast, the following command,

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     --split-output-per-gene-cluster \
                                     --output-file-prefix A_TEST_DIRECTORY/A_NEW_PREFIX
{{ codestop }}

will result in the following files:

```
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000001.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000002.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000003.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000004.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000005.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000006.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000007.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000008.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000009.fa
A_TEST_DIRECTORY/A_NEW_PREFIX_GC_00000010.fa
(...)
```

### Exporting only specific gene clusters

#### Part 1: Choosing gene clusters by collection, bin, or name

You can export only the sequences for a specific %(collection)s or %(bin)s with the parameters `-C` or `-b` respectively. 

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s \
                                     -C %(collection)s 
{{ codestop }}

{:.notice}
You can always display the collections and bins available in your %(pan-db)s by adding `--list-collections` or `--list-bins` flags to your command.

Alternatively, you can export the specific gene clusters by name, either by providing a single gene cluster ID or a file with one gene cluster ID per line. For example: 

{{ codestart }}
anvi-get-sequences-for-gene-clusters -g %(genomes-storage-db)s \
                                     -p %(pan-db)s \
                                     -o %(genes-fasta)s \
                                     --gene-cluster-ids-file gene_clusters.txt
{{ codestop }}

where `gene_clusters.txt` contains the following:

```
GC_00000618
GC_00000643
GC_00000729
```

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
