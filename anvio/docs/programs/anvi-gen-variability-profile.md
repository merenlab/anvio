This program takes the variability data stored within a %(profile-db)s and compiles it from across samples into a single matrix that comprehensively describes your SNVs, SCVs or SAAVs (a %(variability-profile-txt)s).

This program is described on [this blog post](http://merenlab.org/2015/07/20/analyzing-variability/#the-anvio-way), so take a look at that for more details. 

## Let's talk parameters 

Here is a basic run with no bells or whisles: 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s
{{ codestop }}

You can add structural annotations by providing a %(structure-db)s. 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -s %(structure-db)s 
{{ codestop }}

### Focusing on a subset of the input 

You can focus on a specific %(collection)s, %(bin)s, genes (by providing a file or list of caller IDs) or list of splits (in the form of a %(splits-txt)s). 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             --gene-caller-ids GENE_1,GENE_2,GENE_3
{{ codestop }}

When providing a %(structure-db)s, you can also limit your analysis to only genes that have structures in your database. 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -s %(structure-db)s \
                             -C %(collection)s \
                             --only-if-structure
{{ codestop }}

You can also choose to look at only data from specific samples by providing a file with one sample name per line. For example

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             --samples-of-interest my_samples.txt
{{ codestop }}

where `my_samples.txt` looks like this:

    DAY_17A
    DAY_18A
    DAY_22A
    
### SNVs vs. SCVs vs. SAAVs 

Which one you're analyzing depends entirely on the `engine` parameter, which you can set to `NT` (nucleotides), `CDN` (codons), or `AA` (amino acids). The default value is nucleotides. Note that to analyze SCVs or SAAVs, you'll have needed to use the flag `--profile-SCVs` when you ran %(anvi-profile)s or %(anvi-merge)s. 

For example, to analyze SAAVs, run 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -s %(structure-db)s \
                             --engine AA
{{ codestop }}

When analyzing single codon variants, you can choose to skip computing synonymity to save on run time, as so: 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -s %(structure-db)s \
                             --engine CDN \
                             --skip-synonymity
{{ codestop }}

### Filtering the output 

You can filter the output in various ways, so that you can get straight to the variability positions that you're most interested in. Here are some of the filters that you can set:

* The maximum number of variable positions that can come from a single split (e.g. to look at a max of only two random SCVs from each split)
* The maximum and minimum departure from the reference or consensus position
* The minimum coverage value in all samples (if a position is covered less than that value in a even single sample, it will not be reported)

### Adding additional information

You can also set `--quince-mode`, which reports the variability data across all samples for each position reported (even if that position isn't variable in some samples). For example, if nucleotide position 34 of contig 1 was a SNV in one sample, the output would contain the data for nucleotide position 34 for all of your samples. 

You can also ask the program to report the contig names, split names, and gene-level coverage statistics. 

