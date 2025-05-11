
This program takes the variability data stored within a %(profile-db)s and compiles it from across samples into a single matrix that comprehensively describes your SNVs, SCVs or SAAVs (a %(variability-profile-txt)s).  

This program is described on [this blog post](http://merenlab.org/2015/07/20/analyzing-variability/#the-anvio-way), so take a look at that for more details. 

## Let's talk parameters 

Here is a basic run with no bells or whisles: 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \ 
                             -C DEFAULT \
                             -b EVERYTHING
{{ codestop }}

Note that this program requires you to specify a subset of the databases that you want to focus on, so to focus on everything in the databases, run %(anvi-script-add-default-collection)s and use the resulting %(collection)s and %(bin)s, as shown above. 

You can add structural annotations by providing a %(structure-db)s. 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C DEFAULT \
                             -b EVERYTHING \
                             -s %(structure-db)s 
{{ codestop }}

You can also output your %(variability-profile-txt) to a specific location, which can be useful when working with multiple `engine` parameters.

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --output-file /path/to/your/variability.txt
{{ codestop }}

### Focusing on a subset of the input 

Instead of focusing on everything (providing the collection `DEFAULT` and the bin `EVERYTHING`), there are three ways to focus on a subset of the input: 

1. Provide a list of gene caller IDs (as a parameter with the flag `--gene-caller-ids` as shown below, or as a file with the flag `--genes-of-interest`)

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 --gene-caller-ids 1,2,3
    {{ codestop }}

2. Provide a %(splits-txt)s to focus only on a specific set of splits. 

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \
                                 --splits-of-intest %(splits-txt)s
    {{ codestop }}
    
3. Provide some other %(collection)s and %(bin)s. 

    {{ codestart }}
    anvi-gen-variability-profile -p %(profile-db)s \
                                 -c %(contigs-db)s \ 
                                 -C %(collection)s \
                                 -b %(bin)s
    {{ codestop }}

### Additional ways to focus the input 

When providing a %(structure-db)s, you can also limit your analysis to only genes that have structures in your database. 

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -s %(structure-db)s \
                             --only-if-structure
{{ codestop }}

You can also choose to look at only data from specific samples by providing a file with one sample name per line. For example

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --samples-of-interest my_samples.txt
{{ codestop }}

where `my_samples.txt` looks like this:

{{ codestart }}
DAY_17A
DAY_18A
DAY_22A
...
{{ codestop }}

### SNVs vs. SCVs vs. SAAVs 

Which one you're analyzing depends entirely on the `engine` parameter, which you can set to `NT` (nucleotides), `CDN` (codons), or `AA` (amino acids). The default value is nucleotides. Note that to analyze SCVs or SAAVs, you'll have needed to use the flag `--profile-SCVs` when you ran %(anvi-profile)s.

For example, to analyze SAAVs, run

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine AA
{{ codestop }}

To analyze SCVs, run

{{ codestart }}
anvi-gen-variability-profile -p %(profile-db)s \
                             -c %(contigs-db)s \
                             -C %(collection)s \
                             -b %(bin)s \
                             --engine CDN
{{ codestop }}

### Filtering the output 

You can filter the output in various ways, so that you can get straight to the variability positions that you're most interested in. Here are some of the filters that you can set:

* The maximum number of variable positions that can come from a single split (e.g. to look at a max of 100 SCVs from each split, randomly sampled)
* The maximum and minimum departure from the reference or consensus
* The minimum coverage value in all samples (if a position is covered less than that value in _one_ sample, it will not be reported for _all_ samples)


### --quince-mode

You can also set `--quince-mode`, which reports the variability data across all samples for each position reported (even if that position isn't variable in some samples). For example, if nucleotide position 34 of contig 1 was a SNV in one sample, the output would contain data for nucleotide position 34 for all of your samples. 

### --kiefl-mode

The default behavior is to report codon/amino-acid frequencies only at positions where variation was reported during profiling (which by default uses some heuristics to minimize the impact of error-driven variation). Fair enough, but for some diabolical cases, you may want to report _even_ invariant positions. When this flag is used, all positions are reported, regardless of whether they contained variation in any sample. The reference codon for all such entries is given a codon frequency of 1. All other entries (aka those with legitimate variation to be reported) remain unchanged. This flag can only be used with `--engine AA` or `--engine CDN` and is incompatible wth `--quince-mode`.

This flag was added in this [pull request](https://github.com/merenlab/anvio/pull/1794) where you can read about all of the tests that were performed to ensure this mode is behaving properly.

### Adding additional information

You can also ask the program to report the contig names, split names, and gene-level coverage statistics, which appear as additional columns in the output.


