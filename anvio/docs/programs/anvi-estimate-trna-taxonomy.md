This program **uses the taxonomy associates of your tRNA sequences to estimate the taxonomy for genomes, metagenomes, or %(collection)s stored in your %(contigs-db)s**. 

This is the final step in the trna-taxonomy workflow. Before running this program, you'll need to have run %(anvi-run-trna-taxonomy)s on the %(contigs-db)s that you're inputting to this program.

## Input options 

### 1: Running on a single genome

By default, this program will assume that your %(contigs-db)s contains only a single genome and will determine the taxonomy of that single genome.   

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s
{{ codestop }}

This will give you only the best taxonomy hit for your genome based on your tRNA data. If you want to look under the hood and see what results from %(anvi-run-trna-taxonomy)s it's using to get there, add the `--debug` flag. 

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s \
                           --debug 
{{ codestop }}

### 2: Running on a metagenome

In metagenome mode, this program will assume that your %(contigs-db)s contains multiple genomes and will try to give you an overview of the taxa within it.  To do this, anvi'o will determine which anticodon has the most hits in your contigs (for example `GGG`), and then will look at the taxnomy hits for tRNA with that anticodon across your contigs. 

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s \
                           --metagenome-mode 
{{ codestop }}

If instead you want to look at a specific anticodon, you can specify that with the `-S` parameter. For example, to look at `GGT`, just run the following: 

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s \
                           --metagenome-mode \
                           -S GGT
{{ codestop }}

### 3: Running on multiple metagenomes

You can use this program to look at multiple metagenomes by providing a %(metagenomes)s artifact. This is useful to get an overview of what kinds of taxa might be in your metagenomes, and what kinds of taxa they share. 

Running this

{{ codestart }}
anvi-estimate-trna-taxonomy --metagenomes %(metagenomes)s \
                           --output-file-prefix EXAMPLE
{{ codestop }}

will give you an output file containing all taxonomic levels found and their coverages in each of your metagenomes, based on their tRNA. 

### 4: Estimating the taxonomy of bins 

You can use this program to estimate the taxonomy of all of the %(bin)ss in a %(collection)s by providing the %(collection)s and the associated %(profile-db)s. 

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s \
                           --C %(collection)s  \
                           --p %(profile-db)s 
{{ codestop }}

When doing this, you can also put the final results into your %(profile-db)s as a %(misc-data-layers)s with the flag `--update-profile-db-with-taxonomy`

### 5: I don't even have a contigs-db. Just a fasta file. 

This program can run the entire ad hoc sequence search without a %(contigs-db)s involved (just a fasta and number of target sequences as a percent of the total; default: 20 percent), but this is not recommended. However, if you provide other parameters, they will be ignored. 

{{ codestart }}
anvi-estimate-trna-taxonomy --dna-sequence %(fasta)s \
                           --max-num-target-sequences 10
{{ codestop }}

## The Output

Now that you've inputted your desired inputs, you think about whether you want an output and what it will look like. By default, this program won't give you an output (just %(genome-taxonomy)s information in your %(contigs-db)s. However, if you add any of these output options, it will instead produce a %(genome-taxonomy-txt)s. 

### Anticodon Frequencies

If you want to look at the anticodon frequencies before getting taxonomy info at all (for example because you can't decide which anticodon to use for input option 2), add the flag `--report-anticodon-frequencies`. This will report the anticodon frequencies to a tab-delimited file and quit the program. 

### A single output 

To get a single output (a fancy table for your viewing pleasure), just add the output file path. 

In this example, the input will be a single %(contigs-db)s (input option 1), 

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s \
                           -o path/to/output.txt  
{{ codestop }}

This will give you a tab-delimited matrix with all levels of taxonomic information for the genome stored in your %(contigs-db)s. Specifically, the output is a %(genome-taxonomy-txt)s. 

If you want to focus on a single taxonomic level, use the parameter `--taxonomic-level`, like so:

{{ codestart }}
anvi-estimate-trna-taxonomy -c %(contigs-db)s \
                           -o path/to/output.txt  \
                           --taxonomic-level genus 
{{ codestop }}

You can also simplify the taxonomy names in the table with the flag `--simplify-taxonomy-information`

If you're running on a %(profile-db)s, you can also choose to add the anticodon coverage to the output with `--compute-anticodon-coverages`. 

### Multiple outputs

If you have multiple outputs (i.e. you are looking at multiple metagenomes (input option number 3) or you are looking at each anticodon individually with `--per-anticodon-output-file`), you should instead provide a output filename prefix.  

{{ codestart }}
anvi-estimate-trna-taxonomy --metagenomes %(metagenomes)s \
                           --output-file-prefix EXAMPLE
{{ codestop }}

The rest of the options listed for the single output (i.e. focusing on a taxonomic level, simplifying taxonomy information, etc.) still apply. 
