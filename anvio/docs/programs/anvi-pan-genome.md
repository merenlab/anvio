This program **finds, clusters, and organizes the genes** within a %(genomes-storage-db)s to create a %(pan-db)s. 

This is the program that does the brunt of the work when running a pangenomic workflow. Check out [the pangenomic tutorial](http://merenlab.org/2016/11/08/pangenomics-v2) for a more in-depth overview of the contents of this page and the capabilities of a %(pan-db)s. 

### Before running this program

Before running this program, you'll want to make sure your dependencies are all set, since this program requires some aditional dependencies than base anvi'o. If the following command runs without errors, then you're all good. 

{{ codestart }}
anvi-self-test --suite pangenomics
{{ codestop }}

If that command doesn't run smoothly, check out [this page](http://merenlab.org/2016/11/08/pangenomics-v2/#dependencies).

### What this program does

This program finds and organizes your gene clusters to give you all of the data that is displayed when you run %(anvi-pan-genome)s. Almost all of the work described in [this gif that explains the common steps involved in pangenomics](http://merenlab.org/momics/#pangenomics) is done by this program. 

In a little more detail, this program will do three major things for you:

* Calculate the similarity between the all of the gene calls in all of the genomes in your %(genomes-storage-db)s. By default this uses [DIAMOND](https://www.wsi.uni-tuebingen.de/lehrstuehle/algorithms-in-bioinformatics/software/diamond/) to do this, but Meren strongly recommends that you use the `--use-ncbi-blast` flag to use [`blastp`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) instead.  

    *   When doing this, this will look at every genome in your %(genomes-storage-db)s (unless you use `--genome-names`) and will use every gene call, whether or not they are complete (unless you used `--exclude-partial-gene-calls`).   
    
    *   After doing this, it will use the minbit heuristic (originally from ITEP ([Benedict et al., 2014](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-8)) to throw out weak matches. This removes a lot of noise before clustering. 
    
* Use the [MCL](http://micans.org/mcl/) algorithm to identify clusters in your search results.  

* Organize your gene clusters and genomes using their `euclidean` distance and `ward` linkage. 

This program is very smart, and if you're already run it, it will try to use the data that it's already calculated. This way you can change smaller parameters without all of the run time. However, this also means you need to tell it to rerun the process (if that's what you want) with the flag `--overwrite-output-destinations`. 

### Cool. How about some examples and specific parameters?

Who doesn't love a good example? The simplest way to run this is as follows:

{{ codestart }}
anvi-pan-genome -g %(genomes-storage-db)s
{{ codestop }}

But there are many parameters you can alter to your liking. For example, here's a run that specifies that it wants to use NCBI's blastp to find sequence similarities and muscle to align genes and defines its output 

{{ codestart }}
anvi-pan-genomes -g %(genomes-storage-db)s \
                 --align-with muscle \
                 --use-ncbi-blast \ 
                 -n MY_PROJECT_NAME \
                 --description description.txt \
                 -o PATH/TO/%(pan-db)s 
{{ codestop }}

Here's another example that only looks at the complete gene calls within a subset of the genomes, eliminates gene clusters that only have hits in a single genome, and uses DIAMOND but with the sensitive setting enabled:

{{ codestart }}
anvi-pan-genomes -g %(genomes-storage-db)s \
                 -n MY_PROJECT_NAME \
                 --genome-names GENOME_1,GENOME_2,GENOME_3 \
                 --exclude-partial-gene-calls \ 
                 --min-occurance 2 \
                 --sensitive \
                 -o PATH/TO/%(pan-db)s 
{{ codestop }}

Some other parameters available to you allow you to  

- Change the minimum minbit value to elimate weak matches. The default value is 0.5.
- Change the MCL inflation parameter. The default value is 2. 
- Specify a minimum percent identity between two sequences to give that link an edge in MCL clustering. 
- Skip or speed up the calculation of homogeneity values for your clusters
- Enable multithreading with `-T`
