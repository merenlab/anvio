This program **generates a %(genomes-storage-db)s, which stores information about your genomes, primarily for use in pangenomic analysis.** 

Genomes storage databases are to anvi'o's pangenomic workflow what a %(contigs-db)s is to a metagenomic workflow: they store vital information and are passed to most programs you'll want to run. 

Once you've generated a %(genomes-storage-db)s, you can run %(anvi-pan-genome)s, which creates a %(pan-db)s and runs various pangenomic analyses (including calculating similarities between your sequences, identifying gene clusters, and organizing your gene clusters and genomes). After that, you can display your pangenome with %(anvi-display-pan)s. For more information, check out [the pangenomic workflow](http://merenlab.org/2016/11/08/pangenomics-v2/#generating-an-anvio-genomes-storage).

### Inputs: internal and external genomes

You can initialize your genomes storage database with %(internal-genomes)s, %(external-genomes)s, or both. 

%(internal-genomes)s describe genomes that are represented by a %(bin)s within a %(collection)s that already exists within an anvi'o %(profile-db)s. For example, if you had gone through [the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/) and had several metagenome-assembled genomes (MAGs) that you wanted to run pangenomic analyses on. 

{{ codestart }}
anvi-gen-genomes-storage -i %(internal-genomes)s \
                         -o %(genomes-storage-db)s
{{ codestop }}

{:.notice}
The name of your genomes storage database (which follows the `-o` flag) must end with `-GENOMES.db`. This helps differentiate it from other types of anvi'o databases, such as the %(contigs-db)s and %(profile-db)s. 

In contrast, %(external-genomes)s describe genomes that are contained in a %(fasta)s file that you've converted into a %(contigs-db)s (using %(anvi-gen-contigs-database)s). For example, if you had downloaded genomes from [NCBI](https://www.ncbi.nlm.nih.gov/). 

{{ codestart }}
anvi-gen-genomes-storage -e %(external-genomes)s \
                         -o %(genomes-storage-db)s
{{ codestop }}

You can also create a genomes storage database from both types of genomes simultaneously. For example, if you had MAGs from a metagenomic analysis of an environmental sample and wanted to compare them with reference genomes from [NCBI](https://www.ncbi.nlm.nih.gov/). To execute this, simply provide both types of genomes as parameters:

{{ codestart }}
anvi-gen-genomes-storage -i %(internal-genomes)s \
                         -e %(external-genomes)s \
                         -o %(genomes-storage-db)s
{{ codestop }}

### Changing the gene caller

By default, anvi'o will use [Prodigal](https://github.com/hyattpd/Prodigal) and will inform you if you have gene calls identified by other gene callers. However, you are welcome to explicitly use a specific gene caller with the flag `--gene-caller`. 

If you're wondering what gene callers are available in your %(contigs-db)s, you can check by running the program %(anvi-export-gene-calls)s on a specific %(contigs-db)s with the flag `--list-gene-callers`. 
