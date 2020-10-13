This is an Anvi'o database that **stores information about your genomes, primarily for use in pangenomic analyses.**

You can think of it like this: in a way, a genomes-storage-db is to the [the pangenomic workflow](http://merenlab.org/2016/11/08/pangenomics-v2/#generating-an-anvio-genomes-storage) what a %(contigs-db)s is to the [the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). They both describe key information unique to your particular dataset and are required to run the vast majority of programs. 

### What kind of information? 

A genomes storage database contains information about the genomes that you inputted to create it, as well as the genes within them. 

Specifically, there are three tables stored in a genomes storage database: 

* A table describing the information about each of your genomes, such as their name, type (internal or external), GC content, number of contigs, completition, redunduncy, number of genes, etc. 
* A table describing the genes within your genomes. For each gene, this includes its gene caller id, associated genome and position, sequence, length, and whether or not it is partial. 
* A table describing the functions of your genes, including their sources and e-values. 

### Cool. How do I make one? 

You can generate one of these from an %(internal-genomes)s (genomes described in %(bin)ss), %(external-genomes)s (genomes described in %(contigs-db)ss), or both using the program %(anvi-gen-genomes-storage)s. 

### Cool cool. What can I do with one? 

With one of these, you can run %(anvi-pan-genome)s to get a %(pan-db)s. If a genomes storage database is the %(contigs-db)s of pangenomics, then a %(pan-db)s is the %(profile-db)s. It contains lots of information that is vital for analysis, and most programs will require both the %(pan-db)s and its genomes storage database as an input. 
