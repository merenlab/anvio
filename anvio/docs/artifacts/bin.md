A bin is, in its simplest form, **a group of contigs**.  (Think of a literal bin that you're putting data into.)

In Anvi'o, you'll most commonly work with bins both in this form and in the form of %(internal-genomes)s, espeically when you want to work with bins contained in more than one %(profile-db)s. A group of bins is called a %(collection)s.

## What can you use bins for?

### Bins are all over 'omics
In general, you can also use bins and %(collection)ss to limit what you are analyzing downstream. A ton of anvi'o programs are able to take in a bin or a collection (a group of bins) so that you don't have to analyze your entire %(contigs-db)s when you just want to look at one section. This is true of many, many anvi'o analyses you can run, from %(anvi-estimate-metabolism)s to %(anvi-get-sequences-for-hmm-hits)s.

Below are some more specific examples of binning in action!

### Metagenomic binning
A common use of binning is in **metagenomics**. See [this tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/) for details, but essentially, the data you're working with in metagenomics is a sample of all of the genetic material in an environmental sample. It's like pulling a bunch of random DNA fragments out of a bucket of ocean water. If you want to try to rebuild the individual genomes from that mess, one common strategy is to piece together genomes de novo by trying to group the contigs together. This is called genome-resolved metagenomic binning.

Basically, in metagenomic binning, you're trying to group together a bunch of contigs that all belong to the same genome using various metrics like tetranucleotide frequency, differential coverage, completion, etc. You can do this either using various algorithms (for instance, those used by %(anvi-cluster-contigs)s) or manually through the interactive interface (%(anvi-interactive)s). In this type of binning, a lot of your bins will be complete metagenome-assembled genomes, or MAGs; however, if you find an interesting group of contigs (for example, a prophage or a plasmid or even a particular domain), you can also put that into a bin. Then, you can group these bins in different ways using %(collection)ss.

For more information, you might want to watch [this lovely lecture on genome-resolved metagenomics](https://www.youtube.com/watch?v=RjNdHGK4ruo).

### Pangenomic Workflows
You can also use bins to group together gene clusters. This is useful if you want a specific group of contigs to remain together through your entire analysis. Just provide your %(internal-genomes)s file to %(anvi-gen-genomes-storage)s.

Wow, this binning thing seems BINcredible! (not sorry)
