An %(ngrams)s object is a DataFrame that contains count data of synteny patterns collected from a group of similar loci or genomes. It is produced by running %(anvi-analyze-synteny)s when given a %(genomes-storage-db)s and an annotation source.

An `ngram` is a group of neighboring genes that include precisely `n` genes, inspired by the term ngram in [linguistics and natural language processing](https://en.wikipedia.org/wiki/N-gram). This object was inspired by kmer count tables but is inherently different because it is counting adjacent genes and not nucleotides.
