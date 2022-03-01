A contigs database is an anvi'o database that **contains key information associated with your sequences**.

In a way, **an anvi'o contigs database is a modern, more talented form of a FASTA file**, where you can store additional information about your sequences in it and others can query and use it. Information storage and access is primarily done by anvi'o programs, however, it can also be done through the command line interface or programmatically.

The information a contigs database contains about its sequences can include the positions of open reading frames, tetra-nucleotide frequencies, functional and taxonomic annotations, information on individual nucleotide or amino acid positions, and more.

### Another (less computation-heavy) way of thinking about it

When working in anvi'o, you'll need to be able to access previous analysis done on a genome or transcriptome. To do this, anvi'o uses tools like contigs databases instead of regular fasta files. So, you'll want to convert the data that you have into a contigs database to use other anvi'o programs (using %(anvi-gen-contigs-database)s). As seen on the page for %(metagenomes)s, you can then use this contigs database instead of your fasta file for all of your anvi'o needs.

In short, to get the most out of your data in anvi'o, you'll want to use your data (which was probably originally in a %(fasta)s file) to create both a %(contigs-db)s and a %(profile-db)s. That way, anvi'o is able to keep track of many different kinds of analysis and you can easily interact with other anvi'o programs.

## Usage Information

### Creating and populating a contigs database

Contigs databases will be initialized using **%(anvi-gen-contigs-database)s** using a %(contigs-fasta)s. This will compute the k-mer frequencies for each contig, soft-split your contigs, and identify open reading frames. To populate a contigs database with more information, you can then run various other programs.

**Key programs that populate an anvi'o contigs database with essential information** include,

* %(anvi-run-hmms)s (which uses HMMs to annotate your genes against an %(hmm-source)s)
* %(anvi-run-scg-taxonomy)s (which associates its single-copy core gene with taxonomic data)
* %(anvi-scan-trnas)s (which identifies the tRNA genes)
* %(anvi-run-ncbi-cogs)s (which tries to assign functions to your genes using the COGs database)

Once an anvi'o contigs database is generated and populated with information, it is **always a good idea to run %(anvi-display-contigs-stats)s** to see a numerical summary of its contents.

Other programs you can run to populate a contigs database with functions include,

* %(anvi-run-kegg-kofams)s (which annotates the genes in the database with the KEGG KOfam database)

### Analysis on a populated contigs database

Other essential programs that read from a contigs database and yield key information include %(anvi-estimate-genome-completeness)s, %(anvi-get-sequences-for-hmm-hits)s, and %(anvi-estimate-scg-taxonomy)s.

If you wish to run programs like %(anvi-cluster-contigs)s, %(anvi-estimate-metabolism)s, and %(anvi-gen-gene-level-stats-databases)s, or view your database with %(anvi-interactive)s, you'll need to first use your contigs database to create a %(profile-db)s.

## Variants

Contigs databases, like %(profile-db)ss, are allowed have different variants, though the only currently implemented variant, the %(trnaseq-contigs-db)s, is for tRNA transcripts from tRNA-seq experiments. The default variant stored for "standard" contigs databases is `unknown`. Variants should indicate that substantially different information is stored in the database. For instance, open reading frames are applicable to protein-coding genes but not tRNA transcripts, so ORF data is not recorded for the `trnaseq` variant. The $(trnaseq-workflow)s generates %(trnaseq-contigs-db)ss using a very different approach to %(anvi-gen-contigs-database)s.

## For programmers

Tips and use cases for programmers. Send us your questions so we can extend this section with useful examples.

### Get number of approximate number of genomes

You can get the number of genomes once %(anvi-run-hmms)s is run on an contigs database. Here are some examples:

``` python
from anvio.hmmops import NumGenomesEstimator

# the raw data, where each key is one of the HMM collections
# of type `singlecopy` run on the contigs-db
NumGenomesEstimator('CONTIGS.db').estimates_dict
>>> {'Bacteria_71': {'num_genomes': 9, 'domain': 'bacteria'},
     'Archaea_76': {'num_genomes': 1, 'domain': 'archaea'},
     'Protista_83': {'num_genomes': 1, 'domain': 'eukarya'}}

# slightly fancier output with a single integer for
# estimated number of genomes summarized, along with
# domains used
num_genomes, domains_included = NumGenomesEstimator('CONTIGS.db').num_genomes()
print(num_genomes)
>>> 11

print(domains_included)
>>> ['bacteria', 'archaea', 'eukarya']

# limiting the domains
num_genomes, domains_included = NumGenomesEstimator('CONTIGS.db').num_genomes(for_domains=['archaea', 'eukarya'])
print(num_genomes)
>>> 2

print(domains_included)
>>> ['archaea', 'eukarya']
```
