An anvi'o database that **contains key information associated with your sequences**.

In a way, **an anvi'o contigs database is a modern, more talented form of a FASTA file**, where you can store additional information about your sequences in it and others can query and use it. Information storage and access is primarily done by anvi'o programs, however, it can also be done through the command line interface or programmatically. 

The information a contigs database contains about its sequences include the positions of open reading frames, tetra-nucleotide frequencies, functional and taxonomic annotations, information on individual nucleotide or amino acid positions, and more.

Key programs that populate an anvi'o contigs database with essential information include %(anvi-run-hmms)s, %(anvi-run-scg-taxonomy)s, and %(anvi-scan-trnas)s.

Once an anvi'o contigs database is generated and populated with information, it is **always a good idea to run %(anvi-display-contigs-stats)s** to see a numerical summary of its contents.

Other essential programs that read from a contigs database and yield key information include %(anvi-estimate-genome-completeness)s, %(anvi-get-sequences-for-hmm-hits)s, %(anvi-estimate-scg-taxonomy)s.