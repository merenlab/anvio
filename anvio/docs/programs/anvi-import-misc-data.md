This program enables extending anvi'o projects with many kinds of **additional data**. Additional data will extend anvio' %(interactive)s displays, and appear in %(summary)s files, and become accessible to other anvi'o programs thorughout.

This program can add additional data for your items or layers in a %(pan-db)s or %(profile-db)s, or add additional data for your nucleotides or amino acids in a %(contigs-db)s

You also have the option to associate keys with only a specific data group, or transpose the input before processing.

Also see the program %(anvi-show-misc-data)s, %(anvi-export-misc-data)s, and %(anvi-delete-misc-data)s.

## Items Data, Layers Data, and Orders

Please see [this blog post](http://merenlab.org/2017/12/11/additional-data-tables) for a comprehensive documentation on these misc data types.

## Nucleotides, Amino Acids, and Contigs Databases

This feature lets you import additional data about specfic residues or specific base pairs into your %(contigs-db)s. This is especially useful for strucutral analysis (so when running programs like %(anvi-display-structure)s) and will be very relevant to the InteracDome functionality when it's added in anvi'o v7 (curious readers can take a look at [this blog post](http://merenlab.org/2020/07/22/interacdome/)).

When adding additional data, unlike with layers and items, you do not have to provide values for every single nucleotide in your database. With this program, you can easily provide data for only a select few.

Basically, you can add two types of data to your contigs database:

1. %(misc-data-nucleotides)s by providing a %(misc-data-nucleotides-txt)s. This contains information about *specific nucleotides in your database.*

{{ codestart }}
anvi-import-misc-data -c %(contigs-db)s \
                      -t nucleotides \
                      %(misc-data-nucleotides-txt)s
{{ codestop }}

2. %(misc-data-amino-acids)s by providing a %(misc-data-amino-acids-txt)s. This contains information about *specific amino acid residues in your database*

{{ codestart }}
anvi-import-misc-data -c %(contigs-db)s \
                      -t amino_acids \
                      %(misc-data-amino-acids-txt)s
{{ codestop }}
