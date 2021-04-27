This program tells you the completeness and redundency of single-copy gene sources available for your %(contigs-db)s. 

For example, some of the defaults are collections of single-copy core genes named  `Protista_83`, `Archaea_76`, and `Bacteria_71`. This program will give you a rough estimate of how many Protist, Archaeal, and Bacterial genomes are included in your dataset using these single-copy core genes. 

You can use the following run to list available completeness sources in your %(contigs-db)s:

{{ codestart }}
anvi-compute-completeness -c %(contigs-db)s \
                          --list-completeness-sources
{{ codestop }}
                              
Then you can run this program on a specifc source as folows:

{{ codestart }}
anvi-compute-completeness -c %(contigs-db)s \
                          --completeness-source Bacteria_71
{{ codestop }}
                              
You can also provide a %(splits-txt)s to focus on a specific set of splits, or declare a minimum e-value for a gene to count as a hit. The default is `1e-15`.
