This program calculates the completeness and redundancy of single-copy gene collections available in your %(contigs-db)s. 

Single-copy core genes (SCGs) are genes that are expected to occur exactly once in most genomes within a taxonomic group. For example, the default collections include single-copy core gene sets named `Protista_83`, `Archaea_76`, and `Bacteria_71`. This program provides rough estimates of how many Protist, Archaeal, and Bacterial genomes are present in your dataset by analyzing the occurrence patterns of these single-copy core genes. 

To list all available completeness sources in your %(contigs-db)s, use the following command:

{{ codestart }}
anvi-compute-completeness -c %(contigs-db)s \
                          --list-completeness-sources
{{ codestop }}
                              
You can then execute this program on a specific source as follows:

{{ codestart }}
anvi-compute-completeness -c %(contigs-db)s \
                          --completeness-source Bacteria_71
{{ codestop }}
                              
Additional options include providing a %(splits-txt)s file to focus the analysis on a specific subset of splits, or specifying a minimum e-value threshold for a gene to be counted as a valid hit. The default e-value threshold is `1e-15`.
