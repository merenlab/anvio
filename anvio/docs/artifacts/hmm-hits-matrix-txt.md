This file is the output of %(anvi-script-gen-hmm-hits-matrix-across-genomes)s and describes the %(hmm-hits)s across multiple genomes or bins for a single %(hmm-source)s. 

The first column describes each of the genomes (if the input was an %(external-genomes)s) or bins (if the input was an %(internal-genomes)s) that the matrix describes. The following columns describe each of the genes in your %(hmm-source)s. The data within the table describes the number of hits that gene had in that genome or bin. 

For example, if you were to run %(anvi-script-gen-hmm-hits-matrix-across-genomes)s with the `Bacteria_71` %(hmm-source)s on two hypothetical genomes, you would get a file like this:

    genome_or_bin    ADK    AICARFT_IMPCHas    ATP-synt    ATP-synt_A    Adenylsucc_synt    Chorismate_synt    EF_TS    ...
    Genome_1         11     10                 9           9             11                 8                  9        ...
    Genome_2         2      1                  1           2             3                  2                  2        ...
