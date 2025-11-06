This file describes all of the gene calls contained in a %(contigs-db)s from a specified list of sources. It is the output of %(anvi-export-gene-calls)s. 

For each gene identified, this file provides various information, including the caller ID, [start and stop position](https://anvio.org/help/main/artifacts/external-gene-calls/#gene-startstop-positions), direction, whether or not the gene is partial, the [call type](https://anvio.org/help/main/artifacts/external-gene-calls/#call-type), source and version (if available ), and the amino acid sequence. 

{:.notice}
Want more information? This file is in the same format as an %(external-gene-calls)s, so check out that page. 

Here is an example from the Infant Gut Dataset: 

    gene_callers_id    contig              start    stop    direction    partial    call_type    source     version   aa_sequence
    0                  Day17a_QCcontig1    0        186     f            1          1            prodigal   v2.60     GSSPTAGVEQKQKPTWFLLFLFYSLFFDKLEEGTLKTFIRLKGSYRRMNTSNFSYGIMCLL
    1                  Day17a_QCcontig1    214      1219    f            0          1            prodigal   v2.60     MKILLYFEGEKILAKSGIGRALDHQKRALSEVGIEYTLDADCSDYDILHINTYGVNSHRMVRKARKLGKKVIYHAHSTEEDFRNSFIGSNQLAPLVKKYLISLYSKADHLITPTPYSKTLLEGYGIKVPISAISNGIDLSRFYPSEEKEQKFREYFKIDEEKKVIICVGLFFERKGITDFIEVARQLPEYQFIWFGDTPMYSIPKNIRQLVKEDHPENVIFPGYIKGDVIEGAYAAANLFFFPSREETEGIVVLEALASQQQVLVRDIPVYQGWLVANENCYMGHSIEEFKKYIEGLLEGKIPSTREAGYQVAEQRSIKQIGYELKEVYETVLS
    2                  Day17a_QCcontig1    1265     2489    f            0          1            prodigal   v2.60     MKIGFFTDTYFPQVSGVATSIKTLKDELEKHGHEVYIFTTTDPNATDFEEDVIRMPSVPFVSFKDRRVVVRGMWYAYLIAKELELDLIHTHTEFGAGILGKMVGKKMKIPVIHTYHTMYEDYLHYIAKGKVVRPSHVKFFSRVFTNHTTGVVCPSERVIEKLRDYGVTAPMRIIPTGIEIDKFLRPDITEEMIAGMRQQLGIEEQQIMLLSLSRISYEKNIQAIIQGLPQVIEKLPQTRLVIVGNGPYLEDLKELAEELEVSEYVQFTGEVPNEEVAIYYKAADYFVSASTSETQGLTYTEAMAAGVQCVAEGNAYLNNLFDHESLGKTFKTDSDFAPTLIDYIQANIKMDQTILDEKLFEISSTNFGNKMIEFYQDTLIYFDQLQMEKENADSIKKIKVKFTSLRK
    ...

