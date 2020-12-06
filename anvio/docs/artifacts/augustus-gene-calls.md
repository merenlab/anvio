A gene call file from [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/). 

[AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) is a tool to predict genes from a variety of Eurkaryotic genomes. This includes predicting the 5' UTR and 3' UTR, as well as introns. You can search a sequence in the [Augustus web interface](http://bioinf.uni-greifswald.de/augustus/submission.php). After a search, you can export the results as a `.gff` text file.  

{:.notice}
As of now, Anvi'o (specifically %(anvi-script-augustus-output-to-external-gene-calls)s) is only tested with AUGUSTUS v3.3.3. Feel free to be adventurous and try other versions if you feel so inclined. 

You can convert this file into an anvi'o %(external-gene-calls)s file using %(anvi-script-augustus-output-to-external-gene-calls)s. 

Here is an example of a `.gff` file for the [Homo sapiens RNAP III subunit D sequence](https://www.ncbi.nlm.nih.gov/nuccore/NM_001722.3?report=fasta): 

    # This output was generated with AUGUSTUS (version 3.3.3).
    # AUGUSTUS is a gene prediction tool written by M. Stanke (mario.stanke@uni-greifswald.de),
    # O. Keller, S. KÃ¶nig, L. Gerischer, L. Romoth and Katharina Hoff.
    # Please cite: Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008),
    # Using native and syntenically mapped cDNA alignments to improve de novo gene finding
    # Bioinformatics 24: 637-644, doi 10.1093/bioinformatics/btn013
    # No extrinsic information on sequences given.
    # Initializing the parameters using config directory /data/www/augustus/augustus/config/ ...
    # human version. Using default transition matrix.
    # Looks like /data/www/augustus/webservice/data/AUG-707407769/input.fa is in fasta format.
    # We have hints for 0 sequences and for 0 of the sequences in the input set.
    #
    # ----- prediction on sequence number 1 (length = 5336, name = unnamed-1) -----
    #
    # Predicted genes for sequence number 1 on both strands
    # start gene g1
    unnamed-1    AUGUSTUS    gene    57    1253    1    +    .    g1
    unnamed-1    AUGUSTUS    transcript    57    1253    1    +    .    g1.t1
    unnamed-1    AUGUSTUS    start_codon    57    59    .    +    0    transcript_id "g1.t1"; gene_id "g1";
    unnamed-1    AUGUSTUS    single    57    1253    1    +    0    transcript_id "g1.t1"; gene_id "g1";
    unnamed-1    AUGUSTUS    CDS    57    1253    1    +    0    transcript_id "g1.t1"; gene_id "g1";
    unnamed-1    AUGUSTUS    stop_codon    1251    1253    .    +    0    transcript_id "g1.t1"; gene_id "g1";
    # coding sequence = [atgtcggaaggaaacgccgccggcgagcccagcacgccgggagggccccgacctctcctgactggggcccgggggctca
    # tcgggcggcggccggcgcctcccctcacccccggccgccttccctccatccgttccagggacctcaccctcgggggagtcaagaagaaaaccttcacc
    # ccaaatatcatcagtcggaagatcaaggaagagcccaaggaagaagtaactgtcaagaaggagaagcgtgaaagggacagagaccgacaacgagaggg
    # gcatggacgagggcgaggccgtccagaagtgatccagtctcactccatctttgagcagggcccagctgaaatgatgaagaaaaaagggaactgggata
    # agacagtggatgtgtcagacatgggaccttctcatatcatcaacatcaaaaaagagaagagagagacagacgaagaaactaaacagatcttgcgtatg
    # ctggagaaggacgatttcctcgatgaccccggcctgaggaacgacactcgaaatatgcctgtgcagctgccgctggctcactcaggatggctttttaa
    # ggaagaaaatgacgaaccagatgttaaaccttggctggctggccccaaggaagaggacatggaggtggacatacctgctgtgaaagtgaaagaggagc
    # cacgagatgaggaggaagaggccaagatgaaggctcctcccaaagcagccaggaagactccaggcctcccgaaggatgtatctgtggcagagctgctg
    # agggagctgagcctcaccaaggaagaggaactgctgtttctgcagctgccagacaccctccctggccagccacccacccaggacatcaagcctatcaa
    # gacagaggtgcagggcgaggacggacaggtggtgctcatcaagcaggagaaagaccgagaagccaaattggcagagaatgcttgtaccctggctgacc
    # tgacagagggtcaggttggcaagctactcatccgcaagtctggaagggtgcaactcctcttgggcaaggtgactctggacgtgaccatgggaactgcc
    # tgctccttcctgcaggagctggtgtccgtgggccttggagacagtaggacaggggagatgacagtcctgggacacgtgaagcacaaacttgtatgttc
    # ccctgattttgaatccctcttggatcacaaacaccggtaa]
    # protein sequence = [MSEGNAAGEPSTPGGPRPLLTGARGLIGRRPAPPLTPGRLPSIRSRDLTLGGVKKKTFTPNIISRKIKEEPKEEVTVK
    # KEKRERDRDRQREGHGRGRGRPEVIQSHSIFEQGPAEMMKKKGNWDKTVDVSDMGPSHIINIKKEKRETDEETKQILRMLEKDDFLDDPGLRNDTRNM
    # PVQLPLAHSGWLFKEENDEPDVKPWLAGPKEEDMEVDIPAVKVKEEPRDEEEEAKMKAPPKAARKTPGLPKDVSVAELLRELSLTKEEELLFLQLPDT
    # LPGQPPTQDIKPIKTEVQGEDGQVVLIKQEKDREAKLAENACTLADLTEGQVGKLLIRKSGRVQLLLGKVTLDVTMGTACSFLQELVSVGLGDSRTGE
    # MTVLGHVKHKLVCSPDFESLLDHKHR]
    # end gene g1
    ###
    # command line:
    # /data/www/augustus/augustus/bin/augustus --species=human --strand=both --singlestrand=false --genemodel=partial --codingseq=on --sample=100 --keep_viterbi=true --alternatives-from-sampling=true --minexonintronprob=0.2 --minmeanexonintronprob=0.5 --maxtracks=2 /data/www/augustus/webservice/data/AUG-707407769/input.fa --exonnames=on

