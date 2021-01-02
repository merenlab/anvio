A tRNA-seq database is an anvi'o database that **contains tRNA sequences predicted from a %(trnaseq-fasta)s for a sample, and associated information**.

This database is the key output of **%(anvi-trnaseq)s**, which predicts which input sequences are tRNA, clusters them into discrete biological sequences, and predicts positions in the sequences that are the sites of nucleotide modifications. The database therefore contains tables of information produced by each part of the process.

* Unique sequences predicted to be tRNA, including read counts
* Primary sequence and secondary structural features (stems and loops) predicted in each unique tRNA
* Unconserved nucleotides in the primary sequence that differ from expectation
* Unpaired nucleotides in the stems
* "Trimmed" tRNA sequences, formed from unique sequences only differing by 5' nucleotides beyond the acceptor stem and 3' nucleotides of the CCA acceptor region
* "Normalized" tRNA sequences, formed by dereplicating trimmed tRNA sequences that are 3' fragments produced by incomplete reverse transcription and by mapping biological tRNA fragments
* Potentially modified tRNA sequences, formed by clustering normalized tRNA sequences and retaining those clusters that differ by 3-4 nucleotides at potentially modified positions

This database is the key input to **%(anvi-convert-trnaseq-database)s**, which takes one or more databases comprising the samples of an experiment and generates a %(contigs-db)s of tRNA seed sequences and %(profile-db)s. These can then be displayed and manipulated in anvi'o like other 'omics data.