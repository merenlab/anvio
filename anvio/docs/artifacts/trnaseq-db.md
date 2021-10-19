A tRNA-seq database **contains information on tRNA sequences predicted from a single tRNA-seq sample**.

This database is the key output of **%(anvi-trnaseq)s**. That program predicts which reads are tRNA through structural profiling, clusters tRNA reads into discrete biological sequences, and predicts the positions of nucleotide modifications.

The series of steps implemented in %(anvi-trnaseq)s sequentially adds the following information to the database.

* Unique sequences predicted to be tRNA, including read counts
* Primary sequence and secondary structural features (stems and loops) predicted in each profiled tRNA
* Unconserved nucleotides in the primary sequence that differ from expectation
* Unpaired nucleotides in the stems
* "Trimmed" tRNA sequences, formed from unique sequences only differing by 3' nucleotides of the CCA acceptor region and 5' nucleotides beyond the acceptor stem
* "Normalized" tRNA sequences, formed by dereplicating trimmed tRNA sequences that are 3' fragments from incomplete reverse transcription and by mapping biological 5' and interior tRNA fragments
* Potentially modified tRNA sequences, formed by clustering normalized tRNA sequences and retaining those clusters that differ by 3-4 nucleotides at aligned positions

This database is the key input to **%(anvi-merge-trnaseq)s**, which takes one or more databases comprising the samples in an experiment and generates a %(trnaseq-contigs-db)s of tRNA seed sequences and %(trnaseq-profile-db)ss. These tRNA-seq variant contigs and profile databases can then be manipulated and displayed in anvi'o like normal %(contigs-db)ss and %(profile-db)ss.
