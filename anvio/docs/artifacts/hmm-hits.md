The search results for an %(hmm-source)s in a %(contigs-db)s. Essentially, this is the part of a %(contigs-db)s that handles the HMM data. In anvi'o, this is usually functional annotations, such as identifying specfic ribosomal RNAs, various single-copy core genes, and transfer RNAs, though the user can also define their own HMM sources. 

Upon creation, a %(contigs-db)s will not contain any HMM results. In order to populate it, users can run %(anvi-run-hmms)s using any %(hmm-source)s. The program %(anvi-scan-trnas)s also populates a %(contigs-db)s's hmm-hits with potential tranfer RNA hits.
