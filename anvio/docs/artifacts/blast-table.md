This describes the BLAST table that is outputted when you run [Protein BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) from the terminal. 

When given to %(anvi-script-filter-fasta-by-blast)s, which is currently the only program that uses this artifact, it expects output form 6. By default, this incldues the following data columns: 

    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen
    
However, you'll have to provide the columns in your file and their order to the program wirth the flag `--outfmt`. For the program to work properly, your table must at least include the columns `qseqid`, `bitscore`, `length`, `qlen`, and `pident`.
