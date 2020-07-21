A FASTA-formatted file that does not necessarily meet the standards of a %(contigs-fasta)s.

%(anvi-script-reformat-fasta)s can turn a regular fasta into a %(contigs-fasta)s, which anvi'o will be able to utilize better.

### But what is a FASTA file? 

A FASTA file contains sequences (in this case, nucleotide sequences, though they can also describe peptide sequences) that are formatted as follows: 

    >SEQUENCE_ID VARIOUS_SEQUENCE_DATA
    SEQUENCE
    
The `VARIOUS_SEQUENCE_DATA` region can contain data such as the NCBI taxon ID, GI accession number, a text description of the sequence, or the start and end positions if the sequence is a portion of a larger sample. All of this information is optional. 

The sequence itself is written in standard IUPAC format (though it can be written in lower-case letters).  

For a concrete example, you can download sequences from the NCBI database in FASTA format. 
