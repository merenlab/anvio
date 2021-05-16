This program takes a %(contigs-fasta)s and %(blast-table)s and removes sequences without BLAST hits of a certain level of confidence. 

For example, you could use this program to filter out sequences that do not have high-confidence taxonomy assignments before running a phylogenomic analysis. 

To run this program, you'll need to provide the %(contigs-fasta)s that you're planning to filter, the %(blast-table)s, a list of the column headers in your %(blast-table)s (as given to BLAST by `-outfmt`), and a `proper_pident` threshold at which to remove the sequences. This threshold will remove sequences less than the given percent of the query amino acids that were identical to the corresponding matched amino acids. Note that this diffres from the `pident` blast parameter because it doesn't include unaligned regions. 

For example, if you ran 

{{ codestart }}
anvi-script-filter-fasta-by-blast -f %(contigs-fasta)s \
                                  -o path/to/%(contigs-fasta)s \
                                  -b %(blast-table)s \
                                  -s qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
                                  -t 30
{{ codestop }}
        
Then the output file would be a %(contigs-fasta)s that contains only the sequences in your input file that have a hit in your blast table with more than 30 percent of the amino acids aligned. 
