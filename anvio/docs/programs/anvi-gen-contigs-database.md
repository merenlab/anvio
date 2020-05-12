### Create a contigs database from a FASTA file

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s
{{ codestop }}

A FASTA file that contains one or more sequences. These sequences may belong
to a single genome, or could be contigs obtained from an assembly.

### Create a contigs database with external gene calls

{{ codestart }}
anvi-gen-contigs-database -f %(contigs-fasta)s \
                          -o %(contigs-db)s \
                          -e %(external-gene-calls)s
{{ codestop }}

