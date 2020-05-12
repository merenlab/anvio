### Converting a FASTA file to a contigs FASTA

{{ codestart }}
anvi-script-reformat-fasta %(fasta)s \
                           -o %(contigs-fasta)s \
                           --simplify-names
{{ codestop }}

{:.notice}
If you use the flag *--report-file*, it will also create a TAB-delimited file for you to keep track of which defline in the new file corresponds to which defline in the original file.

### Removing short reads from FASTA

Removing short contigs from a FASTA file will improve the performance of the %(contigs-db)s later. Running the same command this way will also remove sequences that are shorter than 1,000 nts:

{{ codestart }}
anvi-script-reformat-fasta %(fasta)s \
                           -o %(contigs-fasta)s \
                           -l 1000 \
                           --simplify-names
{{ codestop }}

